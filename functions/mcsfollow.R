require(raster)
require(rasterVis)
# require(PP,lib.loc="/project/ukmo/rhel6/R")
require(rgdal)
require(sm)
require(plotrix)
require(RColorBrewer)
require(classInt)

mcsfollow <- function(mcs=mcs, precip=precip, myveg=myveg, id="s", threshold=1000, timestep="10min", indatadir=indatadir, dlresultsdir=dlresultsdir, resultsdir=resultsdir, overwrite=F){
    
    load(paste(dlresultsdir,"Tracking/djzx",id,"_10min.RData",sep=""))
    mcs.bylength <- sort(tapply(results$pixcount, INDEX=as.factor(results$ID), FUN=function(x){length(x)}))
    mcs.bypixsum <- sort(tapply(results$pixcount, INDEX=as.factor(results$ID), FUN=function(x){sum(x)}))
    mcs.byduratn <- sort(tapply(results$timestep, INDEX=as.factor(results$ID), FUN=function(x){length(x)}))
    mcs.ids <- as.numeric(dimnames(mcs.bypixsum)[[1]])
    mydates <- getZ(precip)
    
    # Subset only BL forest for plotting
    my.bl <- subset(shift(myveg, x=-360), 1)
    # Set green theme for plotting
    gntheme <- rasterTheme(pch=19, region=brewer.pal(9, 'YlGn'))
    
    # Plot all MCS tracks for the period
    pdf(paste(resultsdir,"MCS_tracks_djzx",id,".pdf",sep=""), width=12, height=8, onefile=T)
    my.plot <- levelplot(my.bl, margin=F, colorkey=list(space='bottom'), par.settings=gntheme, at=seq(0,1,0.1), xlim=c(-13,16), ylim=c(-9,9), main="MCS tracks (16/08 to 18/08)") + layer(sp.polygons(land_simple, col="dark grey"))
    my.lines <- list()
    i <- 1
    for (ii in 1:length(mcs.ids)){
        print(mcs.ids[ii])
        my.xy <- results[results$ID == mcs.ids[ii],c("x","y")]
        if(dim(my.xy)[1] > 1){
            my.lines[[i]] <- Lines(Line(my.xy),ID=i)
            i <- i+1
        }
    }
    my.lines.sp <- SpatialLines(my.lines)
    my.plot <- my.plot + layer(sp.lines(my.lines.sp, col="black"))
    print(my.plot)
    dev.off()
    
    # Plot individual MCS tracks
    idtoplot <- 469
    my.plot <- levelplot(my.bl, margin=F, colorkey=list(space='bottom'), par.settings=gntheme, at=seq(0,1,0.1), xlim=c(-13,16), ylim=c(-9,9), main=paste("MCS tracks (16/08 to 18/08) for MCS",idtoplot)) + layer(sp.polygons(land_simple, col="dark grey"))
    my.xy <- results[results$ID == idtoplot, c("x","y")]
    my.lines <- list()
    my.lines[[1]] <- Lines(Line(my.xy),ID=idtoplot)
    my.lines.sp <- SpatialLines(my.lines)
    my.plot <- my.plot + layer(sp.lines(my.lines.sp, col="black"))
    print(my.plot)
    
    # Plot Sensible and Latent Heat timeseries for a point along the track
    sensible <- brick("/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/layers_nc/3217.nc")
    latent <- brick("/Users/ajh235/Work/DataLocal/ModelData/WAFR/djzxs/layers_nc/3234.nc")
    samp.pts <- SpatialPoints(data.frame("x"=c(9.82868917, 3.07832476, 1.27895533, -0.20784297)+360, "y"=c(3.021809166, -0.857591780, -0.671397601, -1.825403176)))
    print(my.plot + layer(sp.points(samp.pts, pch=16, cex=1)))
    sens.pts <- extract(sensible, samp.pts, method="simple")
    ltnt.pts <- extract(latent, samp.pts, method="simple")
    prec.pts <- extract(precip, samp.pts, method="simple")
    
    # Create data.frame for ggplot ...
    rm(mydf)
    for (x in 1:nrow(sens.pts)){ 
        if (!exists("mydf")){ 
            mydf <- data.frame(time=mydates, sensible=sens.pts[1,], latent=ltnt.pts[1,], precip=prec.pts[1,], pointID=1) 
        } else { 
            mydf <- rbind(mydf, data.frame(time=mydates, sensible=sens.pts[x,], latent=ltnt.pts[x,], precip=prec.pts[x,], pointID=x)) 
        }
    }
    # Now plot with ggplot ...
    # Sensible heat
    h <- ggplot(data=mydf, aes(x=time, y=sensible))
    h + geom_line() + facet_grid(pointID ~ .) + labs(title=paste("Sensible Heat along MCS",idtoplot), x="Time", y="Wm-2")
    # Latent heat
    h <- ggplot(data=mydf, aes(x=time, y=latent))
    h + geom_line() + facet_grid(pointID ~ .) + labs(title=paste("Latent Heat along MCS",idtoplot), x="Time", y="Wm-2")
    # Precipitation
    h <- ggplot(data=mydf, aes(x=time, y=sensible))
    h + geom_line() + facet_grid(pointID ~ .) + labs(title=paste("Precipitation along MCS",idtoplot), x="Time", y="mms-1")
    
    # Check the sensible and latent heat in front of and within selected MCS
    mydates <- getZ(precip)
    for (i in 2:(length(mydates)-1)){
        t.mydatefmt <- format(results[results$ID == idtoplot,"timestep"][i], "%d.%H%M")
        tm1.mydatefmt <- format(results[results$ID == idtoplot,"timestep"][i-1], "%d.%H%M")
        tp1.mydatefmt <- format(results[results$ID == idtoplot,"timestep"][i+1], "%d.%H%M")
        t.mcsids.f <- paste(dlresultsdir,"Tracking/djzx",id,"_",timestep,"/mcs_tracking_",threshold,"km_",t.mydatefmt,".tif",sep="")
        tp1.mcsids.f <- paste(dlresultsdir,"Tracking/djzx",id,"_",timestep,"/mcs_tracking_",threshold,"km_",tp1.mydatefmt,".tif",sep="")
        t.mcsids <- raster(t.mcsids.f)
        tp1.mcsids <- raster(tp1.mcsids.f)
        
        # Load Latent and Sensible heat at t and t+1
        t.mydatefmt <- format(results[results$ID == idtoplot,"timestep"][i], "%d%H.%M")
        tp1.mydatefmt <- format(results[results$ID == idtoplot,"timestep"][i+1], "%d%H.%M")
        t.sens <- raster(paste(indatadir,"djzx",id,"/layers_nc/3217.",t.mydatefmt,".nc",sep=""))
        tp1.sens <- raster(paste(indatadir,"djzx",id,"/layers_nc/3217.",tp1.mydatefmt,".nc",sep=""))
        t.sens <- raster(paste(indatadir,"djzx",id,"/layers_nc/3234.",t.mydatefmt,".nc",sep=""))
        tp1.sens <- raster(paste(indatadir,"djzx",id,"/layers_nc/3234.",tp1.mydatefmt,".nc",sep=""))
        
        # Recode IDS to pull out 
    }
    
    
}

# Load data for a given date ...
loadHt <- function(mydate, id, indatadir){
    
}