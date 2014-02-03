# Need to run first ...
# source("functions/loadVeg.R")
# First 75 lines of loadData_diurnal.R
# Need to adapt script to be called from the above script ...
# 
# mcsIntensity <- function(results, myveg, mcs, inbr){
    
    source("functions/makeLines.R")
    library(ggplot2)
    library(reshape)
    
    adjCoords <- function(r){
        r <- shift(r, x=0.01798)
        res(r) <- 0.036
        r <- shift(r, y=0.018)
        return(r)
    }
    
    pdf("/Users/ajh235/Work/Projects/InternalSabbatical/Results/MCS_intensity_bysize_v3.pdf", width=14, height=9)
    mcsrst.path <- "/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/Tracking/djzxs_10min/"
    gntheme <- rasterTheme(pch=19, region=brewer.pal(9, 'YlGn'))
    
    # Load veg fractions
    if(round(extent(myveg)@xmax + extent(myveg)@xmin) != 0){
        myveg <- shift(myveg, x=-360)
    }
    myveg <- adjCoords(myveg)
    blfrac <- myveg[[1]]; blfrac[is.na(blfrac)] <- 0
    grassfrac <- myveg[[3]] + myveg[[4]]; grassfrac[is.na(grassfrac)] <- 0
    if (!exists("mycl")){
        allveg <- vegPrep(model.nm=models[1], id=id[1], myveg, myorog, mylandfrac, land_simple, sppa, spp.r, plots=F, vegThreshold=0.3, overwrite=F) # return(mycl, mycl.z) and creates pdf plots
        mycl <- allveg[[1]]
    }
    
    # Get unique list of MCS ids
    mcs.ids <- unique(results$ID)
    alllines <- makeLines(results)[[1]]
    
    # Plot lines with BL tree fraction in background
    print(
        levelplot(blfrac, margin=F, auto.key=list(title="Broadleaf tree fraction"), par.settings=gntheme, at=seq(0,1,0.1), xlim=c(-13,16), ylim=c(-9,9), main="MCS results (16/08 to 18/08)") + 
        latticeExtra::layer(sp.polygons(land_simple, col="dark grey")) + 
        latticeExtra::layer(sp.lines(alllines))
    )
    
    # Create a subset of all MCS, either by length of by area
    # mcs.length <- sort(table(results$ID))
    # iii <- as.numeric(dimnames(mcs.length[mcs.length > 100])[[1]]) # Get all MCS with length > 100
    id.pixcount <- aggregate(pixcount ~ ID, data=results, FUN=sum)
    iii <- id.pixcount[which(id.pixcount$pixcount > 20000),"ID"]
    spresults <- makeLines(results[results$ID %in% iii,])
    myspdf <- spresults[[1]] # Spatial Lines
    mytext <- spresults[[2]] # Spatial points of starting point
    
    # Get coordinates for starting points
    coordinates(mytext) <- ~x+y
    
    print(
        levelplot(blfrac, margin=F, par.settings=gntheme) + 
              latticeExtra::layer(sp.polygons(land_simple, col="grey")) + 
              latticeExtra::layer(sp.lines(myspdf)) + 
              latticeExtra::layer(sp.text(loc=coordinates(mytext), txt=mytext$ID)) 
        )
    
    # Loop through each track
    for (ii in iii){
        
        mydates <- results[results$ID == ii,"timestep"]
        if (exists("mydf")){ rm(mydf) }
        
        for (x in 1:length(mydates)){
            print(mydates[x])
            
            mcs.now <- shift(raster(paste(mcsrst.path,"mcs_tracking_1000km_",format(mydates[x], "%d.%H%M"), ".tif",sep="")), x=-360)
            mcs.now <- adjCoords(mcs.now)
            
            mcsiinow <- mcs.now == ii
            mcsiinow[mcsiinow == 0] <- NA
            precip.now <- shift(rb5216.4km.std[[which(getZ(rb5216.4km.std) == mydates[x])]], x=-360)
            precip.now <- adjCoords(precip.now)
            
            # Mask precip that's not in the MCS
            mcsprecip   <- mask(precip.now, mcsiinow)
            precip.mask <- getValues(mcsprecip)
            precip.mask <- precip.mask[!is.na(precip.mask)]
            
            # Mask blfrac that's not in the MCS
            mcstree     <- mask(blfrac, mcsiinow)
            blfrac.mask <- getValues(mcstree)
            blfrac.mask <- blfrac.mask[!is.na(blfrac.mask)]
            
            # Mask grassfrac that's not in the MCS
            mcsgrass       <- mask(grassfrac, mcsiinow)
            grassfrac.mask <- getValues(mcsgrass)
            grassfrac.mask <- grassfrac.mask[!is.na(grassfrac.mask)]
            
            # What is the max precip rate, and what are veg fractions at that location?
            maxprecip <- cellStats(mcsprecip, 'max')*60*60
            mcspoints <- rasterToPoints(mcsprecip, spatial=FALSE)
            maxpoint  <- matrix(mcspoints[which(mcspoints[,3]==maxprecip/3600), 1:2], byrow=T, nrow=1, ncol=2)
            maxgrass  <- extract(grassfrac, maxpoint)
            maxtree   <- extract(shift(blfrac, x=-360), maxpoint)
            
            # Where the rainfall is most intense, what's the rate over grass and the rate over tree?
            # 1. Classify intense rain (>20mm)
            mcsprecip.intense <- calc(mcsprecip, fun=function(x){x[x<(20/3600)]<-NA; return(x)})
            
            # 2. Clip intense rain to forest or non-forest class 
            # forest: (tree >= 0.3) and non-forest: (tree < 0.3)
            browser()
            tree.intprecip <- mcsprecip.intense * 3600
            grass.intprecip <- mcsprecip.intense * 3600
            grass.intprecip[mcstree >= 0.3] <- NA
            tree.intprecip[mcstree < 0.3] <- NA
            
            # 3. Create summary stats for intense rain under forest / non-forest
            tree.int.mean  <- cellStats(tree.intprecip, 'mean', na.rm=T)
            tree.int.sd    <- cellStats(tree.intprecip, 'sd', asSample=F, na.rm=T)
            grass.int.mean <- cellStats(grass.intprecip, 'mean', na.rm=T)
            grass.int.sd   <- cellStats(grass.intprecip, 'sd', asSample=F, na.rm=T)
            
            # 4. Add results to a data.frame
            if (!exists("mydf")){
                mydf <- data.frame(time=mydates[x], blfrac=mean(blfrac.mask) , grassfrac=mean(grassfrac.mask) , precip=mean(precip.mask)*60*60, maxprecip=maxprecip, maxgrass=maxgrass, maxtree=maxtree, treemean=tree.int.mean, treesdmin=tree.int.mean-tree.int.sd, treesdmax=tree.int.mean+tree.int.sd, grassmean=grass.int.mean, grasssdmin=grass.int.mean-grass.int.sd, grasssdmax=grass.int.mean+grass.int.sd, areasqkm=sum(getValues(mcsiinow), na.rm=T)*4)
            } else {
                mydf <- rbind(mydf, data.frame(time=mydates[x], blfrac=mean(blfrac.mask) , grassfrac=mean(grassfrac.mask) , precip=mean(precip.mask)*60*60, maxprecip=maxprecip, maxgrass=maxgrass, maxtree=maxtree, treemean=tree.int.mean, treesdmin=tree.int.mean-tree.int.sd, treesdmax=tree.int.mean+tree.int.sd, grassmean=grass.int.mean, grasssdmin=grass.int.mean-grass.int.sd, grasssdmax=grass.int.mean+grass.int.sd, areasqkm=sum(getValues(mcsiinow), na.rm=T)*4))
            }
            
        }
        
        # 5. Reshape data.frame, and plot using ggplot ...        
        mydf4gg <- melt(mydf, id.vars=1)
        j <- which(mydf4gg$variable %in% c("precip","blfrac","grassfrac","areasqkm"))
        f1 <- ggplot(mydf4gg[j,], aes(time, value, ymin=0, ymax=value)) + scale_colour_identity() + facet_grid(variable ~ ., scales="free") + labs(title = paste("MCS id:",ii))
        p1 <- f1 + geom_line()
        print(p1)
        
        # Plot mean precip intensity over forest and grass
        treedf <- data.frame(mydf[,c("time","treemean","treesdmin","treesdmax")], landcover="tree")
        grassdf <- data.frame(mydf[,c("time","grassmean","grasssdmin","grasssdmax")], landcover="grass")
        colnames(treedf) <- c("time", "mean", "sdmin", "sdmax", "landcover")
        colnames(grassdf) <- c("time", "mean", "sdmin", "sdmax", "landcover")
        ribbondf <- rbind(treedf, grassdf)
        f2 <- ggplot(ribbondf, aes(x=time, group=landcover, fill=landcover))
        p2 <- f2 + geom_ribbon(aes(ymin=sdmin, ymax=sdmax), alpha=0.3) + geom_line(aes(y=mean, colour=landcover)) + scale_colour_manual(values = c("dark green","#fe9929")) + scale_fill_manual(values = c("dark green","#fed98e")) + labs(title=paste("Mean Precipitation Intensity by Land Cover Type\nMCS id:",ii), y=expression(mm~hour^-1), x="Time")
        print(p2)
        
        
        rm(mydf)
        rm(mydf4gg)
        rm(ribbondf)
        
    }
    
    dev.off()
# }
