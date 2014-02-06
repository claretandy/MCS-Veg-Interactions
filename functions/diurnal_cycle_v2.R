library(raster)
library(rasterVis)
library(rgdal)
library(sm)
library(plotrix)
library(RColorBrewer)
library(classInt)
library(reshape)
library(ggplot2)
library(scales)
library(plyr)
library(Hmisc)


# Produce line plots of Precip vs time for different veg types 
# Time: Both 16/08 0Z to 19/08 6Z and averaged per hour
# Precip: Both total precip and precip split by MCS/POP # This has not been done yet!!!
# Type: Specify whether we are plotting intense precip or all precip over veg types
diurnalcycle2 <- function(inbr, type="intense", patch=F, model.nm="avg5216.4km.std", id="s", spp.r=spp.r, sppa=sppa, mycl=mycl, land_simple, overwrite=F){
	print("Plotting diurnal cycle ...")
	
	if (Sys.info()[["sysname"]] == "Darwin"){
		indatadir <- "/Users/ajh235/Work/DataLocal/ModelData/WAFR/"
		resultsdir <- "/Users/ajh235/Work/Projects/InternalSabbatical/Results/"
		scratchdir <- "/Users/ajh235/Work/Scratch/"
	} else {
		indatadir <- "/data/local/hadhy/ModelData/WAFR/"
		resultsdir <- "/home/h02/hadhy/Projects/InternalSabbatical/Results/"
		scratchdir <- "/data/local/hadhy/Scratch/"
		require(PP,lib.loc="/project/ukmo/rhel6/R")
	}
	rasterOptions(tmpdir=scratchdir, todisk=F)
	# Set indir and outpdf
# 	datadir <- paste("/data/local/hadhy/ModelData/WAFR/djzx",id,"/",sep="") 
	outpdf <- paste(resultsdir,"timeseries_byzone-stdveg_",id,"_",type,"_v5.pdf", sep="")
	
	# Load masks etc ...
# 	load("~/Projects/InternalSabbatical/Scripts/myWorkspace_29082012.Rdata")
# 	# Load land boundaries
# 	land_simple <- readOGR(dsn="/data/local/hadhy/WAFR3/ancils", layer="land_rp")
	# List box IDs
	myboxes <- 1:15
	
    # Recode "boundary tree" and "boundary grass" to "boundary"
    if (length(which(freq(mycl)[,1] %in% 5:6)) != 0){
        mycl[mycl==5] <- 4
        mycl[mycl==6] <- 4        
    }

    if (round(extent(mycl)@xmin + (extent(mycl)@xmax - extent(mycl)@xmin)/2) == 360){
        mycl <- shift(mycl, x=-360)
    }

	# Open PDF 
	pdf(outpdf, width=12, height=6)
	for (bx in myboxes){

		print(bx)
#         browser()
		mymsk <- spp.r == bx
		mymsk[mymsk == 0] <- NA
# 		browser()
		e <- extent(sppa[bx,])
		# Plot close-up of grid box vegetation ...
        print("Plot close-up of grid box vegetation ...")
        
        myLUT <- data.frame(ID=c(1,2,3,4,7), Landcover=factor(c("tree", "grass", "sparse", "boundary", "orography"), levels=c("tree", "grass", "sparse", "boundary", "orography")[c(3,2,4,1,5)]), Colours=c("dark green", "yellow", "orange", "brown", "dark grey"), plotOrder=c(4,2,1,3,5))

        mycl.bx <- crop(mycl, e)
        mycl.f <- as.factor(mycl.bx)
        fdat <- myLUT[myLUT$ID %in% levels(mycl.f)[[1]]$ID, ]
        levels(mycl.f) <- fdat #[rev(order(fdat$plotOrder)),]
        print(
            levelplot(mycl.f, att="Landcover", main=paste("Vegetation classes for Zone",bx), xlab=NULL, ylab=NULL, scales=list(draw=FALSE), col.regions=as.character(fdat$Colours)) + 
                latticeExtra::layer(sp.polygons(land_simple))
        )		
		
		# 1) From 16/08 0Z to 19/08 6Z by total precip
        print("pr.bx ...")
        bx.e <- extent(spp[bx,])
		pr.bx <- crop(inbr, bx.e, filename=paste(scratchdir,"pr.bx.tif",sep=""), overwrite=T)
        pr.bx <- calc(pr.bx, fun=function(x){return(x*3600)})
		mycl.bx <- crop(shift(mycl, x=360), bx.e)
		
        if (type == "all"){
            mystats <- zonal(pr.bx, mycl.bx, 'mean')
            colnames(mystats)[-1] <- format(getZ(inbr), "%d%b%H%M")
            mystats.df <- data.frame(mystats)
            mystats.df$Landcover <- myLUT[which(myLUT$ID %in% mystats.df$zone),"Landcover"]
            meltedstats <- melt(mystats.df, id.vars=c("zone","Landcover"), measure.vars=2:475)
            meltedstats$time <- as.POSIXct(paste("2006_",as.character(meltedstats$variable),sep=""), format="%Y_X%d%b%H%M")
            meltedstats$variable <- NULL
            meltedstats$TOD <- format(meltedstats$time, "%H:%M")
            meltedstats$value[meltedstats$value == 0] <- NA
            clrs <- as.character(myLUT[which(myLUT$ID %in% mystats.df$zone),"Colours"])
            ordr <- as.numeric(myLUT[which(myLUT$ID %in% mystats.df$zone),"plotOrder"])
            
            # Plot time series
            print(
                ggplot(data=meltedstats, aes(x=time, y=value, group=Landcover, colour=Landcover)) + 
                    geom_line() + 
                    scale_colour_manual(values = clrs[order(ordr)]) + 
                    labs(title=paste("Mean Precipitation by Land Cover Type\nZone id:",bx), y=expression(mm~hour^-1), x="Time")
                )
            
            # Plot the time of day means
            tmp <- ddply(meltedstats, c("Landcover","TOD"), summarise, TODmean=mean(value, na.rm=T))
            tmp$TOD <- as.POSIXct(paste("2006-08-16",tmp$TOD), "2006-08-16 %H:%M")
            print(
                ggplot(tmp, aes(x=TOD, y=TODmean, group=Landcover, colour=Landcover)) +
                geom_line() +
                scale_x_datetime(breaks=date_breaks("2 hour"), minor_breaks=date_breaks("1 hour"), labels=date_format("%H:%M")) +
                scale_colour_manual(values = clrs[order(ordr)]) +
                labs(title=paste("Mean Precipitation When Raining by Land Cover Type\nZone id:",bx), y=expression(mm~hour^-1), x="Time of Day")

                )

        }

        if (type == "intense"){
            pr.int.bx <- calc(pr.bx, fun=function(x){x[x<10] <- NA; return(x)})
            mystats.int.mean <- zonal(pr.int.bx, mycl.bx, 'mean')
            mystats.int.sd <- zonal(pr.int.bx, mycl.bx, 'sd')
            #            
            colnames(mystats.int.mean)[-1] <- format(getZ(inbr), "%d%b%H%M")
            colnames(mystats.int.sd)[-1] <- format(getZ(inbr), "%d%b%H%M")
            mystats.int <- data.frame(mystats.int.mean)
            mystats.int$Landcover <- myLUT[which(myLUT$ID %in% mystats.int$zone),"Landcover"]
            meltedstats <- melt(mystats.int, id.vars=c("zone","Landcover"), measure.vars=2:475)
            meltedstats$time <- as.POSIXct(paste("2006_",as.character(meltedstats$variable),sep=""), format="%Y_X%d%b%H%M")
            meltedstats$variable <- NULL
            colnames(meltedstats) <- c("zone","Landcover","mean","time")
            meltedstats$TOD <- format(meltedstats$time, "%H:%M")
            clrs <- as.character(myLUT[which(myLUT$ID %in% mystats.int$zone),"Colours"])
            ordr <- as.numeric(myLUT[which(myLUT$ID %in% mystats.int$zone),"plotOrder"])
            
#             browser()
            
            print(
                ggplot(data=meltedstats, aes(x=time, y=mean, group=Landcover, colour=Landcover)) + 
                    geom_line() + 
                    scale_colour_manual(values = clrs[order(ordr)]) + 
                    labs(title=paste("Mean Intense Precipitation (>10mm per hour)\nby Land Cover Type in Zone",bx), y=expression(mm~hour^-1), x="Time")
                )
            
            # Plot time of day means of intense precip over different cover types
            tmp <- ddply(meltedstats, c("Landcover","TOD"), summarise, TODmean=mean(mean, na.rm=T))
            tmp$TOD <- as.POSIXct(paste("2006-08-16",tmp$TOD), "2006-08-16 %H:%M")
            print(
                ggplot(tmp, aes(x=TOD, y=TODmean, group=Landcover, colour=Landcover)) +
                    geom_line() +
                    scale_x_datetime(breaks=date_breaks("2 hour"), minor_breaks=date_breaks("1 hour"), labels=date_format("%H:%M")) +
                    scale_colour_manual(values = clrs[order(ordr)]) +
                    labs(title=paste("Mean Intense Precipitation (>10mm per hour) by Land Cover Type\nand Time of Day in Zone",bx), y=expression(mm~hour^-1), x="Time of Day")
                )
            
        }

	}
	dev.off()

}