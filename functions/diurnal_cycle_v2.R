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
source("functions/getLUT.R")


# Produce line plots of Precip vs time for different veg types 
# Time: Both 16/08 0Z to 19/08 6Z and averaged per hour
# Precip: Both total precip and precip split by MCS/POP 
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
	outpdf <- paste(resultsdir,"timeseries_byzone-stdveg_",id,"_",type,"_v9.pdf", sep="")
	
	# Load masks etc ...
# 	load("~/Projects/InternalSabbatical/Scripts/myWorkspace_29082012.Rdata")
# 	# Load land boundaries
# 	land_simple <- readOGR(dsn="/data/local/hadhy/WAFR3/ancils", layer="land_rp")
	# List box IDs
	myboxes <- 1:15 #7
	
    # Recode "boundary tree" and "boundary grass" to "boundary"
    if (length(which(freq(mycl)[,1] %in% 5:6)) != 0){
        mycl[mycl==5] <- 4
        mycl[mycl==6] <- 4        
    }

    if (xmax(mycl) > 360){
        mycl <- shift(mycl, x=-360)
    }

	# Open PDF 
    pdf(outpdf, width=12, height=6)
    
    # Output dataframe with results in
    if (exists("boxdata")){rm(boxdata)} 

	for (bx in myboxes){

		print(bx)
#         browser()
		mymsk <- spp.r == bx
		mymsk[mymsk == 0] <- NA
# 		browser()
		e <- extent(sppa[bx,])
		# Plot close-up of grid box vegetation ...
        print("Plot close-up of grid box vegetation ...")
        
        myLUT <- getLUT(nBndClass=1)

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
        bx.e <- extent(sppa[bx,])
		pr.bx <- crop(inbr, bx.e, filename=paste(scratchdir,"pr.bx.tif",sep=""), overwrite=T)
        pr.bx <- calc(pr.bx, fun=function(x){return(x*3600)})
		mycl.bx <- crop(mycl, bx.e)
		
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
            # Calculate Zonal mean, sd, and se of intense precipitation
            pr.int.bx <- calc(pr.bx, fun=function(x){x[x<10] <- NA; return(x)})
            mystats.int.mean <- zonal(pr.int.bx, mycl.bx, 'mean')
            mystats.int.sd <- zonal(pr.int.bx, mycl.bx, 'sd')
            mystats.int.se <- zonal(pr.int.bx, mycl.bx, fun=function(x, na.rm=T){x <- x[!is.na(x)]; return(sd(x)/sqrt(length(x)))})
            #            
            
            # Add time information to field headers for later ...
            colnames(mystats.int.mean)[-1] <- format(getZ(inbr), "%d%b%H%M")
            colnames(mystats.int.sd)[-1] <- format(getZ(inbr), "%d%b%H%M")
            colnames(mystats.int.se)[-1] <- format(getZ(inbr), "%d%b%H%M")
            #            
            # Prepare for melting ...
            mystats.int <- data.frame(mystats.int.mean)
            mystats.int.sd <- data.frame(mystats.int.sd)
            mystats.int.se <- data.frame(mystats.int.se)
            mystats.int$Landcover <- myLUT[which(myLUT$ID %in% mystats.int$zone),"Landcover"]
            mystats.int.sd$Landcover <- myLUT[which(myLUT$ID %in% mystats.int.sd$zone),"Landcover"]
            mystats.int.se$Landcover <- myLUT[which(myLUT$ID %in% mystats.int.se$zone),"Landcover"]
            # Do the melt man ...
            meltedstats <- melt(mystats.int, id.vars=c("zone","Landcover"), measure.vars=2:475)
            meltedstats.sd <- melt(mystats.int.sd, id.vars=c("zone","Landcover"), measure.vars=2:475)
            meltedstats.se <- melt(mystats.int.se, id.vars=c("zone","Landcover"), measure.vars=2:475)
            # Tidy up and merge all melts together ...
            meltedstats$time <- as.POSIXct(paste("2006_",as.character(meltedstats$variable),sep=""), format="%Y_X%d%b%H%M")
            meltedstats$variable <- NULL
            meltedstats.sd$time <- as.POSIXct(paste("2006_",as.character(meltedstats.sd$variable),sep=""), format="%Y_X%d%b%H%M")
            meltedstats.sd$variable <- NULL
            meltedstats.se$time <- as.POSIXct(paste("2006_",as.character(meltedstats.se$variable),sep=""), format="%Y_X%d%b%H%M")
            meltedstats.se$variable <- NULL
            
            tmp <- merge(meltedstats, meltedstats.sd, by.x=c("zone","Landcover","time"), by.y=c("zone","Landcover","time"))
            meltedstats <- merge(tmp, meltedstats.se, by.x=c("zone","Landcover","time"), by.y=c("zone","Landcover","time"))
            
            colnames(meltedstats) <- c("zone","Landcover","time","mean","sd","se")
            meltedstats$TOD <- format(meltedstats$time, "%H:%M")
            clrs <- as.character(myLUT[which(myLUT$ID %in% mystats.int$zone),"Colours"])
            ordr <- as.numeric(myLUT[which(myLUT$ID %in% mystats.int$zone),"plotOrder"])
            
#             browser()
            
            print(
                ggplot(data=meltedstats, aes(x=time, y=mean, group=Landcover, fill=Landcover)) + 
                    geom_line(aes(y=mean, colour=Landcover)) + 
                    geom_ribbon(aes(ymin=mean-se, ymax=mean+se), alpha=0.3) +
                    scale_colour_manual(values = clrs[order(ordr)]) + 
                    scale_fill_manual(values = clrs[order(ordr)]) +
                    labs(title=paste("Mean Intensity of Convective Precipitation (>10mm per hour)\nby Land Cover Type in Zone",bx), y=expression(mm~hour^-1), x="Time")
                )
            
            # Put error bars on the time of day plots. 
            # Involves looping through all times of day, and classes
            sefun <- function(x, na.rm=T){x <- x[!is.na(x)]; return(sd(x)/sqrt(length(x)))}
            lnfun <- function(x, na.rm=T){x <- x[!is.na(x)]; return(length(x))}
            alltod <- dimnames(table(format(getZ(inbr), "%H%M")))[[1]]
            if(exists("todstats")){rm(todstats)}
            for (tod in alltod){
                allhrs <- format(getZ(inbr), "%H%M")
                i <- which(allhrs %in% tod)
                pr.int.bxss <- subset(pr.int.bx, i)
                mycl.bxvals <- getValues(mycl.bx)
                pr.int.bxssvals <- getValues(pr.int.bxss)
                
                for (cl in unique(mycl.bxvals)){
                    cl.alldata <- as.vector(pr.int.bxssvals[which(mycl.bxvals %in% cl),])
                    cl.mean <- mean(cl.alldata, na.rm=T)
                    cl.se   <- sefun(cl.alldata)
                    cl.sum <- sum(cl.alldata, na.rm=T)
                    cl.area <- length(which(mycl.bxvals %in% cl))*16 # One grid cell is 16kmsq
                    cl.quan <- quantile(cl.alldata, probs=c(0.25,0.5,0.75), names=F, na.rm=T)
                    if(!exists("todstats")){
                        todstats <- data.frame(TOD=tod, Class=cl, Landcover=myLUT[myLUT$ID == cl,"Landcover"], mean=cl.mean, se=cl.se, pc25=cl.quan[1], pc50=cl.quan[2], pc75=cl.quan[3], sum=cl.sum, area=cl.area)
                    } else {
                        todstats <- rbind(todstats, data.frame(TOD=tod, Class=cl, Landcover=myLUT[myLUT$ID == cl,"Landcover"], mean=cl.mean, se=cl.se, pc25=cl.quan[1], pc50=cl.quan[2], pc75=cl.quan[3], sum=cl.sum, area=cl.area))   
                    }
                }
                
            }

            
            # Plot time of day means of intense precip over different cover types
            # tmp <- ddply(meltedstats, c("Landcover","TOD"), summarise, TODmean=mean(mean, na.rm=T))
            # tmp$TOD <- as.POSIXct(paste("2006-08-16",tmp$TOD), "2006-08-16 %H:%M")
            todstats$TODpos <- as.POSIXct(paste("2006-08-16",as.character(todstats$TOD)), format="%Y-%m-%d %H%M")
            
            # Mean and standard error plot
            print(
                ggplot(todstats, aes(x=TODpos, y=mean, group=Landcover)) + 
                    geom_line(aes(colour=Landcover)) + 
                    scale_colour_manual(values = clrs[order(ordr)]) + 
                    geom_ribbon(aes(ymin=mean-se, ymax=mean+se, fill=Landcover), alpha=0.2) +
                    scale_fill_manual(values = clrs[order(ordr)]) +
                    scale_x_datetime(breaks=date_breaks("2 hour"), minor_breaks=date_breaks("1 hour"), labels=date_format("%H:%M")) +
                    labs(title=paste("Mean Intense Precipitation (>10mm per hour) by Land Cover Type\nand Time of Day in Zone",bx), y=expression(mm~hour^-1), x="Time of Day")
            )
            
            # Median and 25th - 75th percentiles
            print(
                ggplot(todstats, aes(x=TODpos, y=pc50, group=Landcover)) + 
                    geom_line(aes(colour=Landcover)) + 
                    scale_colour_manual(values = clrs[order(ordr)]) + 
                    geom_ribbon(aes(ymin=pc25, ymax=pc75, fill=Landcover), alpha=0.2) +
                    scale_fill_manual(values = clrs[order(ordr)]) +
                    scale_x_datetime(breaks=date_breaks("2 hour"), minor_breaks=date_breaks("1 hour"), labels=date_format("%H:%M")) +
                    labs(title=paste("Median Intense Precipitation (>10mm per hour) by Land Cover Type\nand Time of Day in Zone",bx), y=expression(mm~hour^-1), x="Time of Day")
                )
            
            # Plot total accumulated precip from intense MCS precip
            print(
                ggplot(todstats, aes(x=TODpos, y=sum/area, group=Landcover)) + 
                    geom_line(aes(colour=Landcover)) + 
                    scale_colour_manual(values = clrs[order(ordr)]) + 
                    scale_x_datetime(breaks=date_breaks("2 hour"), minor_breaks=date_breaks("1 hour"), labels=date_format("%H:%M")) +
                    labs(title=paste("Total Intense Precipitation (>10mm per hour) normalised by class area\nby Land Cover Type and Time of Day in Zone",bx), y=expression(mm~hour^-1), x="Time of Day")
            )
            
            if (!exists("boxdata")){
                boxdata <- cbind(zone=bx, todstats)
            } else {
                boxdata <- rbind(boxdata, cbind(zone=bx, todstats))
            }
            
        }

	}
    browser()
    # For ALL Zones, plot total accumulated precip from intense MCS precip
    print(
        ggplot(boxdata, aes(x=TODpos, y=sum/area, group=Landcover)) + 
            facet_wrap(~ zone, scales = "free_y", ncol=2) +
            geom_line(aes(colour=Landcover)) + 
            scale_colour_manual(values = clrs[order(myLUT$plotOrder)]) + 
            scale_x_datetime(breaks=date_breaks("2 hour"), minor_breaks=date_breaks("1 hour"), labels=date_format("%H:%M")) +
            labs(title=paste("Total Intense Precipitation (>10mm per hour) normalised by class area\nby Land Cover Type and Time of Day in Zone",bx), y=expression(mm~hour^-1), x="Time of Day")
    )

	dev.off()

}