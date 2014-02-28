# Need to run first ...
# source("functions/loadVeg.R")
# First 75 lines of loadData_diurnal.R
# Need to adapt script to be called from the above script ...
# 
# mcsIntensity <- function(results, myveg, mcs, inbr){
    
source("functions/makeLines.R")
source("functions/multiplot.R")
library(ggplot2)
library(reshape)
    
inthr <- 10
pdf(paste("/Users/ajh235/Work/Projects/InternalSabbatical/Results/MCS_intensity_bysize_gt",inthr,"mm_v6.pdf",sep=""), width=14, height=9)
mcsrst.path <- "/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/Tracking/djzxs_10min/"
gntheme <- rasterTheme(pch=19, region=brewer.pal(9, 'YlGn'))
myLUT <- data.frame(ID=c(1,2,3,4,7), Landcover=factor(c("tree", "grass", "sparse", "boundary", "orography"), levels=c("tree", "grass", "sparse", "boundary", "orography")[c(3,2,4,1,5)]), Colours=c("dark green", "yellow", "orange", "brown", "dark grey"), plotOrder=c(4,2,1,3,5))

# Load veg fractions
if (round(extent(mycl)@xmin + (extent(mycl)@xmax - extent(mycl)@xmin)/2) == 360){
    mycl <- shift(mycl, x=-360)
}

if (round(extent(myveg)@xmin + (extent(myveg)@xmax - extent(myveg)@xmin)/2) == 360){
    myveg <- shift(myveg, x=-360)
}

blfrac <- myveg[[1]]; blfrac[is.na(blfrac)] <- 0
grassfrac <- myveg[[3]] + myveg[[4]]; grassfrac[is.na(grassfrac)] <- 0
if (!exists("mycl")){
    allveg <- vegPrep(model.nm=models[1], id=id[1], myveg, myorog, mylandfrac, land_simple, sppa, spp.r, plots=F, vegThreshold=0.3, overwrite=F) # return(mycl, mycl.z) and creates pdf plots
    mycl <- allveg[[1]]
    mycl <- shift(mycl, x=-360)
} 

# Recode "boundary tree" and "boundary grass" to "boundary"
mycl[mycl==5] <- 4
mycl[mycl==6] <- 4

# Get unique list of MCS ids
mcs.ids <- unique(results$ID)
alllines <- makeLines(results)[[1]]

# Plot lines with BL tree fraction in background
#     print(
#         levelplot(blfrac, margin=F, auto.key=list(title="Broadleaf tree fraction"), par.settings=gntheme, at=seq(0,1,0.1), xlim=c(-13,16), ylim=c(-9,9), main="MCS tracks (16/08 to 18/08)") + 
#         latticeExtra::layer(sp.polygons(land_simple, col="dark grey")) + 
#         latticeExtra::layer(sp.lines(alllines))
#     )

# Create a subset of all MCS, either by length of by area
# mcs.length <- sort(table(results$ID))
# iii <- as.numeric(dimnames(mcs.length[mcs.length > 100])[[1]]) # Get all MCS with length > 100
# id.pixcount <- aggregate(pixcount ~ ID, data=results, FUN=sum)
#     iii <- id.pixcount[which(id.pixcount$pixcount > 20000),"ID"]
results$TOD <- format(results$timestep, "%H:%M")
mcslist <- sort(table(results[results$TOD > 15 & results$TOD <= 21,"ID"]))
iii <- as.numeric(attributes(mcslist[which(mcslist > 2)])$dimnames[[1]]) # These are all the MCS that occur between 15:00 and 21:00

#     iii <- c(110, 136, 150, 174, 328, 334, 345, 347, 469, 542, 586, 593, 611, 782, 812, 818, 863, 884, 919, 934, 949, 987) # This is a subset of the above > 20,000 pixel threshold
spresults <- makeLines(results[results$ID %in% iii,])
myspdf <- spresults[[1]] # Spatial Lines
mytext <- spresults[[2]] # Spatial points of starting point

# Get coordinates for starting points
coordinates(mytext) <- ~x+y

print(
    levelplot(blfrac, main="Large MCS tracks (16/08 to 18/08)", margin=F, par.settings=gntheme) + 
          latticeExtra::layer(sp.polygons(land_simple, col="grey")) + 
          latticeExtra::layer(sp.lines(myspdf)) + 
          latticeExtra::layer(sp.text(loc=coordinates(mytext), txt=mytext$ID)) 
    )

if (exists("bigdata")){ rm(bigdata) }

# Loop through each track
for (ii in iii){
    
    mydates <- results[results$ID == ii,"timestep"]
    if (exists("mydf")){ rm(mydf) }
    
    print(paste("MCS: ",ii, "; ",length(mydates)," timesteps",sep=" "))
    
    for (x in 1:length(mydates)){
#         print(mydates[x])
        
        mcs.now <- shift(raster(paste(mcsrst.path,"mcs_tracking_1000km_",format(mydates[x], "%d.%H%M"), ".tif",sep="")), x=-360)
        mcs.now <- adjCoords(mcs.now)
        
        mcsiinow <- mcs.now == ii
        mcsiinow[mcsiinow == 0] <- NA
        precip.now <- shift(rb5216.4km.std[[which(getZ(rb5216.4km.std) == mydates[x])]], x=-360)
        
        # Mask precip that's not in the MCS
        mcsprecip   <- raster::mask(precip.now, mcsiinow)*3600

        # Where the rainfall is most intense, what's the rate over grass and the rate over tree?
        # 1. Classify intense rain (>20mm)
        mcsprecip.intense <- calc(mcsprecip, fun=function(x){x[x<inthr]<-NA; return(x)})
        
        # 2. Mask veg classes to intense precip patches
        mycl.mask <- raster::mask(mycl, mcsprecip.intense)
        
        # Check that some cells have data in ...
        haveData <- length(which(is.na(getValues(mycl.mask)))) != ncell(mycl.mask)
        
        if (haveData){
            
            # Write all data to a big dataframe
            loc <- which(!is.na(getValues(mcsprecip.intense)))
            prvals <- getValues(mcsprecip.intense)[loc]
            clvals <- getValues(mycl.mask)[loc]
            if (!exists("bigdata")){
                bigdata <- data.frame("ID"=ii, "time"=mydates[x], "TOD"=format(mydates[x], "%H:%M"), "Precip"=prvals, "Class"=clvals)
            } else {
                bigdata <- rbind(bigdata, data.frame("ID"=ii, "time"=mydates[x], "TOD"=format(mydates[x], "%H:%M"), "Precip"=prvals, "Class"=clvals))
            }
            
            # 3. Zonal mean and std dev 
            cl.mean <- zonal(mcsprecip.intense, mycl, fun='mean')
            cl.std  <- zonal(mcsprecip.intense, mycl, fun='sd')
            cl.len  <- freq(mycl.mask, useNA="no")
            
            zon.stats <- data.frame(cl.mean, std=cl.std[,2], count=NA)
            zon.stats[which(cl.mean[,"zone"] %in% cl.len[,"value"]),"count"] <- cl.len[,"count"]
            zon.stats$landcover <- c("tree","grass","sparse","boundary", "orography")
            
            # 4. Add results to a data.frame
            if (!exists("mydf")){
                mydf <- data.frame(zon.stats[,1:2], sdmin=zon.stats[,2]-zon.stats[,3], sdmax=zon.stats[,2]+zon.stats[,3], zon.stats[,4:5], time=mydates[x])                
            } else {
                mydf <- rbind(mydf, data.frame(zon.stats[,1:2], sdmin=zon.stats[,2]-zon.stats[,3], sdmax=zon.stats[,2]+zon.stats[,3], zon.stats[,4:5], time=mydates[x]))
            }
            
        } else {
            print("We don\'t have any data for this timestep")
        }
        
    }
    
    if (exists("mydf")){
        if (length(unique(mydf$time))>3){
            # Plot results for this MCS
            mydf$landcover <- factor(x=mydf$landcover, levels=c("sparse", "grass","boundary", "tree", "orography"))
            mydf$sdmin[mydf$sdmin < inthr] <- inthr
            
            f2 <- ggplot(mydf, aes(x=time, group=landcover, fill=landcover))
            p2 <- f2 + geom_ribbon(aes(ymin=sdmin, ymax=sdmax), alpha=0.3) + geom_line(aes(y=mean, colour=landcover)) + scale_colour_manual(values = c("orange","yellow","brown","dark green","dark grey")) + scale_fill_manual(values = c("orange","yellow","brown","#31a354","dark grey")) + labs(title=paste("MCS Mean Intense Precipitation (>",inthr,"mm) by Land Cover Type\nMCS id:",ii), y=expression(mm~hour^-1), x="Time")
            
            p3 <- f2 + geom_line(aes(y=count, colour=landcover)) + scale_colour_manual(values = c("orange","yellow","brown","dark green","dark grey")) + labs(title=paste("Grid Cell Count of Intense Precipitation (>",inthr,"mm) per Land Cover Type",sep=""), y="Grid cells", x="Time")
            
            #         multiplot(p2, p3, layout=layout)
            print(p2)
            print(p3)
            
            rm(mydf)
            
        }
    }
    
#     browser()
    
} # End of ii

dev.off()

# For all MCS that occur between 15Z and 21Z, how many intense rainfall points occur over each cover type?
bigdata$Hour <- as.numeric(format(bigdata$time, "%H"))
rm(out)
for (i in 15:21){
    print(i)
    myfreq <- table(bigdata[bigdata$Hour > i & bigdata$Hour <= i+1,"Class"])
    myfreq <- data.frame(myfreq)
    colnames(myfreq)[2] <- paste("Hr", i, sep="_")
    if (!exists("out")){
        out <- cbind(myLUT[which(myLUT$ID %in% myfreq$Var1),], myfreq)
        out$Var1 <- NULL
    } else {
        out$smthg <- NA
        out[which(out$ID %in% myfreq$Var1),"smthg"] <- myfreq[,2]
        z <- length(colnames(out))
        colnames(out)[z] <- paste("Hr", i, sep="_")
    }
    out$Var1 <- NULL
}
barplot(as.matrix(out[,-c(1:4)]), names.arg=15:21, col=as.character(out$Colours))

# Find MCS that occur over tree/boundary/grass at the same time, and compare mean precip rate
for (x in unique(bigdata$ID)){
    print(x)
    alltime <- unique(bigdata[bigdata$ID == x, "time"])
    for (y in 1:length(alltime)){
        indata <- bigdata[bigdata$ID == x & bigdata$time == alltime[y], c("Class", "Precip")]
        # Test if >10% of each class
        infreq <- data.frame(table(indata$Class))
        infreq$Freq <- infreq$Freq/sum(infreq$Freq)
        this.lc <- myLUT[which(myLUT$ID %in% infreq$Var1), "Landcover"]
        len.lc <- length(which(c("tree", "grass", "boundary") %in% this.lc))
        if (len.lc == 3){
            gt10pc <- length(which(infreq[infreq$Var1 %in% c(1,2,4),"Freq"] > c(0.1, 0.1, 0.1)))
            if (gt10pc == 3){
                # Here we are!!!
                if (!exists("smalldata")){
                    smalldata <- bigdata[bigdata$ID == x & bigdata$time == alltime[y], ]
                } else {
                    smalldata <- rbind(smalldata, bigdata[bigdata$ID == x & bigdata$time == alltime[y], ])
                }
            }
        }
    }
}


# Idea: Plot scatterplot of tree precip vs grass precip for only MCS where >10% of each cover type exists
ids <- as.numeric(attributes(which(table(smalldata$ID)>20))$names)
rm(compall)
for (xx in ids){
    try(
        treepr <- aggregate(Precip ~ time, data=smalldata[smalldata$ID == xx & smalldata$Hour > 14 & smalldata$Hour < 21 & smalldata$Class==1,], FUN=mean)[,2]
        , silent=T)
    try(
        grasspr <- aggregate(Precip ~ time, data=smalldata[smalldata$ID == xx & smalldata$Hour > 14 & smalldata$Hour < 21 & smalldata$Class==2,], FUN=mean)[,2]
        , silent=T)
    try(
        bndrypr <- aggregate(Precip ~ time, data=smalldata[smalldata$ID == xx & smalldata$Hour > 14 & smalldata$Hour < 21 & smalldata$Class==4,], FUN=mean)[,2]
        , silent=T)
    
    if(exists("treepr") & exists("grasspr") & exists("bndrypr")){
        pltmin <- floor(min(treepr, grasspr))
        pltmax <- ceiling(max(treepr, grasspr))
        plot(x=treepr, y=grasspr, xlim=c(pltmin,pltmax), ylim=c(pltmin, pltmax), main=xx)
        lines(c(pltmin-5,pltmax+5), c(pltmin-5,pltmax+5))
        hr <- as.numeric(format(unique(smalldata[which(smalldata$ID == xx & smalldata$Hour > 14 & smalldata$Hour < 21 & smalldata$Class==2),"time"]), "%H"))
        if (!exists("compall")){
            compall <- data.frame(ID=xx, Treepr=treepr, Grasspr=grasspr, Bndpr=bndrypr, Hour=hr)
        } else {
            compall <- rbind(compall, data.frame(ID=xx, Treepr=treepr, Grasspr=grasspr, Bndpr=bndrypr, Hour=hr))
        }
        rm(grasspr)
        rm(treepr)
        rm(bndrypr)
    }
    
}

