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
pdf(paste("/Users/ajh235/Work/Projects/InternalSabbatical/Results/MCS_intensity_bysize_gt",inthr,"mm_v5.pdf",sep=""), width=14, height=9)
mcsrst.path <- "/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/Tracking/djzxs_10min/"
gntheme <- rasterTheme(pch=19, region=brewer.pal(9, 'YlGn'))

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
iii <- as.numeric(attributes(sort(table(results[results$TOD > 15 & results$TOD <= 21,"ID"])))$dimnames[[1]]) # These are all the MCS that occur between 15:00 and 21:00
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
        
#             # Mask precip that's not in the MCS
        mcsprecip   <- raster::mask(precip.now, mcsiinow)*3600
#             precip.mask <- getValues(mcsprecip)
#             precip.mask <- precip.mask[!is.na(precip.mask)]
#             
#             # Mask blfrac that's not in the MCS
#             mcstree     <- mask(blfrac, mcsiinow)
#             blfrac.mask <- getValues(mcstree)
#             blfrac.mask <- blfrac.mask[!is.na(blfrac.mask)]
#             
#             # Mask grassfrac that's not in the MCS
#             mcsgrass       <- mask(grassfrac, mcsiinow)
#             grassfrac.mask <- getValues(mcsgrass)
#             grassfrac.mask <- grassfrac.mask[!is.na(grassfrac.mask)]
#             
#             # What is the max precip rate, and what are veg fractions at that location?
#             maxprecip <- cellStats(mcsprecip, 'max')*60*60
#             mcspoints <- rasterToPoints(mcsprecip, spatial=FALSE)
#             maxpoint  <- matrix(mcspoints[which(mcspoints[,3]==maxprecip/3600), 1:2], byrow=T, nrow=1, ncol=2)
#             maxgrass  <- extract(grassfrac, maxpoint)
#             maxtree   <- extract(shift(blfrac, x=-360), maxpoint)
        
        # Where the rainfall is most intense, what's the rate over grass and the rate over tree?
        
        # 1. Classify intense rain (>20mm)
        mcsprecip.intense <- calc(mcsprecip, fun=function(x){x[x<inthr]<-NA; return(x)})
        
        # 2. Mask veg classes to intense precip patches
        mycl.mask <- raster::mask(mycl, mcsprecip.intense)
        browser()

        # Check that some cells have data in ...
        haveData <- length(which(is.na(getValues(mycl.mask)))) != ncell(mycl.mask)
        
        if (haveData){
            
        } else {
            print("We don\'t have any data")
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
        
    }
    
    mydf$landcover <- factor(x=mydf$landcover, levels=c("sparse", "grass","boundary", "tree", "orography"))
    mydf$sdmin[mydf$sdmin < inthr] <- inthr

    f2 <- ggplot(mydf, aes(x=time, group=landcover, fill=landcover))
    p2 <- f2 + geom_ribbon(aes(ymin=sdmin, ymax=sdmax), alpha=0.3) + geom_line(aes(y=mean, colour=landcover)) + scale_colour_manual(values = c("orange","yellow","brown","dark green","dark grey")) + scale_fill_manual(values = c("orange","yellow","brown","#31a354","dark grey")) + labs(title=paste("MCS Mean Intense Precipitation (>",inthr,"mm) by Land Cover Type\nMCS id:",ii), y=expression(mm~hour^-1), x="Time")
#         print(p2)

    p3 <- f2 + geom_line(aes(y=count, colour=landcover)) + scale_colour_manual(values = c("orange","yellow","brown","dark green","dark grey")) + labs(title=paste("Grid Cell Count of Intense Precipitation (>",inthr,"mm) per Land Cover Type",sep=""), y="Grid cells", x="Time")

#         multiplot(p2, p3, layout=layout)
    print(p2)
    print(p3)

    rm(mydf)
#         rm(mydf4gg)
#         rm(ribbondf)
    
}

dev.off()
# }
