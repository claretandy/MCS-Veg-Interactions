# Need to run first ...
# source("functions/loadVeg.R")
# First 75 lines of loadData_diurnal.R
# Need to adapt script to be called from the above script ...
library(ggplot2)
library(reshape)
source("functions/makeLines.R")

# Load results table
load("/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/Tracking/djzxs_10min.RData")

# Load veg fractions
vegfrac <- brick("/Users/ajh235/Work/DataLocal/ModelData/WAFR/ancils/km4/qrparm.veg.frac_4km.std.nc")
vegfrac <- shift(vegfrac, x=-360)
blfrac <- vegfrac[[1]]; blfrac[is.na(blfrac)] <- 0
grassfrac <- vegfrac[[3]] + vegfrac[[4]]; grassfrac[is.na(grassfrac)] <- 0

# Get unique list of MCS ids
mcs.ids <- unique(results$ID)

# Make MCS tracks from results
my.lines.sp <- makeLines(results)[[1]]

# Plot lines with BL tree fraction in background
pdf("/Users/ajh235/Work/Projects/InternalSabbatical/Results/MCS_intensity_bysize_v2.pdf", width=14, height=9)
gntheme <- rasterTheme(pch=19, region=brewer.pal(9, 'YlGn'))
blfrac <- myveg[[1]]
my.plot <- rasterVis::levelplot(shift(blfrac, x=-360), margin=F, auto.key=list(title="Broadleaf tree fraction"), par.settings=gntheme, at=seq(0,1,0.1), xlim=c(-13,16), ylim=c(-9,9), main="MCS tracks (16/08 to 18/08)") + latticeExtra::layer(sp.polygons(land_simple, col="dark grey"))
my.plot <- my.plot + latticeExtra::layer(sp.lines(my.lines.sp))
print(my.plot)

# What is the intensity of rainfall over forest and grass within an MCS?
# load all MCS raster data
# mcs.length <- sort(table(results$ID))
# iii <- as.numeric(dimnames(mcs.length[mcs.length > 100])[[1]]) # Get all MCS with length > 100
id.pixcount <- aggregate(pixcount ~ ID, data=results, FUN=sum)
iii <- id.pixcount[which(id.pixcount$pixcount > 20000),"ID"]
myspdf <- makeLines(results[results$ID %in% iii,])
mytext <- myspdf[[2]]
myspdf <- myspdf[[1]]
# startingpt <- lapply(coordinates(myspdf), FUN=function(x){x[[1]][1,]})
# mytext.check <- cbind(matrix(unlist(startingpt), ncol=2, byrow=T), iii)

mcsrst.path <- "/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/Tracking/djzxs_10min/"

coordinates(mytext) <- ~x+y

print(levelplot(shift(blfrac, x=-360), margin=F, par.settings=gntheme) + latticeExtra::layer(sp.polygons(land_simple, col="grey")) + latticeExtra::layer(sp.lines(myspdf)) + latticeExtra::layer(sp.text(loc=coordinates(mytext), txt=mytext$ID)) )


for (ii in iii){
    
    mydates <- results[results$ID == ii,"timestep"]
    
    for (x in 1:length(mydates)){
        print(mydates[x])
        mcs.now <- shift(raster(paste(mcsrst.path,"mcs_tracking_1000km_",format(mydates[x], "%d.%H%M"), ".tif",sep="")), x=-360)
        
        mcsiinow <- mcs.now == ii
        mcsiinow[mcsiinow == 0] <- NA
        precip.now <- shift(rb5216.4km.std[[which(getZ(rb5216.4km.std) == mydates[x])]], x=-360)
        
        # Mask precip that's not in the MCS
        precip.mask <- getValues(mask(precip.now, mcsiinow))
        precip.mask <- precip.mask[!is.na(precip.mask)]
        
        # Mask blfrac that's not in the MCS
        blfrac.mask <- getValues(mask(shift(blfrac, x=-360), mcsiinow))
        blfrac.mask <- blfrac.mask[!is.na(blfrac.mask)]
        
        # Mask grassfrac that's not in the MCS
        grassfrac.mask <- getValues(mask(grassfrac, mcsiinow))
        grassfrac.mask <- grassfrac.mask[!is.na(grassfrac.mask)]
        
        # What is the max precip rate, and what are veg fractions at that location?
        maxprecip <- cellStats(mask(precip.now, mcsiinow), 'max')*60*60
        mcspoints <- rasterToPoints(mask(precip.now, mcsiinow), spatial=FALSE)
        maxpoint  <- matrix(mcspoints[which(mcspoints[,3]==maxprecip/3600), 1:2], byrow=T, nrow=1, ncol=2)
        maxgrass  <- extract(grassfrac, maxpoint)
        maxtree   <- extract(shift(blfrac, x=-360), maxpoint)
        
        # Where the rainfall is most intense, what's the rate over grass and the rate over tree?
        # Try tree threshold >0.3, but may need to experiment
        browser()
        
        if (!exists("mydf")){
            mydf <- data.frame(time=mydates[x], blfrac=mean(blfrac.mask) , grassfrac=mean(grassfrac.mask) , precip=mean(precip.mask)*60*60, maxprecip=maxprecip, maxgrass=maxgrass, maxtree=maxtree, areasqkm=sum(getValues(mcsiinow), na.rm=T)*4)
        } else {
            mydf <- rbind(mydf, data.frame(time=mydates[x], blfrac=mean(blfrac.mask) , grassfrac=mean(grassfrac.mask) , precip=mean(precip.mask)*60*60, maxprecip=maxprecip, maxgrass=maxgrass, maxtree=maxtree, areasqkm=sum(getValues(mcsiinow), na.rm=T)*4))
        }
    }
    
    mydf4gg <- melt(mydf, id.vars=1)
    f <- ggplot(mydf4gg, aes(time, value, ymin=0, ymax=value)) + scale_colour_identity() + facet_grid(variable ~ ., scales="free") + labs(title = paste("MCS id:",ii))
    p <- f + geom_line(subset=.(variable=="blfrac")) + geom_line(subset=.(variable=="grassfrac")) + geom_line(subset=.(variable=="precip")) + geom_line(subset=.(variable=="maxprecip")) + geom_line(subset=.(variable=="maxtree")) + geom_line(subset=.(variable=="maxgrass")) + geom_line(subset=.(variable=="areasqkm"))
    print(p)
    
    rm(mydf)
    rm(mydf4gg)
}

dev.off()
    
# a) At all times of day
# b) Split by time of day