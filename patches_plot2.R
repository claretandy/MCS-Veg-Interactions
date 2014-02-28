# Rerun patches_plot.R with djzxs simulation
# All in rotated pole projection
library(raster)
library(rasterVis)
library(rgdal)
library(scales)
library(plyr)
library(googleVis)
source("functions/MaskPrecip.R") # Loads MaskPrecip function
source("functions/loadAllAncils.R")
source("functions/patches.R")
source("getMyData.R")
source("functions/adjCoords.R")
source("functions/getLUT.R")
source("functions/MaskPrecip.R")
source("functions/subsetHours.R")
source("functions/plotBars2.R")

if (Sys.info()[["sysname"]] == "Darwin"){
    indatadir <- "/Users/ajh235/Work/DataLocal/ModelData/WAFR/"
    dlresultsdir <- "/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/"
    resultsdir <- "/Users/ajh235/Work/Projects/InternalSabbatical/Results/"
    scratchdir <- "/Users/ajh235/Work/Scratch/"
} else {
    indatadir <- "/data/local/hadhy/ModelData/WAFR/"
    dlresultsdir <- "/data/local/hadhy/Projects/InternalSabbatical/Results/"
    resultsdir <- "/home/h02/hadhy/Projects/InternalSabbatical/Results/"
    scratchdir <- "/data/local/hadhy/Scratch/"
    require(PP,lib.loc="/project/ukmo/rhel6/R")
}


# Set a few variables first
timestep <- "10min" # "avg"
threshold <- 1000
myproj <- "rp" # "ll"
myLUT <- getLUT(nBndClass=1)
overwrite <- F

# Get Veg classes
ancils <- loadAllAncils(myproj="rp", nBndClass=1, model="rb5216.4km.std", overwrite=F)
mycl.f <- ancils[[5]]
mycl.z <- ancils[[10]]
mylandfrac <- ancils[[2]]
land_simple <- ancils[[9]]

# Get precip
mydata <- getMyData(timestep=timestep, var="lsrain", overwrite=F)
rb5216.4km.std <- shift(adjCoords(mydata[[2]]), x=-360)

# Get MCS patches
mcs.infile <- paste(indatadir,"djzxs/patches/",threshold,"km_",timestep,"/allmcs.1000km.vrt", sep="")
mcs <- shift(adjCoords(brick(mcs.infile)), x=-360)

# Get POP patches
pop.infile <- paste(indatadir,"djzxs/patches/",threshold,"km_",timestep,"/allpop.1000km.vrt", sep="")
pop <- shift(adjCoords(brick(pop.infile)), x=-360)

# Sum precip within POP patches
tot.pr.10min.mcs <- MaskPrecip(rb5216.4km.std, mcs, paste(indatadir, "/djzxs/derived/tot.precip.10min.mcs.",threshold,"km.tif",sep=""), overwrite=F)*60*10
tot.pr.10min.pop <- MaskPrecip(rb5216.4km.std, pop, paste(indatadir, "/djzxs/derived/tot.precip.10min.pop.",threshold,"km.tif",sep=""), overwrite=F)*60*10

plotBars2(tot.pr.10min.mcs, mycl.z, mytitle=paste("Accumulated MCS precipitation divided by zonal area of class",sep=""))
plotBars2(tot.pr.10min.pop, mycl.z, mytitle=paste("Accumulated small scale precipitation divided by zonal area of class",sep=""))


# Create a subset for certain times of day
# Could loop through different hours here 
# hours <- c(18, 19, 6)
# plens <- c(2, 4, 6)
# hr=18
# plen=2

pdf("../../Results/Precip_byZone_byClass_byTime_MCS_15-21Z.pdf", width=12, height=8)
# pdf("../../Results/Precip_byZone_byClass_byTime_MCS_1-24Z.pdf", width=12, height=8)
dmcs <- dev.cur()
pdf("../../Results/Precip_byZone_byClass_byTime_POP_15-21Z.pdf", width=12, height=8)
# pdf("../../Results/Precip_byZone_byClass_byTime_POP_1-24Z.pdf", width=12, height=8)
dpop <- dev.cur()

if (exists("myresult.mcs.all")){rm(myresult.mcs.all)}
if (exists("myresult.pop.all")){rm(myresult.pop.all)}

# plen=1
# for (hr in 1:24){

plen=6
for (hr in 21){
    
    print(paste("Hour:",hr," Period:",plen))
    
    myss <- subsetHours(rb5216.4km.std, hr=hr, plen=plen)[[1]]
    phr <- subsetHours(rb5216.4km.std, hr=hr, plen=plen)[[2]]
    pr.tmp <- subset(stack(rb5216.4km.std), myss)
    mcs.tmp <- subset(mcs, myss)
    pop.tmp <- subset(pop, myss)
    
    # Do the masking for the subset of hours
    pr.sum.mcs.f <- paste(indatadir, "/djzxs/derived/pr.sum.mcs.",threshold,"km.",hr,"_",plen,".tif",sep="")
    pr.sum.pop.f <- paste(indatadir, "/djzxs/derived/pr.sum.pop.",threshold,"km.",hr,"_",plen,".tif",sep="")
    if (!file.exists(pr.sum.mcs.f) | overwrite == T){
        pr.sum.mcs <- MaskPrecip(pr.tmp, mcs.tmp, pr.sum.mcs.f, overwrite=T)*60*10        
    } else {
        pr.sum.mcs <- raster(pr.sum.mcs.f)*60*10
    }

    if (!file.exists(pr.sum.pop.f) | overwrite == T){
        pr.sum.pop <- MaskPrecip(pr.tmp, pop.tmp, pr.sum.pop.f, overwrite=T)*60*10        
    } else {
        pr.sum.pop <- raster(pr.sum.pop.f)*60*10
    }
        
    # Do the plotting 
    dev.set(dmcs)
    mytitle1 <- paste("Ratio of MCS Precipitation Fraction : Fractional coverage of Vegetation Class by Zone\n",phr,"Z to ",hr,"Z",sep="")
    mytitle2 <- paste("Ratio of MCS Precipitation Fraction : Fractional coverage of Vegetation Class by Zone\n",phr,"Z to ",hr,"Z",sep="")
    myresult.mcs <- plotBars3(pr.sum.mcs, mycl.z, mytitle1=mytitle1, mytitle2=mytitle2)
    plotBars2(pr.sum.mcs, mycl.z, mytitle=paste("Accumulated MCS precipitation for ",phr,"Z to ",hr,"Z divided by zonal area of class",sep=""))
    
    dev.set(dpop)
    mytitle1 <- paste("Ratio of small scale Precipitation Fraction : Fractional coverage of Vegetation Class by Zone\n",phr,"Z to ",hr,"Z",sep="")
    mytitle2 <- paste("Ratio of small scale Precipitation Fraction : Fractional coverage of Vegetation Class by Zone\n",phr,"Z to ",hr,"Z",sep="")
    myresult.pop <- plotBars3(pr.sum.pop, mycl.z, mytitle1=mytitle1, mytitle2=mytitle2)
    plotBars2(pr.sum.pop, mycl.z, mytitle=paste("Accumulated small scale precipitation for ",phr,"Z to ",hr,"Z divided by zonal area of class",sep=""))
    
    
    if (!exists("myresult.mcs.all")){
        myresult.mcs.all <- cbind(myresult.mcs, Hour=hr)
    } else {
        myresult.mcs.all <- rbind(myresult.mcs.all, cbind(myresult.mcs, Hour=hr))
    }

    if (!exists("myresult.pop.all")){
        myresult.pop.all <- cbind(myresult.pop, Hour=hr)
    } else {
        myresult.pop.all <- rbind(myresult.pop.all, cbind(myresult.pop, Hour=hr))
    }
    
}
dev.off(dmcs) 
dev.off(dpop)

# Save results for 1:24
save(list=c("myresult.mcs.all", "myresult.pop.all"), file="myresult.all.Rdata")

# Load results from previous run of the above ...
load("myresult.all.Rdata")

# googleVis code ...
incols <- c("ID","Hour","Landcover","ClFrac","PrFrac","ratio")
myresult.pop.all$ID <- myresult.pop.all$Class + (as.numeric(myresult.pop.all$Zone) * 10)
myresult.mcs.all$ID <- myresult.mcs.all$Class + (as.numeric(myresult.mcs.all$Zone) * 10)
pop.chart <- gvisMotionChart(data=myresult.pop.all[,incols], idvar="ID", timevar="Hour", colorvar="Landcover", xvar="ClFrac", yvar="PrFrac", sizevar="ratio")
mcs.chart <- gvisMotionChart(data=myresult.mcs.all[,incols], idvar="ID", timevar="Hour", colorvar="Landcover", xvar="ClFrac", yvar="PrFrac", sizevar="ratio")
plot(pop.chart)
plot(mcs.chart)


pdf("../../Results/Actual-Expected.pdf", width=12, height=8)

# total land area per land cover class by zone
ggplot(myresult.pop.all[myresult.pop.all$Hour==10,], aes(x=Zone, y=areaClZone, fill=Landcover)) + geom_bar(stat="identity") + scale_fill_manual(values=as.character(myLUT$Colours)[c(3,2,4,1,5)]) + labs(title="Total Land Area per Land Cover Class by Zone", y="Area sqkm")

# Total amount of precipitation from small scale precipitation per land cover class by zone
tmp <- ddply(myresult.pop.all, .(Zone, Landcover), summarize, totPr=sum(prClZone))
ggplot(tmp, aes(x=Zone, y=totPr, fill=Landcover)) + 
    scale_fill_manual(values=as.character(myLUT$Colours)[c(3,2,4,1,5)]) + 
    geom_bar(stat="identity") + 
    labs(title="Total Accumulated Small-Scale Precipitation by Zone and Land Cover class", y="Precipitation (mm)")

# Total amount of precipitation from MCS precipitation per land cover class by zone
tmp <- ddply(myresult.mcs.all, .(Zone, Landcover), summarize, totPr=sum(prClZone))
ggplot(tmp, aes(x=Zone, y=totPr, fill=Landcover)) + 
    scale_fill_manual(values=as.character(myLUT$Colours)[c(3,2,4,1,5)]) + 
    geom_bar(stat="identity") + 
    labs(title="Total Accumulated MCS Precipitation by Zone and Land Cover class", y="Precipitation (mm)")


# Actual minus Expected precipitation by zone and land cover class
# Small scale precipitation
myresult.pop.all$Hour.POSIXct <- as.POSIXct(paste("2006-08-16 ",myresult.pop.all$Hour,":00",sep=""))
ggplot(data=myresult.pop.all[myresult.pop.all$Zone %in% c(2,4,7,8,9,12,13),], aes(x=Hour.POSIXct, y=prClZone-(ClFrac*prZone))) + # c(2,4,7,8,9,12,13)
    geom_line(aes(colour=Landcover)) + 
    labs(title="Actual minus Expected Small-Scale Precipitation over Land Cover types\nIn Zones 2,4,7,8,9,12 and 13", x="Hour", y="Accumulated Precipitation (mm per zone)") + 
    facet_grid(Zone ~ .) + 
    scale_colour_manual(values=as.character(myLUT$Colours)[c(3,2,4,1,5)]) + 
    scale_x_datetime(breaks = date_breaks("3 hours"), labels = date_format("%H:%M")) + 
    coord_cartesian(xlim = c(as.POSIXct("2006-08-16 00:30"), as.POSIXct("2006-08-17 00:30")))

# MCS precipitation
myresult.mcs.all$Hour.POSIXct <- as.POSIXct(paste("2006-08-16 ",myresult.mcs.all$Hour,":00",sep=""))
ggplot(data=myresult.mcs.all[myresult.mcs.all$Zone %in% c(2,4,7,8,9,12,13),], aes(x=Hour.POSIXct, y=prClZone-(ClFrac*prZone))) + # c(2,4,7,8,9,12,13)
    geom_line(aes(colour=Landcover)) + 
    labs(title="Actual minus Expected MCS Precipitation over Land Cover types\nIn Zones 2,4,7,8,9,12 and 13", x="Hour", y="Accumulated Precipitation (mm per zone)") + 
    facet_grid(Zone ~ .) + 
    scale_colour_manual(values=as.character(myLUT$Colours)[c(3,2,4,1,5)]) + 
    scale_x_datetime(breaks = date_breaks("3 hours"), labels = date_format("%H:%M")) + 
    coord_cartesian(xlim = c(as.POSIXct("2006-08-16 00:30"), as.POSIXct("2006-08-17 00:30")))

dev.off()
