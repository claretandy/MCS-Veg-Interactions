# Plots to go in the paper
# All plots MUST be in Lat Long

library(raster)
library(rasterVis)
library(rgdal)
source("functions/loadAllAncils.R")
source("functions/loadVeg.R") # Returns myveg = fractional veg cover for each pft tile
source("functions/loadOtherAncils.R")
source("functions/makeBoxes.R") 
source("functions/vegPrep.R") # Returns allveg = vegetation classess (1 to 6) AND veg classes intersected with zones (i.e. boxes)
source("functions/patches.R")
source("functions/mcsStats.R")
source("functions/popStats.R")
source("functions/initiations.R")
source("functions/makeLines.R")
source("getMyData.R")
source("trackCheck.R")
# source("mcsIntensity.R")

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

rasterOptions(tmpdir=scratchdir, todisk=F)

timestep <- "10min" # "avg"
threshold <- 1000
myproj <- "ll" # "rp"
models <- c("rb5216.4km.std", "rb5216.4km.50k", "rb5216.4km.300k") # [2:3]
id <- "s" # c("s","w","y")#[2:3]

# Get precip data
mydata <- getMyData(timestep=timestep, var="lsrain", overwrite=F)
rb5216.4km.std <- mydata[[2]]

land_simple <- readOGR(dsn=paste(indatadir,"ancils",sep=""), layer="land_ll") # Lat Long

# Get boxes made in RP, but projected to LatLong
spp <- readOGR(dsn="/Users/ajh235/Work/DataLocal/ModelData/WAFR/ancils", layer="boxes_rp2ll")
spp.r <- rasterize(spp, mylandfrac)

myveg <- loadVeg(model.nm=models[1], proj=myproj, overwrite=F)
mylandfrac <- loadOtherAncils(model.nm=models[1], ancil="landfrac", proj=myproj, overwrite=F)
myorog <- loadOtherAncils(model.nm=models[1], ancil="orog", proj=myproj,overwrite=F)  

allveg <- vegPrep(model.nm=models[1], id=id[1], myproj=myproj, myveg, myorog, mylandfrac, land_simple, spp, spp.r, plots=F, vegThreshold=0.3, overwrite=T) # return(mycl, mycl.z) and creates pdf plots
mycl <- allveg[[1]]

myLUT <- data.frame(ID=c(1,2,3,4,5,6,7), Landcover=factor(c("tree", "grass", "sparse", "boundary", "boundary, tree", "boundary, grass", "orography"), levels=c("tree", "grass", "sparse", "boundary", "boundary, tree", "boundary, grass", "orography")[c(3,2,6,4,5,1,7)]), Colours=c("dark green", "yellow", "orange", "sienna", "yellow green", "gold", "dark grey"), plotOrder=c(4,2,1,3,5,6,7))

mycl.f <- as.factor(mycl)
ftab <- levels(mycl.f)[[1]]
ftab$Name <- myLUT$Landcover
levels(mycl.f) <- ftab

# Plot model domains
# 12km
e12km <- extent(c(xmin=-21.99668292,  xmax=14.06752545, ymin=-0.18893839,  ymax=23.96993351))
land12k.pol <- as(extent(e12km), "SpatialPolygons")
land12k.pol <- SpatialPolygonsDataFrame(land12k.pol, data=data.frame(id=1))
writeOGR(land12k.pol, dsn="/Users/ajh235/Work/DataLocal/ModelData/WAFR/ancils/km12/", layer="extent_12km_ll", driver="ESRI Shapefile", check_exists=T, overwrite_layer=T)

# 4km
e4km <- extent(c(xmin=-20.62620937,  xmax=12.58290222, ymin=1.29328735,  ymax=22.85308582))
land4k.pol <- as(extent(e4km), "SpatialPolygons")
land4k.pol <- SpatialPolygonsDataFrame(land4k.pol, data=data.frame(id=1))
writeOGR(land4k.pol, dsn="/Users/ajh235/Work/DataLocal/ModelData/WAFR/ancils/km4", layer="extent_4km_ll", driver="ESRI Shapefile", check_exists=T, overwrite_layer=T)


# Plot Vegetation classes w/ all boundary classes
png("../../Results/Vegetation_classes2.png", width=1000, height=600)
print(
    levelplot(mycl.f, maxpixels=600000, par.settings=rasterTheme(region=myLUT$Colours), xlab=NULL, ylab=NULL, xlim=c(-12,10), ylim=c(4,18), main="Vegetation classes and zones") + # , scales=list(draw=FALSE),  xlim=c(-24,15), ylim=c(-1,26), 
        latticeExtra::layer(sp.polygons(land_simple, col="black", lty=2)) + 
        latticeExtra::layer(sp.polygons(spp)) +
        latticeExtra::layer(sp.text(loc=coordinates(spp), txt=1:nrow(spp@data), cex=3)) +
        latticeExtra::layer(sp.polygons(land12k.pol)) +
        latticeExtra::layer(sp.polygons(land4k.pol))
)
dev.off()

# Plot Vegetation classes w/ ONE boundary class
myLUT <- data.frame(ID=c(1,2,3,4,7), Landcover=factor(c("tree", "grass", "sparse", "boundary", "orography"), levels=c("tree", "grass", "sparse", "boundary", "orography")[c(3,2,4,1,5)]), Colours=c("dark green", "yellow", "orange", "sienna", "dark grey"), plotOrder=c(4,2,1,3,5))

mycl.1b <- mycl
mycl.1b[mycl.1b == 5] <- 4
mycl.1b[mycl.1b == 6] <- 4
mycl.1bf <- as.factor(mycl.1b)
ftab <- myLUT[myLUT$ID %in% levels(mycl.f)[[1]]$ID, ]
levels(mycl.1bf) <- ftab

png("../../Results/Vegetation_classes_1bnd.png", width=1000, height=600)
print(
    levelplot(mycl.1bf, maxpixels=600000, par.settings=rasterTheme(region=myLUT$Colours), xlab=NULL, ylab=NULL, xlim=c(-12,10), ylim=c(4,18), main="Vegetation classes and zones") + # , scales=list(draw=FALSE),  xlim=c(-24,15), ylim=c(-1,26), 
        latticeExtra::layer(sp.polygons(land_simple, col="black", lty=2)) + 
        latticeExtra::layer(sp.polygons(spp)) +
        latticeExtra::layer(sp.text(loc=coordinates(spp), txt=1:nrow(spp@data), cex=3)) +
        latticeExtra::layer(sp.polygons(land12k.pol)) +
        latticeExtra::layer(sp.polygons(land4k.pol))
)
dev.off()

# Plot afternoon initiations
aftinit <- results[results$class == 'generation' & (as.numeric(format(results$timestep, "%H")) >= 16 & as.numeric(format(results$timestep, "%H")) <= 17),c("x","y")]
initpts_rp <- SpatialPoints(aftinit, CRS("+proj=ob_tran +o_proj=longlat +o_lon_p=175.3000030517578 +o_lat_p=77.4000015258789 +lon_0=180 +ellps=sphere"))
# Reproject afternoon initiation points
initpts_ll <- spTransform(initpts_rp, CRSobj=CRS("+init=epsg:4326"), use_ob_tran=T)
initpts_ll <- initpts_ll[!is.na(extract(mycl.f, initpts_ll)),]

png("../../Results/Vegetation_AfternoonInitiations.png", width=1000, height=600)
print(
    levelplot(mycl.1bf, att="Landcover", maxpixels=600000, main="Afternoon (16-18Z) MCS Initiations Over Vegetation Classes", xlim=c(-18,11), ylim=c(4,20), xlab=NULL, ylab=NULL, col.regions=as.character(fdat$Colours)) + # scales=list(draw=FALSE), xlim=c(-3,10), ylim=c(12,20)
        latticeExtra::layer(sp.polygons(land_rp2ll, lty=2)) +
        latticeExtra::layer(sp.points(initpts_ll, pch="+", cex=4, col="black")) #+
    )
dev.off()

# Plot POP and MCS precipitation statistics
source("patches_plot2.R")

# Get Veg classes in rotated pole
ancils <- loadAllAncils(myproj="rp", nBndClass=1, model="rb5216.4km.std", overwrite=F)
mycl <- ancils[[4]]
mycl.z <- ancils[[10]]
mylandfrac <- ancils[[2]]
land_simple <- ancils[[9]]
spp.r <- ancils[[8]][[1]]
sppa <- ancils[[8]][[3]]

# MCS intense precipitation 
#diurnalcycle2(rb5216.4km.std, type="all", patch=F, model.nm="rb5216.4km.std", id="s", spp.r=spp.r, sppa=sppa, mycl=mycl, land_simple, overwrite=F) # Creates pdf plots
diurnalcycle2(rb5216.4km.std, type="intense", patch=F, model.nm="rb5216.4km.std", id="s", spp.r=spp.r, sppa=sppa, mycl=mycl, land_simple, overwrite=F) # Creates pdf plots

