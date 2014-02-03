# Load and process data from all 4 model variants
# Need to re-work many of old functions so that they read and write to/from the correct locations.

library(raster)
library(rasterVis)
library(rgdal)
source("functions/consec_drywet.R")
source("functions/diurnal_cycle.R")
source("functions/loadVeg.R") # Returns myveg = fractional veg cover for each pft tile
source("functions/loadOtherAncils.R")
source("functions/makeBoxes.R") 
source("functions/vegPrep.R") # Returns allveg = vegetation classess (1 to 6) AND veg classes intersected with zones (i.e. boxes)
source("functions/patches.R")
source("functions/movies.R")
source("functions/tracker.R")
source("functions/tracker_v2.R")
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

adjCoords <- function(r){
    this.res <- res(r)[2]
    if ( this.res <= 0.0136){
        r <- shift(r, x=0.0068, y=0.00675)
    }
    
    if (this.res <= 0.0365 & this.res > 0.0355){
        r <- shift(r, x=0.01798, y=0.018)
        res(r) <- 0.036
    }
    return(r)
}

timestep <- "10min" # "avg"
threshold <- 1000
land_simple <- readOGR(dsn=paste(indatadir,"ancils",sep=""), layer="land_rp")

mydata <- getMyData(timestep=timestep, var="lsrain", overwrite=F)
rb5216.1km.std <- adjCoords(mydata[[1]])
rb5216.4km.std <- adjCoords(mydata[[2]])
rb5216.4km.50k <- adjCoords(mydata[[3]])
rb5216.4km.300k <- adjCoords(mydata[[4]])

# u = 1km std veg
# s = 4km std veg
# w = 4km smoothed 50km veg
# y = 4km smoothed 300km veg

# models <- c("avg5216.4km.std", "avg5216.1km.std", "avg5216.4km.300k", "avg5216.4km.50k")
# models <- c("rb5216.4km.std", "rb5216.1km.std", "rb5216.4km.300k", "rb5216.4km.50k")
# id <- c("s","u","y","w")
models <- c("rb5216.4km.std", "rb5216.4km.50k", "rb5216.4km.300k") # [2:3]
id <- "s" # c("s","w","y")#[2:3]
testing <- FALSE # TRUE

if (testing == FALSE){
    # Load veg, orog and landfrac
    myveg <- adjCoords(loadVeg(model.nm=models[1], overwrite=F))
    myorog <- loadOtherAncils(model.nm=models[1], ancil="orog", overwrite=F)
    mylandfrac <- loadOtherAncils(model.nm=models[1], ancil="landfrac", overwrite=F)
    
    # Prepare boxes 
    myboxes <- makeBoxes2(xmin=354, ymin=-8, xmax=374, ymax=4, xdim=4, ydim=4, temp=myorog)
    spp.r <- myboxes[[1]] # non-adjusted raster
    sppa.r <- myboxes[[2]] # adjusted raster
    spp <- myboxes[[3]] # non-adjusted polygons
    sppa <- myboxes[[4]] # adjusted polygons
    
    # Prepare vegetation data for use in clipping etc
    allveg <- vegPrep(model.nm=models[1], id=id[1], myveg, myorog, mylandfrac, land_simple, sppa, spp.r, plots=F, vegThreshold=0.3, overwrite=F) # return(mycl, mycl.z) and creates pdf plots
    mycl <- allveg[[1]]
}


for (x in 1:length(id)){
	print(models[x])
	# load data
	inbr <- get(models[x])
	
# 	allpatch <- mypatch(threshold=1000, inbr=inbr, id=id[x], land_simple, timestep=timestep, indatadir=indatadir, dlresultsdir=dlresultsdir, spp=spp, sppa=sppa, overwrite=T)
# 	mcs <- allpatch[[1]]
# 	pop <- allpatch[[2]]
	mcs <- adjCoords(brick(paste(indatadir,"djzx",id[x],"/patches/",threshold,"km_",timestep,"/allmcs.1000km.vrt", sep="")))
	
# 	pop.stats(pop=pop, id=id[x], threshold=threshold, timestep="10min", indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=F)
    
# 	mymcsstats <- mcs.stats(mcs=mcs, precip=inbr, id=id[x], threshold=1000, timestep="10min", indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=F)
	
#     initiations(inpop=pop, id=id[x], threshold=1000, timestep="10min", indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=T)
    
	results <- tracker2(threshold=threshold, mcs, inbr, id=id[x], timestep=timestep, indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=F)

    # Check MCS tracks match MCS blocks
#     trackCheck(results, mcs, inbr, id=id[x], threshold=threshold, timestep=timestep)
    
# 	mcsfollow <- mcs.stats(mcs=mcs, precip=inbr, id=id[x], threshold=1000, timestep=timestep, indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=F)
    
	# Create plots of diurnal cycle
# 	diurnalcycle(inbr, type="all", patch=F, model.nm=models[x], id=id[x], spp.r=spp.r, sppa=sppa, mycl=mycl, land_simple, overwrite=T) # Creates pdf plots
# 	diurnalcycle(inbr, type="pop", patch=pop, model.nm=models[x], id=id[x], spp.r=spp.r, sppa=sppa, mycl=mycl, land_simple, overwrite=T) # Creates pdf plots
	
    # Consecutive dry and wet days ...
# 	consec_drywet(inbr, id=id[x], indatadir=indatadir, scratchdir=scratchdir, overwrite=T)
    
    # Within MCS intensity
    source("mcsIntensity.R")
#     mcsIntensity(results, myveg, mcs, inbr)
}

# Plot time series of avg.5216 for each model ...
# makeMovie(rb5216.4km.std, rb5216.4km.std, rb5216.4km.50k, rb5216.4km.300k, land_simple, sppa, dlresultsdir=dlresultsdir, models=models, boxes=F, type="precip", overwrite=T)
