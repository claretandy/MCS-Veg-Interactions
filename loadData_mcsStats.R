# Load and process data from all 4 model variants
# Need to re-work many of old functions so that they read and write to/from the correct locations.

library(raster)
library(rasterVis)
library(rgdal)
source("functions/consec_drywet.R")
source("functions/diurnal_cycle.R")
source("functions/loadVeg.R")
source("functions/loadOtherAncils.R")
source("functions/makeBoxes.R") 
source("functions/vegPrep.R")
source("functions/patches.R")
source("functions/movies.R")
source("functions/tracker.R")
source("functions/mcsStats.R")
source("functions/popStats.R")
source("functions/initiations.R")
source("getMyData.R")

if (Sys.info()[["sysname"]] == "Darwin"){
	indatadir <- "/Volumes/MYBOOK/WAFR/"
	dlresultsdir <- "/Users/Andy/Work/InternalSabbatical/Results/"
	resultsdir <- "/Users/Andy/Work/InternalSabbatical/Results/"
	scratchdir <- "/Users/Andy/Work/Scratch/"
} else {
	indatadir <- "/data/local/hadhy/ModelData/WAFR/"
	dlresultsdir <- "/data/local/hadhy/Projects/InternalSabbatical/Results/"
	resultsdir <- "/home/h02/hadhy/Projects/InternalSabbatical/Results/"
	scratchdir <- "/data/local/hadhy/Scratch/"
	require(PP,lib.loc="/project/ukmo/rhel6/R")
}

setOptions(tmpdir=scratchdir, todisk=F)

land_simple <- readOGR(dsn=paste(indatadir,"ancils",sep=""), layer="land_rp")

mydata <- getMyData(timestep="10min", overwrite=F)
rb5216.1km.std <- mydata[[1]]
rb5216.4km.std <- mydata[[2]]
rb5216.4km.50k <- mydata[[3]]
rb5216.4km.300k <- mydata[[4]]

plot.hr <- function(in1k, in4k, in4k.50, in4k.300){
	in1k.std <- setMinMax(avg5216.1km.std[[18]]*60*60*24)
	in4k.std <- setMinMax(avg5216.4km.std[[18]]*60*60*24)
	in4k.50k <- setMinMax(avg5216.4km.50k[[18]]*60*60*24)
	in4k.300k <- setMinMax(avg5216.4km.300k[[18]]*60*60*24)
	
	maxprecip <- c(maxValue(in1k.std), maxValue(in4k.std), maxValue(in4k.50k), maxValue(in4k.300k))
	print(maxprecip)
	par(mfrow=c(2,2))
	plot(shift(in1k.std, x=-360), main="1km Standard Vegetation"); plot(land_simple, add=T)
	plot(shift(in4k.std, x=-360), main="4km Standard Vegetation"); plot(land_simple, add=T)
	plot(shift(in4k.50k, x=-360), main="4km Smoothed Vegetation (50km)"); plot(land_simple, add=T)
	plot(shift(in4k.300k, x=-360), main="4km Smoothed Vegetation (300km)"); plot(land_simple, add=T)
}


# plot.hr(avg5216.1km.std[[18]], avg5216.4km.std[[18]], avg5216.4km.std[[18]], avg5216.4km.std[[18]])

# u = 1km std veg
# s = 4km std veg
# w = 4km smoothed 50km veg
# y = 4km smoothed 300km veg

# models <- c("avg5216.4km.std", "avg5216.1km.std", "avg5216.4km.300k", "avg5216.4km.50k")
# models <- c("rb5216.4km.std", "rb5216.1km.std", "rb5216.4km.300k", "rb5216.4km.50k")
# id <- c("s","u","y","w")
models <- c("rb5216.4km.std", "rb5216.4km.50k", "rb5216.4km.300k")
id <- c("s","w","y")

for (x in 1:length(id)){
	print(models[x])
	# load data
	inbr <- get(models[x])
	
	# Load veg, orog and landfrac
# 	myveg <- loadVeg(model.nm=models[1], overwrite=F)
	myorog <- loadOtherAncils(model.nm=models[1], ancil="orog", overwrite=F)
# 	mylandfrac <- loadOtherAncils(model.nm=models[1], ancil="landfrac", overwrite=F)
	
	# Create zones ...
# 	myboxes <- makeBoxes(4, myorog) # DON'T USE! 0.5, 1, 2, 4 degree box sizes
	myboxes <- makeBoxes2(xmin=354, ymin=-8, xmax=374, ymax=4, xdim=4, ydim=4, temp=myorog)
	spp.r <- myboxes[[1]] # non-adjusted raster
	sppa.r <- myboxes[[2]] # adjusted raster
	spp <- myboxes[[3]] # non-adjusted polygons
	sppa <- myboxes[[4]] # adjusted polygons

	allpatch <- mypatch(threshold=1000, inbr=inbr, id=id[x], land_simple, timestep="10min", indatadir=indatadir, dlresultsdir=dlresultsdir, spp=spp, sppa=sppa, overwrite=F)
	
	mcs <- allpatch[[1]]
	pop <- allpatch[[2]]
	
# 	pop.stats(pop=pop, id=id[x], threshold=1000, timestep="10min", indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=F)
    
	mymcsstats <- mcs.stats(mcs=mcs, precip=inbr, id=id[x], threshold=1000, timestep="10min", indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=T)
	
#     initiations(inpop=pop, id=id[x], threshold=1000, timestep="10min", indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=T)
    
	#tracker(threshold=1000, mcs, inbr, id=id[x], timestep="10min", indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=T)
	# Prepare vegetation data for use in clipping etc
# 	allveg <- vegPrep(model.nm=models[1], id=id[1], myveg, myorog, mylandfrac, land_simple, sppa, spp.r, plots=F, overwrite=F) # return(mycl, mycl.z) and creates pdf plots
	
	# Create plots of diurnal cycle
# 	diurnalcycle(inbr, model.nm=models[x], id=id[x], spp.r, sppa, allveg[[1]], land_simple, overwrite=T) # Creates pdf plots
	# Consecutive dry and wet days ...
# 	consec_drywet(inbr, id=id[x], indatadir=indatadir, scratchdir=scratchdir, overwrite=T)
}

# Plot time series of avg.5216 for each model ...
# makeMovie(rb5216.4km.std, rb5216.4km.std, rb5216.4km.50k, rb5216.4km.300k, land_simple, sppa, dlresultsdir=dlresultsdir, models=models, boxes=F, type="precip", overwrite=T)
