require(raster)
require(rasterVis)
# require(PP,lib.loc="/project/ukmo/rhel6/R")
require(rgdal)
require(sm)
require(plotrix)
require(RColorBrewer)
require(classInt)

pop.stats <- function(pop=pop, id="s", threshold=1000, timestep="10min", indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=T){
	
	print("Running popStats ...")
	
	# Set indir and outdir
	datadir <- paste(indatadir,"djzx",id,"/",sep="") 
	outdir <- paste(indatadir,"djzx",id,"/popStats/", threshold, "km_", timestep, "/", sep="")
	# resultsdir <- paste(dlresultsdir,"Clumps/djzx",id,"_",timestep,"/",sep="")
	if (!file.exists(outdir)){ dir.create(outdir, recursive=T) }
	# if (!file.exists(resultsdir)){ dir.create(resultsdir, recursive=T) }
	outfile <- paste(outdir, "allpoptot.tif", sep="")
# browser()
	# How many pop and how many POP in each grid cell?
	if (!file.exists(outfile) | overwrite==T){
		poptot <- pop[[1]]; poptot[is.na(poptot)] <- 0
		writeRaster(poptot, filename=paste(outdir, "poptot_",format(getZ(pop)[1], "%d%H.%M"), ".tif", sep=""), overwrite=T)
		for (x in 2:nlayers(pop)){
			mydate <- getZ(pop)[x]
			thispop.f <- paste(outdir, "poptot_",format(mydate, "%d%H.%M"),".tif", sep="")
			if(!file.exists(thispop.f) | overwrite==T){
				print(mydate)
				thispop <- pop[[x]]
				thispop[is.na(thispop)] <- 0
				poptot <- overlay(poptot, thispop, fun=sum, filename=thispop.f, format="GTiff", overwrite=T)
			} else {
				poptot <- raster(thispop.f)
			}
		}
		writeRaster(poptot, filename=outfile, format="GTiff", overwrite=T)
	} else {
		poptot <- raster(outfile)
	}
	
	return(poptot)	
}