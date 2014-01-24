require(raster)
require(rasterVis)
# require(PP,lib.loc="/project/ukmo/rhel6/R")
require(rgdal)
require(sm)
require(plotrix)
require(RColorBrewer)
require(classInt)

mcs.stats <- function(mcs=mcs, precip=precip, id="s", threshold=1000, timestep="10min", indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=T){
	
	print("Running mcsStats ...")
	
	# Set indir and outdir
	datadir <- paste(indatadir,"djzx",id,"/",sep="") 
	outdir <- paste(indatadir,"djzx",id,"/mcsStats/", threshold, "km_", timestep, "/", sep="")
	# resultsdir <- paste(dlresultsdir,"Clumps/djzx",id,"_",timestep,"/",sep="")
	if (!file.exists(outdir)){ dir.create(outdir, recursive=T) }
	# if (!file.exists(resultsdir)){ dir.create(resultsdir, recursive=T) }
	outfile <- paste(outdir, "allmcstot.tif", sep="")
# browser()
    
    # Make sure that we have dates in mcs brick ...
    if (is.null(getZ(mcs))){
        mcs <- setZ(mcs, getZ(precip))
    }
    
	# How many MCS and how many POP in each grid cell?
	if (!file.exists(outfile) | overwrite==T){
		mcstot <- mcs[[1]]; mcstot[is.na(mcstot)] <- 0
		writeRaster(mcstot, filename=paste(outdir, "mcstot_",format(getZ(mcs)[1], "%d%H.%M"), ".tif", sep=""), overwrite=T)
		for (x in 2:nlayers(mcs)){
			mydate <- getZ(mcs)[x]
			thisMCS.f <- paste(outdir, "mcstot_",format(mydate, "%d%H.%M"),".tif", sep="")
			if(!file.exists(thisMCS.f) | overwrite==T){
				print(mydate)
				thisMCS <- mcs[[x]]
				thisMCS[is.na(thisMCS)] <- 0
				mcstot <- overlay(mcstot, thisMCS, fun=sum, filename=thisMCS.f, format="GTiff", overwrite=T)
			} else {
				mcstot <- raster(thisMCS.f)
			}
		}
		writeRaster(mcstot, filename=outfile, format="GTiff", overwrite=T)
	} else {
		mcstot <- raster(outfile)
	}
	
	# How much precipitation per MCS in each grid cell?
	# Set indir and outdir
    print("MCS precip statistics ...")
#     overwrite=T
	outdir <- paste(indatadir,"djzx",id,"/mcsStats/", threshold, "km_", timestep, "/", sep="")
	if (!file.exists(outdir)){ dir.create(outdir, recursive=T) }
	outfile <- paste(outdir, "allmcs.totpreciprate.tif", sep="")
	
    mydate <- getZ(mcs)[1]
    mcstot <- mcs[[1]]; mcstot[is.na(mcstot)] <- 0
    mcstotpr <- overlay(mcstot, precip[[1]], fun=function(x,y){x*y}, filename=paste(outdir, "mcs.precip.mmsec_",format(mydate, "%d%H.%M"),".tif", sep=""), overwrite=T) # This outfile is correct
	if (!file.exists(outfile) | overwrite==T){
    for (x in 1:nlayers(mcs)){
        
        mydate <- getZ(mcs)[x]
        thisMCS.f <- paste(outdir, "mcs.precip.mmsec_",format(mydate, "%d%H.%M"),".tif", sep="")
#                 thisMCStotpr.f <- paste(outdir, "mcs.accumprecip.mmsec_",format(mydate, "%d%H.%M"),".tif", sep="")
#         if (!file.exists(thisMCS.f) | overwrite==T){
        thisPrecip <- precip[[x]]
        if(!file.exists(thisMCS.f) | overwrite==T){
            print(mydate)
            thisMCS <- mcs[[x]]
            thisMCS[is.na(thisMCS)] <- 0
            # Overlay MCS and precip at timestep x
            thisMCSprecip <- overlay(thisMCS, thisPrecip, fun=function(x,y){x[is.na(x)] <- 0; z <- x*y}, filename=thisMCS.f, format="GTiff", overwrite=T) # Precip rate per grid cell

        } else {
            thisMCSprecip <- raster(thisMCS.f)
# 	            mcstotpr <- raster(thisMCStotpr.f)
        }
        # Add MCS precip at this timestep to tot rate 
        mcstotpr <- overlay(thisMCSprecip, mcstotpr, fun=function(x,y){return(x+y)}, filename=outfile, format="GTiff", overwrite=T) # Precip rate per grid cell
    }
# 	    writeRaster(mcstotpr, filename=outfile, format="GTiff", overwrite=T)
	} else {
	    mcstotpr <- raster(outfile)
	}
    
    # Put both parts together ...
    mcs.avpr <- overlay(mcstot, mcstotpr, fun=function(x,y){return(y/x)}, filename=paste(outdir, "allmcs.meanpreciprate.tif", sep=""), overwrite=T)
	return(list(mcstot, mcstotpr, mcs.avpr))
}