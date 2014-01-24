loadOtherAncils <- function(model.nm="avg5216.4km.std", ancil="orog", overwrite=T){
	if (Sys.info()[["sysname"]] == "Darwin"){
	    indatadir <- "/Users/ajh235/Work/DataLocal/ModelData/WAFR/"
	    resultsdir <- "/Users/ajh235/Work/Projects/InternalSabbatical/Results/"
	} else {
		indatadir <- "/data/local/hadhy/ModelData/WAFR/"
		resultsdir <- "/home/h02/hadhy/Projects/InternalSabbatical/Results/"
		require(PP,lib.loc="/project/ukmo/rhel6/R")
	}
	
	print(paste("Loading ",ancil," ancillaries ...", sep=""))

	km <- unlist(strsplit(model.nm, ".", fixed=T))[2]
	km.num <- unlist(strsplit(km, ""))[1]
	veg <- unlist(strsplit(model.nm, ".", fixed=T))[3]
	file <- paste(indatadir,"ancils/km", km.num, "/qrparm.",ancil,"_",km,".pp", sep="")
	
	if (!file.exists(file) | overwrite==T){
		myanc <- readPP(file, outfile=gsub(".pp", ".tif", x=file))
		myanc <- subset(myanc, 1)
	} else {
		myanc <- raster(gsub(".pp", ".tif", x=file))
	}

	return(myanc)
}