loadOtherAncils <- function(model.nm="avg5216.4km.std", ancil="orog", proj="ll", overwrite=T){
    source("~/Scripts/R/readPP.R")
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
	file_ll <- gsub(".pp", "_ll.nc", x=file)
	file_rp <- gsub(".pp", "_rp.nc", x=file)
	
	if (proj == 'll'){
	    if (!file.exists(file_ll) | overwrite==T){
	        system(paste('source $HOME/scitools/bin/activate ; python repro.py "',file, '"', sep=""), intern=T)
	    }
	    myanc <- raster(file_ll)
	}
	
	if (proj == 'rp'){
	    if (!file.exists(file_rp) | overwrite==T){
	        myanc <- readPP(file)
	    } else {
	        myanc <- raster(file_rp)
	    }
	}

	return(myanc)
}