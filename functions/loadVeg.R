# Loads veg fraction from ancillary files for each model run
loadVeg <- function(model.nm="avg5216.4km.std", proj='ll', overwrite=F, doplot=F){
    
	if (Sys.info()[["sysname"]] == "Darwin"){
		indatadir <- "/Users/ajh235/Work/DataLocal/ModelData/WAFR/"
		resultsdir <- "/Users/ajh235/Work/Projects/InternalSabbatical/Results/"
	} else {
		indatadir <- "/data/local/hadhy/ModelData/WAFR/"
		resultsdir <- "/home/h02/hadhy/Projects/InternalSabbatical/Results/"
		require(PP,lib.loc="/project/ukmo/rhel6/R")
	}

	print("Loading vegetation fractions ...")
	
	require(raster)
	require(rasterVis)
	source("~/Scripts/R/readPP.R")
    
	km <- unlist(strsplit(model.nm, ".", fixed=T))[2]
	km.num <- unlist(strsplit(km, ""))[1]
	veg <- unlist(strsplit(model.nm, ".", fixed=T))[3]
	file <- paste(indatadir,"ancils/km", km.num, "/qrparm.veg.frac_",km,".",veg, ".pp", sep="")
    file_ll <- gsub(".pp", "_ll.nc", x=file)
    file_rp <- gsub(".pp", ".nc", x=file)
    
    if (proj == 'll'){
        if (!file.exists(file_ll) | overwrite==T){
            system(paste('activatescitools ; python2.7 repro.py ',file, sep=""), intern=T)
        }
        myveg <- brick(file_ll)
    }
    
    if (proj == 'rp'){
        if (!file.exists(file_ll) | overwrite==T){
            myveg <- readPP(file, outfile=file_rp)
        } else {
            myveg <- brick(file_rp)
        }
    }
    
    if (doplot == T){
        png(paste(resultsdir,"vegfrac_maps_",km,".",veg,".png", sep=""), width=1000, height=1000)
        gntheme <- rasterTheme(pch=19, region=brewer.pal(9, 'YlGn'))
        names(myveg) <- c("BL", "NL", "C3 grass", "C4 grass", "Shrub", "Urban", "Water", "Bare", "Ice")
        print(levelplot(myveg, margins=F, par.settings=gntheme))
        dev.off()        
    }
	return(myveg)
}
