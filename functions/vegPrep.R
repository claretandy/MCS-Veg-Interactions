# Prepare land cover masks etc ...

vegPrep <- function(model.nm, id="s", myveg, myorog, mylandfrac, land_simple, sppa, spp.r, plots=F, overwrite=T){
	
	print("Running vegPrep ...")
	if (Sys.info()[["sysname"]] == "Darwin"){
		indatadir <- "/Users/ajh235/Work/DataLocal/ModelData/WAFR/"
		resultsdir <- "/Users/ajh235/Work/Projects/InternalSabbatical/Results/"
		scratchdir <- "/Users/ajh235/Work/Scratch/"
    } else {
		indatadir <- "/data/local/hadhy/ModelData/WAFR/"
		resultsdir <- "/home/h02/hadhy/Projects/InternalSabbatical/Results/"
		scratchdir <- "/data/local/hadhy/Scratch/"
		require(PP,lib.loc="/project/ukmo/rhel6/R")
	}

	
	outdir <- paste(indatadir,"djzx",id,"/derived/",sep="")
	if (!file.exists(outdir)){ dir.create(outdir) }

	km <- unlist(strsplit(model.nm, ".", fixed=T))[2]
	km.num <- unlist(strsplit(km, ""))[1]
	veg <- unlist(strsplit(model.nm, ".", fixed=T))[3]

	# Print some plots to pdf 
    if (plots==T){
        pdf(paste(resultsdir,"vegPrep_plots.",km,".",veg,".pdf", sep=""), width=12, height=8, onefile=T)
    }
	
	
	# Calculate tree and grass fractions ...
	tree <- calc(subset(myveg, c(1,2,5)), fun=sum)
	grass <- calc(subset(myveg, c(3,4)), fun=sum)

    # Creating orog mask...
	orogmsk <- myorog > 500
	orogmsk[orogmsk == 1] <- NA
    
	if (plots==T){
    	print("Plotting orography")
    	plot(shift(mask(myorog,orogmsk), x=-360), col=terrain.colors(22), main="Orography <= 500m", breaks=seq(-25,500,25), xlim=c(-10,16), ylim=c(-10,6))
    	plot(land_simple, add=T)
    	plot(sppa, add=T)
    	centroids <- coordinates(sppa)
    	text(centroids, labels=1:nrow(centroids), cex=3)
	}
    
	# Choose 50% threshold of vegetation (may need to change this) ...
	print("Calculating tree/grass fractions ...")
	tr.cl <- tree > 0.5
	gr.cl <- grass > 0.5
	mix.cl <- (grass > 0.2 & grass <= 0.5 & tree > 0.2 & tree <= 0.5)
	all.cl <- stack(tr.cl, gr.cl, mix.cl)
	names(all.cl) <- c("Tree", "Grass", "Mix")
	mylandfrac[(mylandfrac < 1)] <- NA
	mycl <- raster(tr.cl)
	mycl[tr.cl == 1] <- 1
	mycl[gr.cl == 1] <- 2
	mycl[mix.cl == 1]<- 3

	# Mask out orography from zones ...
	spp.r <- mask(spp.r, orogmsk)
	if (plots==T){	    
    	plot(shift(spp.r, x=-360), xlim=c(-10,16), ylim=c(-10,6))
    	plot(land_simple, add=T)
    	plot(sppa, add=T)
    	text(centroids, labels=1:nrow(centroids), cex=3)
	}
    
	# Fractions of each class per box (using land and orography mask)
	if (plots==T){
        print("Barplots of fractions per zone ...")
    	barplot(t(zonal(all.cl, spp.r, na.rm=T, stat='sum')[,2:4] / zonal(mylandfrac, spp.r, na.rm=T, stat='sum')[,2]), beside=T, names.arg=1:nrow(centroids), col=c("dark green","yellow","brown"), xlab="Zones")
    	legend(x=59, y=0.9, xpd=T, c("tree","grass","mix"), fill=c("dark green","yellow","brown"))
	}
    
	# Where are the tree:grass boundaries?
	print("Detecting boundaries ...")
	if (!file.exists(paste(outdir,"mycl_masked.tif",sep="")) | overwrite==T){
		mybnds <- focal(mycl, w=3, fun=function(x){ifelse(table(x)[1]>2 & table(x)[2]>2,3,0)}, filename=paste(outdir,"mycl_focal.tif",sep=""), format="GTiff", overwrite=T)
		mybnds.buf <- focal(mybnds, w=3, fun=max, na.rm=T, filename=paste(outdir,"mybnds_focal.tif",sep=""), format="GTiff", overwrite=T)
		mycl[mybnds == 3] <- 4
		mycl[mybnds.buf == 3 & mycl == 1] <- 5 # Tree close to boundary
		mycl[mybnds.buf == 3 & mycl == 2] <- 6 # Grass close to boundary
		mycl <- mask(mycl, orogmsk, filename=paste(outdir,"mycl_masked.tif",sep=""), format="GTiff", overwrite=T)
	} else {
		mycl <- raster(paste(outdir,"mycl_masked.tif",sep=""))
	}
	
	if (plots==T){
    	print("Maps of vegetation zones and boxes ...")
    	plot(shift(mycl, x=-360), col=c("green", "yellow", "grey", "brown","red","blue"), main="Vegetation boundaries masked by orography > 500m", legend=F, xlim=c(-10,16), ylim=c(-10,6))
    	legend(x=17, y=5, c("tree","grass","mix","boundary","boundary tree","boundary grass"), fill=c("green", "yellow", "grey", "brown","red","blue"), xpd=T)
    	plot(land_simple, add=T)
    	plot(sppa, add=T)
    	text(centroids, labels=1:nrow(centroids), cex=3)
    	
    	dev.off()
	}
    
	# Combine veg/bnd classes with zones
	print("Creating mycl.z")
	if (!file.exists(paste(outdir,"mycl.z.tif",sep="")) | overwrite==T){
		mycl.z <- overlay(mycl,spp.r*10, fun=sum, filename=paste(outdir,"mycl.z.tif",sep=""), overwrite=T, format="GTiff")
	} else {
		mycl.z <- raster(paste(outdir,"mycl.z.tif",sep=""))
	}

	return(list(mycl, mycl.z))
}