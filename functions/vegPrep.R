# Prepare land cover masks etc ...

vegPrep <- function(model.nm, id="s", myproj="ll", myveg, myorog, mylandfrac, land_simple, sppa, spp.r, plots=F, vegThreshold=0.3, nBndClass=1, overwrite=T){
	
    require(rgdal)
    require(raster)
    
	print("Running vegPrep ...")
	if (Sys.info()[["sysname"]] == "Darwin"){
		indatadir <- "/Users/ajh235/Work/DataLocal/ModelData/WAFR/"
		resultsdir <- "/Users/ajh235/Work/myprojects/InternalSabbatical/Results/"
		scratchdir <- "/Users/ajh235/Work/Scratch/"
    } else {
		indatadir <- "/data/local/hadhy/ModelData/WAFR/"
		resultsdir <- "/home/h02/hadhy/myprojects/InternalSabbatical/Results/"
		scratchdir <- "/data/local/hadhy/Scratch/"
		require(PP,lib.loc="/myproject/ukmo/rhel6/R")
	}
    
    myLUT <- data.frame(ID=c(1,2,3,4,5,6,7), Landcover=factor(c("tree", "grass", "sparse", "boundary", "boundary, tree", "boundary, grass", "orography"), levels=c("tree", "grass", "sparse", "boundary", "boundary, tree", "boundary, grass", "orography")[c(3,2,6,4,5,1,7)]), Colours=c("dark green", "yellow", "orange", "sienna", "red", "blue", "dark grey"), plotOrder=c(4,2,1,3,5,6,7))

	outdir <- paste(indatadir,"djzx",id,"/derived/",sep="")
	if (!file.exists(outdir)){ dir.create(outdir) }

	km <- unlist(strsplit(model.nm, ".", fixed=T))[2]
	km.num <- unlist(strsplit(km, ""))[1]
	veg <- unlist(strsplit(model.nm, ".", fixed=T))[3]
    
    if (myproj == "ll"){
        require(maptools)
        # Very simple boundaries
        # land_simple <- get(data(wrld_simpl))
        # Complex boundaries ...
        land_simple <- getData("countries")
    } 
    if (myproj == "rp"){
        land_simple <- readOGR(dsn=paste(indatadir,"ancils",sep=""), layer="land_rp")
    }
    
	# Print some plots to pdf 
    if (plots==T){
        pdf(paste(resultsdir,"vegPrep_plots_",myproj,".",km,".",veg,".pdf", sep=""), width=12, height=8, onefile=T)
    }
	    
	# Calculate tree and grass fractions ...
	tree <- calc(subset(myveg, c(1,2,5)), fun=sum)
	grass <- calc(subset(myveg, c(3,4)), fun=sum)
    bare <- myveg[[8]]

    # Creating orog mask...
	orogmsk <- myorog > 500
	orogmsk[orogmsk == 1] <- NA
    
	if (plots==T){
    	print("Plotting orography")
    	plot(mask(myorog,orogmsk), col=terrain.colors(22), main="Orography <= 500m", breaks=seq(-25,500,25), xlim=c(-12,10), ylim=c(4,18))
    	plot(land_simple, add=T)
    	plot(sppa, add=T)
    	centroids <- coordinates(spp)
    	text(centroids, labels=1:nrow(centroids), cex=3)
	}
    
	# Choose 30% threshold of vegetation (may need to change this) ...
	print("Calculating tree/grass fractions ...")
	vt <- vegThreshold
    tr.cl <- tree > vt
	gr.cl <- (tr.cl != 1 & grass > vt)
	mix.cl <- (grass > 0 & grass <= vt & tree > 0 & tree <= vt & bare > vt) # Actually sparse veg when threshold set to 0.3
# 	mix.cl <- (grass > 0 & grass <= 0.3 & tree > 0 & tree <= 0.3)
	all.cl <- stack(tr.cl, gr.cl, mix.cl)
	names(all.cl) <- c("Tree", "Grass", "Mix")
	mylandfrac[(mylandfrac < 1)] <- NA
	mycl <- raster(tr.cl)
	mycl[tr.cl == 1] <- 1
	mycl[gr.cl == 1] <- 2
	mycl[mix.cl == 1]<- 3 # Sparse

	# Mask out orography from zones ...
# 	spp.r <- raster::mask(spp.r, orogmsk)
	if (plots==T){	    
    	plot(spp.r, xlim=c(-12,10), ylim=c(4,18))
    	plot(land_simple, add=T)
    	plot(sppa, add=T)
    	text(centroids, labels=1:nrow(centroids), cex=3)
	}
    
	# Fractions of each class per box (using land and orography mask)
	if (plots==T){
        print("Barplots of fractions per zone ...")
    	barplot(t(zonal(all.cl, spp.r, na.rm=T, fun='sum')[,2:4] / zonal(mylandfrac, spp.r, na.rm=T, fun='sum')[,2]), beside=T, names.arg=1:nrow(centroids), col=c("dark green","yellow","brown"), xlab="Zones")
    	legend(x=59, y=0.9, xpd=T, c("tree","grass","mix"), fill=c("dark green","yellow","brown"))
	}
    
	# Where are the tree:grass boundaries?
	print("Detecting boundaries ...")
    myclfile <- paste(outdir,"mycl_masked_",myproj,"_",vt,".tif",sep="")
	if (!file.exists(myclfile) | overwrite==T){
        # We're only interested in forest-grass boundaries, so recode everything to forest or grass
        focalin <- mycl
        focalin[focalin == 3] <- 2
		mybnds <- focal(focalin, w=matrix(1, nrow=3, ncol=3), fun=function(x){z <- table(x); ifelse(z[1]>2 & z[2]>2,3,0)}, filename=paste(outdir,"mycl_focal_",myproj,"_",vt,".tif",sep=""), format="GTiff", overwrite=T)
		mybnds.buf <- focal(mybnds, w=matrix(1, nrow=3, ncol=3), fun=max, na.rm=T, filename=paste(outdir,"mybnds_focal_",myproj,"_",vt,".tif",sep=""), format="GTiff", overwrite=T)
		mycl[mybnds == 3] <- 4
		mycl[mybnds.buf == 3 & mycl == 1] <- 5 # Tree close to boundary
		mycl[mybnds.buf == 3 & mycl == 2] <- 6 # Grass close to boundary
        mycl[myorog > 500]  <- 7 # Orography
		# Mask out sparse veg in zones 6 to 15
		mycl[mycl==3 & spp.r > 5] <- NA
        # Mask out orography in zone 4
		mycl[mycl==7 & spp.r == 4] <- NA
        mycl <- writeRaster(mycl, filename=myclfile, format="GTiff", overwrite=T)
# 		mycl <- mask(mycl, orogmsk, filename=myclfile, format="GTiff", overwrite=T)
	} else {
		mycl <- raster(myclfile)
	}

	if (plots==T){
    	print("Maps of vegetation zones and boxes ...")
#     	plot(shift(mycl, x=-360), col=c("green", "yellow", "orange", "brown","red","blue", "grey"), main="Vegetation boundaries masked by orography > 500m", legend=F, xlim=c(-10,16), ylim=c(-10,6))
#     	legend(x=17, y=5, c("tree","grass","sparse","boundary","boundary tree","boundary grass", "orogrpahy"), fill=c("green", "yellow", "orange", "brown","red","blue", "grey"), xpd=T)
#     	plot(land_simple, add=T)
#     	plot(sppa, add=T)
#     	text(centroids, labels=1:nrow(centroids), cex=3)
#     	
#         browser()
        mycl.f <- as.factor(mycl)
        ftab <- levels(mycl.f)[[1]]
        ftab$Name <- myLUT$Landcover
        levels(mycl.f) <- ftab
        print(
            levelplot(mycl.f, par.settings=rasterTheme(region=myLUT$Colours), xlab=NULL, ylab=NULL, xlim=c(-12,10), ylim=c(4,18), main="Vegetation classes and zones") + # , scales=list(draw=FALSE)
                latticeExtra::layer(sp.polygons(land_simple, col="black", lty=2)) + 
                latticeExtra::layer(sp.polygons(spp)) +
                latticeExtra::layer(sp.text(loc=coordinates(centroids), txt=1:nrow(centroids), cex=3)) 
            )
    	
        
    	dev.off()
	}
    
	# Combine veg/bnd classes with zones
	print("Creating mycl.z")
    mycl.zfile <- paste(outdir,"mycl.z_",myproj,"_",vt,".tif",sep="")
    if (nBndClass==1){
        mycl[mycl==5] <- 4
        mycl[mycl==6] <- 4
    }

	if (!file.exists(mycl.zfile) | overwrite==T){
		mycl.z <- overlay(mycl,spp.r*10, fun=sum, filename=paste(outdir,"mycl.z.tif",sep=""), overwrite=T, format="GTiff")
	} else {
		mycl.z <- raster(mycl.zfile)
	}

	return(list(mycl, mycl.z))
}