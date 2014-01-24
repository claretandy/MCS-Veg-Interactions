# First attempt at identifying patches ...
require(raster)
require(rasterVis)
# require(PP,lib.loc="/project/ukmo/rhel6/R")
require(rgdal)
require(sm)
require(plotrix)
require(RColorBrewer)
require(classInt)

mypatch <- function(threshold=100, inbr=inbr, id="s", land_simple, timestep="avg", indatadir="/data/local/hadhy/ModelData/WAFR/", dlresultsdir="/data/local/hadhy/Projects/InternalSabbatical/Results/", spp=spp, sppa=sppa, overwrite=T){
	
	print("Running patch detection ...")
	
	# Set indir and outdir
	datadir <- paste(indatadir,"djzx",id,"/",sep="") 
	outdir <- paste(indatadir,"djzx",id,"/patches/", threshold, "km_", timestep, "/", sep="")
	resultsdir <- paste(dlresultsdir,"Clumps/djzx",id,"_",timestep,"/",sep="")
	if (!file.exists(outdir)){ dir.create(outdir, recursive=T) }
	if (!file.exists(resultsdir)){ dir.create(resultsdir, recursive=T) }

	for (f in 1:nlayers(inbr)){
		mydate <- getZ(inbr)[f]
		print(mydate)
		# Set output filenames
		mcs.f <- paste(outdir, "mcs.",threshold,"km_", format(mydate, "%d%H.%M"),".tif",sep="")
		pop.f <- paste(outdir, "pop.",threshold,"km_", format(mydate, "%d%H.%M"),".tif",sep="")
		
		if (!file.exists(mcs.f) | overwrite==T){
		
			if (timestep=="10min"){
				# Convert precip from kg/m2/s to total mm in last 10 minutes
				r <- calc(inbr[[f]], fun=function(x){x*60*10}, filename=paste(scratchdir,"myrain.tif",sep=""), format="GTiff", overwrite=T)
#                 browser()
				# Do the clump. 10min precip > 0.16667 ~ 1mm per hour (i.e. 10 * 1/60)
				myclump <- clump(r > 0.16667, directions=4, gaps=F)
			}
			if (timestep=="avg"){
				# avg precip is average rainfall rate in the last 1 hour, expressed in kg/m2/second
	# 			r <- inbr[[f]]
				r <- calc(inbr[[f]], fun=function(x){x*60*60}, filename=paste(scratchdir,"myrain.tif",sep=""), format="GTiff", overwrite=T)
				#if(f == 13){ browser() }
				# Do the clump. avg precip > 1mm per hour 
				r.gt1 <- r>1
				
				myclump <- clump(r.gt1, directions=4, gaps=F)
				# browser()
			}
			
			# Get the approximate number of kmsq per clump. 
			# 1 grid cell ~ 4*4km = 16 sqkm
			# 1 grid cell ~ 1.5*1.5km = 2.25sqkm
			sqkm <- ifelse(id=="u", 2.25, 16) 
			mytab <- as.data.frame(table(getValues(myclump)))
			mytab$Freq <- mytab$Freq * sqkm
			# Create empty raster
			myr <- raster(myclump)
			# Match raster values to table values, and return Freq
			myrasterdata <- mytab[match(getValues(myclump), as.integer(mytab$Var1)),"Freq"]
			# Set values to the empty raster
			myr <- setValues(myr, myrasterdata)
			# Clumps > 1000kmsq
			mcs <- myr > threshold
			pop <- myr <= threshold
			writeRaster(mcs, filename=mcs.f, format="GTiff", overwrite=T)
			writeRaster(pop, filename=pop.f, format="GTiff", overwrite=T)
			
			# Write jpeg out
			jpeg(paste(resultsdir,"/clumplot_",threshold,"km_",format(mydate, "%d%H.%M"),".jpg",sep=""), quality=100, width=900, height=600)
	
			# Set plot bounds
			xbnd <- c(extent(spp)@xmin,extent(spp)@xmax)-360
			ybnd <- c(extent(spp)@ymin,extent(spp)@ymax)
			
			# Plot mcs (pop are 0 areas)
# 			plot(shift(mcs, x=-360), main=format(mydate, "%b%d %H:%M"), xlim=xbnd, ylim=ybnd, col=c("green","orange"), legend=F )
# 			legend(x=0, y=-9, horiz=T, c("MCS","POP"), fill=rev(c("green","orange")), xpd=T)
# 			plot(land_simple, add=T)
# 			plot(sppa, add=T)
			mcs.rat <- ratify(mcs)
            rat <- levels(mcs.rat)[[1]]
            rat$code <- c("Pop","MCS")
            levels(mcs) <- rat
            levelplot(shift(mcs, x=-360), margin=F, col.regions=c("green","orange"), main=format(mydate, "%b%d %H:%M")) + layer(sp.polygons(land_simple))
            
			# browser()
			dev.off()
		
		} # End of if(file.exists) etc...
	} # End of loop through layers of inbr

	allmcs.f <- paste(outdir,"allmcs.",threshold,"km.vrt",sep="")
	allpop.f <- paste(outdir,"allpop.",threshold,"km.vrt",sep="")
		
	if (!file.exists(allmcs.f)){
		system(paste("/Library/Frameworks/GDAL.framework/Versions/Current/Programs//gdalbuildvrt -separate ",allmcs.f," ",outdir,"mcs.",threshold,"km*",sep="")) }
	if (!file.exists(allpop.f)){
		system(paste("/Library/Frameworks/GDAL.framework/Versions/Current/Programs//gdalbuildvrt -separate ",allpop.f," ",outdir,"pop.",threshold,"km*",sep="")) }
	# browser()
	allmcs <- brick(allmcs.f); allmcs <- setZ(allmcs, getZ(inbr))
	allpop <- brick(allpop.f); allpop <- setZ(allpop, getZ(inbr))
	
	return(list(allmcs,allpop))
}
