# First attempt at identifying patches ...
require(raster)
require(rasterVis)
# require(PP,lib.loc="/project/ukmo/rhel6/R")
require(rgdal)
require(sm)
require(plotrix)
require(RColorBrewer)
require(classInt)

initiations <- function(threshold=1000, inpop=pop, id="s", timestep="10min", indatadir="/data/local/hadhy/ModelData/WAFR/", dlresultsdir="/data/local/hadhy/Projects/InternalSabbatical/Results/", overwrite=T){
	
	print("Running the POP initiation tracker ...")
	
	# Set indir and outdir
	datadir <- paste(indatadir,"djzx",id,"/",sep="") 
	outdir <- paste(indatadir,"djzx",id,"/initiations/", threshold, "km_", timestep, "/", sep="")
	resultsdir <- paste(dlresultsdir,"Clumps/djzx",id,"_",timestep,"/",sep="")
	if (!file.exists(outdir)){ dir.create(outdir, recursive=T) }
	if (!file.exists(resultsdir)){ dir.create(resultsdir, recursive=T) }
	
	mydates <- getZ(inpop)
	
	# Start new ID list
	uniqID <- vector("numeric")
	t1 <- raster(inpop, layer=1)
#     browser()
	t1.clump <- clump(t1, directions=8, gaps=F, filename=paste(outdir,"IDs.", format(mydates[1], "%d%H.%M"), ".tif", sep=""), overwrite=overwrite)
    
	myfreq <- freq(t1.clump)
	myfreq <- myfreq[which(myfreq[,1] != "NA"),]
# 	results <- data.frame("ID"=myfreq[,1], "pixcount"=myfreq[,2], "timestep"=rep(mydates[1], nrow(myfreq)),"class"=rep("regularTracking", nrow(myfreq)) )
# 	results$class <- as.character(results$class) # Force it to be character rather than factor
	uniqID <- myfreq[,1]
	prev.clust <- t1.clump
	prev.clust[is.na(prev.clust)] <- 0

	for (f in 2:(nlayers(inpop)-1)){
	# for (f in 2:3){
		fromID <- vector("numeric")
		toID <- vector("numeric")
		mydate <- getZ(inpop)[f]
		print(mydate)
		# Set output filenames
		clus.f <- paste(outdir,"initiations.", format(mydates[f], "%d%H.%M"), ".tif", sep="")
		
		# Read minus1, current and plus1 time steps, and clump each 
		rt <- raster(inpop, layer=f) ; rt.clump <- clump(rt, directions=8, gaps=F); rt.clump[is.na(rt.clump)] <- 0
		rtm1.clump <- prev.clust; rtm1.clump[is.na(rtm1.clump)] <- 0
		
		# browser()
		# Crosstab backwards and forwards
		tm1xtab <- crosstab(rt.clump, rtm1.clump)
		
		# Remove zero overlaps
		tm1xtab.nozero <- tm1xtab[which(row.names(tm1xtab) != "0"),which(colnames(tm1xtab) != "0")]

		# Get pixel count for each class in t
		myfreq <- as.data.frame(freq(rt.clump))
		tm1freq <- as.data.frame(freq(rtm1.clump))
		
		# Classification of cases
		###################################################################
		### Generation: MCS in t that are not in t-1
		#### Create new IDs
		print("Generation ...")
		generation <- as.numeric(attributes(which(tm1xtab[,which(colnames(tm1xtab)=="0")] == apply(tm1xtab, 1, sum)))$names)
		if (length(generation)!=0){
            
			this.fromID <- generation
			this.toID <- seq(max(uniqID)+1, max(uniqID)+length(generation))
			fromID <- c(fromID, this.fromID)
			toID <- c(toID, this.toID)
			uniqID <- c(uniqID, this.toID)
			
#             browser()

            rectab<- data.frame("from"=fromID, "to"=toID)
			rt.clump.rc <- subs(rt.clump, rectab, by="from", which="to", filename=clus.f, overwrite=T)
			mypolys <- rasterToPolygons(shift(rt.clump.rc, x=-360), dissolve=T, na.rm=T, fun=function(x){x>0})
			mypts <- SpatialPointsDataFrame(coords=coordinates(mypolys), data=mypolys@data)
			writeOGR(mypts, dsn=outdir, layer=paste("centre_pts.", format(mydates[f], "%d%H.%M"), ".shp", sep=""), driver="ESRI Shapefile", check_exists=T, overwrite_layer=T)
		}
		
		prev.clust <- raster(clus.f)
	}
}