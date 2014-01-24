# First attempt at identifying patches ...
require(raster)
require(rasterVis)
# require(PP,lib.loc="/project/ukmo/rhel6/R")
require(rgdal)
require(sm)
require(plotrix)
require(RColorBrewer)
require(classInt)

tracker <- function(threshold=1000, inmcs, inprecip, id="s", timestep="10min", indatadir="/data/local/hadhy/ModelData/WAFR/", dlresultsdir="/data/local/hadhy/Projects/InternalSabbatical/Results/", overwrite=T){
	
	print("Running the MCS tracker ...")
	
	# Set indir and outdir
	datadir <- paste(indatadir,"djzx",id,"/",sep="") 
	outdir <- paste(indatadir,"djzx",id,"/tracker/", threshold, "km_", timestep, "/", sep="")
	resultsdir <- paste(dlresultsdir,"Clumps/djzx",id,"_",timestep,"/",sep="")
	if (!file.exists(outdir)){ dir.create(outdir, recursive=T) }
	if (!file.exists(resultsdir)){ dir.create(resultsdir, recursive=T) }
	
	mydates <- getZ(inprecip)
	
	# Start new ID list
	uniqID <- vector("numeric")
	t1 <- raster(inmcs, layer=1)
	t1.clump <- clump(t1, directions=8, gaps=F, filename=paste(outdir,"IDs.", format(mydates[1], "%d%H.%M"), ".tif", sep=""), overwrite=overwrite)
	myfreq <- freq(t1.clump)
	myfreq <- myfreq[which(myfreq[,1] != "NA"),]
	results <- data.frame("ID"=myfreq[,1], "pixcount"=myfreq[,2], "timestep"=rep(mydates[1], nrow(myfreq)),"class"=rep("regularTracking", nrow(myfreq)) )
	results$class <- as.character(results$class) # Force it to be character rather than factor
	uniqID <- myfreq[,1]

	for (f in 2:(nlayers(inmcs)-1)){
		fromID <- vector("numeric")
		toID <- vector("numeric")
		mydate <- getZ(inprecip)[f]
		print(mydate)
		# Set output filenames
		clus.f <- paste(outdir,"IDs.", format(mydates[f], "%d%H.%M"), ".tif", sep="")
		
		# Read minus1, current and plus1 time steps, and clump each 
		rt <- raster(inmcs, layer=f) ; rt.clump <- clump(rt, directions=8, gaps=F); rt.clump[is.na(rt.clump)] <- 0
		rtm1 <- raster(inmcs, layer=f-1) ; rtm1.clump <- clump(rtm1, directions=8, gaps=F); rtm1.clump[is.na(rtm1.clump)] <- 0
		rtp1 <- raster(inmcs, layer=f+1) ; rtp1.clump <- clump(rtp1, directions=8, gaps=F); rtp1.clump[is.na(rtp1.clump)] <- 0
		
		# Read inprecip for current time step
		pr <- raster(inprecip, layer=f)
		
		# Crosstab backwards and forwards
		tm1xtab <- crosstab(rt.clump, rtm1.clump)
		tp1xtab <- crosstab(rt.clump, rtp1.clump)
		# tm1xtab <- crosstab(rt.clump, rtm1.clump)
		# tp1xtab <- crosstab(rt.clump, rtp1.clump)
		
		# For each clump in r. clump, which overlaps the most (backwards and forwards)
		# tm1.overlap <- apply(tm1xtab, 1, which.max)
		# tp1.overlap <- apply(tp1xtab, 1, which.max)
		
		# Remove zero overlaps
		tm1xtab.nozero <- tm1xtab[-1,-1]

		# Get pixel count for each class in t
		myfreq <- freq(rt.clump)
		tm1freq <- freq(rtm1.clump)
		
		# Classification of cases
		###################################################################
		### Generation: MCS in t that are not in t-1
		#### Create new IDs
		print("Generation ...")
		generation <- as.numeric(attributes(which(tm1xtab[,1] == apply(tm1xtab, 1, sum)))$names)
		this.fromID <- generation
		this.toID <- seq(max(uniqID)+1, max(uniqID)+length(generation))
		fromID <- c(fromID, this.fromID)
		toID <- c(toID, this.toID)
		uniqID <- c(uniqID, this.toID)
		results <- rbind(results, data.frame("ID"=this.toID, "pixcount"=myfreq[which(myfreq[,"value"] %in% this.fromID), "count"], "timestep"=rep(mydate, length(this.fromID)),"class"=rep("generation", length(this.fromID)) ))
		
		###################################################################
		### Dissipation: MCS from t-1 that are not in t
		#### Change record in results for t-1 to dissipation
		print("Dissipation ...")
		dissipation <- as.numeric(attributes(which(tm1xtab[1,] == apply(tm1xtab, 2, sum)))$names)
		results[which(results$ID %in% dissipation),"class"] <- rep("dissipation", length(dissipation))
		
		
		###################################################################
		### Splitting: 1 MCS in t-1 becomes 2 or more in t. Largest MCS in t keeps ID from t-1
		print("Splitting 1 ...")
		splitting <- as.numeric(attributes(which(apply(tm1xtab.nozero, 2, FUN=function(x){x <- x[x>0]; length(x)}) > 1))$names) # IDs refer to t-1
		# These ones keep ID from t-1 (treat same as regulartrk)
		rowIDofsplits <- as.numeric(row.names(tm1xtab.nozero)[apply(tm1xtab.nozero[,splitting], 2, which.max)]) # Returns the ID of the row in which the max for each column exists ???? Checked: is it row ID?
		# If splitting.regular pulls out 2 clusters that are the same, for different clusters in t-1, then use the class from t-1 that has the highest number of pixels
		browser()
		if (anyDuplicated(splitting.regular) != 0){
			dupes <- anyDuplicated(splitting.regular)
			dupes.tab <- tm1freq[which(tm1freq[,"value"] %in% splitting[which(splitting.regular %in% dupes)]), ]
			not2keep <- dupes.tab[!which.max(dupes.tab[,"count"]),"value"]
			splitting <- splitting[!(splitting %in% not2keep)]
			splitting.regular <- as.numeric(row.names(tm1xtab.nozero)[as.numeric(apply(tm1xtab.nozero, 2, FUN=function(x){which.max(x)})[which(as.numeric(colnames(tm1xtab.nozero)) %in% splitting)])])
		}
		this.fromID <- splitting.regular
		this.toID <- splitting
		fromID <- c(fromID, this.fromID)
		toID <- c(toID, this.toID)
		browser()
		results <- rbind(results, data.frame("ID"=this.toID, "pixcount"=myfreq[which(myfreq[,"value"] %in% this.fromID), "count"], "timestep"=rep(mydate, length(this.fromID)),"class"=rep("regularTracking", length(this.fromID)) ))
		# These are true splits, and get new IDs
		# browser()
		print("Splitting 2 ...")
		splitting.true <- as.numeric(row.names(tm1xtab.nozero[as.numeric(which(apply(tm1xtab.nozero[,which(colnames(tm1xtab.nozero) %in% splitting)], 1, FUN=function(x){x <- x[x>0]; length(x)}) > 0)),])) # IDs refer to t
		# Remove IDs 
		splitting.true <- splitting.true[which(!(splitting.true %in% splitting.regular))]
		this.fromID <- splitting.true
		this.toID <- seq(max(uniqID)+1, max(uniqID)+length(this.fromID))
		uniqID <- c(uniqID, this.toID)
		fromID <- c(fromID, this.fromID)
		toID <- c(toID, this.toID)
		results <- rbind(results, data.frame("ID"=this.toID, "pixcount"=myfreq[which(myfreq[,"value"] %in% this.fromID), "count"], "timestep"=rep(mydate, length(this.fromID)),"class"=rep("splitting", length(this.fromID)) ))
		
		###################################################################
		### Merging: 2 or more MCS in t-1 becomes 1 in t. Largest MCS in t-1 gives ID to MCS in t
		#### Identify clusters in t ...
		print("Merging ...")
		merging <- as.numeric(attributes(which(apply(tm1xtab.nozero, 1, FUN=function(x){x <- x[x>0]; length(x)}) > 1))$names) # IDs refer to t
		merging.regular <- as.numeric(colnames(tm1xtab.nozero)[as.numeric(apply(tm1xtab.nozero, 1, FUN=function(x){which.max(x)})[which(myfreq[-1,"value"] %in% merging)])]) # Returns the ID of the cluster in t-1 that has the greatest contribution to the merger
		this.fromID <- merging
		this.toID <- merging.regular
		#### If the cluster(s) already exist in the from/to list, remove them
		fromID <- c(fromID[-which(fromID %in% merging)], this.fromID)
		toID <- c(toID[-which(fromID %in% merging)], this.toID)
		results <- rbind(results[-which(fromID %in% merging),], data.frame("ID"=this.toID, "pixcount"=myfreq[which(myfreq[,"value"] %in% this.fromID), "count"], "timestep"=rep(mydate, length(this.fromID)),"class"=rep("splitting", length(this.fromID)) ))
		
		###################################################################
		### Regular tracking
		print("Regular tracking ...")
		mycols <- as.numeric(attributes(which(apply(tm1xtab.nozero, 2, FUN=function(x){x <- x[x>0]; length(x)}) == 1))$names)
		regulartrk <- as.numeric(attributes(which(apply(tm1xtab.nozero[,which(colnames(tm1xtab.nozero) %in% mycols)], 1, FUN=function(x){x <- x[x>0]; length(x)}) == 1))$names)
		#regulartrk <- unique(c(regulartrk, splitting.regular))
		this.fromID <- regulartrk[-which(regulartrk %in% fromID)]
		this.toID <- as.numeric(attributes(which(apply(tm1xtab.nozero[which(row.names(tm1xtab.nozero) %in% this.fromID),], 2, FUN=function(x){x <- x[x>0]; return(length(x))}) > 0))$names)
		fromID <- c(fromID, this.fromID)
		toID <- c(toID, this.toID)
		browser()
		results <- rbind(results, data.frame("ID"=this.toID, "pixcount"=myfreq[which(myfreq[,"value"] %in% this.fromID), "count"], "timestep"=rep(mydate, length(this.fromID)),"class"=rep("regularTracking", length(this.fromID)) ))
		###################################################################
		
browser()
	mypolys <- rasterToPolygons(rt.clump, dissolve=T, na.rm=T, fun=function(x){x>0})
	mypolys.tm1 <- rasterToPolygons(rtm1.clump, dissolve=T, na.rm=T, fun=function(x){x>0})
	e <- extent(mypolys.tm1[4,])
	plot(rtm1.clump, xlim=c(e@xmin, e@xmax), ylim=c(e@ymin, e@ymax))
	plot(mypolys, add=T)
	text(SpatialPoints(coordinates(mypolys)), 1:16, col="red")


		
		# Remove zeros
		# tm1.overlap <- tm1.overlap[which(attributes(tm1.overlap)$names != "0")]
		# tp1.overlap <- tp1.overlap[which(attributes(tp1.overlap)$names != "0")]
		
		## Reclassify clumps according to overlaps with t-1
		# Append both t-1 and t together
		t.all <- data.frame("tm1"=as.integer(tm1.overlap), "t"=as.numeric(attributes(tm1.overlap)$names)) 
		# Which IDs in t-1 have only 1 overlap in t?
		tm1.ids <- as.numeric(attributes(which(table(tm1.overlap) == 1))$names) 
		# Which IDs in t do these correspond to?
		t.ids <- tm1.overlap[tm1.overlap %in% tm1.ids]
		rectab <- data.frame("from"=as.numeric(attributes(t.ids)$names), "to"=as.integer(t.ids))
		# Add new IDs to rectab to include new initiations
		newids <- seq(max(uniqID)+1, max(uniqID)+length(t.all[!(t.all$t %in% rectab$from),1]), by=1)
		# Add old and new IDs together
		rectab <- rbind(rectab, data.frame("from"=t.all$t[!(t.all$t %in% rectab$from)],"to"=newids))
		
		# Convert clumps to polygons and get central point
		# rt.clump.poly <- rasterToPolygons(rt.clump, dissolve=T, na.rm=T, fun=function(x){x>0})
		# rt.clump.cntr < SpatialPoints(coordinates(rt.clump.poly))
		
		## Classify each potential case
		# 1) Regular tracking - MCS overlaps with a cluster both t-1 and t+1 
		# 	- overlaps in forward and backward images are unique
		# 	- recode overlapping clusters with the ID from t-1
		uniqID <- c(uniqID, newids)

		rt.clump.rc <- subs(rt.clump, rectab, filename=clus.f, overwrite=T)
		
		browser()
		
		results <- rbind(results, data.frame())
		
		# 2) Identify centre point of each cluster
		cntr <- SpatialPoints(coordinates(mypolys))
		
		# 3) Identify point of most intense precipitation in cluster
		myzones <- zonal(pr, rt.clump, stat="max")
		rt.clump.max <- subs(rt.clump, as.data.frame(myzones), by=1, which=2, subsWithNA=T)
		cntr.max <- rasterToPoints(overlay(pr, rt.clump.max, fun=function(x,y){x==y}), spatial=T, fun=function(x){x>0})
		
	}

	
}