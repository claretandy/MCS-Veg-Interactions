require(raster)
require(rasterVis)
# require(PP,lib.loc="/project/ukmo/rhel6/R")
require(rgdal)
require(sm)
require(plotrix)
require(RColorBrewer)
require(classInt)

tracker2 <- function(threshold=1000, inmcs, inprecip, id="s", timestep="10min", indatadir="/data/local/hadhy/ModelData/WAFR/", dlresultsdir="/data/local/hadhy/Projects/InternalSabbatical/Results/", overwrite=T){
    
    print("Running the MCS tracker ...")
    
    start=1 # Timestep to start at for testing
    
    # Set indir and outdir
    datadir <- paste(indatadir,"djzx",id,"/",sep="") 
    outdir <- paste(indatadir,"djzx",id,"/tracker/", threshold, "km_", timestep, "/", sep="")
    resultsdir <- paste(dlresultsdir,"Tracking/djzx",id,"_",timestep,"/", sep="")
    resultsfile <- paste(dlresultsdir,"Tracking/djzx",id,"_",timestep,".RData", sep="")
    
    if (!file.exists(outdir)){ dir.create(outdir, recursive=T) }
    if (!file.exists(resultsdir)){ dir.create(resultsdir, recursive=T) }
    
    if (overwrite == T | !file.exists(resultsfile)){
            
        mydates <- getZ(inprecip)
        mydate <- mydates[start]
        
        # Start new ID list
        uniqID <- vector("numeric")
        t1 <- raster(inmcs, layer=start)
        
        # Create (or read) clumps file for first time step. These IDs are used to recode clumps in the next time step.
        myfile <- paste(resultsdir,"IDs.", format(mydates[start], "%d%H.%M"), ".tif", sep="")
        if (overwrite==T){
            if (file.exists(myfile)){
                file.remove(myfile)
            }
            t1.clump <- clump(t1, directions=8, gaps=F, filename=myfile)
        } else {
            if (file.exists(myfile)){
                t1.clump <- raster(myfile)
            } else {
                t1.clump <- clump(t1, directions=8, gaps=F, filename=myfile)   
            }
        }
    
        # Create a results table to hold classification of all MCS for all dates
        myfreq <- freq(t1.clump, useNA="no")
        t1.poly <- rasterToPolygons(shift(t1.clump, x=-360), dissolve=T)
        mypts <- coordinates(t1.poly)
        results <- data.frame("ID"=myfreq[,1], "pixcount"=myfreq[,2], "timestep"=rep(mydates[start], nrow(myfreq)),"class"=rep("regularTracking", nrow(myfreq)), "x"=mypts[,1], "y"=mypts[,2] )
        results$class <- as.character(results$class) # Force it to be character rather than factor
        
        # Write outfile and set up for looping
        outfile <- paste(resultsdir,"mcs_tracking_",threshold,"km_",format(mydate, "%d.%H%M"),".tif",sep="")
        if (overwrite==T){
            tnew <- writeRaster(t1.clump, filename=outfile, overwrite=overwrite)
        } else {
            tnew <- t1.clump
        }
        prev.clust <- tnew
        prev.clust[is.na(prev.clust)] <- 0
        uniqID <- myfreq[,1]
        
        # Loop through all timesteps 
        # Compare current timestep to the previous timestep
        for (f in 2:nlayers(inmcs)){
        #for (f in (start+1):(start+20)){
            mydate <- mydates[f]
            mydate.tm1 <- mydates[f-1]
            print(mydate)
            # Set output filenames
            clus.f <- paste(outdir,"IDs.", format(mydates[f], "%d%H.%M"), ".tif", sep="")
            
            # Create clumps for t and read clumps for t-1
            rt <- raster(inmcs, layer=f) # Subset inmcs
            rt.clump <- clump(rt, directions=8, gaps=F) # Create clumps for current timestep
            rt.clump <- rt.clump * 1000 # Makes IDs that are very different from t-1, helping error identification
            rt.clump[is.na(rt.clump)] <- 0 # Set NA to 0 for crosstab
            rtm1.clump <- prev.clust # Get clusters from t-1
            rtm1.clump[is.na(rtm1.clump)] <- 0
            
    #         # Creat clumps for t+1
    #         rtp1 <- raster(inmcs, layer=f+1)
    #         rtp1.clump <- clump(rtp1, directions=8, gaps=F)
    #         rtp1.clump[is.na(rtp1.clump)] <- 0
            
            # Read inprecip for current time step
    #         pr <- raster(inprecip, layer=f)
            # Crosstab backwards
            myxtab <- mycrosstab(rt.clump, rtm1.clump)
    #         myxtab.p1 <- mycrosstab(rt.clump, rtp1.clump)
            tnew <- rt.clump
            tnew[tnew == 0] <- NA
            
            # Run series of tracing algorithms
            # Inputs: tm1, t, tnew, myxtab, results
            # Outputs: return(list(tnew, tremain, results, uniqID))
            
            trk <- regular_tracking(rtm1.clump, rt.clump, tnew, myxtab, results, uniqID, mydate)
    #         print(unique(trk[[1]]))
            trk <- dissipation(rtm1.clump, trk[[2]], trk[[1]], myxtab, trk[[3]], trk[[4]], mydate, mydate.tm1) 
    #         print(unique(trk[[1]]))
            trk <- generation(rtm1.clump, trk[[2]], trk[[1]], myxtab, trk[[3]], trk[[4]], mydate) 
    #         print(unique(trk[[1]]))
            trk <- merging(rtm1.clump, trk[[2]], trk[[1]], myxtab, trk[[3]], trk[[4]], mydate, mydate.tm1) 
    #         print(unique(trk[[1]]))
            trk <- splitting(rtm1.clump, trk[[2]], trk[[1]], myxtab, trk[[3]], trk[[4]], mydate)
            trk <- remaining(rtm1.clump, trk[[2]], trk[[1]], myxtab, trk[[3]], trk[[4]], mydate)
    
            # Write outfile
            outfile <- paste(resultsdir,"mcs_tracking_",threshold,"km_",format(mydate, "%d.%H%M"),".tif",sep="")
            tnew <- writeRaster(trk[[1]], filename=outfile, overwrite=overwrite)
            # End: set clumps for t as prev.clust
            prev.clust <- tnew # tnew
            results <- trk[[3]]
            uniqID <- trk[[4]]
            
            resultsfile <- paste(dlresultsdir,"Tracking/djzx",id,"_",timestep,"_temp.RData", sep="")
            if(file.exists(resultsfile)){ file.remove(resultsfile) }
            save(results, file=resultsfile)
        }
        resultsfile <- paste(dlresultsdir,"Tracking/djzx",id,"_",timestep,".RData", sep="")
        if(file.exists(resultsfile)){ file.remove(resultsfile) }
        file.remove(paste(dlresultsdir,"Tracking/djzx",id,"_",timestep,"_temp.RData", sep=""))
        save(results, file=resultsfile)
    
    } else {
        load(resultsfile)
    }
    
    return(results)
}

my10pc.test <- function(x){
    mypc <- 100*(x / sum(x)); 
    pos <- which(mypc > 1)
    gt10pc <- as.numeric(attributes(x[pos])$names)
    gt10pc <- gt10pc[gt10pc > 0]
    if(length(gt10pc)==1){
        return(gt10pc)
    } else {
        return(NA)
    }
}

regular_tracking <- function(tm1, t, tnew, myxtab, results, uniqID, mydate){
    ### Regular tracking: MCS in t that ARE in t-1
    ### Classify cases that have a 1:1 relationship
    ### Uses IDs from t-1
    
    print("Regular tracking ...")
    
    # First, remove overlaps with zero
#     myxtab.nz <- myxtab[rownames(myxtab) != 0, colnames(myxtab) != 0]
    
    # How many clusters in t-1 overlap with clusters in t?
    # Clusters in t-1 (t) may overlap clusters in t (t-1) by as much as 10% before they are called splits or merges
#     no.ovrcl.tm1.t <- apply(myxtab.nz, 2, FUN=function(x){z <- which(x>0); return(length(z))})
    
    keep.t <- apply(myxtab[which(rownames(myxtab)>0),which(colnames(myxtab)>0)], 1, FUN=function(x){x <- x[x>0]; return(length(x))})
    keep.t <- as.numeric(attributes(keep.t[keep.t == 1])$names)
    
    keep.tm1 <- apply(myxtab[which(rownames(myxtab) %in% keep.t),which(colnames(myxtab)>0)], 2, FUN=function(x){x<-x[as.numeric(attributes(x)$names)>0]; x <- x[x>0]; return(length(x))})
    keep.tm1 <- as.numeric(attributes(keep.tm1[keep.tm1 == 1])$names)

    myxtab.regtrk <- myxtab[which(as.numeric(dimnames(myxtab)[[1]]) %in% keep.t), which(as.numeric(dimnames(myxtab)[[2]]) %in% keep.tm1)]
    myxtab.regtrk <- myxtab.regtrk[which(rowSums(myxtab.regtrk) != 0),]
    
#     if (mydate == as.POSIXct("2006-08-16 13:10:00 BST")){browser()}
#     
#     # Which clusters in t-1 have 1 overlap with clusters in t? 
#     no.ovrcl.tm1.t <- apply(myxtab.regtrk, 2, my10pc.test)
#     # Which clusters in t have 1 overlap with clusters in t-1?
#     no.ovrcl.t.tm1 <- apply(myxtab.regtrk, 1, my10pc.test)
#     
#     keepall <- no.ovrcl.tm1.t[!(no.ovrcl.tm1.t %in% as.numeric(attributes(no.ovrcl.t.tm1[is.na(no.ovrcl.t.tm1)])$names)) & !is.na(no.ovrcl.tm1.t)]
    
#     # Which clusters are they in t-1?
#     tm1keep <- as.numeric(attributes(which(!is.na(no.ovrcl.tm1.t)))$names)
#     # Which clusters are they in t-1?
#     tkeep <- as.numeric(attributes(which(!is.na(no.ovrcl.t.tm1)))$names)
#     
#     # Do another apply to remove cases of merges or splits
#     browser()
#     myxtab.mands <- myxtab.nz[which(rownames(myxtab.nz) %in% tkeep), which(colnames(myxtab.nz) %in% tm1keep)]
#     myxtab.tm1.nz <- apply(myxtab.mands, 2, FUN=function(x){z <- which(x>0); return(length(z))})
#     myxtab.t.nz <- apply(myxtab.mands, 1, FUN=function(x){z <- which(x>0); return(length(z))})
#     tm1keep <- as.numeric(attributes(which(myxtab.tm1.nz == 1))$names)
#     tkeep <- as.numeric(attributes(which(myxtab.t.nz == 1))$names)
#     
#     # Subset myxtab including only clusters that have a 1 to 1 match
#     myxtab.1to1 <- myxtab.nz[which(rownames(myxtab.nz) %in% tkeep), which(colnames(myxtab.nz) %in% tm1keep)]
#     
    # Set up recoding data.frame
    myrecodes <- apply(myxtab.regtrk, 2, FUN=function(x){cn <- colnames(myxtab.regtrk); ii <- which(x>0); as.numeric(attributes(ii)$names)} )
    this.fromID <- as.numeric(myrecodes)
    this.toID <- as.numeric(attributes(myrecodes)$names)
    
    # Make sure they are both the same length
    if (length(this.fromID) != length(this.toID)){
        # This means there are either splits or merges
        print("this.fromID != this.toID")
        browser()
    }
    
    # Set up recoding data.frame
    rctab <- data.frame("from"=this.fromID, "to"=this.toID)
    t.rc <- subs(t, rctab, by="from", which="to", subWithNA=T) # Substitutes from with to, and sets all other cells to NA. 
    
    # Add new clumps into tnew raster
    tnew <- cover(tnew, t.rc)
    
    # Which clumps still need to be classified from t?
    tremain <- overlay(t,tnew, fun=function(x,y){z <- x > 0 & is.na(y); return(z*x) } )
    
    # A few more fields for the results table
    myfreq <- freq(t, useNA="no")
    trc.poly <- rasterToPolygons(shift(t.rc, x=-360), dissolve=T)
    mypts <- coordinates(trc.poly)
    
    # Append to results data.frame
    pixcount <- as.numeric(myfreq[which(myfreq[,"value"] %in% this.fromID), "count"])
    newresults <- data.frame("ID"=this.toID, "pixcount"=pixcount, "timestep"=rep(mydate, length(this.fromID)),"class"=rep("regularTracking", length(this.fromID)), "x"=mypts[,1], "y"=mypts[,2] )
    results <- rbind(results, newresults)
    
    # NB: uniqID stays the same because all the new IDs in tnew are in t-1
    return(list(tnew, tremain, results, uniqID))
}

generation <- function(tm1, t, tnew, myxtab, results, uniqID, mydate){
    ### Generation: MCS in t that are not in t-1
    ### Creates new IDs
    print("Generation ...")
    
    # Which IDs in t have no overlapping cluster at t-1?
    mygens <- as.numeric(attributes(which(myxtab[,which(colnames(myxtab)=="0")] == apply(myxtab, 1, sum)))$names)
#     browser()

    # If some generations exist, give them new IDs following on from uniqID
    if (length(mygens)>0){
#         print(paste(length(mygens), "generations"))
        # Set up recoding data.frame
        this.fromID <- mygens
        this.toID <- seq(max(uniqID)+1, max(uniqID)+length(mygens))
        rctab <- data.frame("from"=this.fromID, "to"=this.toID)
        t.rc <- subs(t, rctab, by="from", which="to", subWithNA=T)
        
        # Add new clumps into tnew raster
        tnew <- cover(tnew, t.rc)
        
        # Which clumps still need to be classified from t?
        tremain <- overlay(t,tnew, fun=function(x,y){z <- x > 0 & is.na(y); return(z*x) } )
#         print("check that these generations are being recorded and removed from tremain")
#         browser()
        # A few more fields for the results table
        myfreq <- freq(t, useNA="no")
        trc.poly <- rasterToPolygons(shift(t.rc, x=-360), dissolve=T)
        mypts <- coordinates(trc.poly)
        
        # Add new clumps into results data.frame
        uniqID <- c(uniqID, this.toID)
        pixcount <- as.numeric(myfreq[which(myfreq[,"value"] %in% this.fromID),"count"])
        newresults <- data.frame("ID"=this.toID, "pixcount"=pixcount, "timestep"=rep(mydate, length(this.fromID)),"class"=rep("generation", length(this.fromID)), "x"=mypts[,1], "y"=mypts[,2] )
        results <- rbind(results, newresults)
        
    } else {
        tremain <- t
    }
    return(list(tnew, tremain, results, uniqID))
}

dissipation <- function(tm1, t, tnew, myxtab, results, uniqID, mydate, mydate.tm1){
    ### Dissipation: MCS in t that are not in t+1
    ### Looks backwards (t-1), and recodes "regular tracking" MCS that do NOT exist in t 
    print("Dissipation ...")
    
    # Crosstab tnew (i.e. MCS just classified as regular tacking) with tp1
#     myxtab.new <- crosstab(tm1, t)
#     browser()
    
    # Which IDs in t have no overlapping cluster at t+1?
    ## First, what's the total number of pixels in each cluster at t
    t.colsum <- apply(myxtab, 2, sum)
    ## Second, how many of these cells are 0 in t+1?
    tm1.row0 <- myxtab[which(rownames(myxtab) == 0),]
    
    # If First == Second, then none of the cells in t overlap with cells in t+1, meaning that the MCS has dissipated. Following tests for each MCS in t if first == second
    z <- which(t.colsum == tm1.row0)
    if (length(z) > 0){
#         print("Some MCS have dissipated, so now I have to write some more code ...")
#         browser()   
        # Check this code works when it occurs for the first time
        # 1. Pullout correct ID from t.rowsum
        diss.ids <- as.numeric(attributes(t.colsum[z])$names)
        
        # 2. Recode entry in results table from "regular tracking" to "dissipation"
        oldclass <- results[which(results$timestep == mydate.tm1 & results$ID %in% diss.ids),"class"]
        results[which(results$timestep == mydate.tm1 & results$ID %in% diss.ids),"class"] <- paste(oldclass,"/ dissipation")
    } 
    # tremain stays the same because we are not classifying any new clusters
    tremain <- t
    return(list(tnew, tremain, results, uniqID))
}

notmax <- function(x){
    mypc <- 100*(x / sum(x))
    pos <- which(mypc > 0 & x != max(x)) # Removed 10% threshold pos <- which(mypc > 10 & x != max(x))
    notmaxtm1 <- as.numeric(attributes(pos)$names)
    return(notmaxtm1[notmaxtm1 > 0])
}

merging <- function(tm1, t, tnew, myxtab, results, uniqID, mydate, mydate.tm1){
    ### Merging: 2 or more MCS in t-1 becomes 1 in t. 
    #### Largest MCS in t-1 gives ID to MCS in t
    print("Merging ...")
    tm1[is.na(tm1)] <- 0
    myxtab.remain <- mycrosstab(t, tm1)
    # How many MCS in t-1 overlap with remaining clusters in t?
    no.overlaps.gt10perc <- apply(myxtab.remain, 1, FUN=function(x){y<-100*(x/sum(x)); z <- which(y>0 & as.numeric(attributes(x)$names))!=0; length(z)}) # Removed 10% threshold z <- which(y>10 & as.numeric(attributes(x)$names))!=0;
    # Exclude zero in t
    no.overlaps.gt10perc.nz <- no.overlaps.gt10perc[attributes(no.overlaps.gt10perc)$names > 0]
    # How many clusters in t are merges?
    merging <- as.numeric(attributes(which(no.overlaps.gt10perc.nz > 1))$names)
    
    if (length(merging)==0){
        
        print("No More clusters to classify")
        tremain <- t
    } else {
        remain.merges <- myxtab.remain[rownames(myxtab.remain) %in% merging,]
#         if (mydate == as.POSIXct("2006-08-16 17:30:00 BST")){browser()}
        
        if (length(merging)==1){
            # merging = cluster(s) in t that need to be recoded to largest overlap in t-1
            # maxtm1 = value they will be recoded to
            maxtm1 <- as.numeric(attributes(which.max(remain.merges[attributes(remain.merges)$names > 0]))$names)
            # Which clusters in t-1 need to be recoded in results to "merges next timestep"?
            notmaxtm1 <- notmax(remain.merges) 
            this.fromID <- merging
            this.toID <- maxtm1
        } else {
            # What is the main overlapping cluster from t-1? Use this for recoding.
            maxtm1 <- apply(remain.merges[,which(colnames(remain.merges) != 0)], 1, FUN=function(x){as.numeric(attributes(x[which.max(x)])$names)})
            # What are the minor overlapping clusters from t-1? Use this for recoding results in t-1
            # Minor overlaps must be > 10% of the cluster in t
            notmaxtm1 <- apply(remain.merges, 1, notmax)
            
            this.fromID <- as.numeric(attributes(maxtm1)$names)
            this.toID <- as.numeric(maxtm1)
        }
#         if (mydate == as.POSIXct("2006-08-16 19:40:00 BST")){browser()}
        rctab <- data.frame("from"=this.fromID, "to"=this.toID)
        t.rc <- subs(t, rctab, by="from", which="to", subWithNA=T)
        
        # Add new clumps into tnew raster
        tnew <- cover(tnew, t.rc)
        
        # Which clumps still need to be classified from t?
        tremain <- overlay(t,tnew, fun=function(x,y){z <- x > 0 & is.na(y); return(z*x) } ) # Returns 1*t if t>0 and tnew is NA
        
        # A few more fields for the results table
        myfreq <- freq(t, useNA="no")
        trc.poly <- rasterToPolygons(shift(t.rc, x=-360), dissolve=T)
        mypts <- coordinates(trc.poly)
        
        # Add new clumps into results data.frame
#         browser()
        pixcount <- as.numeric(myfreq[which(myfreq[,"value"] %in% this.fromID),"count"])
        newresults <- data.frame("ID"=this.toID, "pixcount"=pixcount, "timestep"=rep(mydate, length(this.fromID)),"class"=rep("merging", length(this.fromID)), "x"=mypts[rank(this.toID),1], "y"=mypts[rank(this.toID),2] )
        results <- rbind(results, newresults)
        
        # Edit results for IDs in t-1 that contributed to the merge, but lost their ID in t
        if (class(notmaxtm1) == "numeric" | class(notmaxtm1) == "matrix"){
            tm1.merged <- unique(as.numeric(notmaxtm1))
        } else {
            tm1.merged <- unique(as.numeric(unlist(notmaxtm1)))
        }
#         if (mydate == as.POSIXct("2006-08-17 08:20:00 BST")){ browser() }
        results[which(results$timestep == mydate.tm1 & results$ID %in% tm1.merged),"class"] <- rep("merges in next timestep", length(tm1.merged))
        

    }
#     if (mydate == as.POSIXct("2006-08-16 19:40:00 BST")){browser()}
    return(list(tnew, tremain, results, uniqID))
}

splitting <- function(tm1, t, tnew, myxtab, results, uniqID, mydate){
    
    print("Splitting ...")
    # Cross-tabulate remaining clusters in t with t-1
    myxtab.remain <- mycrosstab(t, tm1)
    # How many MCS in t-1 overlap with remaining clusters in t?
    no.overlaps.gt10perc <- apply(myxtab.remain, 2, FUN=function(x){y<-100*(x/sum(x)); z <- which(y>0 & as.numeric(attributes(x)$names))!=0; length(z)}) # Removed 10% threshold z <- which(y>10 & as.numeric(attributes(x)$names))!=0;
    # Exclude zero
    no.overlaps.gt10perc.nz <- no.overlaps.gt10perc[attributes(no.overlaps.gt10perc)$names > 0]
    # Which clusters in t-1 are splits in t?
    splits <- as.numeric(attributes(which(no.overlaps.gt10perc.nz > 1))$names)
    
    if (length(splits)==0){
        print("No More clusters to classify")
        tremain <- t
    } else {
        print("Found a split!")
        splitting.remain <- myxtab.remain[which(rownames(myxtab.remain)>0),which(colnames(myxtab.remain) %in% splits)]
        if (mydate == as.POSIXct("2006-08-16 11:30:00 BST") & id == "w"){ browser() }
        if (length(splits)==1){
            # Check that the cluster in t-1 does not split evenly between 2 or more clusters in t
            if (length(unique(splitting.remain)) != length(splitting.remain)){
                # In this case, maxt is the cluster with the most pixels (ie including 0)
                maxt <- as.numeric(
                    attributes(
                        which.max(
                            apply(myxtab.remain[rownames(myxtab.remain)!=0,], 1, sum)))$names)
                allt <- as.numeric(attributes(myxtab.remain[rownames(myxtab.remain)!=0,colnames(myxtab.remain)==splits])$names)
                notmaxt <- allt[!(allt %in% maxt)]
            } else {
                maxt <- as.numeric(attributes(which.max(splitting.remain[attributes(splitting.remain)$names > 0]))$names)
                # Which clusters in t-1 need to be recoded in results to "minor split"?
                notmaxt <- notmax(splitting.remain) 
            }
            this.fromID <- maxt
            this.toID <- splits

        } else {
            # If there are > 1 splits ...
            # What is the main overlapping cluster from t-1? Use this for recoding.
            maxt <- apply(splitting.remain, 2, FUN=function(x){y <- x[sum(x)>0]; as.numeric(attributes(which.max(y))$names)})
            # Which clusters in t-1 need to be recoded in results to "minor split"?
            notmaxt <- apply(splitting.remain, 2, notmax)
            
            this.fromID <- as.numeric(maxt)
            this.toID <- as.numeric(attributes(maxt)$names)
        }
        ##################################
        # Remove zeros
        maxt <- maxt[maxt > 0]
        
        # Recode the larger parts of the split
        rctab <- data.frame("from"=this.fromID, "to"=this.toID)
        t.rc <- subs(t, rctab, by="from", which="to", subWithNA=T)
        
        # Add new clumps into tnew raster
        tnew <- cover(tnew, t.rc)
        
        # Which clumps still need to be classified from t?
        tremain <- overlay(t,tnew, fun=function(x,y){z <- x > 0 & is.na(y); return(z*x) } ) # Returns 1*t if t>0 and tnew is NA
#         if (mydate == as.POSIXct("2006-08-16 17:10:00 BST")){browser()}
        # A few more fields for the results table
        myfreq <- freq(t, useNA="no")
        trc.poly <- rasterToPolygons(shift(t.rc, x=-360), dissolve=T)
        mypts <- coordinates(trc.poly)
        
        # Add new clumps into results data.frame
        pixcount <- as.numeric(myfreq[which(myfreq[,"value"] %in% this.fromID),"count"])
        newresults <- data.frame("ID"=this.toID, "pixcount"=pixcount, "timestep"=rep(mydate, length(this.fromID)),"class"=rep("splitting - largest part", length(this.fromID)), "x"=mypts[,1], "y"=mypts[,2] )
        results <- rbind(results, newresults)
#         browser()
        # Now treat the smaller splits as generations
        # What are the minor overlapping clusters from t-1? Use this for recoding results in t-1
        split.remain <- myxtab.remain[,colnames(myxtab.remain) %in% splits] # rownames(myxtab.remain)!=0
        if (!exists("notmaxt")){
            if (class(split.remain) == "integer"){
                notmaxt <- notmax(split.remain)
            } else {
                notmaxt <- apply(split.remain, 2, notmax)
            }            
        }
         
        # Result of the above (notmaxt) could be either:
        #   - numeric (if there's only 1 split, and 1 minor overlap)
        #   - a matrix (if the number of minor overlaps is the same for all splits), or 
        #   - a list (if the number of minor overlaps is different)
        # Remove zeros
        if (class(notmaxt) == "numeric"){
            this.fromID <- notmaxt
            this.toID <- seq(max(uniqID)+1, max(uniqID)+length(this.fromID))
            uniqID <- c(uniqID, this.toID)
            # Recode the minor parts of the split
            rctab <- data.frame("from"=this.fromID, "to"=this.toID)
            t.rc <- subs(tremain, rctab, by="from", which="to", subWithNA=T)
            
            # Add new clumps into tnew raster
            tnew <- cover(tnew, t.rc)
            
            # Which clumps still need to be classified from t?
            tremain <- overlay(tremain,tnew, fun=function(x,y){z <- x > 0 & is.na(y); return(z*x) } ) # If tremain>0 and tnew is NA Returns 1*tremain
            
            # A few more fields for the results table
            myfreq <- freq(t, useNA="no")
            trc.poly <- rasterToPolygons(shift(t.rc, x=-360), dissolve=T)
            mypts <- coordinates(trc.poly)
            
            # Add new clumps into results data.frame
            pixcount <- as.numeric(myfreq[which(myfreq[,"value"] %in% this.fromID),"count"])
            newresults <- data.frame("ID"=this.toID, "pixcount"=pixcount, "timestep"=rep(mydate, length(this.fromID)),"class"=rep("splitting - small parts", length(this.fromID)), "x"=mypts[,1], "y"=mypts[,2] )
            results <- rbind(results, newresults)
        } else {
            
            if (class(notmaxt) == "list"){
                notmaxt <- lapply(notmaxt, FUN=function(x){x[which(x!=0)]})
                notmaxt.length <- length(notmaxt)
            }
            
            if (class(notmaxt) == "matrix"){
                notmaxt.length <- dim(notmaxt)[2]
            }
            
            # Recode smaller parts of the split as new IDs
            for (i in 1:notmaxt.length){
                if (class(notmaxt) == "list"){
                    this.fromID <- notmaxt[[i]]
                }
                if (class(notmaxt) == "matrix"){
                    this.fromID <- notmaxt[notmaxt[,i] > 0,i]
                }
                if (class(notmaxt) == "numeric"){
                    this.fromID <- notmaxt[which(notmaxt[i] > 0)]
                }
                
                this.toID <- seq(max(uniqID)+1, max(uniqID)+length(this.fromID))
                uniqID <- c(uniqID, this.toID)
                
                #                 browser()
                # Recode the larger parts of the split
                rctab <- data.frame("from"=this.fromID, "to"=this.toID)
                t.rc <- subs(tremain, rctab, by="from", which="to", subWithNA=T)
                
                # Add new clumps into tnew raster
                tnew <- cover(tnew, t.rc)
                
                # Which clumps still need to be classified from t?
                tremain <- overlay(tremain,tnew, fun=function(x,y){z <- x > 0 & is.na(y); return(z*x) } ) # If tremain>0 and tnew is NA Returns 1*tremain
#                 browser()
                # A few more fields for the results table
                myfreq <- freq(t, useNA="no")
                trc.poly <- rasterToPolygons(shift(t.rc, x=-360), dissolve=T)
                mypts <- coordinates(trc.poly)
                
                # Add new clumps into results data.frame
                pixcount <- as.numeric(myfreq[which(myfreq[,"value"] %in% this.fromID),"count"])
                newresults <- data.frame("ID"=this.toID, "pixcount"=pixcount, "timestep"=rep(mydate, length(this.fromID)),"class"=rep("splitting - small parts", length(this.fromID)), "x"=mypts[,1], "y"=mypts[,2] )
                results <- rbind(results, newresults)
            }
        }
        
        
    }
        
    return(list(tnew, tremain, results, uniqID))
}

remaining <- function(tm1, t, tnew, myxtab, results, uniqID, mydate){
    
    if (cellStats(t, "max") == 0){
        tremain <- t
    } else {
        print("More cases to classify!")
        
        # Cases that are flagged up here :
        #   - a split that also contains a merge
        
        myxtab.remain <- mycrosstab(t, tm1)
        no.ovrcl.tm1.t <- apply(myxtab.remain, 2, my10pc.test)
        # Which clusters in t have 1 overlap with clusters in t-1?
        no.ovrcl.t.tm1 <- apply(myxtab.remain, 1, my10pc.test)
        
        keepall <- no.ovrcl.tm1.t[!(no.ovrcl.tm1.t %in% as.numeric(attributes(no.ovrcl.t.tm1[is.na(no.ovrcl.t.tm1)])$names)) & !is.na(no.ovrcl.tm1.t)]
        
        # Set up recoding data.frame
        this.fromID <- as.numeric(keepall)
        this.toID <- as.numeric(attributes(keepall)$names)
        rctab <- data.frame("from"=this.fromID, "to"=this.toID)
        t.rc <- subs(t, rctab, by="from", which="to", subWithNA=T) # Substitutes from with to, and sets all other cells to NA. 
        
        # Add new clumps into tnew raster
        tnew <- cover(tnew, t.rc)
        
        # Which clumps still need to be classified from t?
        tremain <- overlay(t,tnew, fun=function(x,y){z <- x > 0 & is.na(y); return(z*x) } )
        
        # A few more fields for the results table
        myfreq <- freq(t, useNA="no")
        trc.poly <- rasterToPolygons(shift(t.rc, x=-360), dissolve=T)
        mypts <- coordinates(trc.poly)
        
        # Append to results data.frame
        pixcount <- as.numeric(myfreq[which(myfreq[,"value"] %in% this.fromID), "count"])
        newresults <- data.frame("ID"=this.toID, "pixcount"=pixcount, "timestep"=rep(mydate, length(this.fromID)),"class"=rep("regularTracking", length(this.fromID)), "x"=mypts[,1], "y"=mypts[,2] )
        results <- rbind(results, newresults)
        
    }
    if (cellStats(tremain, "max") > 0){ browser() }
    
    return(list(tnew, tremain, results, uniqID))
}

test.tracking <- function(inmcs, results){
    r <- raster(paste(dlresultsdir,"Tracking/djzxw_10min/mcs_tracking_1000km_16.2020.tif", sep=""))
    load(paste(dlresultsdir,"Tracking/djzxw_10min.Rdata",sep=""))
    print(plot(shift(r, x=-360))); print(plot(SpatialLines( list( Lines(list(Line(results[results$ID==1,c("x","y")])),ID=1) ) ), add=T))
    for (ii in 2:150){ 
        print(plot(SpatialLines( list( Lines(list(Line(results[results$ID==ii,c("x","y")])),ID=ii) ) ), add=T)); print(ii); browser() 
    }
}

mycrosstab <- function(r1, r2){
    
    myxtab <- crosstab(r1, r2)
    myxtab2 <- reshape(myxtab, v.names="Freq", idvar="Var1", timevar="Var2", direction="wide")
    colnames(myxtab2) <- sub("Freq.","",colnames(myxtab2))
    myxtab2 <- myxtab2[,colnames(myxtab2)!="NA"]
    myxtab2 <- myxtab2[!is.na(myxtab2$Var1),] # Drop the row with NA
    rownames(myxtab2) <- myxtab2$Var1
    myxtab2 <- myxtab2[,-1]

    return(myxtab2)
}
