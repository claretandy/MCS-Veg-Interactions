# Check MCS line tracks match MCS blocks by plotting each time step and each line segment

trackCheck <- function(tracks, mcs, inbr, id="s", threshold="1000", timestep="10min"){
    
    require(animation)
    require(colorRamps)
    
    # Make colourful legend for MCS classes ...
    mycols<- rep(rev(matlab.like(50)), length.out=length(unique(tracks$ID)))
#     mycols<- rep(rev(rainbow(50)), length.out=length(unique(tracks$ID)))
        
    # Loop through dates
    mydates <- getZ(inbr)
    
    # Save to animated gif
#     saveGIF({
    
    for (d in 100:110){ #1:length(mydates)){
        
        print(mydates[d])
        # Load MCS image for current date
        r <- raster(paste(dlresultsdir,"Tracking/djzx",id,"_",timestep,"/mcs_tracking_",threshold,"km_",format(mydates[d], "%d"),".",format(mydates[d], "%H%M"),".tif", sep=""))
        
        # Load tracks up to (and including) the current date
        ## subset tracks
        trackss <- tracks[tracks$timestep<= mydates[d],c("x","y","ID")]
        trackLines <- makeLines(trackss)
        curLines <- trackLines[trackLines$MCSid %in% unique(mcs.now),]
        oldLines <- trackLines[!trackLines$MCSid %in% unique(mcs.now),]

        # Plot image and lines
        par(mfrow=c(2,1), mar=c(1,1,1,1))
        rain <- calc(shift(inbr[[d]], x=-360), function(x){x*3600})
        rain[rain<1] <- NA
        
        # Convert mcs.now to factor
        mcs.now <- shift(r, x=-360)
        mcs.nowf <- as.factor(mcs.now)[[1]]
        flev <- levels(f)[[1]]
        flev[["MCSid"]] <- flev$ID
        levels(mcs.nowf) <- flev

        p1 <- levelplot(rain, at=c(0,10,20,30,40,60,100,200), margin=F, ylim=c(-9, 10),
                        col.regions=rev(terrain.colors(10)), auto.key=FALSE, 
                        scales=list(draw=FALSE), ylab=NULL, xlab=NULL) + 
            layer(sp.polygons(land_simple, col='grey'))
        
        p2 <- levelplot(mcs.nowf, margin=F, ylim=c(-9, 10), 
                        main=format(mydates[d], "%d %b %H:%M"),
                        col.regions=mycols[unique(mcs.now)], auto.key=FALSE, 
                        scales=list(draw=FALSE), ylab=NULL, xlab=NULL) + 
            layer(sp.polygons(land_simple, col='grey')) +
            layer(sp.lines(curLines, col="black")) +
            layer(sp.lines(oldLines, col="#E6E6E6"))
        
        print(p1, split=c(1, 1, 1, 2), more=TRUE)
        print(p2, split=c(1, 2, 1, 2))

    }
    
#     }, movie.name="mcs_tracks_movie.gif", interval=0.5, ani.width=1000, ani.height=600) # End of saveGIF
    
}

makeLines <- function(tracks){
    
    iii <- unique(tracks$ID)
    li <- 1
    mylines <- list()
    for (ii in iii){
        if(!exists("mytext")){
            mytext <- data.frame(tracks[tracks$ID == ii,c("x","y","ID")][1,])
        } else {
            mytext <- rbind(mytext, tracks[tracks$ID == ii,c("x","y","ID")][1,])
        }
        myline <- Line(tracks[tracks$ID == ii,c("x","y")])
        mylines[[li]] <- Lines(list(myline), ID=li) # or ID=li
#         assign(x=paste("myline",li,sep=""), value=mylines)
        li <- li+1
    }
    myspl <- SpatialLines(mylines)
    myspdf <- SpatialLinesDataFrame(sl=myspl, data=data.frame(lineid=1:length(iii), MCSid=iii))
    
    return(myspdf)
}

