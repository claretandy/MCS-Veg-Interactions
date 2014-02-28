loadAllAncils <- function(myproj="rp", nBndClass=1, model="rb5216.4km.std", overwrite=F){
    
    require(raster)
    require(rgdal)
    source("functions/loadVeg.R")
    source("functions/loadOtherAncils.R")
    source("functions/makeBoxes.R")
    source("functions/vegPrep.R")
    source("functions/getLUT.R")
    
    land_simple <- readOGR(dsn=paste(indatadir,"ancils",sep=""), layer=paste("land_",myproj,sep=""))
    
    myveg <- loadVeg(model.nm=model, proj=myproj, overwrite=overwrite)
    mylandfrac <- loadOtherAncils(model.nm=model, ancil="landfrac", proj=myproj, overwrite=overwrite)
    myorog <- loadOtherAncils(model.nm=model, ancil="orog", proj=myproj, overwrite=overwrite)  
    
    if (myproj == "rp"){
        myveg <- shift(adjCoords(myveg), x=-360)
        mylandfrac <- shift(adjCoords(mylandfrac), x=-360)
        myorog <- shift(adjCoords(myorog), x=-360)
        # Prepare boxes (rp) ...
        myboxes <- makeBoxes2(xmin=-6, ymin=-8, xmax=14, ymax=4, xdim=4, ydim=4, temp=myorog)
    } else {
        # Prepare boxes 
        myboxes <- makeBoxes2(xmin=-11, ymin=5, xmax=9, ymax=17, xdim=4, ydim=4, temp=myorog)
    }
    
    spp.r <- myboxes[[1]] # non-adjusted raster
    spp <- myboxes[[3]] # non-adjusted polygons
        
    # Prepare vegetation data for use in clipping etc
    allveg <- vegPrep(model.nm="rb5216.4km.std", id="s", myproj=myproj, myveg, myorog, mylandfrac, land_simple, spp, spp.r, plots=F, vegThreshold=0.3, nBndClass=1, overwrite=overwrite) # return(mycl, mycl.z) and creates pdf plots
    mycl <- allveg[[1]]
    mycl.z <- allveg[[2]]
    
    if (myproj == "rp" & raster::xmax(mycl) > 360){
        mycl <- shift(mycl, x=-360)
        mycl.z <- shift(mycl.z, x=-360)
    }
    
    myLUT <- getLUT(nBndClass=nBndClass)
        
    mycl.f <- as.factor(mycl)
    fdat <- myLUT[myLUT$ID %in% levels(mycl.f)[[1]]$ID, ]
    levels(mycl.f) <- fdat 
    mycl.f <- as.factor(mycl)
    fdat <- myLUT[myLUT$ID %in% levels(mycl.f)[[1]]$ID, ]
    levels(mycl.f) <- fdat
    
    return(list(myveg, mylandfrac, myorog, mycl, mycl.f, spp, spp.r, myboxes, land_simple, mycl.z))
}