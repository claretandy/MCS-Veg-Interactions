
library(animation)

mydates <- getZ(rb5216.4km.std)

mcs_tracks_all <- shift(brick(paste(mcsrst.path,"mcs_tracking.vrt",sep="")), x=-360)
precip.all <- shift(stack(rb5216.4km.std), x=-360)

# Mask only MCS precip
precip.mcsall <- mask(precip.all, mcs_tracks_all)*3600
# precip.mcsall.int <- calc(precip.mcsall, fun=function(x){x[x<20 | is.na(x)]<-0; return(x)})

iii <- c(110, 136, 150, 174, 328, 334, 345, 347, 469, 542, 586, 593, 611, 782, 812, 818, 863, 884, 919, 934, 949, 987) # This is a subset of the above > 20,000 pixel threshold

for (ii in iii){
    print(ii)
    
    # Get the correct time slices
    mcs.dates <- unique(results[results$ID == ii,"timestep"])    
    ss <- which(mydates %in% mcs.dates)

#     this.mcs <- subset(mcs_tracks_all, which(getZ(inbr) %in% mydates))
#     this.precip <- subset(precip.all, which(getZ(inbr) %in% mydates))
    
    # Mask out other MCS and precip < 20mm per hour 
    this.mcs <- calc(mcs_tracks_all[[ss]], fun=function(x){x[x!=ii]<-NA; return(x/ii)}) #  y <- ifelse(x==ii,1,NA); return(y)}
#     precip.int <- calc(this.precip, fun=function(x){x[x<(20/3600)]<-NA; return(x)})
#     mcs.inttrack <- mask(this.mcs, precip.int) precip.mcsall.int
#     mcs.inttrack <- mask(this.mcs, precip.mcsall.int[[ss]]) 
    
    # Combine all MCS patches together to get the extent of the MCS
    mcs.int.sum <- calc(this.mcs, fun=function(x){x[is.na(x)]<-0; return(sum(x))})
    mcs.int.sum <- trim(mcs.int.sum, value=0, padding=1)
    e <- extent(mcs.int.sum)
    
    # Calculate all combined intense precip for all timesteps of the MCS
    this.precip.mcs <- mask(crop(precip.mcsall[[ss]], e), crop(this.mcs, e))
    precip.int.sum <- calc(this.precip.mcs, fun=function(x){
        x[x<20 | is.na(x)] <- 0;
        return(sum(x))})
    
#     browser()
    # Plot movie zoomed in on the MCS
#     saveGIF({ 
    saveHTML({
        ani.options(interval = 0.25)
        for (x in ss){
            print(mydates[x]); 
            
            spresults <- makeLines(results[which(results$ID %in% ii & results$timestep <= mydates[x]), ])[[1]]
            
            print(plot(crop(blfrac, e), col=colorRampPalette(brewer.pal(5,"Greys"), space="Lab")(10), zlim=c(0,1), breaks=seq(0,1,0.1), legend=F, axes=F, cex.main=1.2, main=paste("MCS ",ii," on ",format(mydates[x],"%d %b %H:%M"),sep="") ) )
            
            print(plot(land_simple, add=T, border="white", axes=F))
            
            print(plot(precip.mcsall[[x]], col=c(rev(terrain.colors(10)),rep("#00A600FF",10)), add=T, legend=F, axes=F)) # alpha=0.8, 
            
            print(plot(spresults, add=T, axes=F, border="black"))
            
            # Plots Precip in the extent of the current MCS, but including all surrounding MCS
            print(plot(precip.mcsall[[x]], zlim=c(0,300), legend.only=T, col=c(rev(terrain.colors(10)),rep("#00A600FF",10)), legend.shrink=0.65, legend.width=1.25, horizontal=T, # smallplot=c(0.2, 0.75, 0.17, 0.2),
                       axis.args=list(at=seq(0,300, 50),
                                      labels=seq(0, 300, 50), 
                                      cex.axis=0.7),
                       legend.args=list(text=expression(Precipitation~'( mm'~hour^-1~')'), side=1, font=2, line=2, cex=0.8)))
            
            # Plots Precip ONLY for the current MCS
#             print(plot(this.precip.mcs[[x]], zlim=c(0,300), legend.only=T, col=c(rev(terrain.colors(10)),rep("#00A600FF",10)), legend.shrink=0.65, legend.width=1.25, horizontal=T, # smallplot=c(0.2, 0.75, 0.17, 0.2),
#                        axis.args=list(at=seq(0,300, 50),
#                                       labels=seq(0, 300, 50), 
#                                       cex.axis=0.7),
#                        legend.args=list(text=expression(Precipitation~'( mm'~hour^-1~')'), side=1, font=2, line=2, cex=0.8)))
            
            print(plot(crop(blfrac, e), zlim=c(0,1), breaks=seq(0,1,0.1), legend.only=T, col=colorRampPalette(brewer.pal(5,"Greys"), space="Lab")(10), legend.shrink=0.75, legend.width=1.75, horizontal=F, 
                       axis.args=list(at=seq(0,1,0.2),
                                      labels=seq(0,1,0.2), 
                                      cex.axis=0.7),
                       legend.args=list(text='Tree Fraction', side=4, font=2, line=2, cex=0.8)))
            ani.pause(0.1)
        } 
    }, img.name=paste("MCS_id",ii,sep=""), outdir="/Users/ajh235/Work/Projects/InternalSabbatical/Results/MCS_movies/html", htmlfile="index.html", ani.width=1000, ani.height=500, verbose=F)
#         }, movie.name=paste("MCS_",ii,".gif",sep=""), interval=0.1, outdir="/Users/ajh235/Work/Projects/InternalSabbatical/Results/MCS_movies", ani.width=1000, ani.height=500)
    # imgdir="/Users/ajh235/Work/Projects/InternalSabbatical/Results/MCS_movies", htmlfile='index.html',
    
    # Contour of number of 10 minute time steps with heavy rain per grid cell
    spresults <- makeLines(results[results$ID %in% ii,])
    myspdf <- spresults[[1]] # Spatial Lines
    
    png(paste("../../Results/TotalIntensity_overBLfrac/MCS_",ii,".png",sep=""), width=1000, height=500)
    print(
        levelplot(crop(blfrac, e), margin=F, par.settings=gntheme, main=paste("Totals from intense precipitation for MCS",ii), at=seq(0,1,0.1), xlab=NULL, ylab=NULL, scales=list(draw=FALSE)) + 
            latticeExtra::layer(sp.polygons(land_simple, col="white")) + 
            contourplot(precip.int.sum, at=c(1,500,1000), labels=list(cex=0.7, labels=c("",500,1000))) + 
            latticeExtra::layer(sp.lines(myspdf))
    )
    dev.off()
        
}
