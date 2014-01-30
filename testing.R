
library(animation)

mcs_tracks_all <- shift(brick(paste(mcsrst.path,"mcs_tracking.vrt",sep="")), x=-360)
precip.all <- shift(stack(rb5216.4km.std), x=-360)

# Mask only MCS precip
precip.mcsall <- mask(precip.all, mcs_tracks_all)
precip.mcsall.int <- calc(precip.mcsall, fun=function(x){x[x<(20/3600)]<-NA; return(x)})

for (ii in iii[c(11,16)]){
    print(ii)
    
    # Get the correct time slices
    mydates <- unique(results[results$ID == ii,"timestep"])
    this.mcs <- subset(mcs_tracks_all, which(getZ(inbr) %in% mydates))
    this.precip <- subset(precip.all, which(getZ(inbr) %in% mydates))
    
    # Mask out other MCS and precip < 20mm per hour 
    this.mcs <- calc(this.mcs, fun=function(x){y <- ifelse(x==ii,1,NA); return(y)})
    precip.int <- calc(this.precip, fun=function(x){x[x<(20/3600)]<-NA; return(x)})
    mcs.inttrack <- mask(this.mcs, precip.int)
    
    # Combine all MCS intense patches together to show where the MCS intense precip has fallen
    mcs.int.sum <- calc(mcs.inttrack, fun=function(x){x[is.na(x)]<-0; return(sum(x))})
    mcs.int.sum <- trim(mcs.int.sum, value=0, padding=1)
    e <- extent(mcs.int.sum)
    
    # Calculate all combined intense precip for all timesteps of the MCS
    this.precip.mcs <- mask(crop(this.precip, e), crop(this.mcs, e))*3600
    precip.int.sum <- calc(this.precip.mcs, fun=function(x){
        x[x<20 | is.na(x)] <- 0;
        return(sum(x))})
    
#     browser()
    # Plot movie zoomed in on the MCS
#     saveGIF({ 
    saveHTML({
        ani.options(interval = 0.25)
        for (x in 1:length(mydates)){
            print(mydates[x]); 
            
            spresults <- makeLines(results[which(results$ID %in% ii & results$timestep <= mydates[x]), ])[[1]]
            
            print(plot(crop(blfrac, e), col=colorRampPalette(brewer.pal(5,"Greys"), space="Lab")(10), zlim=c(0,1), breaks=seq(0,1,0.1), legend=F, axes=F, cex.main=1.2, main=paste("MCS ",ii," on ",format(mydates[x],"%d %b %H:%M"),sep="") ) )
            print(plot(land_simple, add=T, border="white", axes=F))
            print(plot(this.precip.mcs[[x]], col=c(rev(terrain.colors(10)),rep("#00A600FF",10)), add=T, legend=F, axes=F)) # alpha=0.8, 
            print(spresults, add=T)
            print(plot(this.precip.mcs[[x]], zlim=c(0,300), legend.only=T, col=c(rev(terrain.colors(10)),rep("#00A600FF",10)), legend.shrink=0.65, legend.width=1.25, horizontal=T, # smallplot=c(0.2, 0.75, 0.17, 0.2),
                       axis.args=list(at=seq(0,300, 50),
                                      labels=seq(0, 300, 50), 
                                      cex.axis=0.7),
                       legend.args=list(text=expression(Precipitation~'( mm'~hour^-1~')'), side=1, font=2, line=2, cex=0.8)))
            print(plot(crop(blfrac, e), zlim=c(0,1), breaks=seq(0,1,0.1), legend.only=T, col=colorRampPalette(brewer.pal(5,"Greys"), space="Lab")(10), legend.shrink=0.75, legend.width=1.75, horizontal=F, 
                       axis.args=list(at=seq(0,1,0.2),
                                      labels=seq(0,1,0.2), 
                                      cex.axis=0.7),
                       legend.args=list(text='Tree Fraction', side=4, font=2, line=2, cex=0.8)))
            ani.pause(0.1)
        } 
    }, img.name=paste("MCS_id",ii,sep=""), outdir="/Users/ajh235/Work/Projects/InternalSabbatical/Results/MCS_movies/html", htmlfile="myindex.html", ani.width=1000, ani.height=500, verbose=F)
#         }, movie.name=paste("MCS_",ii,".gif",sep=""), interval=0.1, outdir="/Users/ajh235/Work/Projects/InternalSabbatical/Results/MCS_movies", ani.width=1000, ani.height=500)
    # imgdir="/Users/ajh235/Work/Projects/InternalSabbatical/Results/MCS_movies", htmlfile='index.html',
    
    # Contour of number of 10 minute time steps with heavy rain per grid cell
    f <- levelplot(crop(blfrac, e), par.settings=gntheme, margin=F, main=paste("MCS",ii), at=seq(0,1,0.1)) +
            contourplot(precip.int.sum, at=c(1,500,1000), labels=list(cex=0.7, labels=c("",500,1000)))
    
        
    # Create contours of intense rainfall at different dates
#     for (di in 1:length(mydates)){
#         ifstate <- (format(mydates[di], "%H") == "00" | format(mydates[di], "%H") == "03" | format(mydates[di], "%H") == "06" | format(mydates[di], "%H") == "09" | format(mydates[di], "%H") == "12" | format(mydates[di], "%H") == "15" | format(mydates[di], "%H") == "18" | format(mydates[di], "%H") == "21")
#         if (format(mydates[di], "%M") == "00"){ #  & ifstate
#             print(mydates[di])
# #             browser()
#             precip.di <- mask(crop(this.precip[[di]], e), crop(this.mcs[[di]], e))*3600
#             precip.di[is.na(precip.di)] <- 0
#             if (format(mydates[di], "%H") != "00"){
#                 f <- f + contourplot(precip.di, at=c(1,20), labels=list(cex=0.6, labels=c("one",format(mydates[di], "%H"))), margin=F, pretty=F)                
#             } else {
#                 f <- f + contourplot(precip.di, at=c(1,20), labels=list(cex=0.6, labels=c("one",format(mydates[di], "%d/%m"))), margin=F, pretty=F)
#             }
#         }
#     }
    
    # Get lines for MCS ii, and plot over mcs.int.sum
    
    spresults <- makeLines(results[results$ID %in% ii,])[[1]]
    f + latticeExtra::layer(sp.lines(spresults))
    print(f)
    
}
