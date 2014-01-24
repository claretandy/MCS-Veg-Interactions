makeMovie <- function(inbr1, inbr2, inbr3, inbr4, land_simple, sppa, dlresultsdir, models, boxes=F, type="mcs", overwrite=T){
	moddates <- getZ(inbr1)

	resultsdir <- paste(dlresultsdir, "movie_4runs/",type,"/", sep="")
	if (!file.exists(resultsdir)){ dir.create(resultsdir) }
	# for (d in 1:2){
	for (d in 1:length(moddates)){
		print(format(moddates[d], "%b-%d %H:%M"))
		jpeg(paste(resultsdir, type, ".",format(moddates[d],"%d%H.%M"),".jpg", sep=""), width=1100,height=800, quality=100)
		par(mfrow=c(2,2))
		par(mar=c(1,1,1,1))
		for (m in 1:4){
			r <- get(paste("inbr",m,sep=""))[[d]]
			r <- calc(r, function(x){x*3600})
			r[r<1] <- NA
            if (type=="precip"){
                plot(shift(r, x=-360), breaks=c(0,10,20,30,40,60,100,200), col=rev(terrain.colors(7)), axes=F, ext=extent(sppa), main=models[m], zlim=c(0,200))
            }
			if (type=="mcs"){
			    plot(shift(r, x=-360), col=c("orange","green"), axes=F, ext=extent(sppa), main=models[m], zlim=c(0,200))
			}
			
			plot(land_simple, border="grey", add=T)
            if (boxes==T){
                plot(sppa, add=T)
            }
			
			# browser()
			text(x=4, y=-9, label=format(moddates[d], "%b-%d %H:%M"), cex=2)
		}
		
		dev.off()
#         browser()
	}

}