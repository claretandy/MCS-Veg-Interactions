# Calculate variables, and do the 'killer' plot

library(raster)
library(ncdf)
library(RColorBrewer)
require(rasterVis)

ems <- c("s", "w", "y") # "u", 

windDir <- function(u, v) {
	return((180 / pi) * atan(u/v) + ifelse(v>=0,180,ifelse(u>=0,360,0)))
}

if (Sys.info()[["sysname"]] == "Darwin"){
    indatadir <- "/Users/ajh235/Work/DataLocal/ModelData/WAFR/"
    dlresultsdir <- "/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/"
    resultsdir <- "/Users/ajh235/Work/Projects/InternalSabbatical/Results/"
    scratchdir <- "/Users/ajh235/Work/Scratch/"
} else {
    indatadir <- "/data/local/hadhy/ModelData/WAFR/"
    dlresultsdir <- "/data/local/hadhy/Projects/InternalSabbatical/Results/"
    resultsdir <- "/home/h02/hadhy/Projects/InternalSabbatical/Results/"
    scratchdir <- "/data/local/hadhy/Scratch/"
    require(PP,lib.loc="/project/ukmo/rhel6/R")
}


# Loop through each ensemble member 
for (em in ems){
	print(em)
	# Load data ...
	mypath <- paste(indatadir,"djzx",em,"/layers_nc/",sep="")
	myx <- brick(paste(mypath,"3209.nc",sep=""), varname="x-wind")
	myy <- brick(paste(mypath,"3209.nc",sep=""), varname="y-wind")
	
	print("Crop to smaller area ...")
	myx <- crop(myx, extent(355,365,-5,5))
	myy <- crop(myy, extent(355,365,-5,5))
	
	myscale <- ifelse(em=="s","std",ifelse(em=="w","50k","300k"))
	vegfrac <- brick(paste(indatadir,"ancils/km4/qrparm.veg.frac_4km.",myscale,".tif", sep=""))
	vegfrac <- crop(vegfrac, extent(355,365,-5,5))
	tree <- calc(subset(vegfrac, c(1,2,5)), fun=sum, filename=paste(mypath,"tree_frac.tif", sep=""), format="GTiff", overwrite=T)
# 	grass <- calc(subset(vegfrac, c(3,4)), fun=sum)
	
	# Filter tree in x and y directions
	tree.u <- focal(tree, w=matrix(c(0,-1,0,0,0,0,0,1,0), nrow=3), pad=T, na.rm=F, fun=sum, filename=paste(mypath,"tree_u.tif", sep=""), format="GTiff", overwrite=T) # x direction
	tree.v <- focal(tree, w=matrix(c(0,0,0,1,0,-1,0,0,0), nrow=3), pad=T, na.rm=F, fun=sum, filename=paste(mypath,"tree_v.tif", sep=""), format="GTiff", overwrite=T) # y direction
	
	print("Resampling x.wind ...")
	myx.reg <- resample(myx, vegfrac, method="bilinear", filename=paste(mypath,"xwind.resampled.tif", sep=""), format="GTiff", overwrite=T) # subset(myx, 110:120)
	print("Resampling y.wind ...")
	myy.reg <- resample(myy, vegfrac, method="bilinear", filename=paste(mypath,"ywind.resampled.tif", sep=""), format="GTiff", overwrite=T) # subset(myy, 110:120)
	
	print("Cropping rain ...")
	rain <- brick(paste(mypath,"5216.nc", sep=""), varname="lsrain")
	rain <- crop(rain, extent(355,365,-5,5), filename=paste(mypath,"rain.cropped.tif", sep=""), format="GTiff", overwrite=T)
	
# 	ws <- overlay(myx.reg, myy.reg, fun=function(x,y){sqrt((x^2 + y^2))}, filename=paste(mypath,"windspeed.tif", sep=""), format="GTiff", overwrite=T)
# 	wdir <- overlay(myx.reg, myy.reg, fun=function(x,y){windDir(x,y)}, filename=paste(mypath,"windDirection.tif", sep=""), format="GTiff", overwrite=T)
	
# 	browser()
	print("Calculating windVeg ...")
	windveg <- overlay(myx.reg, myy.reg, tree.u, tree.v, fun=function(u,v,tu,tv){return((u*tu) + (v*tv))}, recycle=T, filename=paste(mypath,"windVeg.tif", sep=""), format="GTiff", overwrite=T)

}

# Plot some things ....

# Get time for temporal subsetting later on ..
myx <- brick(paste(mypath,"3209.nc",sep=""), varname="x-wind")
mytimes <- round(as.POSIXct(getZ(myx)*86400, origin="2006-08-15 06:00:00"))

for (em in ems){                                                                                                                   
	print(em)
	mypath <- paste(indatadir,"djzx",em,"/layers_nc/",sep="")
	myscale <- ifelse(em=="s","std",ifelse(em=="w","50k","300k"))
# 	vegfrac <- brick(paste("/data/local/hadhy/ModelData/WAFR/ancils/km4/qrparm.veg.frac_4km.",myscale,".tif", sep=""))
# 	vegfrac <- crop(vegfrac, extent(355,365,-5,5))
	tree <- raster(paste(mypath,"tree_frac.tif", sep=""))

	tree.u <- raster(paste(mypath,"tree_u.tif", sep="")) # x direction
	tree.v <- raster(paste(mypath,"tree_v.tif", sep="")) # y direction

	myx.reg <- brick(paste(mypath,"xwind.resampled.tif", sep=""))
	myy.reg <- brick(paste(mypath,"ywind.resampled.tif", sep=""))
	windveg <- brick(paste(mypath,"windVeg.tif", sep=""))
	rain <- brick(paste(mypath,"rain.cropped.tif", sep=""))

	assign(paste("tree.",em,sep=""), tree)
	assign(paste("windveg.",em,sep=""), windveg)
	assign(paste("windx.",em,sep=""), myx.reg)
	assign(paste("windy.",em,sep=""), myy.reg)
	assign(paste("rain.",em,sep=""), rain)
	assign(paste("tree.u.",em,sep=""),tree.u)
	assign(paste("tree.v.",em,sep=""),tree.v)
} 

pdf(paste(resultsdir,"tree_frac.pdf",sep=""), width=10, height=15)
print(levelplot(stack(tree.s, tree.w, tree.y), col.regions=rev(terrain.colors(20)), names.attr=c("Std veg", "50km smooth", "300km smooth"), main="Tree Fractions"))
# print(levelplot(stack(tree.u.s, tree.v.s, tree.u.w, tree.v.w, tree.u.y, tree.v.y), names.attr=c("i+1 - i-1 Std veg","j+1 - j-1 Std veg", "i+1 - i-1 50km smooth","j+1 - j-1 50km smooth", "i+1 - i-1 300km smooth", "j+1 - j-1 300km smooth"), main="Vegetation Gradients", col.regions=colorRampPalette(brewer.pal(n=9, "RdYlGn"))(20)))
print(levelplot(stack(tree.u.s, tree.v.s), names.attr=c("i+1 - i-1 Std veg","j+1 - j-1 Std veg"), main="Vegetation Gradients", col.regions=colorRampPalette(brewer.pal(n=9, "RdYlGn"))(20)))
print(levelplot(stack(tree.u.w, tree.v.w), names.attr=c("i+1 - i-1 50km smooth","j+1 - j-1 50km smooth"), main="Vegetation Gradients", col.regions=colorRampPalette(brewer.pal(n=9, "RdYlGn"))(20)))
print(levelplot(stack(tree.u.y, tree.v.y), names.attr=c("i+1 - i-1 300km smooth","j+1 - j-1 300km smooth"), main="Vegetation Gradients", col.regions=colorRampPalette(brewer.pal(n=9, "RdYlGn"))(20)))

# Summarise by time of day ...
ss1 <- which(as.numeric(format(mytimes, "%H")) >= 0 & as.numeric(format(mytimes, "%H")) <6 )
ss2 <- which(as.numeric(format(mytimes, "%H")) >= 6 & as.numeric(format(mytimes, "%H")) < 12 )
ss3 <- which(as.numeric(format(mytimes, "%H")) >= 12 & as.numeric(format(mytimes, "%H")) < 18 )
ss4 <- which(as.numeric(format(mytimes, "%H")) >= 18 & as.numeric(format(mytimes, "%H")) < 24 )
ss5 <- which(as.numeric(format(mytimes, "%H")) >= 15 & as.numeric(format(mytimes, "%H")) < 21 )

dev.off()


pdf(paste(resultsdir,"wind-tree_gradients.pdf",sep=""), width=17, height=11)

# Sample the more extreme tree fractions
treeGrad <- tree.u.s
treeGrad[treeGrad > -0.3 & treeGrad < 0.3] <- NA
mypts <- sampleStratified(treeGrad*10, size=2, exp=1, na.rm=T, xy=F, sp=T)

tree.u.s.pts <- extract(tree.u.s, mypts)
tree.v.s.pts <- extract(tree.v.s, mypts)

rasterVis::levelplot(stack(tree.u.s, tree.v.s), col.regions=rev(terrain.colors(20))) + latticeExtra::layer(sp.points(mypts, cex=1.5, pch="+")) + latticeExtra::layer(sp.text(loc=coordinates(mypts), txt=paste(row.names(mypts), "\n(u: ", round(tree.u.s.pts, digits=2), ")", "\n(v: ", round(tree.v.s.pts, digits=2), ")" , sep=""), pos=3, cex=0.8))

mydata <- extract(windveg.s, mypts)
rm(mydf)
for (x in 1:nrow(mydata)){ 
    if (!exists("mydf")){ 
        mydf <- data.frame(time=mytimes, values=mydata[1,], pointID=1) 
    } else { 
        mydf <- rbind(mydf, data.frame(time=mytimes, values=mydata[x,], pointID=x)) 
    }
}
h <- ggplot(data=mydf, aes(x=time, y=values))
h + geom_line() + facet_grid(pointID ~ .) + labs(title="Wind:Vegetation Gradient", x="Time", y="Wind:Veg Gradient")

myrain <- extract(rain.s, mypts)*60*60
rm(myrr)
for (x in 1:nrow(myrain)){ 
    if (!exists("myrr")){ 
        myrr <- data.frame(time=mytimes, values=myrain[1,], pointID=1) 
    } else { 
        myrr <- rbind(myrr, data.frame(time=mytimes, values=myrain[x,], pointID=x)) 
    }
}
h <- ggplot(data=myrr, aes(x=time, y=values))
h + geom_line(colour="blue") + facet_grid(pointID ~ .) + labs(title="Instantaneous Precipitation", x="Time", y="Precipitation mm/hour")

##############################################
# Now, do the same thing, but on nearby points with no gradient ...
mypts.r <- rasterize(mypts, tree.u.s, field="cell")
mypts.r.buff <- buffer(mypts.r, width=(0.036*5))
treeGrad <- tree.u.s
treeGrad[treeGrad < -0.05 | treeGrad > 0.05] <- NA
treeGrad.pts <- mask(treeGrad, mypts.r.buff)

mypts.nograd <- sampleStratified(treeGrad.pts*10, size=16, exp=1, na.rm=T, xy=F, sp=T)
tree.u.s.ptsng <- extract(tree.u.s, mypts.nograd)
tree.v.s.ptsng <- extract(tree.v.s, mypts.nograd)

rasterVis::levelplot(stack(tree.u.s, tree.v.s), col.regions=rev(terrain.colors(20))) + latticeExtra::layer(sp.points(mypts.nograd, cex=1.5, pch="+")) + latticeExtra::layer(sp.text(loc=coordinates(mypts.nograd), txt=paste(row.names(mypts.nograd), "\n(u: ", round(tree.u.s.ptsng, digits=2), ")", "\n(v: ", round(tree.v.s.ptsng, digits=2), ")" , sep=""), pos=3, cex=0.8))

mydata <- extract(windveg.s, mypts.nograd)
rm(mydf)
for (x in 1:nrow(mydata)){ 
    if (!exists("mydf")){ 
        mydf <- data.frame(time=mytimes, values=mydata[1,], pointID=1) 
    } else { 
        mydf <- rbind(mydf, data.frame(time=mytimes, values=mydata[x,], pointID=x)) 
    }
}
h <- ggplot(data=mydf, aes(x=time, y=values))
h + geom_line() + facet_grid(pointID ~ .) + labs(title="Wind:Vegetation Gradient on NO Tree Gradient", x="Time", y="Wind:Veg Gradient")

myrain <- extract(rain.s, mypts.nograd)*60*60
rm(myrr)
for (x in 1:nrow(myrain)){ 
    if (!exists("myrr")){ 
        myrr <- data.frame(time=mytimes, values=myrain[1,], pointID=1) 
    } else { 
        myrr <- rbind(myrr, data.frame(time=mytimes, values=myrain[x,], pointID=x)) 
    }
}
h <- ggplot(data=myrr, aes(x=time, y=values))
h + geom_line(colour="blue") + facet_grid(pointID ~ .) + labs(title="Instantaneous Precipitation on NO Tree Gradient", x="Time", y="Precipitation mm/hour")


# Plot example at time of initiation ...


dev.off()

