library(raster)
library(rasterVis)
library(rgdal)
source("functions/consec_drywet.R")
source("functions/diurnal_cycle.R")
source("functions/loadVeg.R")
source("functions/loadOtherAncils.R")
source("functions/makeBoxes.R") 
source("functions/vegPrep.R")
source("functions/patches.R")
source("functions/movies.R")
source("functions/tracker.R")
source("functions/mcsStats.R")
source("functions/initiations.R")
source("getMyData.R")

if (Sys.info()[["sysname"]] == "Darwin"){
    indatadir <- "/Volumes/MYBOOK/WAFR/"
    dlresultsdir <- "/Users/Andy/Work/InternalSabbatical/Results/"
    resultsdir <- "/Users/Andy/Work/InternalSabbatical/Results/"
    scratchdir <- "/Users/Andy/Work/Scratch/"
} else {
    indatadir <- "/data/local/hadhy/ModelData/WAFR/"
    dlresultsdir <- "/data/local/hadhy/Projects/InternalSabbatical/Results/"
    resultsdir <- "/home/h02/hadhy/Projects/InternalSabbatical/Results/"
    scratchdir <- "/data/local/hadhy/Scratch/"
    require(PP,lib.loc="/project/ukmo/rhel6/R")
}

setOptions(tmpdir=scratchdir, todisk=F)

land_simple <- readOGR(dsn=paste(indatadir,"ancils",sep=""), layer="land_rp")

mydata <- getMyData(timestep="10min", overwrite=F)
rb5216.1km.std <- mydata[[1]]
rb5216.4km.std <- mydata[[2]]
rb5216.4km.50k <- mydata[[3]]
rb5216.4km.300k <- mydata[[4]]

models <- c("rb5216.4km.std", "rb5216.4km.50k", "rb5216.4km.300k")
id <- c("s","w","y")

for (x in 1:length(id)){
    print(models[x])
    assign(paste(id[x],".mcs",sep=""), raster(paste(indatadir, "djzx",id[x],"/mcsStats/1000km_10min/allmcstot.tif", sep="")))
    assign(paste(id[x],".pop",sep=""), raster(paste(indatadir, "djzx",id[x],"/popStats/1000km_10min/allpoptot.tif", sep="")))
    
    # MCS: Total precip rate
    assign(paste(id[x],".mcsprtot",sep=""), raster(paste(indatadir, "djzx",id[x],"/mcsStats/1000km_10min/allmcs.totpreciprate.tif", sep="")))
    # MCS: Mean precip rate: broken, but fixed below ...
    #assign(paste(id[x],".mcsprmean",sep=""), raster(paste(indatadir, "djzx",id[x],"/mcsStats/1000km_10min/allmcs.meanpreciprate.tif", sep="")))
  
}
pdf(paste(resultsdir,"MCSanalysis.pdf",sep=""), width=8.27, height=11.69, onefile=T)
# Is there a difference in the total number of MCS between std veg and 50k or 300k smoothing?
mydiffs <- stack(shift(s.mcs - w.mcs, x=-360),shift(s.mcs - y.mcs, x=-360))
myplot <- levelplot(mydiffs, par.settings=RdBuTheme, at=seq(-200,200,25)) + layer(sp.polygons(land_simple))
myplot <- update(myplot, strip = strip.custom(factor.levels = c("MCS count: Std Veg - 50km smoothing", "MCS count: Std Veg - 300km smoothing")))
print(myplot)

# Is there a difference in the total number of POP between std veg and 50k or 300k smoothing?
mydiffs <- stack(shift(s.pop - w.pop, x=-360),shift(s.pop - y.pop, x=-360))
myplot <- levelplot(mydiffs, par.settings=RdBuTheme, at=c(-200,-150,seq(-100,100,10),150,200)) + layer(sp.polygons(land_simple))
myplot <- update(myplot, strip = strip.custom(factor.levels = c("POP count: Std Veg - 50km smoothing", "POP count: Std Veg - 300km smoothing")))
print(myplot)

# Is there a difference in the total "within MCS" precip rate between std veg, 50k and 300k?
mydiffs <- stack(shift(s.mcsprtot - w.mcsprtot, x=-360),shift(s.mcsprtot - y.mcsprtot, x=-360))
myplot <- levelplot(mydiffs, par.settings=RdBuTheme, at=seq(-1,1,0.05)) + layer(sp.polygons(land_simple)) # 
myplot <- update(myplot, strip = strip.custom(factor.levels = c("MCS tot pr: Std Veg - 50km smoothing", "MCS tot pr: Std Veg - 300km smoothing")))
print(myplot)

# Is there a difference in the mean "within MCS" precip rate between std veg, 50k and 300k?
# Something wrong with the following ...
# mydiffs <- stack(shift(s.mcsprmean - w.mcsprmean, x=-360),shift(s.mcsprmean - y.mcsprmean, x=-360))
# myplot <- levelplot(mydiffs, par.settings=RdBuTheme) + layer(sp.polygons(land_simple)) # at=seq(-1,1,0.05)
# myplot <- update(myplot, strip = strip.custom(factor.levels = c("MCS mean pr: Std Veg - 50km smoothing", "MCS mean pr: Std Veg - 300km smoothing")))
# print(myplot)

# Following replicates the above - mean "within MCS" precip rate for std veg, 50k and 300k
mymeans <- stack(shift(s.mcsprtot / s.mcs, x=-360), shift(w.mcsprtot / w.mcs, x=-360), shift(y.mcsprtot / y.mcs, x=-360))*3600
levelplot(mymeans, main="MCS Mean precipitation rate (mm/hr)", names.attr=c("Std veg", "50km smoothing", "300km smoothing"), par.settings=rasterTheme(region=rev(terrain.colors(52))), at=seq(0,50,1))+ layer(sp.polygons(land_simple))

# What is the POP count, and how does it relate to vegetation boundaries
plot(shift(s.pop, x=-360), zlim=c(0,30), xlim=c(-10,16), ylim=c(-10,6), main="POP count @ 10min intervals, with 30% tree contours"); contour(shift(myveg[[1]], x=-360), levels=c(0.3), drawlabels=F, add=T); plot(land_simple, border="grey", add=T)

# Mean "within MCS" precip rate, with (50km) veg boundaries overlaid 
plot(mymeans[[1]], zlim=c(0,30), xlim=c(-10,16), ylim=c(-10,6), main="MCS mean precip rate, with 30% tree contours"); contour(shift(myveg[[1]], x=-360), levels=c(0.3), drawlabels=F, add=T); plot(land_simple, border="grey", add=T)

# Difference in the mean rate between std veg and 50k or 300km smoothing ...
plot(mymeans[[1]] - mymeans[[2]], zlim=c(0,30), xlim=c(-10,16), ylim=c(-10,6), main="MCS mean precip rate (std veg - 50k smooth), with 30% tree contours"); contour(shift(myveg[[1]], x=-360), levels=c(0.3), drawlabels=F, add=T); plot(land_simple, border="grey", add=T)

plot(mymeans[[1]] - mymeans[[3]], zlim=c(0,30), xlim=c(-10,16), ylim=c(-10,6), main="MCS mean precip rate (std veg - 50k smooth), with 30% tree contours"); contour(shift(myveg[[1]], x=-360), levels=c(0.3), drawlabels=F, add=T); plot(land_simple, border="grey", add=T)


dev.off()

