# Statistics on initiations over boundaries
# Need to load results from loadData_diurnal.R ....

source("functions/loadAllAncils.R")
source("functions/getLUT.R")
source("functions/adjCoords.R")

if (Sys.info()[["sysname"]] == "Darwin"){
    indatadir <- "/Users/ajh235/Work/DataLocal/ModelData/WAFR/"
    resultsdir <- "/Users/ajh235/Work/Projects/InternalSabbatical/Results/"
    dlresultsdir <- "/Users/ajh235/Work/DataLocal/Projects/InternalSabbatical/Results/"
    scratchdir <- "/Users/ajh235/Work/Scratch/"
} else {
    indatadir <- "/data/local/hadhy/ModelData/WAFR/"
    resultsdir <- "/home/h02/hadhy/Projects/InternalSabbatical/Results/"
    scratchdir <- "/data/local/hadhy/Scratch/"
    require(PP,lib.loc="/project/ukmo/rhel6/R")
}

outpdf <- paste(resultsdir,"initiations.pdf", sep="")
pdf(outpdf, width=12, height=9)

myLUT <- getLUT(nBndClass=1)

# Get results data.frame
id="s"
x=1
threshold=1000
timestep="10min"
mcs.infile <- paste(indatadir,"djzx",id[x],"/patches/",threshold,"km_",timestep,"/allmcs.1000km.vrt", sep="")
mcs <- adjCoords(brick(mcs.infile))  
mydata <- getMyData(timestep=timestep, var="lsrain", overwrite=F)
inbr <- adjCoords(mydata[[2]])
results <- tracker2(threshold=threshold, mcs, inbr, id=id[x], timestep=timestep, indatadir=indatadir, dlresultsdir=dlresultsdir, overwrite=F)

# First, get data in Rotated Pole
ancils <- loadAllAncils(myproj="rp", nBndClass=1, model="rb5216.4km.std")
mycl.f <- ancils[[5]]
land_simple <- ancils[[9]]

# Initiations in the afternoon / evening 
initiations <- results[results$class == 'generation' & (as.numeric(format(results$timestep, "%H")) >= 15 & as.numeric(format(results$timestep, "%H")) <= 21),c("x","y")]
coordinates(initiations) <- ~x+y

# Total Land Area per veg class
areasqkm <- area(mycl)
mycl.f.area <- zonal(areasqkm, mycl.f, fun='sum')
sum(mycl.f.area) # Total land area = 4,411,886 km2
mycl.f.area <- zonal(areasqkm, mycl.f, fun='sum') # Total land per class
cbind(myLUT, FractionalArea=100*mycl.f.area[,"sum"] / sum(mycl.f.area)) # Fractional land per class

# How many initiations altogether over land? A: 427
allinit <- results[results$class == 'generation',]
coordinates(allinit) <- ~x+y
sum(table(extract(mycl.f, allinit)))

# How many initiations over boundaries at all times of day? A: 117 (27.4%)
allinit.count <- table(extract(mycl.f, allinit)) # Count per class
allinit.perc  <- 100 * table(extract(mycl.f, allinit)) / sum(table(extract(mycl.f, allinit))) # %
cbind(myLUT, InitCount=as.vector(allinit.count), InitPercent=as.vector(allinit.perc))

# What is the percentage of coverage of each surface type? 
# Answer: Tree (21.6%); Grass (25.3%); Sparse (5.3%); Boundary (27.5%); Orography (20.2%)
data.frame(Class=freq(mycl.f)[!is.na(freq(mycl.f)[,"value"]),"value"] , Coverage=100 * freq(mycl.f)[!is.na(freq(mycl.f)[,"value"]),"count"] / sum(freq(mycl.f)[!is.na(freq(mycl.f)[,"value"]),"count"]) )

# How many initiations at different times of day?
initbytime <- results[results$class == 'generation',c("x","y","timestep")]
initbytime$Hour <- as.numeric(format(initbytime$timestep, "%H"))
coordinates(initbytime) <- ~x+y
initbytime$class <- extract(mycl.f, initbytime)
freqdata <- data.frame(table(initbytime@data[!is.na(initbytime$class),"Hour"]))
ggplot(data=freqdata, aes(x=Var1, y=Freq)) + geom_bar(fill="white", colour="black") + labs(title="Frequency of MCS initiations by time of day", x="Time of day", y="Frequency")
# hist(initbytime@data[!is.na(initbytime$class),"Hour"], xlab="Time of Day", main="Frequency of MCS initiations by time of day", breaks=0:24)
table(initbytime@data[!is.na(initbytime$class),"Hour"])
#  0  1  2  3  4  5  6  7  8 10 11 12 13 14 15 16 17 18 19 20 21 22 23 
# 11 10 10  9 14  9 12  3  2  4  5  7 20 17 45 56 59 42 26 34 16  8  8 

# How many afternoon initiations over boundaries?
aftinit <- results[results$class == 'generation' & (as.numeric(format(results$timestep, "%H")) >= 16 & as.numeric(format(results$timestep, "%H")) <= 17),c("x","y")]
coordinates(aftinit) <- ~x+y
aftinit.count <- data.frame(table(extract(mycl.f, aftinit)))[,2]
aftinit.perc  <- data.frame(100* table(extract(mycl.f, aftinit)) / sum(table(extract(mycl.f, aftinit))))[,2]
cbind(myLUT, AftInitCount=aftinit.count, AftInitPerc=aftinit.perc)
initpts_rp <- SpatialPoints(aftinit[!is.na(extract(mycl.f, aftinit)),], CRS("+proj=ob_tran +o_proj=longlat +o_lon_p=175.3000030517578 +o_lat_p=77.4000015258789 +lon_0=180 +ellps=sphere"))

# Answer: (land coverage)
# 1:tree (21.6%) 2:grass (25.3%) 3:sparse (5.3%) 4:boundary (27.5%) 7:orography (20.2%)
# 21.656051      21.656051       5.732484        35.031847          15.923567 

# Where are afternoon initiations over boundaries?
# First, get data in Lat/Long
ancils <- loadAllAncils(myproj="ll", nBndClass=1, model="rb5216.4km.std")
mycl.f <- ancils[[5]]
land_simple <- ancils[[9]]

# Check that spTransform reprojects the country boundaries correctly ...
land_rp <- readOGR(dsn=paste(indatadir,"ancils",sep=""), layer="land_rp")
projection(land_rp) <- CRS("+proj=ob_tran +o_proj=longlat +o_lon_p=175.3000030517578 +o_lat_p=77.4000015258789 +lon_0=180 +ellps=sphere")
land_rp2ll <- spTransform(land_rp, CRSobj=CRS("+init=epsg:4326"), use_ob_tran=T)

# Reproject afternoon initiation points
initpts_ll <- spTransform(initpts_rp, CRSobj=CRS("+init=epsg:4326"), use_ob_tran=T)

print(
    levelplot(mycl.f, att="Landcover", maxpixels=600000, main="Afternoon (16-18Z) MCS Initiations Over Vegetation Classes", xlim=c(-12,10), ylim=c(4,18), xlab=NULL, ylab=NULL, col.regions=as.character(myLUT$Colours)) + # scales=list(draw=FALSE), xlim=c(-3,10), ylim=c(12,20)
        latticeExtra::layer(sp.polygons(land_simple, lty=2)) +
        latticeExtra::layer(sp.points(initpts_ll, pch="+", cex=2, col="black")) #+
        #latticeExtra::layer(sp.polygons(spp))
)

dev.off()
