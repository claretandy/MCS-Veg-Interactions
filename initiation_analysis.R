# Statistics on initiations over boundaries

if (Sys.info()[["sysname"]] == "Darwin"){
    indatadir <- "/Users/ajh235/Work/DataLocal/ModelData/WAFR/"
    resultsdir <- "/Users/ajh235/Work/Projects/InternalSabbatical/Results/"
    scratchdir <- "/Users/ajh235/Work/Scratch/"
} else {
    indatadir <- "/data/local/hadhy/ModelData/WAFR/"
    resultsdir <- "/home/h02/hadhy/Projects/InternalSabbatical/Results/"
    scratchdir <- "/data/local/hadhy/Scratch/"
    require(PP,lib.loc="/project/ukmo/rhel6/R")
}

mycl[mycl==5] <- 4
mycl[mycl==6] <- 4

outpdf <- paste(resultsdir,"initiations.pdf", sep="")
pdf(outpdf, width=12, height=9)

myLUT <- data.frame(ID=c(1,2,3,4,7), Landcover=factor(c("tree", "grass", "sparse", "boundary", "orography"), levels=c("tree", "grass", "sparse", "boundary", "orography")[c(3,2,4,1,5)]), Colours=c("dark green", "yellow", "orange", "brown", "dark grey"), plotOrder=c(4,2,1,3,5))

mycl.f <- as.factor(mycl)
fdat <- myLUT[myLUT$ID %in% levels(mycl.f)[[1]]$ID, ]
levels(mycl.f) <- fdat #[rev(order(fdat$plotOrder)),]


initiations <- results[results$class == 'generation' & (as.numeric(format(results$timestep, "%H")) > 15 & as.numeric(format(results$timestep, "%H")) <= 21),c("x","y")]
coordinates(initiations) <- ~x+y

# Total Land Area per veg class
mycl.f.area <- zonal(areasqkm, mycl.f, fun='sum')
sum(mycl.f.area) # Total land area
mycl.f.area <- zonal(areasqkm, mycl.f, fun='sum') # Total land per class
cbind(myLUT, FractionalArea=100*mycl.f.area[,"sum"] / sum(mycl.f.area)) # Fractional land per class

# How many initiations altogether over land? A: 427
allinit <- results[results$class == 'generation',]
coordinates(allinit) <- ~x+y
sum(table(extract(mycl.f, allinit)))

# How many initiations over boundaries at all times of day? A: 117 (27.4%)
allinit.count <- table(extract(mycl.f, allinit)) # Count per class
allinit.perc  <- 100 * table(extract(mycl.f, allinit)) / sum(table(extract(mycl.f, allinit))) # %
cbind(myLUT, InitCount=allinit.count, InitPercent=allinit.perc)

# What is the percentage of coverage of each surface type? 
# Answer: Tree (21.6%); Grass (25.3%); Sparse (5.3%); Boundary (27.5%); Orography (20.2%)
data.frame(Class=freq(mycl.f)[!is.na(freq(mycl.f)[,"value"]),"value"] , Coverage=100 * freq(mycl.f)[!is.na(freq(mycl.f)[,"value"]),"count"] / sum(freq(mycl.f)[!is.na(freq(mycl.f)[,"value"]),"count"]) )

# How many initiations at different times of day?
initbytime <- results[results$class == 'generation',c("x","y","timestep")]
initbytime$Hour <- as.numeric(format(initbytime$timestep, "%H"))
coordinates(initbytime) <- ~x+y
initbytime$class <- extract(mycl.f, initbytime)
hist(initbytime@data[!is.na(initbytime$class),"Hour"], xlab="Time of Day", main="Frequency of MCS initiations by time of day", breaks=24)
table(initbytime@data[!is.na(initbytime$class),"Hour"])
#  0  1  2  3  4  5  6  7  8 10 11 12 13 14 15 16 17 18 19 20 21 22 23 
# 11 10 10  9 14  9 12  3  2  4  5  7 20 17 45 56 59 42 26 34 16  8  8 

# How many afternoon initiations over boundaries?
aftinit <- results[results$class == 'generation' & (as.numeric(format(results$timestep, "%H")) > 15 & as.numeric(format(results$timestep, "%H")) <= 18),c("x","y")]
coordinates(aftinit) <- ~x+y
aftinit.count <- data.frame(table(extract(mycl.f, aftinit)))[,2]
aftinit.perc  <- data.frame(100* table(extract(mycl.f, aftinit)) / sum(table(extract(mycl.f, aftinit))))[,2]
cbind(myLUT, AftInitCount=aftinit.count, AftInitPerc=aftinit.perc)
# Answer: (land coverage)
# 1:tree (21.6%) 2:grass (25.3%) 3:sparse (5.3%) 4:boundary (27.5%) 7:orography (20.2%)
# 21.656051      21.656051       5.732484        35.031847          15.923567 

# Where are afternoon initiations over boundaries?
png()
print(
    levelplot(mycl.f, att="Landcover", main="Afternoon MCS Initiations Over Vegetation Classes", xlab=NULL, ylab=NULL, scales=list(draw=FALSE), col.regions=as.character(fdat$Colours)) + 
        latticeExtra::layer(sp.polygons(land_simple)) +
        latticeExtra::layer(sp.points(aftinit[!is.na(extract(mycl.f, aftinit)),], pch=1, col="black"))
)    	

dev.off()
