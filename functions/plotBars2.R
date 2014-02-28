plotBars2 <- function(inpr, mycl.z, mytitle){
    require(ggplot2)
    require(raster)
    source("functions/getLUT.R")
    
    tmpresult <- data.frame(zonal(x=inpr, z=mycl.z, fun=function(x, na.rm=T){return(mean(x))}))
    myresult <- data.frame(zone=floor(tmpresult$zone/10), class=tmpresult$zone - (10*floor(tmpresult$zone/10)), mean=tmpresult$value)
    myLUT <- getLUT(nBndClass=1)
    
    myresult$se <- data.frame(zonal(x=inpr, z=mycl.z, fun=function(x, na.rm=T){return(sd(x)/sqrt(length(x)))}))[,2]
    myresult$sd <- data.frame(zonal(x=inpr, z=mycl.z, fun=function(x, na.rm=T){return(sd(x))}))[,2]
    myresult$ci <- 1.96 * myresult$sd
    
    myresult <- merge(myresult, myLUT, by.x="class", by.y="ID")
    
    myresult[myresult$Landcover == "sparse" & myresult$zone %in% 6:15, c(3:6)] <- NA
    myresult$zone <- as.factor(myresult$zone)
    myLUT <- getLUT(nBndClass=1)
    print(
        ggplot(myresult, aes(x=zone, y=mean, fill=Landcover)) + 
            geom_bar(position=position_dodge(), stat="identity") +
            scale_fill_manual(values=as.character(myLUT$Colours)[c(3,2,4,1,5)]) +
            labs(title=mytitle, y="Precipitation (mm/hour)", x="Zone") +
            geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                          width=.2,                    # Width of the error bars
                          position=position_dodge(.9))    
    )
}

plotBars3 <- function(inpr, mycl.z, mytitle1, mytitle2){
    require(ggplot2)
    require(raster)
    source("functions/getLUT.R")
    
    spp.r <- floor(mycl.z/10)
    
    # Get total precip per class per zone
    tmpresult <- data.frame(zonal(x=inpr, z=mycl.z, fun=function(x, na.rm=T){return(sum(x))}))
    myresult <- data.frame(zone=floor(tmpresult$zone/10), class=tmpresult$zone - (10*floor(tmpresult$zone/10)), prsum=tmpresult$value)
    
    # Total precip per zone
    zone_prec <- zonal(inpr, spp.r, fun=sum)
    myresult <- merge(myresult, zone_prec, by=c("zone"), all=T)
    
    # Total area per class per zone 
    tmpresult <- data.frame(zonal(area(spp.r), mycl.z, fun=sum))
    clas_area <- data.frame(zone=floor(tmpresult$zone/10), class=tmpresult$zone - (10*floor(tmpresult$zone/10)), clArea=tmpresult$value)
    myresult <- merge(myresult, clas_area, by=c("zone","class"), all=T)
    
    # Total area per zone 
    zone_area <- zonal(area(spp.r), spp.r, fun=sum)
    myresult <- merge(myresult, zone_area, by=c("zone"), all=T)

    colnames(myresult) <- c("Zone", "Class", "prClZone","prZone","areaClZone","areaZone")
    
    # Precip fraction
    myresult$PrFrac <- myresult$prClZone / myresult$prZone
    
    # Landcover fraction
    myresult$ClFrac <- myresult$areaClZone / myresult$areaZone
    
    # Ratio TotPrecipOverClass : FractionalClassCover
    myresult$ratio <- myresult$PrFrac / myresult$ClFrac

    # Plot the results
    
    myLUT <- getLUT(nBndClass=1)
    
    myresult <- merge(myresult, myLUT, by.x="Class", by.y="ID")
    
#     myresult[myresult$Landcover == "sparse" & myresult$Zone %in% 6:15, c(3:9)] <- NA
    myresult$Zone <- as.factor(myresult$Zone)
    myLUT <- getLUT(nBndClass=1)
#     browser()
    print(
        ggplot(myresult, aes(x=Zone, y=ratio, fill=Landcover)) + 
            geom_bar(position=position_dodge(), stat="identity") +
            ylim(0,3) +
            scale_fill_manual(values=as.character(myLUT$Colours)[c(3,2,4,1,5)]) +
            labs(title=mytitle1, y="Precipitation (mm/hour)", x="Zone") 
    )
    
    print(
        ggplot(myresult, aes(x=ClFrac, y=PrFrac, colour=Landcover)) + 
#             geom_point(aes(shape=Zone)) +
            geom_text(aes(label=Zone)) + 
            scale_colour_manual(values=as.character(myLUT$Colours)[c(3,2,4,1,5)]) +
            coord_fixed(ratio=1/1) +
            ylim(0,1) + xlim(0,1) +
            labs(title=mytitle2, x="Fraction of vegetation class in zone", y="Fraction of total precipitation in zone") +
            geom_abline(intercept = 0, slope = 1)
        )
    return(myresult)
}
