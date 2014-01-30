makeLines <- function(mcstracks){

    require(sp)
    
    jjj <- unique(mcstracks$ID)
    li <- 1
    mylines <- list()
#     if(exists("mytext")){rm(mytext)}
    for (jj in jjj){
        if(!exists("this.mytext")){
            this.mytext <- data.frame(mcstracks[mcstracks$ID == jj,c("x","y","ID")][1,])
        } else {
            this.mytext <- rbind(this.mytext, mcstracks[mcstracks$ID == jj,c("x","y","ID")][1,])
        }
        myline <- Line(mcstracks[mcstracks$ID == jj,c("x","y")])
        mylines[[li]] <- Lines(list(myline), ID=li) # or ID=li
        #         assign(x=paste("myline",li,sep=""), value=mylines)
        li <- li+1
    }
    myspl <- SpatialLines(mylines, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    mylines <- SpatialLinesDataFrame(sl=myspl, data=data.frame(lineid=1:length(jjj), MCSid=jjj))
    
    return(list(mylines, this.mytext))
}
