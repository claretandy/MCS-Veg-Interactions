makeLines <- function(tracks){

    require(sp)
    
    iii <- unique(tracks$ID)
    li <- 1
    mylines <- list()
#     if(exists("mytext")){rm(mytext)}
    for (ii in iii){
        if(!exists("this.mytext")){
            this.mytext <- data.frame(tracks[tracks$ID == ii,c("x","y","ID")][1,])
        } else {
            this.mytext <- rbind(this.mytext, tracks[tracks$ID == ii,c("x","y","ID")][1,])
        }
        myline <- Line(tracks[tracks$ID == ii,c("x","y")])
        mylines[[li]] <- Lines(list(myline), ID=li) # or ID=li
        #         assign(x=paste("myline",li,sep=""), value=mylines)
        li <- li+1
    }
    myspl <- SpatialLines(mylines)
    myspdf <- SpatialLinesDataFrame(sl=myspl, data=data.frame(lineid=1:length(iii), MCSid=iii))
    
    return(list(myspdf, this.mytext))
}
