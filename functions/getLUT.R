getLUT <- function(nBndClass=1){
    if (nBndClass == 1){
        myLUT <- data.frame(ID=c(1,2,3,4,7), Landcover=factor(c("tree", "grass", "sparse", "boundary", "orography"), levels=c("tree", "grass", "sparse", "boundary", "orography")[c(3,2,4,1,5)]), Colours=c("dark green", "yellow", "orange", "sienna", "dark grey"), plotOrder=c(4,2,1,3,5))
    } else {
        # Full LUT with all boundary classes
        myLUT <- data.frame(ID=c(1,2,3,4,5,6,7), Landcover=factor(c("tree", "grass", "sparse", "boundary", "boundary, tree", "boundary, grass", "orography"), levels=c("tree", "grass", "sparse", "boundary", "boundary, tree", "boundary, grass", "orography")[c(3,2,6,4,5,1,7)]), Colours=c("dark green", "yellow", "orange", "sienna", "red", "blue", "dark grey"), plotOrder=c(6,2,1,4,5,3,7))
    }
    
    return(myLUT)
}