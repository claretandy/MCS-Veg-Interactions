subsetHours <- function(inpr, hr=15, plen=3){
    phr <- hr-plen # Gets the 6 hour period before the given hour
    phr <- ifelse(phr < 0, phr + 24, phr)
    if(phr > hr){
        myss <- which(as.integer(format(getZ(inpr), "%H")) >= phr | as.integer(format(getZ(inpr), "%H")) < hr)
    } else {
        myss <- which(as.integer(format(getZ(inpr), "%H")) >= phr & as.integer(format(getZ(inpr), "%H")) < hr)
    }
    return(list(myss, phr))
}