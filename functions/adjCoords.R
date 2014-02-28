adjCoords <- function(r){
    this.res <- res(r)[2]
    if ( this.res <= 0.0136){
        r <- shift(r, x=0.0068, y=0.00675)
    }
    
    if (this.res <= 0.0365 & this.res > 0.0355){
        r <- shift(r, x=0.01798, y=0.018)
        res(r) <- 0.036
    }
    return(r)
}
