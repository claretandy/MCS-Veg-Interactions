require(raster)

# Sums the total precip occuring within a mask at all time steps
# myprecip and mymask must have the same number of layers

MaskPrecip <- function(myprecip, mymask, outfile, overwrite=T){
	if (file.exists(outfile) & overwrite == F){
		return(brick(outfile))
	} else {
		for (i in 1:nlayers(mymask)){
# 			print(i)
			tmp <- overlay(myprecip[[i]], mymask[[i]], fun=function(x,y){return(x*y)})
			tmp[is.na(tmp)] <- 0
			if (i == 1){
				tmpout <- raster(myprecip)
				tmpout <- tmp
			} else {
				tmpout <- tmpout + tmp
			}
		}
		return(writeRaster(tmpout, filename=outfile, overwrite=T, format="GTiff"))
	}
}