# Make different size sample boxes for precipitation analysis
# r <- c(0.5,1,2,4)
# temp is a raster object of the standard domain
makeBoxes <- function(r, temp){
	
	print("Running makeBoxes ...")
	
	require(raster)
	require(sp)
	
	# 192 1 degree grid cells in the domain 16 X 12
	# 240 1 degree grid cells in the domain 20 X 12
	b <- 240 / r^2 

	boxes <- raster(matrix(1:b, nrow=12/r, ncol=20/r, byrow=T), xmn=354, xmx=374, ymn=-8, ymx=4)
	boxes.poly <- rasterToPolygons(boxes)
	boxes.r <- resample(boxes, temp, method="ngb")
	boxesa<- raster(matrix(1:b, nrow=12/r, ncol=20/r, byrow=T), xmn=-6, xmx=14, ymn=-8, ymx=4)
	boxesa.poly <- rasterToPolygons(boxesa)
	boxesa.r <- resample(boxesa, shift(temp, x=-360), method="ngb")
	
	return(list(boxes.r, boxesa.r, boxes.poly, boxesa.poly))
}

makeBoxes2 <- function(xmin, ymin, xmax, ymax, xdim, ydim, temp){
	
	print("Running makeBoxes2 ...")
	
	require(raster)
	require(sp)
	
# 	b <- 240 / r^2 
	nr <- ceiling((ymax - ymin) / ydim)
	nc <- ceiling((xmax - xmin) / xdim)

	boxes <- raster(matrix(1:(nr*nc), nrow=nr, ncol=nc, byrow=T), xmn=xmin, xmx=xmin+(nc*xdim), ymn=ymin, ymx=ymin+(nr*ydim))
	boxes.poly <- rasterToPolygons(boxes)
	boxes.r <- resample(boxes, temp, method="ngb")
	if (xmax(boxes.r) > 360){
	    boxesa <- shift(boxes, x=-360)
	    boxesa.poly <- rasterToPolygons(boxesa)
	    boxesa.r <- resample(boxesa, shift(temp, x=-360), method="ngb")
	    return(list(boxes.r, boxesa.r, boxes.poly, boxesa.poly))
	}
	
	return(list(boxes.r, boxes.r, boxes.poly, boxes.poly))
}