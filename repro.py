__author__ = 'ajh235'

import iris
import iris.analysis.cartography as carto
import iris.experimental.regrid
import iris.coords
import numpy as np
import iris.coord_systems as cs
import sys
import cartopy.crs as ccrs

def main(infile):

    thiscube = iris.load_cube(infile)

    xmin_rp = min(thiscube.coord('grid_longitude').points) - (0.036/2)
    xmax_rp = max(thiscube.coord('grid_longitude').points) + (0.036/2)
    ymin_rp = min(thiscube.coord('grid_latitude').points) - (0.036/2)
    ymax_rp = max(thiscube.coord('grid_latitude').points) + (0.036/2)

    pole_lat = thiscube.coord("grid_latitude").coord_system.grid_north_pole_latitude
    pole_lon = thiscube.coord("grid_latitude").coord_system.grid_north_pole_longitude

    bbox_ll = carto.unrotate_pole(np.array([xmin_rp,xmax_rp]), np.array([ymin_rp, ymax_rp]), pole_lon, pole_lat)

    latpts = np.arange(bbox_ll[1][0], bbox_ll[1][1], float(0.036))
    lonpts = np.arange(bbox_ll[0][0], bbox_ll[0][1], float(0.036))

    latitude = iris.coords.DimCoord(latpts, standard_name='latitude', units='degrees')
    longitude = iris.coords.DimCoord(lonpts, standard_name='longitude', units='degrees')

    ll = cs.GeogCS(semi_major_axis=6378137, semi_minor_axis=6356752.314245)
    llcube = iris.cube.Cube(np.zeros((latpts.size, lonpts.size), np.float32), dim_coords_and_dims=[(latitude, 0), (longitude, 1)])
    llcube.coord_system = ll

    llcube.coord(axis='x').coord_system = ll
    llcube.coord(axis='y').coord_system = ll

    thiscube_ll = iris.experimental.regrid.regrid_bilinear_rectilinear_src_and_grid(thiscube, llcube)

    outfile = infile.replace('.pp','_ll.nc')

    iris.save(thiscube_ll, outfile)

    return(thiscube_ll)

def maintest(infile, nl):
    print nl
    thiscube = iris.load(infile)[nl]
    print thiscube

def rp2ll(x, y, rpcube):
    ll = ccrs.Geodetic()
    rp = rpcube.coord('grid_latitude').coord_system.as_cartopy_crs()
    outcoords = ll.transform_point(x, y, rp)
    return(outcoords)

if __name__ == '__main__':
    # infile = '/Users/ajh235/Work/DataLocal/ModelData/WAFR/ancils/km4/qrparm.veg.frac_4km.std.pp'
    infile = sys.argv[1]
    # nl = sys.argv[2]
    print infile
    main(infile)
