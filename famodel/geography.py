""" "Geo" module to hold geography-specific functions for Project class """

import os
import numpy as np
import matplotlib.pyplot as plt

import moorpy as mp
from moorpy.helpers import set_axes_equal
import yaml

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon, LineString

import famodel.seabed_tools as sbt

from pyproj import CRS
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info


# =============================================================================
# COORDINATE REFERENCE SYSTEMS & TRANSFORMATIONS
# =============================================================================

def getLatLongCRS(epsg_code=4326):
    '''Returns a coordinate reference system (CRS) object from the pyproj package of a 'wordly' CRS with units of latitude and longitude
    
    Parameters
    ----------
    epsg_code : int, optional
        The EPSG code for the lat-long CRS to use.
        4326 is the most common (equates to WGS-84)

    Returns
    -------
    latlong_crs : pyproj.CRS object
        Coordinate reference system of the epsg_code provided
    '''

    latlong_crs = CRS.from_epsg(epsg_code)

    return latlong_crs


def getTargetCRS(longitudes, latitudes):
    '''Get a coordinate reference system (CRS) object of a UTM CRS based on an area of interest defined by a series of latitudes and longitudes

    Parameters
    ----------
    longitudes : float, array
        A list of longitude degree coordinates paired to latitude coordinates
    latitudes : float, array
        A list of latitude degree coordinates paired to longitude coordinates

    Returns
    -------
    target_crs : pyproj.CRS object
        Coordinate reference system of the area of interest provided
    '''

    # make sure inputs are in numpy arrays just in case
    longitudes = np.array(longitudes)
    latitudes = np.array(latitudes)

    # set the extents of the area of interest based on lats/longs
    west_long = np.min(longitudes)
    east_long = np.max(longitudes)
    north_lat = np.max(latitudes)
    south_lat = np.min(latitudes)

    # get the list of UTM CRS's that are included in the lat/long extents
    utm_crs_list = query_utm_crs_info(
        datum_name="WGS 84",
        area_of_interest=AreaOfInterest(
            west_lon_degree=west_long,
            east_lon_degree=east_long,
            north_lat_degree=north_lat,
            south_lat_degree=south_lat,
        ),
    )
    # extract the relevant CRS code from the list
    result = CRS.from_epsg(utm_crs_list[0].code)    # get the CRS
    # put it into an EPSG code
    epsg_code = result.srs.split(":")[1]

    # get a pyproj CRS object from that code
    target_crs = CRS.from_epsg(epsg_code)

    return target_crs


def getCustomCRS(long, lat):
    '''Seemingly way too simple of a method to create a pyproj CRS centered around a custom geographical point

    Parameters
    ----------
    long : float
        A longitude degree coordinate
    lat : float
        A latitude degree coordinate

    Returns
    -------
    custom_crs : string (but can be used as input like a pyproj.CRS object)
        Coordinate reference system in meters relative to the input lat/long
    '''
    
    custom_crs = f'+proj=tmerc +lat_0={lat} +lon_0={long} +ellps=WGS84 +units=m +no_defs'

    return custom_crs



def getLeaseCoords(lease_name):

    # read in the BOEM shapefile that contains all Wind Energy Lease Areas (can use other shapefiles for aliquots)
    if 'Humboldt' in lease_name or 'MorroBay' in lease_name:
        lease_areas = gpd.read_file(os.path.dirname(__file__)+'/../geography/Wind_Lease_Outlines_2_2023.shp')
    elif 'GulfofMaine' in lease_name:
        lease_areas = gpd.read_file(os.path.dirname(__file__)+'/../geography/BOEM_Wind_Planning_Area_Outlines_04_29_2024.shp')

    # extract the lease area of interest
    if lease_name=='Humboldt_NE':
        lease_area = lease_areas.loc[lease_areas['LEASE_NUMB']=='OCS-P0561 - Provisional']
    elif lease_name=='Humboldt_SW':
        lease_area = lease_areas.loc[lease_areas['LEASE_NUMB']=='OCS-P0562 - Provisional']
    elif lease_name=='MorroBay_W':
        lease_area = lease_areas.loc[lease_areas['LEASE_NUMB']=='OCS-P0563 - Provisional']
    elif lease_name=='MorroBay_C':
        lease_area = lease_areas.loc[lease_areas['LEASE_NUMB']=='OCS-P0564 - Provisional']
    elif lease_name=='MorroBay_E':
        lease_area = lease_areas.loc[lease_areas['LEASE_NUMB']=='OCS-P0565 - Provisional']
    elif lease_name=='GulfofMaine_ResearchArray':
        lease_area = lease_areas.loc[lease_areas['ADDITIONAL']=='Marine Research Array Requested Lease']
    else:
        raise ValueError(f"The lease area name '{lease_area}' is not supported yet")
    
    # extract the longitude and latitude coordinates of the lease area
    area_longs, area_lats = lease_area.geometry.unary_union.exterior.coords.xy

    # calculate the centroid of the lease area
    centroid = ( lease_area.geometry.centroid.values.x[0], lease_area.geometry.centroid.values.y[0] )

    return area_longs, area_lats, centroid
        




def convertLatLong2Meters(longs, lats, centroid, latlong_crs, target_crs, return_centroid=False):
    '''input longs/lats need to be in EPSG:4326 CRS
    Longs and Lats need to be in pairs, i.e., the first entry to longs and 
    the first entry of lats needs to correspond to a specific point
    
    Parameters
    ----------
    longs : array/list
        array of longitude coordinates (x positions)
    lats : array/list
        array of latitude coordinates (y positions)
    centroid : tuple
        A tuple or list of two values, x and y, of the project 
        reference point or centroid

    Returns
    -------
    xs : array
        x values relative to centroid (from longitudes) [m]
    ys : array
        y values relative to centroid (from latitudes) [m]
    '''

    if len(longs) != len(lats):
        raise ValueError('The list of longs needs to be the same length as the list of lats')
    
    # organize all the long/lat points into a shapely Polygon
    points = [Point(longs[i],lats[i]) for i in range(len(longs))]

    # put list of longs/lats and the centroid into a GeoDataFrame in EPSG:4326
    gdf = gpd.GeoDataFrame({'type':'shape','geometry':points}, crs=latlong_crs)
    gdf_centroid = gpd.GeoDataFrame({'type':'centroid','geometry': [Point(centroid)]}, crs=latlong_crs)
    # combine the geodataframes into one (adding rows to a gdf)
    gdf = pd.concat([gdf, gdf_centroid])

    # convert the GeoDataFrame to UTM coordinates (i.e., meters)
    gdf_utm = gdf.to_crs(target_crs)

    # extract the centroid in x and y [m]
    xcentroid = gdf_utm[gdf_utm['type']=='centroid'].geometry.x.values[0]
    ycentroid = gdf_utm[gdf_utm['type']=='centroid'].geometry.y.values[0]

    # calculate the long/lat distances from the centroid, in units of meters
    xs = np.zeros(len(longs))
    ys = np.zeros(len(lats))
    for i,x in enumerate(gdf_utm[gdf_utm['type']=='shape'].geometry.get_coordinates().x.values):
        xs[i] = np.array(x - xcentroid)
    for i,y in enumerate(gdf_utm[gdf_utm['type']=='shape'].geometry.get_coordinates().y.values):
        ys[i] = np.array(y - ycentroid)
    
    if return_centroid:
        return xs, ys, (xcentroid, ycentroid)
    else:
        return xs, ys


def convertMeters2LatLong(xs, ys, centroid, latlong_crs, target_crs, mesh=False):
    '''Input xs and ys need to be in the target CRS
    Xs and Ys need to be in pairs, i.e. the first entry to xs and the 
    first entry to ys needs to correspond to a specific point
    
    Parameters
    ----------
    xs : array/list
        array of x coordinates
    ys : array/list
        array of y coordinates
    centroid : tuple
        A tuple or list of two values, x and y, of the project 
        reference point or centroid

    Returns
    -------
    longs : array
        longitudes relative to centroid [m]
    lats : array
        latitudes relative to centroid [m]
    '''

    if len(xs) != len(ys):
        raise ValueError('The list of xs needs to be the same length as the list of ys')

    if mesh:
        points = [Point(xs[i,j] + centroid[0], ys[i,j] + centroid[1]) for i in range(len(xs)) for j in range(len(xs[0]))]
    else:
        # organize all the long/lat points into a shapely Polygon
        points = [Point(xs[i] + centroid[0], ys[i] + centroid[1]) for i in range(len(xs))]

    # input the Polygon of longs/lats and the centroid into a GeoDataFrame in EPSG:4326
    gdf = gpd.GeoDataFrame({'type':'shape','geometry':points}, crs=target_crs)
    # convert back to lat/long coordinates
    gdf_wgs = gdf.to_crs(latlong_crs)

    # calculate all the long/lat points distances from the centroid
    longs = np.array(gdf_wgs[gdf_wgs['type']=='shape'].geometry.get_coordinates().x.values)
    lats = np.array(gdf_wgs[gdf_wgs['type']=='shape'].geometry.get_coordinates().y.values)
    
    return longs, lats

# =============================================================================
# BATHYMETRY & GEOSPATIAL DATA PROCESSING
# =============================================================================

def getMapBathymetry(gebcofilename):

    # load the GEBCO bathymetry file
    depths = -np.loadtxt(gebcofilename, skiprows=6)
    # we want the depths matrix to mimic lats/longs on a map, so we need to flip it because the top left value from GEBCO is the first (W) longitude and "first" (S) latitude, when it should be the first (N) latitude
    depths = np.flipud(depths)
    # organize the rest of the variables in the GEBCO file
    nrows = len(depths)
    ncols = len(depths[0])
    newlinedata = []
    with open(gebcofilename) as f:
        lines = f.readlines()
    for i,line in enumerate(lines):
        if i < 6: newlinedata.append(line.split())
        if i==2: xllcorner = float(line.split()[1])
        if i==3: yllcorner = float(line.split()[1])
        if i==4: cellsize = float(line.split()[1])

    # create an array of latitudes and longitudes of the GEBCO bathymetry 
    longs = np.linspace(xllcorner, xllcorner+ncols*cellsize, ncols)
    lats  = np.linspace(yllcorner, yllcorner+nrows*cellsize, nrows)

    return longs, lats, depths, ncols, nrows


def convertBathymetry2Meters(longs, lats, depths, centroid, centroid_utm, 
                             latlong_crs, target_crs, ncols, nrows,
                             xs=[], ys=[]):
    '''
    Convert...
    
    Parameters
    ----------
    ...
    xs : list (optional)
        List of desired x grid values (overrides nrows if provided).
    ys : list (optional)
        List of desired y grid values (overrides nrows if provided).
    '''
    # If converted straight to meters, these longs/lats would not equate to a perfect square due to the distortion of projecting a spherical map into 2D
    # So, what we do is take the four corners of the long/lat rectangle and convert those 4 points to UTM [m] x/y coordinates, relative to a centroid, 
    # which will not line up to a perfect rectangle (will be distorted). We will then make our own arbitrary rectangle in the UTM CRS and 
    # interpolate the depths at each of those new points

    # convert the long/lat coordinates of the four corners of the bathymetry grid to meters, using the target_crs 
    # Default set to [bottom left, bottom right, top right, top left]
    corner_xs, corner_ys = convertLatLong2Meters(
        np.array([longs[0], longs[-1], longs[-1], longs[0]]),
        np.array([lats[0], lats[0], lats[-1], lats[-1]]),
        centroid, latlong_crs, target_crs)

    # now create a new, arbitrary square in the UTM coordinate system based on minimums of corner points (find the rectangle inside the irregular polygon); 
    grid_xs = np.linspace(np.max([corner_xs[0], corner_xs[-1]]), np.min([corner_xs[1], corner_xs[-2]]), ncols)
    grid_ys = np.linspace(np.max([corner_ys[0], corner_ys[1]]), np.min([corner_ys[-2], corner_ys[-1]]), nrows)
    # optional manual choice of x or y grid values (quick add-on for now)
    if len(xs) > 0: grid_xs = np.array(xs)
    if len(ys) > 0: grid_ys = np.array(ys)

    # create a mesh of these new x/y points
    X, Y = np.meshgrid(grid_xs, grid_ys)

    # convert each x/y point in the new bathymetry grid back to lat/long (this will result in another irregular polygon but now in lat/long space)
    grid_mesh_longs_list, grid_mesh_lats_list = convertMeters2LatLong(X, Y, centroid_utm, latlong_crs, target_crs, mesh=True)

    # interpolate the depth at each lat/long pair from the original GEBCO data
    bath_depths = np.zeros([len(grid_ys), len(grid_xs)])
    for i in range(len(grid_mesh_longs_list)):
        depth, _ = sbt.getDepthFromBathymetry(grid_mesh_longs_list[i], grid_mesh_lats_list[i], longs, lats, depths)
        iy = int(i/len(grid_xs))
        ix = i - int(i/len(grid_xs))*len(grid_xs)
        bath_depths[iy,ix] = depth

    bathXs = grid_xs
    bathYs = grid_ys
    bath_depths = bath_depths

    return bathXs, bathYs, bath_depths

        
        
def writeBathymetryFile(moorpy_bathymetry_filename, bathXs, bathYs, bath_depths, soil=False):
    '''Write a MoorDyn/MoorPy-style bathymetry text file based on provided
    x and y grid line values and a 2D array of depth values.'''

    # open the file
    f = open(os.path.join(os.getcwd(), moorpy_bathymetry_filename), 'w')

    f.write('--- MoorPy Bathymetry Input File ---\n')
    f.write(f'nGridX {len(bathXs)}\n')
    f.write(f'nGridY {len(bathYs)}\n')
    f.write(f'      ')

    # x-coordinates
    for ix in range(len(bathXs)):
        f.write(f'{bathXs[ix]:.2f} ')
    f.write('\n')
    #f.write(" ".join(map(str, grid_x)) + "\n")  # old version from seabed_tools

    # y-coordintes
    for iy in range(len(bathYs)):
        f.write(f'{bathYs[iy]:.2f} ')
        for id in range(len(bath_depths[iy])):
            if soil:
                f.write(f'{bath_depths[iy,id]} ')
            else:
                f.write(f'{bath_depths[iy,id]:8.3f} ')
        f.write('\n')
    #for i, y in enumerate(grid_y):         # alternative writing version
        #row = [y] + list(grid_depth[i, :])
        #f.write(" ".join(map(str, row)) + "\n")
    
    # close the file
    f.close()



def getLeaseAndBathymetryInfo(lease_name, bathymetry_file, bath_ncols=100, bath_nrows=100, write_bathymetry=True):

    # initialize the conventional lat/long CRS
    latlong_crs = getLatLongCRS()

    # get lease area coordinates based on BOEM shapefile
    lease_longs, lease_lats, centroid = getLeaseCoords(lease_name)
    
    # get the CRS about the centroid of the lease area of interest
    custom_crs = getCustomCRS(centroid[0], centroid[1])
    
    # convert the lease boundary to meters
    lease_xs, lease_ys, centroid_utm = convertLatLong2Meters(lease_longs, lease_lats, centroid, latlong_crs, custom_crs, return_centroid=True)

    if write_bathymetry:
        # get bathymetry information from a GEBCO file (or other)
        bath_longs, bath_lats, bath_depths, ncols, nrows = getMapBathymetry(bathymetry_file)
        # convert bathymetry to meters
        bath_xs, bath_ys, bath_depths = convertBathymetry2Meters(bath_longs, bath_lats, bath_depths, centroid, centroid_utm, latlong_crs, custom_crs, bath_ncols, bath_nrows)
        # export to MoorPy-readable file
        bathymetry_file = f'bathymetry_{bath_ncols}x{bath_nrows}.txt'
        writeBathymetryFile(bathymetry_file, bath_xs, bath_ys, bath_depths)
    
    ms = mp.System(bathymetry=bathymetry_file)

    info = {}
    info['lease_longs'] = lease_longs
    info['lease_lats'] = lease_lats
    info['lease_centroid'] = centroid
    info['centroid_utm'] = centroid_utm
    info['lease_xs'] = lease_xs
    info['lease_ys'] = lease_ys
    if write_bathymetry:
        info['bath_longs'] = bath_longs
        info['bath_lats'] = bath_lats
        info['bath_xs'] = bath_xs
        info['bath_ys'] = bath_ys
        info['bath_depths'] = bath_depths
    else:
        info['bath_xs'] = ms.bathGrid_Xs
        info['bath_ys'] = ms.bathGrid_Ys
        info['bath_depths'] = ms.bathGrid


    return info



# =============================================================================
# FILE I/O & DATA FORMAT CONVERSION
# =============================================================================

def processGeotiff(filename, lat, lon, outfilename="processGeotiff.txt", **kwargs):
    '''Process a geotiff file containing bathymetry (or other info)
    and convert into a rectangular bathymetry grid in units of m relative to 
    the project reference point.
    
    Parameters
    ----------
    filename : string
        Path and name of geotiff file to load.
    lat : float
        lattitude of reference point to use for array y grid
    long : float
        lattitude of reference point to use for array x grid
    outfilename : string, optional
        If provided, writes a MoorDyn/MoorPy style bathymetry file
    kwargs : dict
        Optional extra arguments that will be relayed to convertBathymetry2Meters,
        can be used to specify desired x and y grid coordinates.

    Returns
    -------
    Xs : array
        x values of grid points [m]
    Ys : array
        y values of grid points [m]
    depths  : 2D array
        water depth grid (positive down) [m]
    '''

    import rasterio
    import rasterio.plot

    tiff = rasterio.open(filename)  # load the geotiff file
    
    #rasterio.plot.show(tiff)  # plot it to see that it works
    
    # note: a CRS is stored with the geotiff, accessible with tiff.crs
    
    # Get lattitude and longitude grid values
    #_, longs = rasterio.transform.xy(tiff.transform, range(tiff.height),0)
    #lats, _  = rasterio.transform.xy(tiff.transform, 0, range(tiff.width-1,-1,-1))
    height, width = tiff.shape
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    longs_mesh, lats_mesh = rasterio.transform.xy(tiff.transform, rows, cols)
    longs = np.array(longs_mesh)[0,:]
    lats = np.flip(np.array(lats_mesh)[:,0])
    # lats data provided from top left corner, i.e., latitudes are descending. It seems that the following interpolation functions (getDepthFromBathymetry)
    # can only work if latitudes start small and increase, meaning that the first latitude entry has to be the bottom left corner
    
    # Depth values in numpy array form
    depths = -tiff.read(1)
    depths = np.flipud(depths)
    # because the interpolation functions require the latitude array to go from small to big (i.e., the bottom left corner), we need to flip the depth matrix to align
    # it will all get sorted out later to what it should be geographically when plotting in MoorPy
    
    
    # extract the coordinate reference systems needed (pyproj CRS objects)
    latlong_crs = tiff.crs
    #target_crs = getTargetCRS(lon, lat)     # should be UTM 10N for Humboldt/California coast
    target_crs = getCustomCRS(lon, lat)     # get custom CRS centered around the lat/long point you want
    
    # get the centroid/reference location in lat/long coordinates
    centroid = (lon, lat)

    # get the centroid/reference location in target_crs coordinates
    #centroid_utm = (lon, lat)
    _, _, centroid_utm = convertLatLong2Meters([centroid[0]], [centroid[1]], centroid, latlong_crs, target_crs, return_centroid=True)
    
    # set the number of rows and columns to use in the MoorPy bathymetry file
    ncols = 100
    nrows = 100
    
    # convert bathymetry to meters
    bath_xs, bath_ys, bath_depths = convertBathymetry2Meters(longs, lats, depths, 
                                                             centroid, centroid_utm, 
                                                             latlong_crs, target_crs, 
                                                             ncols, nrows, **kwargs)
    # export to MoorPy-readable file
    writeBathymetryFile(outfilename, bath_xs, bath_ys, bath_depths)

    return bath_xs, bath_ys, bath_depths



def getSoilType(x, y, centroid, latlong_crs, custom_crs, soil_file):
    """Function to return the name of the soil below a specific x/y coordinate by creating shapely polygons based on the shapefile data.
    It loops through all polygons in the shapefile and if the x/y position is contained in that polygon, it returns the soil of that polygon."""
    
    # create a GeoDataFrame of the shapefile
    soil_gdf = gpd.read_file(soil_file)

    # list of soil names for each geometry in the GeoDataFrame
    soil_names = [soil for soil in soil_gdf['V4_Lith1']]

    # lists of the polygons in the shapefile, and the lat/long coordinates of each polygon
    soil_shapes = [shape for shape in soil_gdf['geometry']]
    soil_longs = [ [x for x,y in shape.exterior.coords] for shape in soil_shapes ]
    soil_lats = [ [y for x,y in shape.exterior.coords] for shape in soil_shapes ]

    # convert the lat/long coordinates to x/y coordinates
    soil_xs = [ convertLatLong2Meters(soil_longs[i], soil_lats[i], centroid, latlong_crs, custom_crs, return_centroid=False)[0] for i in range(len(soil_shapes)) ]
    soil_ys = [ convertLatLong2Meters(soil_longs[i], soil_lats[i], centroid, latlong_crs, custom_crs, return_centroid=False)[1] for i in range(len(soil_shapes)) ]

    # make new coordinates 
    #soil_coords = [ (soil_xs[i], soil_ys[i]) for i in range(len(soil_xs)) ]
    soil_coords = [ [ [soil_xs[i][j], soil_ys[i][j]] for j in range(len(soil_xs[i])) ] for i in range(len(soil_xs)) ]

    # make new shapely polygons based on the new x/y coordinates
    soil_polygons = [ Polygon(soil_coords[i]) for i in range(len(soil_coords)) ]

    # move the large "mud" polygon to the back of the list
    mud_polygon = soil_polygons[0]
    soil_polygons[0] = soil_polygons[-1]
    soil_polygons[-1] = mud_polygon
    # move the mud name to the back of the list
    soil_names[0] = 'hard'
    soil_names[-1] = 'mud'

    # loop through all the new polygons and if the x/y coordinate is within a polygon, return the soil name of that polygon
    for i,polygon in enumerate(soil_polygons):
        if polygon.contains(Point(x,y)):
            soil_type = soil_names[i]
            break

    return soil_type



def getSoilGrid(centroid, latlong_crs, custom_crs, soil_file, nrows=100, ncols=100, xbound=None, ybound=None):
    """Note: can make the outer shapely shape have 'holes' of the inner shapely shapes"""
    
    # create a GeoDataFrame of the shapefile
    soil_gdf = gpd.read_file(soil_file)

    # list of soil names for each geometry in the GeoDataFrame
    soil_names = [soil for soil in soil_gdf['V4_Lith1']]

    # lists of the polygons in the shapefile, and the lat/long coordinates of each polygon
    soil_shapes = [shape for shape in soil_gdf['geometry']]
    soil_longs = [ [x for x,y in shape.exterior.coords] for shape in soil_shapes ]
    soil_lats = [ [y for x,y in shape.exterior.coords] for shape in soil_shapes ]

    # convert the lat/long coordinates to x/y coordinates
    soil_xs = [ convertLatLong2Meters(soil_longs[i], soil_lats[i], centroid, latlong_crs, custom_crs, return_centroid=False)[0] for i in range(len(soil_shapes)) ]
    soil_ys = [ convertLatLong2Meters(soil_longs[i], soil_lats[i], centroid, latlong_crs, custom_crs, return_centroid=False)[1] for i in range(len(soil_shapes)) ]

    # store the length of each soil shape in a separate list (used for reorganization later)
    shape_lengths = [len(soil_xs[i]) for i in range(len(soil_shapes))]
    # sort the shape coordinate lengths from smallest to largest (to eventually loop through from smallest to largest)
    sorted_indices = np.argsort(np.array(shape_lengths))
    # make a new list of soil names based on the length of the shapes
    soil_names = [soil_names[ix] for ix in sorted_indices]

    # make new coordinates 
    soil_coords  = [ [ [soil_xs[i][j], soil_ys[i][j]] for j in range(len(soil_xs[i])) ] for i in range(len(soil_xs)) ]

    # organize coordinates by length of coordinates (so that the smaller shapes are first in the list) - independent of above sorting (but should be equivalent)
    soil_coords.sort(key=len)

    # make new shapely polygons based on the new x/y coordinates
    soil_polygons = [ Polygon(soil_coords[i]) for i in range(len(soil_coords)) ]

    # set up custom grid to store soil data
    #nrows = 100
    #ncols = 100
    if xbound:
        xs = np.linspace(xbound[0], xbound[-1], ncols)
    else:
        xs = np.linspace( np.min([np.min(soil_xs[i]) for i in range(len(soil_shapes))]),  np.max([np.max(soil_xs[i]) for i in range(len(soil_shapes))]),  ncols)
    if ybound:
        ys = np.linspace(ybound[0], ybound[-1], nrows)
    else:
        ys = np.linspace( np.min([np.min(soil_ys[i]) for i in range(len(soil_shapes))]),  np.max([np.max(soil_ys[i]) for i in range(len(soil_shapes))]),  nrows)
    soil_grid = np.zeros([len(ys), len(xs)])

    # for each manmade grid point, loop through all the polygons and determine whether that grid point is within the shape or not
    soil_grid_list = []
    for j in range(len(ys)):
        rowlist = []
        for i in range(len(xs)):
            contains = False
            for p,polygon in enumerate(soil_polygons):
                if polygon.contains(Point(xs[i], ys[j])):
                    rowlist.append(soil_names[p])
                    #soil_grid[j, i] = soil_names[p]
                    contains=True
                    break
            if not contains:
                rowlist.append(0)
                
        soil_grid_list.append(rowlist)
    soil_grid = np.array(soil_grid_list)        # saving to list and then changing to np.array because I couldn't figure out how else to do it with strings

    return xs, ys, soil_grid


        




if __name__ == '__main__':

    """
    # initialize the conventional lat/long CRS
    latlong_crs = getLatLongCRS()

    # get lease area coordinates based on BOEM shapefile
    lease_name = 'Humboldt_SW'
    lease_longs, lease_lats, centroid = getLeaseCoords(lease_name)

    #centroid = (-124.71, 40.93)
    #lease_longs = [centroid[0]]
    #lease_lats = [centroid[1]]

    # based on the lease area, find the target UTM CRS (in m)
    target_crs = getTargetCRS(lease_longs, lease_lats)
    
    custom_crs = getCustomCRS(centroid[0], centroid[1])

    '''
    ##### CRS Tests #####
    # convert lat/long to meters
    longs = [-125.0, -124.75, -124.5]
    lats = [40.8, 40.9, 41.0]
    custom_crs = getCustomCRS(longs[1], lats[1])
    coords = [Point(longs[0], lats[-1]), Point(longs[1], lats[-1]), Point(longs[-1], lats[-1]), Point(longs[0], lats[1]), Point(longs[1], lats[1]), Point(longs[-1], lats[1]), Point(longs[0], lats[0]), Point(longs[1], lats[0]), Point(longs[-1], lats[0])]
    xs, ys = convertLatLong2Meters([point.x for point in coords], [point.y for point in coords], (longs[1], lats[1]), latlong_crs, custom_crs)
    
    # convert meters to lat/long
    x=10000
    turbineList = [Point(-x,x), Point(0,x), Point(x,x), Point(-x,0), Point(0,0), Point(x,0), Point(-x,-x), Point(0,-x), Point(x,-x)]
    longs, lats = convertMeters2LatLong([point.x for point in turbineList], [point.y for point in turbineList], (0,0), latlong_crs, custom_crs)
    ####################
    '''
    
    # convert the lease boundary to meters
    lease_xs, lease_ys, centroid_utm = convertLatLong2Meters(lease_longs, lease_lats, centroid, latlong_crs, custom_crs, return_centroid=True)

    # get bathymetry information from a GEBCO file (or other)
    bath_longs, bath_lats, bath_depths, ncols, nrows = getMapBathymetry('bathymetry/gebco_2023_n41.3196_s40.3857_w-125.2881_e-123.9642.asc')
    # convert bathymetry to meters
    ncols = 500
    nrows = 500
    bath_xs, bath_ys, bath_depths = convertBathymetry2Meters(bath_longs, bath_lats, bath_depths, centroid, centroid_utm, latlong_crs, custom_crs, ncols, nrows)
    # export to MoorPy-readable file
    bathymetryfile = 'bathymetry_large_gebco.txt'
    writeBathymetryFile(bathymetryfile, bath_xs, bath_ys, bath_depths)
    """

    # plot everything
    #plot3d(lease_xs, lease_ys, bathymetryfile, area_on_bath=True, args_bath={'zlim':[-6000, 500], 'cmap': 'gist_earth'})


    # run example
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    lease_name = 'GulfofMaine_ResearchArray'
    gebco_file = __location__+'\\..\\geography\\gebco_2024_n44.1458_s41.4761_w-70.9497_e-66.2146.asc'
    info = getLeaseAndBathymetryInfo(lease_name, bathymetry_file)


    plt.show()

    a = 2






