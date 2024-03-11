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




import famodel.seabed.seabed_tools as sbt




from pyproj import CRS
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info



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
    lease_areas = gpd.read_file('shapefiles/Wind_Lease_Outlines_2_2023.shp')

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

        
        
def writeBathymetryFile(moorpy_bathymetry_filename, bathXs, bathYs, bath_depths):
    '''Write a MoorDyn/MoorPy-style bathymetry text file based on provided
    x and y grid line values and a 2D array of depth values.'''

    f = open(os.path.join(os.getcwd(), moorpy_bathymetry_filename), 'w')
    f.write('--- MoorPy Bathymetry Input File ---\n')
    f.write(f'nGridX {len(bathXs)}\n')
    f.write(f'nGridY {len(bathYs)}\n')
    f.write(f'      ')
    for ix in range(len(bathXs)):
        f.write(f'{bathXs[ix]:.2f} ')
    f.write('\n')
    for iy in range(len(bathYs)):
        f.write(f'{bathYs[iy]:.2f} ')
        for id in range(len(bath_depths[iy])):
            f.write(f'{bath_depths[iy,id]:8.3f} ')
        f.write('\n')
    f.close()






        
def plot3d(lease_xs, lease_ys, bathymetryfilename, area_on_bath=False, args_bath={}):
    '''Plot aspects of the Project object in matplotlib in 3D'''

    # organize the bathymetry arguments
    if len(args_bath)==0:
        args_bath = {'zlim':[-3200,500], 'cmap':'gist_earth'}
 
    fig = plt.figure(figsize=(6,4))
    ax = plt.axes(projection='3d')

    # plot the lease area in a red color, if desired
    ax.plot(lease_xs, lease_ys, np.zeros(len(lease_xs)), color='r', zorder=100)
    
    # plot the bathymetry in matplotlib using a plot_surface

    # !!!! include option to plot entire bathymetry file or not

    if isinstance(bathymetryfilename, str):
        bathGrid_Xs, bathGrid_Ys, bathGrid = sbt.readBathymetryFile(bathymetryfilename)         # parse through the MoorDyn/MoorPy-formatted bathymetry file
        X, Y = np.meshgrid(bathGrid_Xs, bathGrid_Ys)                                    # create a 2D mesh of the x and y values
        bath = ax.plot_surface(X, Y, -bathGrid, rstride=1, cstride=1,
                               vmin=args_bath['zlim'][0], vmax=args_bath['zlim'][1], 
                               cmap=args_bath['cmap'])
    
    '''
    # plot the project boundary
    if boundary:
        ax.plot(self.boundaryXs, self.boundaryYs, np.zeros(len(self.boundaryXs)), color='b', zorder=100, alpha=0.5)
    '''
        
    # plot the projection of the lease area bounds on the seabed, if desired
    if area_on_bath:
        lease_zs = projectAlongSeabed(lease_xs, lease_ys, bathGrid_Xs, bathGrid_Ys, bathGrid)
        ax.plot(lease_xs, lease_ys, -lease_zs, color='k', zorder=10, alpha=0.5)


    set_axes_equal(ax)
    ax.axis('off')



def projectAlongSeabed(x, y, bathXs, bathYs, bath_depths):
    '''Project a set of x-y coordinates along a seabed surface (grid),
    returning the corresponding z coordinates.'''
    
    if len(x) == len(y):
        n = len(x)        
        z = np.zeros(n)   # z coordinate of each point [m]
        a = np.zeros(n)   # could also do slope (dz/dh)
        for i in range(n):
            z[i], nvec = sbt.getDepthFromBathymetry(x[i], y[i], bathXs, bathYs, bath_depths)
    
    else:
        z = np.zeros([len(y), len(x)])
        for i in range(len(y)):
            for j in range(len(x)):
                z[i,j], nvec = sbt.getDepthFromBathymetry(x[j], y[i], bathXs_mesh, bathYs_mesh, bath_depths)
        
    return z








        


    """
        if self.lat0 != 0 and self.lon0 != 0:
                # set the centroid of the project and save in a GeoDataFrame
                self.centroid = (self.lon0, self.lat0)
                self.gdf = gpd.GeoDataFrame({'type':'centroid', 'geometry': [Point(self.centroid)]}, crs=self.latlong_crs)
                # set the target coordinate reference system (CRS) that will switch between regular lat/long system and "meters from centroid" system, based on centroid location
                information
                            # extract the numeric code
                self.  # save the target CRS (UTM 10N = 32610)
        

        gdf_leases = gpd.GeoDataFrame({'type': 'lease_area', 'geometry': lease_area.geometry}, crs=getLatLongCRS() )

        gdf = add2geodataframe(gdf0, gdf_leases)

    # convert lat/long boundary/lease area points to meters away from the centroid
            self.lease_xs, self.lease_ys = self.convertLatLong2Meters(self.area_longs, self.area_lats, self.centroid)


            if which_centroid=='lease_area':
                # make a blank copy of the gdf to switch to the target CRS to get the accurate centroid
                gdf_utm = self.gdf.copy().to_crs(self.target_crs)
                centroid_utm = (gdf_utm.geometry.centroid.values.x[0], gdf_utm.geometry.centroid.values.y[0])
                gdf_centroid = gpd.GeoDataFrame({'type':'centroid', 'geometry': [Point(centroid_utm)]}, crs=self.target_crs)
            else:
                # assuming the input centroid is in a long/lat pair
                gdf_centroid = gpd.GeoDataFrame({'type':'centroid', 'geometry': [Point(which_centroid)]}, crs=self.latlong_crs)
                gdf_centroid.to_crs(self.target_crs)
            
            self.centroid_utm = (gdf_centroid.geometry.values.x[0], gdf_centroid.geometry.values.y[0])
            self.centroid = (gdf_centroid.to_crs(self.latlong_crs).geometry.values.x[0], gdf_centroid.to_crs(self.latlong_crs).geometry.values.y[0]) # assume centroid is focal point of Project

            def initialize_geodataframe(crs, columns=['type','geometry']):

        gdf = gpd.GeoDataFrame(columns=columns, geometry='initialization', crs=crs)

        return gdf

    def add2geodataframe(gdf_to_add_to, gdf_to_add):
        
        # check to make sure they have the same columns and CRS

        gdf_new = pd.concat([gdf_to_add_to, gdf_to_add])

        return gdf_new





        # TODO

        # reference entire CA bathymetry file (maybe)
        # plot2d method and a plotGDF method (with bathymetry in the geodataframe using contours)

        # add the coastline
        #xcoast, ycoast = sbt.getCoast(self.Xs, self.Ys, self.depths)
        #ax.plot(xcoast, ycoast, np.zeros(len(self.Ys)), color='k', zorder=100)

        # need to fix up bounds
        #xbmin, xbmax = sbt.getPlotBounds(self.longs_bath, self.centroid, long=True)
        #ybmin, ybmax = sbt.getPlotBounds(self.lats_bath, self.centroid, long=False)

        plt.show()








    # METHODS USED SPECIFICALLY FOR DEEPFARM LCOE ANALYSIS
    
    def addMap2GDF(self, filename='', states=None):
        '''function to include a shapefile of a provided map'''

        # read in the provided filename to add to the geodataframe
        usa = gpd.read_file(filename)
        # internal list of states, in order, to go with the default U.S. states shapefile
        statenamelist = ['Maryland','Iowa','Delaware','Ohio','Pennsylvania','Nebraska','Washington','Puerto Rico','Alabama','Arkansas','New Mexico',            # 0-10
                         'Texas','California','Kentucky','Georgia','Wisconsin','Oregon','Missouri','Virginia','Tennessee','Louisiana','New York',               # 11-21
                         'Michigan','Idaho','Florida','Alaska','Illinois','Montana','Minnesota','Indiana','Massachusetts','Kansas','Nevada','Vermont',          # 22-33
                         'Connecticut','New Jersey','Washington D.C.','North Carolina','Utah','North Dakota','South Carolina','Mississippi','Colorado',         # 34-42
                         'South Dakota','Oklahoma','Wyoming','West Virginia','Maine','Hawaii','New Hampshire','Arizona','Rhode Island']                         # 43-51
        # insert names of the states into the new gdf
        usa.insert(0, 'type', statenamelist)
        # set the CRS of the USA pdf to the right CRS
        usa.set_crs(crs="EPSG:4326", inplace=True)
        self.usa = usa
        
        for state in states:
            state_gs = usa.loc[usa['type']==state]
            self.gdf = pd.concat([self.gdf, state_gs])
    


    def setFarmLayout(self, style='grid', nrows=10, ncols=10, turbine_spacing=2000, nOSS=2):
        
        if style=='grid':
            # for now, this is very custom code specific to the DeepFarm project
            farmxspacing = (nrows-1)*turbine_spacing
            farmyspacing = (ncols-1)*turbine_spacing

            turbine_distances_from_centroid = []
            oss_distances_from_centroid = []
            for j in reversed(range(ncols)):
                for i in range(nrows):
                    xpos = -(farmxspacing/2)+(i*turbine_spacing)
                    ypos = -(farmyspacing/2)+(j*turbine_spacing)
                    turbine_distances_from_centroid.append((xpos, ypos))
            
            # add positions of two offshore substations (OSSs)
            oss_distances_from_centroid.append((11000.0, 5000.0))
            oss_distances_from_centroid.append((11000.0, -5000.0))
        
        if style=='shared':
            turbine_xspacing = np.sqrt(2000**2-1000**2)
            turbine_yspacing = 2000
            farmxspacing = turbine_xspacing*9
            farmyspacing = turbine_yspacing*9

            turbine_distances_from_centroid = []
            oss_distances_from_centroid = []
            for j in reversed(range(ncols)):
                for i in range(nrows):
                    xpos = -(farmxspacing/2)+(i*turbine_xspacing)
                    ypos = -(farmyspacing/2)+(j*turbine_yspacing) - 1000*np.sin(np.radians(30)) + 1000*(i%2)
                    turbine_distances_from_centroid.append((xpos, ypos))
            
            # add positions of two offshore substations (OSSs)
            oss_distances_from_centroid.append((5.5*turbine_xspacing, 2.0*turbine_yspacing+1000*np.sin(np.radians(30))))
            oss_distances_from_centroid.append((5.5*turbine_xspacing, -2.5*turbine_yspacing-1000*np.sin(np.radians(30))))
        
        if style=='small-shared':
            turbine_xspacing = np.sqrt(2000**2-1000**2)
            turbine_yspacing = 2000
            farmxspacing = turbine_xspacing*1
            farmyspacing = turbine_yspacing*2

            turbine_distances_from_centroid = []
            oss_distances_from_centroid = []
            for j in reversed(range(3)):
                for i in range(2):
                    xpos = -(farmxspacing/2)+(i*turbine_xspacing)
                    ypos = -(farmyspacing/2)+(j*turbine_yspacing) - 1000*np.sin(np.radians(30)) + 1000*(i%2)
                    turbine_distances_from_centroid.append((xpos, ypos))
            
            # add positions of two offshore substations (OSSs)
            oss_distances_from_centroid.append((-0.5*turbine_xspacing, 2.0*turbine_yspacing-1000*np.sin(np.radians(30))))
        

        # create a copy of the global gdf and transform it into the easting/northing coordinate reference system
        gdf_utm = self.gdf.copy().to_crs(self.target_crs)
        xcentroid = gdf_utm.loc[gdf_utm['type']=='centroid'].centroid.x[0]
        ycentroid = gdf_utm.loc[gdf_utm['type']=='centroid'].centroid.y[0]

        # create shapely Point objects of the turbine positions relative to the centroid, in the UTM CRS
        turbine_geoms = []
        for i,(x,y) in enumerate(turbine_distances_from_centroid):
            turbine_geoms.append( Point(xcentroid + x, ycentroid + y) )
        
        oss_geoms = []
        for i,(x,y) in enumerate(oss_distances_from_centroid):
            oss_geoms.append( Point(xcentroid + x, ycentroid + y) )
        
        # make a new gdf to put the turbine data together
        turbine_gdf = gpd.GeoDataFrame({'type': ['turbine']*len(turbine_geoms), 'geometry': turbine_geoms}, crs=self.target_crs)
        # make a new gdf to put the substation data together
        oss_gdf = gpd.GeoDataFrame({'type': 'substation', 'geometry': oss_geoms}, crs=self.target_crs)
        # merge these two geodataframes together into one (best way I can find to "add" rows to a dataframe; ignoring index makes all indices a different number)
        turbine_gdf = pd.concat([turbine_gdf, oss_gdf], ignore_index=True)

        # convert the turbine/oss coordinates back to regular latitude/longitude (EPSG: 4326)
        turbine_gdf = turbine_gdf.to_crs('EPSG:4326')

        # add the turbine gdf to the global gdf
        self.gdf = pd.concat([self.gdf, turbine_gdf], ignore_index=True)

        # add local variables in this method to the turbine_gdf to be used later (but don't need for the global gdf)
        turbine_gdf.insert(2, 'easting_northing_geometry', turbine_geoms + oss_geoms)
        turbine_gdf.insert(3, 'meters_from_centroid', turbine_distances_from_centroid + oss_distances_from_centroid)

        # save this new turbine_gdf for future use
        self.turbine_gdf = turbine_gdf

        # make a layout CSV (used for WHaLE/WAVES)
        self.makeLayoutCSV()

    

    def makeLayoutCSV(self, filename='layout_test.csv'):

        turbine_longs = [point.coords[0][0] for point in self.turbine_gdf.geometry]
        turbine_lats = [point.coords[0][1] for point in self.turbine_gdf.geometry]

        self.turbine_gdf.insert(2, 'longitude', turbine_longs)
        self.turbine_gdf.insert(3, 'latitude', turbine_lats)

        turbine_eastings = [point.coords[0][0] for point in self.turbine_gdf['easting_northing_geometry']]
        turbine_northings = [point.coords[0][1] for point in self.turbine_gdf['easting_northing_geometry']]

        #self.turbine_gdf.insert(5, 'easting', turbine_eastings)
        #self.turbine_gdf.insert(6, 'northing', turbine_northings)

        turbine_x_from_centroid = [point[0] for point in self.turbine_gdf['meters_from_centroid']]
        turbine_y_from_centroid = [point[1] for point in self.turbine_gdf['meters_from_centroid']]

        self.turbine_gdf.insert(5, 'easting', turbine_x_from_centroid)
        self.turbine_gdf.insert(6, 'northing', turbine_y_from_centroid)

        self.turbine_gdf.insert(8, 'floris_x', turbine_x_from_centroid)
        self.turbine_gdf.insert(9, 'floris_y', turbine_y_from_centroid)

        columns = ['type', 'longitude', 'latitude', 'easting', 'northing', 'floris_x', 'floris_y']
        df = pd.DataFrame(self.turbine_gdf)
        df.to_csv(filename, columns=columns)


    def plotGDF(self, kwargs):
        '''2D map-like plot'''
        
        if 'centroid' in kwargs:
            centroid_settings = kwargs['centroid']
            if 'label' in centroid_settings:
                centroid_label = 'centroid'
        if 'map' in kwargs:
            map_settings = kwargs['map']
        if 'farm' in kwargs:
            farm_settings = kwargs['farm']

        fig, ax = plt.subplots(1,1)

        if 'centroid' in kwargs:
            self.gdf.loc[self.gdf['type']=='centroid'].plot(ax=ax, color=centroid_settings['color'], label=centroid_label)

        if 'boundary' in kwargs:
            map_boundary = self.gdf.loc[self.gdf['type']=='California'].boundary
            map_boundary.plot(ax=ax, color=map_settings['color'])
        
        if 'farm' in kwargs:
            self.gdf.loc[self.gdf['type']=='turbine'].plot(ax=ax, color=farm_settings['turbine']['color'], label='turbine')
            self.gdf.loc[self.gdf['type']=='substation'].plot(ax=ax, color=farm_settings['oss']['color'], label='substation')
        
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.legend()

        ax.set_xlim([-124.875, -124.55])
        ax.set_ylim([40.025, 40.25])

        fig.tight_layout()
        
        # Some GeoPandas Help
        # to plot just one entry of a geoseries: gdf.loc[[0],'geometry'].plot()
        # to get the columns of a gdf: gdf.columns
        # merging gdf's
        # adding columns to gdf's

        return fig, ax



    def addPoints(self, ax, pointlist=[], kwargs={}):

        point_settings = kwargs['pointlist']

        points = gpd.GeoDataFrame({'type':['nrel_channel','nrel_humboldt','nrel_crescent_city','hawaii'],
                                   'geometry': pointlist}, crs='EPSG:4326')
        
        points.plot(ax=ax, color=point_settings['color'], marker=point_settings['marker'], label=point_settings['label'])
    

    def addState(self, ax, states=[], kwargs={}):
        
        for state in states:
            state_settings=kwargs[state]

            state_geom = self.usa.loc[self.usa['type']==state]
            if 'boundary' in state_settings:
                state_geom = state_geom.boundary

            newstate = gpd.GeoDataFrame({'type':state, 'geometry':state_geom}, crs='EPSG:4326')

            newstate.plot(ax=ax, color=state_settings['color'])




    

"""



if __name__ == '__main__':

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

    # plot everything
    plot3d(lease_xs, lease_ys, bathymetryfile, area_on_bath=True, args_bath={'zlim':[-6000, 500], 'cmap': 'gist_earth'})




    # store everything in a GDF (for other plotting)

    # plot everything in a matplotlib



    plt.show()

    a = 2






