"""Project class for FAModel, containing information and key methods for
the site information and design information that make up a project."""

import os
import numpy as np
import matplotlib.pyplot as plt

import moorpy as mp
from moorpy.helpers import set_axes_equal
import yaml
# import raft

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon, LineString
from pyproj import CRS
from pyproj.aoi import AreaOfInterest
from pyproj.database import query_utm_crs_info

from anchors.anchor_capacity import anchorCapacity
import seabed.seabed_tools as sbt


class Project():
    '''
    The overall object that defines a floating array for analysis and design 
    purposes. Its main model component is the RAFT model but it also includes
    more general information such as seabed info, metocean info, and...
    
    Ideally, this class can function even if RAFT and MoorPy are not loaded,
    at least for some high-level site processing functions.
    
    '''
    
    def __init__(self, lon=0, lat=0, file=None, depth=100):
        '''Initialize a Project. If input data is not provided, it will
        be empty and can be filled in later.
        
        Parameters
        ----------
        file : string or dict, optional
            Name of YAML file, or a python dictionary, containing input
            information describing the project, following the ontology.
        '''
        
        
        # ----- design information -----
        
        # higher-level design data structures
        self.nPtfm  = 0  # number of floating platforms
        self.nAnch = 0   # number of anchors        
        self.coords = np.zeros([self.nPtfm+self.nAnch, 2]) # x-y coordinate table of platforms and anchors
        
        # more detailed design data structures for submodels
        self.array = None  # RAFT Array
        self.cables = None  # CableSystem
        
        
        # ----- site information -----
        self.latlong_crs = CRS.from_epsg(4326)

        self.lat0  = lat  # lattitude of site reference point [deg]
        self.lon0  = lon  # longitude of site reference point [deg]

        if self.lat0 != 0 and self.lon0 != 0:
            # set the centroid of the project and save in a GeoDataFrame
            self.centroid = (self.lon0, self.lat0)
            self.gdf = gpd.GeoDataFrame({'type':'centroid', 'geometry': [Point(self.centroid)]}, crs=self.latlong_crs)
            # set the target coordinate reference system (CRS) that will switch between regular lat/long system and "meters from centroid" system, based on centroid location
            utm_crs_list = query_utm_crs_info(
                datum_name="WGS 84",
                area_of_interest=AreaOfInterest(
                    west_lon_degree=self.centroid[0],
                    east_lon_degree=self.centroid[0],
                    north_lat_degree=self.centroid[1],
                    south_lat_degree=self.centroid[1],
                ),
            )
            result = CRS.from_epsg(utm_crs_list[0].code)    # get the CRS information
            epsg_code = result.srs.split(":")[1]            # extract the numeric code
            self.target_crs = CRS.from_epsg(epsg_code)  # save the target CRS (UTM 10N = 32610)


        # project boundaries
        self.boundary_Xs = np.zeros(0)
        self.boundary_Ys = np.zeros(0)
        
        self.grid_x      = np.array([0])
        self.grid_y      = np.array([0])
        self.grid_depth  = np.array([[depth]])  # depth at each grid point
        
        self.seabed_type = 'clay'  # switch of which soil property set to use ('clay', 'sand', or 'rock')
        
        # soil parameters at each grid point
        self.soil_mode  = 0                # soil/anchor model level to use (0: none; 1: simple categories; 2: quantitative)
        self.soil_class = [["none"]]       # soil classification name ('clay', 'sand', or 'rock' with optional modifiers)
        self.soil_gamma = np.zeros((1,1))  # soil effective unit weight [kPa] (all soils)
        self.soil_Su0   = np.zeros((1,1))  # undrained shear strength at mudline [kPa] (clay soils)
        self.soil_K     = np.zeros((1,1))  # undrained shear strength gradient [kPa/m] (clay soils)
        self.soil_alpha = np.zeros((1,1))  # soil skin friction coefficient [-] (clay soils)
        self.soil_phi   = np.zeros((1,1))  # angle of internal friction [deg] (sand soils)

        # ----- if an input file has been passed, load it -----
        if file:
            self.load(file)
    


    def load(self, project_yaml):
        '''
        Load a full set of project information from a dictionary or 
        YAML file. This calls other methods for each part of it.
        
        Parameters
        ----------
        input : dict or filename
            Dictionary or YAML filename containing project info.
        '''
        
        # standard function to load dict if input is yaml
        with open(project_yaml) as file:
            project = yaml.load(file, Loader=yaml.FullLoader)
        
        # look for site section
        # call load site method
        self.loadSite(project['site'])
        
        # look for design section
        # call load design method
        self.loadDesign(project)
    

    # ----- Design loading/processing methods -----
    
    def loadDesign(self, d):
        '''Load design information from a dictionary or YAML file
        (specified by input). This should be the design portion of
        the floating wind array ontology.'''
        
        # standard function to load dict if input is yaml
        #d = 
        
        # ===== load FAM-specific model parts =====
        
        # cable types
        
        # dynamic cable basic properties (details are later via MoorPy)
        
        # ----- table of cables -----
        if 'array_cables' in d:
        
            cableInfo = [dict(zip( d['array_cables']['keys'], row))
                         for row in d['array_cables']['data']]
            
            for ci in cableInfo:
                ...
                
                self.cables.addCable(...)
        
        # ----- cables one-by-one -----
        if 'cables' in d:
        
            for ci in d['cables']:
                ...
                
                self.cables.addCable(...)
        
        
        # ===== load RAFT model parts =====
        

    # ----- Site conditions processing functions -----

    def loadSite(self, site):
        '''Load site information from a dictionary or YAML file
        (specified by input). This should be the site portion of
        the floating wind array ontology.'''
        # standard function to load dict if input is yaml
        
        # load general information
        self.depth = getFromDict(site['general'], 'depth', default=100)
        self.rho_water = getFromDict(site['general'], 'rho_water', default=1025.0)
        self.rho_air = getFromDict(site['general'], 'rho_air', default=1.225)
        self.mu_air = getFromDict(site['general'], 'mu_air', default=1.81e-5)

        # load geographical information, if provided
        self.centroid = getFromDict(site['location'], 'centroid', default='lease_area')             # >>>>>>>> might need another if statement or two to sort out the centroid
        self.lease_area_name = getFromDict(site['location'], 'lease_area', dtype=str, default=None)
        if self.lease_area_name is not None:
            self.loadBoundary(self.lease_area_name, which_centroid=self.centroid)
        
        # load bathymetry information, if provided
        self.bathymetry_gebco = getFromDict(site['bathymetry'], 'gebco_file', dtype=str, default='')
        self.bathymetry_moorpy = getFromDict(site['bathymetry'], 'moorpy_file', dtype=str, default='')
        if len(self.bathymetry_gebco) > 0:
            self.loadBathymetry(self.bathymetry_gebco, self.bathymetry_moorpy)

        # load project boundary/grid information
        self.boundary_type = getFromDict(site['boundary'], 'type', dtype=str, default='default')

        if self.boundary_type=='bathymetry':
            if 'gebco_file' not in site['bathymetry']:      # check to make sure you've run loadBathymetry
                raise ValueError("Need to include a bathymetry input file")

            # set the x and y coordinates to be used for a grid based on the extent of the bathymetry file
            bathX_min = np.max(self.bathXs_mesh[:,0])
            bathX_max = np.min(self.bathXs_mesh[:,-1])
            bathY_min = np.max(self.bathYs_mesh[-1,:])
            bathY_max = np.min(self.bathYs_mesh[0,:])
            # initialize the discretization of the grid
            dbathX = np.abs(self.bathXs_mesh[0,1] - self.bathXs_mesh[0,0])*10
            dbathY = np.abs(self.bathYs_mesh[1,0] - self.bathYs_mesh[0,0])*10
            # create new grid x and y arrays however you want (default based on bathymetry)
            xs = np.arange(bathX_min, bathX_max, dbathX)
            ys = np.arange(bathY_min, bathY_max, dbathY)

            #self.extent = [self.bathXs_mesh[:,0], self.bathXs_mesh[:,-1], self.bathYs_mesh[0,:], self.bathYs_mesh[-1,:]]
            print("WARNING: The process to set the grid depth values for a large bathymetry grid can take a long time. Make sure you use a small grid or discretization")
        
        elif self.boundary_type=='lease_area':
            if self.lease_area_name==None:
                raise ValueError('Need to provide the name of a valid lease area to use this boundary type')
            #self.extent = [np.min(self.lease_xs), np.max(self.lease_xs), np.max(self.lease_ys), np.min(self.lease_ys)]
            xs = self.lease_xs
            ys = self.lease_ys
        
        elif self.boundary_type=='extent':      # in long/lat units
            self.extent = getFromDict(site['boundary'], 'data')
        
        elif self.boundary_type=='default':     # in long/lat units
            self.extent = [-1000, 1000, 1000, -1000]
        
        else:
            raise ValueError("Not a valid boundary name")
        
        # and set the project boundary/grid based on the loaded information
        self.setGrid(xs, ys)


        
        # load seabed portions
        
        # load lease area portions
        
        # load metocean portions




    def setGrid(self, xs, ys):
        '''
        Set up the rectangular grid over which site or seabed
        data will be saved and worked with. Directions x and y are 
        generally assumed to be aligned with the East and North 
        directions, respectively, at the array reference point.
        
        Parameters
        ----------        
        xs : float array
            x coordinates relative to array reference point [m]
        ys : float array
            y coordinates relative to array reference point [m]
        '''
        
        xs = np.array(xs)
        ys = np.array(ys)

        dx = np.abs(np.linalg.norm([xs[1], ys[1]]) - np.linalg.norm([xs[0], ys[0]]))

        self.grid_x = np.arange(np.min(xs), np.max(xs)+dx, dx)
        self.grid_y = np.arange(np.min(ys), np.max(ys)+dx, dx)

        self.boundaryXs = np.hstack([self.grid_x, np.ones(len(self.grid_y))*self.grid_x[-1], np.flip(self.grid_x), np.ones(len(self.grid_y))*self.grid_x[0]])
        self.boundaryYs = np.hstack([np.ones(len(self.grid_x))*self.grid_y[0], self.grid_y, np.ones(len(self.grid_x))*self.grid_y[-1], np.flip(self.grid_y)])

        # create a new depth matrix that uses interpolated values from the bathymetry data (be careful: this process can take a long time)
        self.grid_depths = np.zeros([len(self.grid_y), len(self.grid_x)])
        for i in range(len(self.grid_y)):
            for j in range(len(self.grid_x)):
                print(i, j, len(self.grid_x), len(self.grid_y))
                self.grid_depths[i,j] = sbt.getDepthFromBathymetryMesh(self.grid_x[j], self.grid_y[i], self.bathXs_mesh, self.bathYs_mesh, self.bath_depths)
        
        
        #TODO: add check for existing seabed data. If present, convert or raise warning <<<
    



    def loadBoundary(self, lease_name, which_centroid='lease_area'):
        
        # read in the BOEM Shapefile that contains all Wind Energy Lease Areas
        lease_areas = gpd.read_file('shapefiles/Wind_Lease_Outlines_2_2023.shp')    # can use the other shapefile for aliquots

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

        # create a GeoDataFrame of that lease area
        self.gdf = gpd.GeoDataFrame({'type': 'lease_area', 'geometry': lease_area.geometry}, crs=self.latlong_crs)      # gdf to only be used with 2D plotting/visualization

        # set the target coordinate reference system (CRS) that will switch between regular lat/long system and "meters from centroid" system, based on lease area location
        utm_crs_list = query_utm_crs_info(
            datum_name="WGS 84",
            area_of_interest=AreaOfInterest(
                west_lon_degree=self.gdf.geometry.total_bounds[0],
                east_lon_degree=self.gdf.geometry.total_bounds[2],
                north_lat_degree=self.gdf.geometry.total_bounds[3],
                south_lat_degree=self.gdf.geometry.total_bounds[1],
            ),
        )
        result = CRS.from_epsg(utm_crs_list[0].code)    # get the CRS information
        epsg_code = result.srs.split(":")[1]            # extract the numeric code
        self.target_crs = CRS.from_epsg(epsg_code)  # save the target CRS (UTM 10N = 32610)

        if which_centroid=='lease_area':
            # make a blank copy of the gdf to switch to the target CRS to get the accurate centroid
            gdf_utm = self.gdf.copy().to_crs(self.target_crs)
            centroid_utm = (gdf_utm.geometry.centroid.values.x[0], gdf_utm.geometry.centroid.values.y[0])
            gdf_centroid = gpd.GeoDataFrame({'type':'centroid', 'geometry': [Point(centroid_utm)]}, crs=self.target_crs)
            self.centroid = (gdf_centroid.to_crs(self.latlong_crs).geometry.values.x[0], gdf_centroid.to_crs(self.latlong_crs).geometry.values.y[0]) # assume centroid is focal point of Project
        else:
            # assuming the input centroid is in a long/lat pair
            gdf_centroid = gpd.GeoDataFrame({'type':'centroid', 'geometry': [Point(which_centroid)]}, crs=self.latlong_crs)
            gdf_centroid.to_crs(self.target_crs)
            self.centroid = (gdf_centroid.to_crs(self.latlong_crs).geometry.values.x[0], gdf_centroid.to_crs(self.latlong_crs).geometry.values.y[0]) # assume centroid is focal point of Project
 

        # extract lease area bounds from gdf in units of lat/long
        self.area_longs, self.area_lats = self.gdf.geometry.unary_union.exterior.coords.xy
        # convert lat/long boundary/lease area points to meters away from the centroid
        self.lease_xs, self.lease_ys = self.convertLatLong2Meters(self.area_longs, self.area_lats, self.centroid)
        
    



    def loadBathymetry(self, gebcofilename, moorpy_bathymetry_filename='', dbath='gebco'):
        '''Loads a GEBCO .asc bathymetry file into easy-to-use variables (helps with plotting)
        
        Parameters
        ----------
        gebcofilename : string/path
            path name to the GEBCO .asc file
        latlong_extent : array/list
            array of size 4 for the rectangular extent of bathymetry to use
            Convention is [west longitude, east longitude, north latitude, south latitude]
        moorpy_bathymetry_filename : string/path
            path or string to the name of the file to be created, formatted to input
            bathymetry to MoorPy or MoorDyn
        dbath : float
            discretization spacing of bathymetry grid

        Returns
        -------

        '''
        # load the GEBCO bathymetry file
        depths = -np.loadtxt(gebcofilename, skiprows=6)
        # flip the depths matrix upside down because the first value ([0,0]) in the top left corner corresponds to the first longitude, but last latitude
        depths = np.flipud(depths)
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
        
        # create an array of latitudes and longitudes of the GEBCO bathymetry matrix, using above inputs
        longs = np.linspace(xllcorner, xllcorner+ncols*cellsize, ncols)
        lats  = np.linspace(yllcorner, yllcorner+nrows*cellsize, nrows)

        # save the long/lat points in a mesh grid (to be more easily converted into meters)
        longs_mesh, lats_mesh = np.meshgrid(longs, lats)

        # convert long/lat mesh grid point into meters, relative to the centroid (hstack used for ease of use)
        bathXs_hstack, bathYs_hstack = self.convertLatLong2Meters(np.hstack(longs_mesh), np.hstack(lats_mesh), self.centroid)

        # convert back to mesh grid form
        self.bathXs_mesh = np.zeros([len(lats), len(longs)])
        self.bathYs_mesh = np.zeros([len(lats), len(longs)])
        for i in range(len(bathXs_hstack)):
            iy = i - int(i/len(longs))*len(longs)       # extract the row index from the hstack array
            ix = (len(lats)-1) - int(i/len(longs))      # extract the column index from the hstack array (the addition of len(lats)-1 flips the matrix upside down to what we want)
            self.bathXs_mesh[ix,iy] = bathXs_hstack[i]  # size [len(longs), len(lats)] matrix of x values, using long/lat convention of [0,0] in NW corner of the matrix
            self.bathYs_mesh[ix,iy] = bathYs_hstack[i]  # same comment as above for y values
        self.bath_depths = depths                       # save the depths matrix that came from the GEBCO data

        self.bathYs_mesh = np.flipud(self.bathYs_mesh)

        # save a MoorDyn/MoorPy-style bathymetry, file if desired
        if len(moorpy_bathymetry_filename) > 0:

            f = open(os.path.join(os.getcwd(), moorpy_bathymetry_filename), 'w')
            f.write('--- MoorPy Bathymetry Input File ---\n')
            f.write(f'nGridX {ncols}\n')
            f.write(f'nGridY {nrows}\n')
            f.write(f'      ')
            for ix in range(len(self.bathXs_mesh[0,:])):        # different array of x's depending on the latitude (y) - defaults to top row
                f.write(f'{self.bathXs_mesh[0,ix]:.2f} ')
            f.write('\n')
            for iy in range(len(self.bathYs_mesh)):
                f.write(f'{self.bathYs_mesh[iy,0]:.2f} ')         # different array of y's depending on the longitude (x) - defaults to left column
                for id in range(len(self.bath_depths[iy])):
                    f.write(f'{self.bath_depths[iy,id]} ')
                f.write('\n')
            f.close()
            # results in a skewed picture of the bathymetry grid, since the Xs and Ys converted from long/lat are not perfectly rectangular
        
        # NOTE:
        # the data from GEBCO assumes a square long/lat grid, with longs on x-axis and lats on y-axis
        # however, when you convert that to a meters coordinate system, it distorts the square grid (due to curvature of the Earth)
        # for example, the GEBCO grid assumes the world can be mapped on to a rectangaular shape and discretized by longs/lats with horizontal and vertical lines
        # but, the earth is actually curved, so the long/lat grid is actually slanted (in units of meters) and the square grid becomes distorted
        # to use the bathymetry grid how we want it, we will need to make our own grid in the regular coordinate system, and interpolate the bathymetry depths from the GEBCO data

        





    def plot3d(self, ax=None, figsize=(6,4), area=None, bathymetry=None, boundary=None, area_on_bath=None, args_bath={}):
        '''Plot aspects of the Project object in matplotlib in 3D'''

        # organize the bathymetry arguments
        if len(args_bath)==0:
            args_bath = {'zlim':[-3200,500], 'cmap':'gist_earth'}

        # if axes not passed in, make a new figure
        if ax == None:    
            fig = plt.figure(figsize=figsize)
            ax = plt.axes(projection='3d')
        else:
            fig = ax.get_figure()

        # plot the lease area in a red color, if desired
        if area:
            ax.plot(self.lease_xs, self.lease_ys, np.zeros(len(self.lease_xs)), color='r', zorder=100)
        
        # plot the bathymetry in matplotlib using a plot_surface

        # !!!! include option to plot entire bathymetry file or not

        if isinstance(bathymetry, str):
            bathGrid_Xs, bathGrid_Ys, bathGrid = sbt.readBathymetryFile(bathymetry)         # parse through the MoorDyn/MoorPy-formatted bathymetry file
            X, Y = np.meshgrid(bathGrid_Xs, bathGrid_Ys)                                    # create a 2D mesh of the x and y values
            bath = ax.plot_surface(X, Y, -bathGrid, vmin=args_bath['zlim'][0], vmax=args_bath['zlim'][1], cmap=args_bath['cmap'])
        else:           # >>>>>>>> fix this up
            self.default_depth = 600
            xextent = np.array([-1e5, 1e5])
            yextent = np.array([-1e5, 1e5])
            X, Y = np.meshgrid(xextent, yextent)
            depth = np.ones_like(X)*-self.default_depth
            bath = ax.plot_surface(X, Y, depth, color='b', alpha=0.25)
        
        # plot the project boundary
        if boundary:
            ax.plot(self.boundaryXs, self.boundaryYs, np.zeros(len(self.boundaryXs)), color='k', zorder=100)
            
        # plot the projection of the lease area bounds on the seabed, if desired
        if area_on_bath:
            self.lease_zs = self.projectAlongSeabed(self.lease_xs, self.lease_ys)
            ax.plot(self.lease_xs, self.lease_ys, -self.lease_zs, color='k', zorder=10, alpha=0.5)


            

        set_axes_equal(ax)
        ax.axis('off')
        


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






    # Helper functions

    def projectAlongSeabed(self, x, y):
        '''Project a set of x-y coordinates along a seabed surface (grid),
        returning the corresponding z coordinates.'''
        
        if len(x) == len(y):
            n = len(x)        
            z = np.zeros(n)   # z coordinate of each point [m]
            a = np.zeros(n)   # could also do slope (dz/dh)
            for i in range(n):
                z[i], nvec = sbt.getDepthFromBathymetry(x[i], y[i], self.bathXs, self.bathYs, self.bath_depths)
        
        else:
            z = np.zeros([len(y), len(x)])
            for i in range(len(y)):
                for j in range(len(x)):
                    z[i,j], nvec = sbt.getDepthFromBathymetry(x[j], y[i], self.bathXs_mesh, self.bathYs_mesh, self.bath_depths)
            
        return z

    def convertLatLong2Meters(self, longs, lats, centroid):
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
        '''
        points = np.zeros([len(longs), len(lats)])
        for i in range(len(longs)):
            for j in range(len(lats)):
                points[i,j] = Point(longs[i], lats[j])
        '''
        # organize all the long/lat points into a shapely Polygon
        points = [Point(longs[i],lats[i]) for i in range(len(longs))]

        # input the Polygon of longs/lats and the centroid into a GeoDataFrame in EPSG:4326
        gdf = gpd.GeoDataFrame({'type':'shape','geometry':points}, crs=self.latlong_crs)
        gdf_centroid = gpd.GeoDataFrame({'type':'centroid','geometry': [Point(centroid)]}, crs=self.latlong_crs)
        gdf = pd.concat([gdf, gdf_centroid])        # combine the geodataframes into one (adding rows to a gdf)
        gdf_utm = gdf.to_crs(self.target_crs)       # convert the GeoDataFrame to UTM coordinates (i.e., meters)
        xcentroid = gdf_utm[gdf_utm['type']=='centroid'].geometry.x.values[0]   # extract the centroid in x [m]
        ycentroid = gdf_utm[gdf_utm['type']=='centroid'].geometry.y.values[0]   # extract the centroid in y [m]

        # calculate all the long/lat points distances from the centroid, in units of meters
        xs = np.zeros(len(longs))
        ys = np.zeros(len(lats))
        for i,x in enumerate(gdf_utm[gdf_utm['type']=='shape'].geometry.get_coordinates().x.values):
            xs[i] = np.array(x - xcentroid)
        for i,y in enumerate(gdf_utm[gdf_utm['type']=='shape'].geometry.get_coordinates().y.values):
            ys[i] = np.array(y - ycentroid)
        
        return xs, ys



    




    











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





    # METHODS TO USE WITH ANCHOR TOOLS

    def loadSoil(self, filename):
        '''
        Load geoetechnical information from an input file (format TBD), convert to
        a rectangular grid, and save the grid to the floating array object (TBD).
        
        The input file should provide rows with the following entries:
        - x coordinate
        - y coordinate
        - class  - soil classification name ('clay', 'sand', or 'rock' with optional modifiers)
        - gamma* - soil effective unit weight [kPa] (all soils)
        - Su0*   - undrained shear strength at mudline [kPa] (clay 
        - K*     - undrained shear strength gradient [kPa/m] (clay 
        - alpha* - soil skin friction coefficient [-] (clay soils)
        - phi*   - angle of internal friction [deg] (sand soils)
        
        Some (*) parameters are optional depending on the soil class and mode.   

        Irregular sampling points will be supported and interpolated to a 
        rectangular grid.
        
        Paramaters
        ----------
        filename : path
            path/name of file containing soil data
        '''
        
        # load data from file
        
        # interpolate onto grid defined by grid_x, grid_y
        
        # save
        '''
        self.soil_class
        self.soil_gamma
        self.soil_Su0  
        self.soil_K    
        self.soil_alpha
        self.soil_phi  
        '''
        pass
        

    def getSoilAtLocation(self, x, y):
        '''
        Interpolate soil properties at specified location from the soil
        properties grid and return a dictionary of soil properties that
        can be used in anchor capacity calculations.
        
        Parameters
        ----------        
        x : float
            x coordinate in array reference frame [m].        
        y : float
            y coordinate in array reference frame [m].

        Returns
        -------            
        soilProps : dictionary
            Dictionary of standard MoorPy soil properties.
        '''
        
        soilProps = {}
        

        if self.seabed_type == 'clay':
            
            soilProps['class'] = 'clay'
            soilProps['gamma'] = interp2d(x, y, self.seabed_x, self.seabed_y, self.soil_gamma)
            soilProps['Su0'  ] = interp2d(x, y, self.seabed_x, self.seabed_y, self.soil_Su0  )
            soilProps['k'    ] = interp2d(x, y, self.seabed_x, self.seabed_y, self.soil_k    )
            soilProps['alpha'] = interp2d(x, y, self.seabed_x, self.seabed_y, self.soil_alpha)
            soilProps['phi'  ] = None
        
        elif self.seabed_type == 'sand':
            soilProps['class'] = 'sand'
            soilProps['gamma'] = interp2d(x, y, self.seabed_x, self.seabed_y, self.soil_gamma)
            soilProps['Su0'  ] = None
            soilProps['k'    ] = None
            soilProps['alpha'] = None
            soilProps['phi'  ] = interp2d(x, y, self.seabed_x, self.seabed_y, self.soil_phi  )
            
            # note: for sand, can assume homogeneous angle of internal fricton
        else:
            raise ValueError(f"Unsupported seabed type '{self.seabed_type}'.")
            
        return soilProps


    # ----- Anchor capacity calculation functions -----


    def calcAnchorCapacity(self, anchor):
        '''Compute holding capacity of a given anchor based on the soil
        info at its position. The anchor object's anchor properties and
        location will be used to determine the holding capacity, which
        will be saved to the anchor object.
        
        Parameters
        ----------
        anchor : MoorPy Anchor object (derived from Point)
            The anchor object in question.
        '''

        # interpolate soil properties/class based on anchor position
        anchor.soilProps = getSoilAtLocation(anchor.r[0], anchor.r[1])
        
        # fill in generic anchor properties if anchor info not provided
        if not type(anchor.anchorProps) == dict:
            anchor.anchorProps = dict(type='suction', diameter=6, length=12)
        
        # apply anchor capacity model
        capacity, info = anchorCapacity(anchorProps, soilProps)
        
        # save all information to the anchor (attributes of the Point)
        anchor.soilProps = soilProps
        anchor.anchorCapacity = capacity
        anchor.anchorInfo = info
        
        # also return it
        return capacity

    
    def setCableLayout(self):

        # 2-D

        # 3-D
        pass
    
    # ----- general design-related calculations -----
    
    def makeDistanceMatrix(self):
        '''Compute the distance matrix for an array of turbines. This matrix
        is filled with the horizontal distance between every turbine's 
        undisplaced position.
        '''
        
        dists = np.zeros([self.nPtfm, self.nPtfm])  # distance matrix
        
        for i in range(self.nPtfm):
            for j in range(self.nPtfm):
                delta = self.coords[i] - self.coords[j]
                dists[i,j] = np.diag(delta)
                dists[j,i] = dists[i,j]
                
        return dmat
    
    
    # ----- cable calculation methods -----
    
    def calcCableLength(self, cable):
        '''Calculates a cable's length based on its routing.
        '''

        # select cable
        
        # figure out cable length considering end coordinates and path
        
        return length
    
    
    def checkCableExclusions(self, cable):
        '''Checks whether a cable crosses over any exclusions
        or other out of bounds areas.
        '''

        # select cable
        
        # check its path against any exclusion areas or boundaries
        
        # make a list of any exclusion/nodes that it is too close to
        
        return score, list_of_violations
    



    """
    def loadBoundary(self, filename):
        '''
        Load a lease area boundary for the project from an input file.
        
        Parameters
        ----------

        filename : path
            path/name of file containing bathymetry data (format TBD)
        '''
        
        # load data from file
        Xs, Ys = sbt.processBoundary(filename, self.lat0, self.lon0)
        
        # check compatibility with project grid size
        
        # save as project boundaries
        self.boundary_Xs = Xs
        self.boundary_Ys = Ys
        
        # figure out masking to exclude grid data outside the project boundary
        
        
    def loadBathymetry(self, filename):
        '''
        Load bathymetry information from an input file (format TBD), convert to
        a rectangular grid, and save the grid to the floating array object (TBD).
        
        Paramaters
        ----------
        filename : path
            path/name of file containing bathymetry data (format TBD)
        '''
        
        # load data from file
        Xs, Ys, Zs = sbt.processASC(filename, self.lat0, self.lon0)
        
        # ----- map to existing grid -----
        # if no grid, just use the bathymetry grid
        if len(self.grid_x) == 0:
            self.grid_x = np.array(Xs)
            self.grid_y = np.array(Ys)
            self.depths = np.array(Zs)
        # interpolate onto grid defined by grid_x, grid_y
        else:
            for i, x in enumerate(self.grid_x):
                for j, y in enumerate(self.grid_y):
                    self.depths[i,j], _ = sbt.getDepthFromBathymetry(x, y, Xs, Ys, Zs)
        
        
        # also save in RAFT, in its MoorPy System(s)
    """



def getFromDict(dict, key, shape=0, dtype=float, default=None, index=None):
    '''
    Function to streamline getting values from design dictionary from YAML file, including error checking.

    Parameters
    ----------
    dict : dict
        the dictionary
    key : string
        the key in the dictionary
    shape : list, optional
        The desired shape of the output. If not provided, assuming scalar output. If -1, any input shape is used.
    dtype : type
        Must be a python type than can serve as a function to format the input value to the right type.
    default : number or list, optional
        The default value to fill in if the item isn't in the dictionary. 
        Otherwise will raise error if the key doesn't exist. It may be a list
        (to be tiled shape times if shape > 1) but may not be a numpy array.
    '''
    # in future could support nested keys   if type(key)==list: ...

    if key in dict:
        val = dict[key]                                      # get the value from the dictionary
        if shape==0:                                         # scalar input expected
            if np.isscalar(val):
                return dtype(val)
            else:
                raise ValueError(f"Value for key '{key}' is expected to be a scalar but instead is: {val}")
        elif shape==-1:                                      # any input shape accepted
            if np.isscalar(val):
                return dtype(val)
            else:
                return np.array(val, dtype=dtype)
        else:
            if np.isscalar(val):                             # if a scalar value is provided and we need to produce an array (of any shape)
                return np.tile(dtype(val), shape)

            elif np.isscalar(shape):                         # if expecting a 1D array (or if wanting the result to have the same length as the input)
                if len(val) == shape:                        # throw an error if the input is not the same length as the shape, meaning the user is missing data
                    if index == None:
                        return np.array([dtype(v) for v in val])    # if no index is provided, do normally and return the array input
                    else:
                        keyshape = np.array(val).shape              # otherwise, use the index to create the output arrays desired
                        if len(keyshape) == 1:                      # if the input is 1D, shape=n, and index!=None, then tile the indexed value of length shape
                            if index in range(keyshape[0]):
                                return np.tile(val[index], shape)
                            else:
                                raise ValueError(f"Value for index '{index}' is not within the size of {val} (len={keyshape[0]})")
                        else:                                               # if the input is 2D, len(val)=shape, and index!=None
                            if index in range(keyshape[1]):
                                return np.array([v[index] for v in val])    # then pull the indexed value out of each row of 2D input
                            else:
                                raise ValueError(f"Value for index '{index}' is not within the size of {val} (len={keyshape[0]})")
                else:
                    raise ValueError(f"Value for key '{key}' is not the expected size of {shape} and is instead: {val}")

            else:                                            # must be expecting a multi-D array
                vala = np.array(val, dtype=dtype)            # make array

                if list(vala.shape) == shape:                      # if provided with the right shape
                    return vala
                elif len(shape) > 2:
                    raise ValueError("Function getFromDict isn't set up for shapes larger than 2 dimensions")
                elif vala.ndim==1 and len(vala)==shape[1]:   # if we expect an MxN array, and an array of size N is provided, tile it M times
                    return np.tile(vala, [shape[0], 1] )
                else:
                    raise ValueError(f"Value for key '{key}' is not a compatible size for target size of {shape} and is instead: {val}")

    else:
        if default == None:
            raise ValueError(f"Key '{key}' not found in input file...")
        else:
            if shape==0 or shape==-1:
                return default
            else:
                if np.isscalar(default):
                    return np.tile(default, shape)
                else:
                    return np.tile(default, [shape, 1])
    

'''
Other future items:
Cost calc functions
System Reliability/failure analysis functions
Full scenario visualization functions
Load case setup and constraint eval?
'''
