"""Project class for FAModel, containing information and key methods for
the site information and design information that make up a project."""

import numpy as np
import matplotlib.pyplot as plt

import moorpy as mp
# import raft

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
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
    
    def __init__(self, centroid, lat=0, lon=0, depth=100):
        

        # need to check if there is an input centroid or a shapefile (for a lease area)
        self.centroid = centroid
        self.gdf = gpd.GeoDataFrame({'name':'centroid', 'geometry': [Point(self.centroid)]}, crs='EPSG:4326')

        # set the target coordinate reference system (CRS) that will switch between regular lat/long system and "meters from centroid" system
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


        self.lat0  = lat  # lattitude of site reference point [deg]
        self.lon0  = lon  # longitude of site reference point [deg]

        # project boundaries
        self.boundary_Xs = np.zeros(0)
        self.boundary_Ys = np.zeros(0)
        
        #self.array = RAFT_model

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
        usa.insert(0, 'name', statenamelist)
        # set the CRS of the USA pdf to the right CRS
        usa.set_crs(crs="EPSG:4326", inplace=True)
        self.usa = usa
        
        for state in states:
            state_gs = usa.loc[usa['name']==state]
            self.gdf = pd.concat([self.gdf, state_gs])
    


    def setFarmLayout(self, style='grid', nrows=10, ncols=10, turbine_spacing=2000, nOSS=2):

        if style=='grid':
            # for now, this is very custom code specific to the DeepFarm project
            farmxspacing = (nrows-1)*turbine_spacing
            farmyspacing = (ncols-1)*turbine_spacing
            turbine_distances_from_centroid = []
            for j in reversed(range(ncols)):
                for i in range(nrows):
                    xpos = -(farmxspacing/2)+(i*turbine_spacing)
                    ypos = -(farmyspacing/2)+(j*turbine_spacing)
                    turbine_distances_from_centroid.append((xpos, ypos))
            # add positions of two offshore substations (OSSs)
            turbine_distances_from_centroid.append((11000.0, 5000.0))
            turbine_distances_from_centroid.append((11000.0, -5000.0))

        # create a copy of the global gdf and transform it into the easting/northing coordinate reference system
        gdf_utm = self.gdf.copy().to_crs(self.target_crs)
        xcentroid = gdf_utm.loc[gdf_utm['name']=='centroid'].centroid.x[0]
        ycentroid = gdf_utm.loc[gdf_utm['name']=='centroid'].centroid.y[0]

        # create shapely Point objects of the turbine positions relative to the centroid, in the UTM CRS
        turbine_geoms = []
        for i,(x,y) in enumerate(turbine_distances_from_centroid):
            turbine_geoms.append( Point(xcentroid + x, ycentroid + y) )
        
        # make a new gdf to put this data together
        turbine_gdf = gpd.GeoDataFrame({'name': ['turbine_pos']*len(turbine_geoms), 'geometry': turbine_geoms}, crs=self.target_crs)

        # convert the turbine coordinates back to regular latitude/longitude (EPSG: 4326)
        turbine_gdf = turbine_gdf.to_crs('EPSG:4326')

        # add the turbine gdf to the global gdf
        self.gdf = pd.concat([self.gdf, turbine_gdf])


        

    def plot(self, kwargs):

        if 'centroid' in kwargs:
            centroid_settings = kwargs['centroid']
        if 'map' in kwargs:
            map_settings = kwargs['map']
        if 'farm' in kwargs:
            farm_settings = kwargs['farm']

        fig, ax = plt.subplots(1,1)
        self.gdf.loc[self.gdf['name']=='centroid'].plot(ax=ax, color=centroid_settings['color'])

        if 'boundary' in map_settings:
            map_boundary = self.gdf.loc[self.gdf['name']=='California'].boundary
            map_boundary.plot(ax=ax, color=map_settings['color'])
        
        if 'farm' in kwargs:
            self.gdf.loc[self.gdf['name']=='turbine_pos'].plot(ax=ax, color=farm_settings['color'])
        
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        
        # Some GeoPandas Help
        # to plot just one entry of a geoseries: gdf.loc[[0],'geometry'].plot()
        # to get the columns of a gdf: gdf.columns
        # merging gdf's
        # adding columns to gdf's

        return fig, ax


    def addPoints(self, ax, pointlist=[], kwargs={}):

        point_settings = kwargs['pointlist']

        points = gpd.GeoDataFrame({'name':['nrel_channel','nrel_humboldt','nrel_crescent_city','hawaii'],
                                   'geometry': pointlist}, crs='EPSG:4326')
        
        points.plot(ax=ax, color=point_settings['color'], marker=point_settings['marker'])
    

    def addState(self, ax, states=[], kwargs={}):
        
        for state in states:
            state_settings=kwargs[state]

            state_geom = self.usa.loc[self.usa['name']==state]
            if 'boundary' in state_settings:
                state_geom = state_geom.boundary

            newstate = gpd.GeoDataFrame({'name':state, 'geometry':state_geom}, crs='EPSG:4326')

            newstate.plot(ax=ax, color=state_settings['color'])





    # ----- Site conditions processing functions -----

    def setGrid(self, xCoords, yCoords):
        '''
        Set up the rectangular grid over which site or seabed
        data will be saved and worked with. Directions x and y are 
        generally assumed to be aligned with the East and North 
        directions, respectively, at the array reference point.
        
        Parameters
        ----------        
        xCoords : float array
            x coordinates relative to array reference point [m]
        yCoords : float array
            x coordinates relative to array reference point [m]
        '''
        
        self.grid_x = np.array(xCoords)
        self.grid_y = np.array(yCoords)
        
        #TODO: add check for existing seabed data. If present, convert or raise warning <<<
    
    
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
        Xs, Yx, Zs = sbt.processASC(filename, self.lat0, self.lon0)
        
        # interpolate onto grid defined by grid_x, grid_y
        
        # save in object
        
        # also save in RAFT, in its MoorPy System(s)
        
        
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
    

'''
Other future items:
Cost calc functions
System Reliability/failure analysis functions
Full scenario visualization functions
Load case setup and constraint eval?
'''
