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
    
    def __init__(self, centroid, lat=0, lon=0, file=None, depth=100):
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

        # need to check if there is an input centroid or a shapefile (for a lease area)
        self.centroid = centroid
        self.gdf = gpd.GeoDataFrame({'type':'centroid', 'geometry': [Point(self.centroid)]}, crs='EPSG:4326')

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


    def plot(self, kwargs):
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


    def plot3d(self):
        '''3d plot'''

        # >>> to be adjusted to integrate >>>
        xbmin, xbmax = sbt.getPlotBounds(longs, zerozero, long=True)
        ybmin, ybmax = sbt.getPlotBounds(lats, zerozero, long=False)

        ms = mp.System(file='example/sample_deep.txt', bathymetry=bathymetryfilename)
        ms.initialize()

        xcoast, ycoast = sbt.getCoast(Xbath, Ybath, depths)

        fig, ax = ms.plot(hidebox=True, shadow=False, bathymetry=bathymetryfilename, cmap='gist_earth', rang=(-3200, 500), xbounds=(xbmin, xbmax), ybounds=(ybmin, ybmax))
        ax.plot(Xbnds_ne, Ybnds_ne, np.zeros(len(Xbnds_ne)), color='r', zorder=100)
        ax.plot(Xbnds_sw, Ybnds_sw, np.zeros(len(Xbnds_sw)), color='r', zorder=100)
        ax.plot(Xbnds_ne, Ybnds_ne, -Dbnds_ne, color='k', zorder=10, alpha=0.5)
        ax.plot(Xbnds_sw, Ybnds_sw, -Dbnds_sw, color='k', zorder=10, alpha=0.5)
        ax.plot(xcoast, ycoast, np.zeros(len(Ybath)), color='k', zorder=100)


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





    def load(self, input):
        '''
        Load a full set of project information from a dictionary or 
        YAML file. This calls other methods for each part of it.
        
        Parameters
        ----------
        input : dict or filename
            Dictionary or YAML filename containing project info.
        '''
        
        # standard function to load dict if input is yaml
        
        # look for site section
        # call load site method
        self.loadSite(d['site'])
        
        # look for design section
        # call load design method
        self.loadDesign(d)
    

    # ----- Design loading/processing methods -----
    
    def loadDesign(input):
        '''Load design information from a dictionary or YAML file
        (specified by input). This should be the design portion of
        the floating wind array ontology.'''
        
        # standard function to load dict if input is yaml
        d = 
        
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

    def loadSite(input):
        '''Load site information from a dictionary or YAML file
        (specified by input). This should be the site portion of
        the floating wind array ontology.'''
        
        # standard function to load dict if input is yaml
        
        # load seabed portions
        
        # load lease area portions
        
        # load metocean portions
        
        
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
        
    
    def projectAlongSeabed(self, x, y):
        '''Project a set of x-y coordinates along a seabed surface (grid),
        returning the corresponding z coordinates.'''
        
        if len(x) == len(y):
            n = len(x)
        else:
            raise Exception('x and y inputs must be same length.')
        
        z = np.zeros(n)   # z coordinate of each point [m]
        a = np.zeros(n)   # could also do slope (dz/dh)
        
        for i in range(n):
            z[i], nvec = sbt.getDepthFromBathymetry(Xbnds_ne[i], Ybnds_ne[i], Xbath, Ybath, depths_flipped)
            
            # a = ...f(nvec)
            
        return z


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
    

'''
Other future items:
Cost calc functions
System Reliability/failure analysis functions
Full scenario visualization functions
Load case setup and constraint eval?
'''
