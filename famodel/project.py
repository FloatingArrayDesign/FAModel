"""Project class for FAModel, containing information and key methods for
the site information and design information that make up a project."""

import os
import numpy as np
import matplotlib.pyplot as plt

import moorpy as mp
from moorpy.helpers import set_axes_equal
import yaml

#from shapely.geometry import Point, Polygon, LineString

from .anchors.anchor_capacity import anchorCapacity
from .seabed import seabed_tools as sbt
from .mooring.mooring import Mooring


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
        self.turbineList = []
        self.mooringList = []  # A list of Mooring objects
        self.anchorList  = []
        self.cables = None  # CableSystem
        
        
        # ----- site information -----
        self.lat0  = lat  # lattitude of site reference point [deg]
        self.lon0  = lon  # longitude of site reference point [deg]

        # Project boundary (list of x,y coordinate pairs [m])
        self.boundary = np.zeros([0,2])
        
        # Seabed grid
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
        
        
        # ----- mooring line section types -----
        # add support for both list- and table-based formats 
        lineTypes = {}
        
        if 'mooring_line_types' in d:
            
            for k, v in d['mooring_line_types'].items():
                
                lineTypes[k] = v
        
        # ----- mooring connectors -----
        # to be added
        
        # ----- anchor types -----
        # to be added
        
        # ----- mooring line configurations -----
        lineTypes = {}
        
        if 'mooring_line_configs' in d:
            
            for k, v in d['mooring_line_configs'].items():
                pass
                #set up mooring config... (k, v)
                
                # I think we want to just store it as a dictionary
                # validate it,
                # then have it available for use when making Mooring objects and subsystems
        
        
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

        # load bathymetry information, if provided
        if 'bathymetry' in site:
            if 'file' in site['bathymetry']:
                self.loadBathymetry(site['bathymetry']['file'])
            elif 'x_y_z' in site['bathymetry']:
                xyz = site['bathymetry']['x_y_z']
                xs = np.unique( [point[0] for point in xyz] )
                ys = np.unique( [point[1] for point in xyz] )
                depths = np.zeros( [len(ys), len(xs)] )
                for iy in range(len(depths)):
                    for ix in range(len(depths[0])):
                        x = xs[ix]; y = ys[iy]
                        for point in xyz:
                            if point[0]==x and point[1]==y:
                                depths[iy,ix] = point[2]
                depths = np.flipud(depths)      # flip upside down to equate to N-S on local coordinate system

        # Load project boundary, if provided
        if 'boundaries' in site:
            if 'file' in site['boundaries']:  # load boundary data from file
                self.loadBoundary(site['boundaries']['file'])
            elif 'x_y' in site['boundaries']:  # process list of buondary x,y vertices
                xy = site['boundaries']['x_y']
                self.boundary = np.zeros([len(xy),2])
                for i in range(len(xy)):
                    self.boundary[i,0] = xy[i][0]
                    self.boundary[i,1] = xy[i][1]

        # and set the project boundary/grid based on the loaded information
        # TBD, may not be necessary in the short term. self.setGrid(xs, ys)
        
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
        
        # Create a new depth matrix with interpolated values from the original
        depths = np.zeros([len(ys), len(xs)])
        for i in range(len(ys)):
            for j in range(len(xs)):
                depths[i,j], nvec = sbt.getDepthFromBathymetry(xs[j], ys[i], 
                                    self.grid_x, self.grid_y, self.grid_depth)
        
        # Replace the grid data with the updated values
        self.grid_x = np.array(xs)
        self.grid_y = np.array(ys)
        self.grid_depth = depths


    # Helper functions

    def getDepthAtLocation(self, x, y, return_n=False):
        '''Compute the depth at a specified x-y location based on the
        bathymetry grid stored in the project. Setting return_n to 
        True will return depth and the seabed normal vector.
        '''
        
        z, n = sbt.getDepthFromBathymetry(x, y, self.grid_x, self.grid_y, 
                                          self.grid_depth)
        
        if return_n:
            return z, n
        else:
            return z


    def seabedIntersect(self, r, u):
        '''
        Compute location at which a ray (defined by coordinate i and direction
        vector u) crosses the seabed.
        
        Paramaters
        ----------
        r : float, array
            Absolute xyz coordinate of a point along the ray [m].
        u : float, array
            3D direction vector of the ray.

        Returns
        -------
        r_i: float, array
            xyz coordinate of the seabed intersection [m].
        '''
        
        # Initialize 
        r_i = np.array(r)  # calculated seabed intersect starts at the r provided
        ix_last = -1  # index of intersected grid panel (from previous iteration)
        iy_last = -1
        
        for i in range(10):  # iterate up to 10 times (only needed for fine grids)
            
            # Get depth of seabed beneath r_i, and index of grid panel
            depth, nvec, ix, iy = sbt.getDepthFromBathymetry(r_i[0], r_i[1], self.grid_x, 
                                       self.grid_y, self.grid_depth, index=True)
            
            #print(f" {i}  {depth:6.0f}   {ix}, {iy}   {r_i[0]:5.0f},{r_i[1]:5.0f},{r_i[2]:5.0f}")
            
            # Convergence check (if we're on the same grid panel, we're good)
            if ix==ix_last and iy==iy_last:
                break

            # Save some values for convergence checking
            ix_last = ix
            iy_last = iy
            r_last  = np.array(r_i)
            
            # vectors from r_i to a point on the bathymetry panel beneath it
            # w = np.array([r_anch[0], r_anch[1], -depth]) - r_i
            w = np.array([0, 0, -depth - r_i[2]]) # same!! as above
            
            # fraction along u where it crosses the seabed (can be greater than 1)
            fac = np.dot(nvec, w) / np.dot(nvec, u)
            
            r_i = r_i + u*fac  # updated seabed crossing estimate 

            if any(np.isnan(r_i)):
                breakpoint()
        
        return r_i


    def projectAlongSeabed(self, x, y):
        '''Project a set of x-y coordinates along a seabed surface (grid),
        returning the corresponding z coordinates.'''
        
        if len(x) == len(y):
            n = len(x)        
            z = np.zeros(n)   # z coordinate of each point [m]
            a = np.zeros(n)   # could also do slope (dz/dh)
            for i in range(n):
                z[i] = self.getDepthAtLocation(x[i], y[i])
        
        else:
            z = np.zeros([len(y), len(x)])
            for i in range(len(y)):
                for j in range(len(x)):
                    z[i,j] = getDepthAtLocation(x[j], y[i])
            
        return z



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
        
        # SIMPLE HACK FOR NOW
        Xs, Ys, Rockys = sbt.readBathymetryFile(filename)  # read MoorDyn-style file
        
        self.soil_x = np.array(Xs)
        self.soil_y = np.array(Ys)
        self.soil_rocky = np.array(Rockys)
        
        
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
        
        # SIMPLE HACK FOR NOW        
        rocky, _,_,_,_ = sbt.interpFromGrid(x, y, self.soil_x, self.soil_y, self.soil_rocky)
        
        return rocky
        
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
        '''

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
        self.boundary = np.vstack([Xs, Ys])
        
        # figure out masking to exclude grid data outside the project boundary
        
        
    def loadBathymetry(self, filename, interpolate=False):
        '''
        Load bathymetry information from an input file (format TBD), convert to
        a rectangular grid, and save the grid to the floating array object (TBD).
        
        Paramaters
        ----------
        filename : path
            path/name of file containing bathymetry data (format TBD)
        '''
        
        # load data from file
        Xs, Ys, Zs = sbt.readBathymetryFile(filename)  # read MoorDyn-style file
        # Xs, Ys, Zs = sbt.processASC(filename, self.lat0, self.lon0)
        
        # ----- map to existing grid -----
        # if no grid, just use the bathymetry grid
        if not interpolate: #len(self.grid_x) == 0: 
            self.grid_x = np.array(Xs)
            self.grid_y = np.array(Ys)
            self.grid_depth = np.array(Zs)
            
        else:
        # interpolate onto grid defined by grid_x, grid_y
            for i, x in enumerate(self.grid_x):
                for j, y in enumerate(self.grid_y):
                    self.grid_depth[i,j], _ = sbt.getDepthFromBathymetry(x, y, Xs, Ys, Zs)
        
        
        # also save in RAFT, in its MoorPy System(s)
    


    def plot3d(self, ax=None, figsize=(10,8), fowt=None, save=False,
               draw_boundary=True, boundary_on_bath=True, args_bath={}, draw_axes=True):
        '''Plot aspects of the Project object in matplotlib in 3D.
        
        TODO - harmonize a lot of the seabed stuff with MoorPy System.plot...
        
        Parameters
        ----------
        ...
        '''
        
        # color map for soil plotting
        import matplotlib.cm as cm
        from matplotlib.colors import Normalize
        cmap = cm.cividis_r
        norm = Normalize(vmin=-0.5, vmax=1.5)
        #print(cmap(norm(np.array([0,1]))))
        

        # if axes not passed in, make a new figure
        if ax == None:    
            fig = plt.figure(figsize=figsize)
            ax = plt.axes(projection='3d')
        else:
            fig = ax.get_figure()

        # try icnraesing grid density
        xs = np.arange(-1000,8000,500)
        ys = np.arange(-1000,9500,500)
        self.setGrid(xs, ys)

        # plot the bathymetry in matplotlib using a plot_surface
        X, Y = np.meshgrid(self.grid_x, self.grid_y)  # 2D mesh of seabed grid
        '''
        # interpolate soil rockyness factor onto this grid
        xs = self.grid_x
        ys = self.grid_y
        rocky = np.zeros([len(ys), len(xs)])
        for i in range(len(ys)):
            for j in range(len(xs)):
                rocky[i,j], _,_,_,_ = sbt.interpFromGrid(xs[j], ys[i], 
                           self.soil_x, self.soil_y, self.soil_rocky)
        # apply colormap
        rc = cmap(norm(rocky))
        bath = ax.plot_surface(X, Y, -self.grid_depth, facecolors=rc, **args_bath)
        '''
        bath = ax.plot_surface(X, Y, -self.grid_depth, **args_bath)
        
        
        # also if there are rocky bits... (TEMPORARY)
        X, Y = np.meshgrid(self.soil_x, self.soil_y)
        ax.scatter(X, Y, c=self.soil_rocky, s=4, cmap='cividis_r', vmin=-0.5, vmax=1.5)
        
        # plot the project boundary
        if draw_boundary:
            boundary = np.vstack([self.boundary, self.boundary[0,:]])
            ax.plot(boundary[:,0], boundary[:,1], np.zeros(boundary.shape[0]), 
                    'b--', zorder=100, lw=1, alpha=0.5)
            
        # plot the projection of the boundary on the seabed, if desired
        if boundary_on_bath:
            boundary_z = self.projectAlongSeabed(boundary[:,0], boundary[:,1])
            ax.plot(boundary[:,0], boundary[:,1], -boundary_z, 'k--', zorder=10, lw=1, alpha=0.7)

        # plot the Moorings
        for mooring in self.mooringList:
            #mooring.subsystem.plot(ax = ax, draw_seabed=False)
            if mooring.subsystem:
                mooring.subsystem.drawLine(0, ax)
        
        # plot the FOWTs using a RAFT FOWT if one is passed in (TEMPORARY)
        if fowt:
            for i in range(self.nt):
                xy = self.turb_coords[i,:]
                fowt.setPosition([xy[0], xy[1], 0,0,0,0])
                fowt.plot(ax, zorder=20)
        
        # Show full depth range
        ax.set_zlim([-np.max(self.grid_depth), 0])

        set_axes_equal(ax)
        if not draw_axes:
            ax.axis('off')
        
        ax.view_init(20, -130)
        ax.dist -= 3
        fig.tight_layout()
        
        # ----- Save plot with an incremented number if it already exists
        if save:
            counter = 1
            output_filename = f'wind farm 3d_{counter}.png'
            while os.path.exists(output_filename):
                counter += 1
                output_filename = f'wind farm 3d_{counter}.png'
            
            # Increase the resolution when saving the plot
            plt.savefig(output_filename, dpi=300, bbox_inches='tight')  # Adjust the dpi as needed
            


    def plot2d(self, ax=None, plot_seabed=True, plot_boundary=True):
        '''Plot aspects of the Project object in matplotlib in 3D.
        
        TODO - harmonize a lot of the seabed stuff with MoorPy System.plot...
        
        Parameters
        ----------
        Placeholder
        '''
     
        # if axes not passed in, make a new figure
        if ax == None:
            fig, ax = plt.subplots(1,1, figsize=figsize)
        else:
            fig = ax.get_figure()
       
        

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
