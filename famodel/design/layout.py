import os
import moorpy as mp
from moorpy.helpers import getFromDict
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LogNorm
from matplotlib.collections import PolyCollection
from matplotlib.ticker import AutoLocator, AutoMinorLocator
from matplotlib.ticker import FuncFormatter
from scipy.interpolate import NearestNDInterpolator
from scipy import interpolate, optimize
from scipy.optimize import minimize, differential_evolution, NonlinearConstraint
from scipy.spatial.distance import cdist, pdist, squareform
from scipy.spatial import distance
from scipy import optimize
from sklearn.cluster import SpectralClustering # KMeans

import random
import csv
# from moorpy.helpers import getFromDict
# from shapely import Point, Polygon
import shapely as sh
from shapely.geometry import Point, LineString, MultiLineString, Polygon, MultiPolygon
from shapely.ops import unary_union, nearest_points

import shapely.geometry
import shapely.affinity
from shapely.affinity import translate
#import shapely.affinity as sa
import networkx as nx

import yaml
#import raft
from copy import deepcopy

from moorpy.helpers import set_axes_equal

import fadesign

from famodel.project import Project
from famodel.mooring.mooring import Mooring
from famodel.anchors.anchor import Anchor
from famodel.platform.platform import Platform
from famodel.cables.cable import Cable
from famodel.cables.cable_properties import loadCableProps, getCableProps
from famodel.substation.substation import Substation

from fadesign.layout_helpers import getLower, makeMooringListN
from fadesign.CableLayout_functions import getCableLayout

import floris
from floris import FlorisModel

# from floris.turbine_library import TurbineInterface, TurbineLibrary

from pyswarm import pso

# Import PySwarms
#import pyswarms as pso
#from pyswarms.utils.functions import single_obj as fx
#from pyswarms.utils.plotters import (plot_cost_history, plot_contour, plot_surface)
#from pyswarms.utils.plotters.formatters import Mesher


# SPYDER
# Interactive plots on 
#%matplotlib qt 
# Interactive plots off
#%matplotlib inline


class Layout(Project):
    '''A class to store and work with the layout information of a wind farm.'''

    def __init__(self, X, Xu, Xdb =[], wind_rose = [], ss = None, mooringAdjuster = None, **kwargs):

        '''Create a new Layout object that can be used for optimization.
        The Layout class allows storage of various data for the layout design
        problem, such as the boundaries, the seabed bathymetry, and wind rose.
        This initialization function sets those data. Many of the data inputs
        are optional and default to not being used.
        
        For FREE LAYOUT OPTIMIZATION: X = [x,y,phi], Xu = []
        For UNIFORM GRID LAYOUT OPTIMIZATION: X = [], 
        Xu = [ spacing_x, spacing_y, trans_x, trans_y, rotang, skew]
   

        Parameters
        ----------
        X : 1D array
            Design vector considered within the optimization [x,y,phi]
            x, y : turbine positions in (m)
            phi : turbine heading in (deg)
        Xu : 1D array
            Design vector considered within the uniform grid optimization 
            [grid_spacing_x,grid_spacing_y,grid_trans_x, grid_trans_y, grid_rotang, grid_skew]
            grid_spacing_x, y : x,y turbine spacing in (m)
            grid_trans_x, y : x,y translation of entire grid in (m)
            grid_rotang : rotation of grid around centroid of lease area (deg) 
            grid_skew : skew angle of grid (deg)
        Xdb :
            ???
        nTurbines : int
            Number of turbines to work with.
        boundary_coords : 2D array
            List of x coordinates of lease area boundary vertices (m).
            List of y coordinates of lease area boundary vertices (m).

        grid_x : 1D array
            List of x coordinates of bathymetry grid (km).
        grid_y : 1D array
            List of y coordinates of bathymetry grid (km).
        grid_depth : 2D array
            Matrix of depth values corresponding to x,y coordinates (m).

        wind_rose  : FLORIS wind rose
            A wind rose of wind speeds, direction, frequency and TI
        ss : MoorPy Subsystem, optional
            A MoorPy Subsystem to adapt for a 3D representation of each mooring line.
        mode : string
            'LCOE', 'AEP' or 'CAPEX'.
        rotation_mode : Bool
            True : considering rotation as part of the design vector as design variable 
            False: not considering rotation as design variable
        
        turb_minrad=360
            Radius of turbine buffer zone.
        moor_minrad=50
            Radius of mooring buffer zone.
        moorOpt_mode : string
            'basic' : Basic mooring layout, without MoorPy
            'advanced' : Mooring layout considers MoorPy input

        
        
        '''
        # Initialize Project aspects to start with
        super().__init__()
        
        self.display = getFromDict(kwargs, 'display', default=0)
        
        # add seabed bathymetry (based on file for now)
        self.bathymetry_file = getFromDict(kwargs, 'bathymetry_file', dtype=str, default = '')
        self.loadBathymetry(self.bathymetry_file)
        self.soil_file =  getFromDict(kwargs, 'soil_file', dtype=str, default = '')
        self.loadSoil(self.soil_file)
        self.cable_mode= getFromDict(kwargs, 'cable_mode', default = True)
        
        
        # ----- Optimization modes -----
        self.optimizer =  getFromDict(kwargs, 'optimizer', dtype=str, default = '')          # Optimizer
        self.obj_penalty =  getFromDict(kwargs, 'obj_penalty', default = True)      # Use penalty factor in objective function yes (1) or no (0)
        self.mode =  getFromDict(kwargs, 'mode', dtype=str, default = 'LCOE')                    # Optimization mode
        self.rotation_mode =  getFromDict(kwargs, 'rotation_mode', default = True)  # Rotation included as Design Variable or not
        self.alternate_rows = getFromDict(kwargs, 'alternate_rows', default = False )
        self.log = dict(x=[], f=[], g=[])  # initialize a log dict with empty values  
        self.iter = -1  # iteration number of a given optimization run (incremented by updateDesign)
        self.parallel =  getFromDict(kwargs, 'parallel', default = False)
        self.infeasible_obj_update = getFromDict(kwargs, 'infeasible_obj_update', default = False) # set True to update objective function even when layout violates constraints
        
        
        # -----  Turbine quantity and positions -----
        self.nt =  int(getFromDict(kwargs, 'n_turbines', default = 67))                # Number of turbines
        self.turb_coords= np.zeros((self.nt,2)) # Turbine positions (x,y) [m]
        self.turb_depth = np.zeros(self.nt)
        
        self.turb_minrad = getFromDict(kwargs, 'turb_minrad', default = 200)
        self.moor_minrad = getFromDict(kwargs, 'moor_minrad', default = 20)
        self.anchor_minrad = getFromDict(kwargs, 'anchor_minrad', default = 50)
        
        self.turb_mindist_m = np.zeros(self.nt)  # currently inactive
        self.con_turb_turb = np.zeros(self.nt)  # distance to turbine
        # distance to boundary - considering anchor radius
        self.con_turb_boundary = np.zeros(self.nt)
        # distance to boundary - center of WTG
        self.turb_dist_tb2_m = np.zeros(self.nt)
        
        if np.size(Xu) != 0:
            self.Xlast = np.zeros((self.nt)*2) # X vector containting X and Y coordinates
        else:
            self.Xlast = np.zeros_like(X)
        
        self.obj_value = 0
        self.turb_rating_MW = getFromDict(kwargs, 'turb_rating_MW', default = 15)  # Rating of each turbine in MW  
        # IAC System parameters
        self.iac_voltage_kV = getFromDict(kwargs, 'iac_voltage_kV', default = 66)  # Voltage level in kV
        self.iac_type = 'dynamic_cable_66' # Cable type, as defined in cable properties YAML
        
        # Turbine electrical current in ampere A
        self.turb_I = (self.turb_rating_MW  * 1e6) / (self.iac_voltage_kV * 1e3) 
        
        # Cable conductor sizes for 66 kV transmission system
        # List to be further specified
        
        #dir = os.path.dirname(os.path.realpath(__file__))
        #with open(os.path.join(dir,"CableProps_default.yaml")) as file:
        #    source = yaml.load(file, Loader=yaml.FullLoader)    
        #As = source['conductor_size']['size_A_df']  
        self.iac_typical_conductor = getFromDict(kwargs, 'iac_typical_conductor',shape=-1, default = [0])
        if len(self.iac_typical_conductor)==1 and self.iac_typical_conductor[0] == 0:
            self.iac_typical_conductor = np.array([ 70, 95, 120, 150, 185, 240, 300, 400, 500, 630,
            800, 1000, 1200,1400,1600,200,2500,3000,4000])

        
        # -----  Offshore Substation -----
        self.noss = int(getFromDict(kwargs,'noss', default = 1))  
        self.oss_coords_initial = getFromDict(kwargs, 'oss_coords', shape=-1, default = np.zeros((self.noss,2))) # initial OSS coordinates
        # adjust to a nested list [[x,y]] if given oss coords as [x,y] for compatibility with situations with multiple oss
        if self.oss_coords_initial.shape == (2,):
            self.oss_coords_initial = np.array([self.oss_coords_initial])
        # for now we'll set oss_coords = oss_coords_initial but this could change in generateGridPoints if using uniform grid layout
        self.oss_coords = deepcopy(self.oss_coords_initial)
        self.oss_minrad = getFromDict(kwargs, 'oss_minrad', default = self.turb_minrad*2)
        self.static_substations = getFromDict(kwargs, 'static_substations', dtype=bool, default = False)
        
        # create substation platform object
        for oo in range(self.noss):
            r = [self.oss_coords[oo][0], self.oss_coords[oo][1], 0]
            self.platformList[self.nt+oo] = Platform(id=self.nt+oo,r=r,rFair=ss.rad_fair,zFair=ss.z_fair)
            self.platformList[self.nt+oo].entity = 'Substation'
            
            
        
        
        # -----  Turbine Cluster -----
        self.n_cluster = int(getFromDict(kwargs, 'n_cluster', default = 9))# Amount of turbine cluster for cable routing
        #self.n_tcmax = (np.ceil(self.nt/self.n_cluster)).astype(int)
        self.n_tcmax = round(self.nt/(self.n_cluster*self.noss))
        
        # ----- set default obj. values -----
        self.aep = 0
        self.obj_value = None
        
        
        print("setting up areas/geometry")
        # ----- Lease area boundary polygon -----
        #self.boundary = boundary_coords#list(zip(boundary_x, boundary_y))
        self.boundary_coords = getFromDict(kwargs, 'boundary_coords', shape = -1, default = np.array([(0, 0),(10000, 0),(10000, 10000), (0,10000) ]))
        self.setBoundary(self.boundary_coords[:,0], self.boundary_coords[:,1])
        self.boundary_sh = sh.Polygon(self.boundary)
        
        # set up any interior sub boundaries (useful for mulitple separate uniform grids)
        self.sub_boundary_coords = getFromDict(kwargs, 'sub_boundary_coords', shape=-1, 
                                               default = [])
        self.sub_boundary = []
        self.sub_boundary_sh = []
        self.sub_boundary_centroid = []
        self.sub_boundary_centroid_x = []
        self.sub_boundary_centroid_y = []
        for subb in self.sub_boundary_coords:
            subb = np.array(subb)
            # save as project sub boundaries
            self.sub_boundary.append(np.vstack([[subb[i,0],subb[i,1]] for i in range(len(subb))]))
            
            # if the boundary doesn't repeat the first vertex at the end, add it
            if not all(subb[0,:] == subb[-1,:]):
                self.sub_boundary[-1] = np.vstack([self.sub_boundary[-1], subb[0,:]])
            # create sub boundary shapely polygon    
            self.sub_boundary_sh.append(sh.Polygon(self.sub_boundary[-1]))
            # create sub boundary centroid and store centroid coords
            self.sub_boundary_centroid.append(self.sub_boundary_sh[-1].centroid)
            self.sub_boundary_centroid_x.append(self.sub_boundary_centroid[-1].x)
            self.sub_boundary_centroid_y.append(self.sub_boundary_centroid[-1].y)
            
        # trim the bathymetry grid to avoid excess
        self.trim_grids = getFromDict(kwargs,'trimGrids',default=True)
        if self.trim_grids:
            self.trimGrids()
        
        # Get centroid of lease area
        self.boundary_centroid = self.boundary_sh.centroid
        self.boundary_centroid_x, self.boundary_centroid_y = self.boundary_centroid.x, self.boundary_centroid.y       
        
        # Safety margin
        self.boundary_margin =  getFromDict(kwargs, 'boundary_margin', default = 0) #margin applied to exterior boundary
        # Calculate the buffered polygon with safety margin
        # internal: safety margin, to ensure that there is enough space for mooring system without crossing lease area boundaries
        # idea: Safety margin dependent on water depth?
        self.boundary_sh_int = self.boundary_sh.buffer(-self.boundary_margin)
        #self.grid_x_min, self.grid_y_min, self.grid_x_max, self.grid_y_max = self.boundary_sh_ext.bounds
        
        '''
        # Calculate the total area of the buffered polygon
        self.total_area_ext = self.boundary_sh_ext.area
        self.total_area_int = self.boundary_sh_int.area
        self.a_f=round(self.total_area_ext/self.total_area_int)
 
        # Maximum x and y distance in boundary shape    
        def max_distance(x_values):
          if len(x_values) < 2:
              return 0
          x_max = max(x_values)
          x_min = min(x_values)
          return abs(x_max - x_min)

        self.bound_dist_x = max_distance(self.boundary[:,0])  
        self.bound_dist_y = max_distance(self.boundary[:,1])
        '''
            
        # INTIAL: Parse the design vector and store updated positions internally
        #x_pos, y_pos = X[:len(X)//2], X[len(X)//2:]
        # ONLY FOR FREE OPTIMIZATION
        if np.size(Xu) == 0 and np.size(Xdb) == 0:
            if self.rotation_mode:
                x_pos, y_pos, rot_rad = X[:self.nt], X[self.nt:2*self.nt], X[2*self.nt:]
                #self.turb_rot_deg = rot_deg
                self.turb_rot= rot_rad#np.radians(rot_deg) 
            else:
                x_pos, y_pos = X[:len(X)//2], X[len(X)//2:]

                self.turb_rot = getFromDict(kwargs, 'turb_rot', shape = self.nt, default = np.zeros(self.nt))#rot_rad#np.radians(turb_rot)
    
            # Turbine positons: INPUT in [m]
            self.turb_coords[:,0]= x_pos
            self.turb_coords[:,1]= y_pos
        # UNIFORM GRID LAYOUT OPTIMIZATION
        else:
            self.turb_rot = getFromDict(kwargs, 'turb_rot', shape = self.nt, default = np.zeros(self.nt))
            self.turb_coords = np.zeros((self.nt,2))
            
           # X = self.generateGridPoints(Xu) 
            #x_pos, y_pos = X[:len(X)//2], X[len(X)//2:]
            #self.turb_coords[:,0]= x_pos
            #self.turb_coords[:,1]= y_pos
       
        # ----- Exclusion zone polygons -----
        # Turbine distances to exclusion zones (if any) 
        self.exclusion = getFromDict(kwargs, 'exclusion_coords', shape = -1, default = []) 
        self.turb_dist_tez1_m = np.zeros((self.nt*len(self.exclusion)))
        self.turb_dist_tez2_m = np.zeros((self.nt*len(self.exclusion)))
        self.exclusion_polygons_sh = []  # List to store polygons   
        
        # Create exclusion polygons 
        for ie in range(len(self.exclusion)):
            exclusion_polygon = sh.Polygon(self.exclusion[ie])
            self.exclusion_polygons_sh.append(exclusion_polygon)
        
        
        # ----- Wind data -----
        self.wind_rose = wind_rose
        
        
        # ----- Mooring system variables -----
        print("setting up mooringList")   
        self.mooringList = makeMooringListN(ss, 3*self.nt)  # make Moorings

        for mooring in self.mooringList.values():   # hackily set them up
            mooring.dd['sections'] = []
            mooring.dd['connectors'] = []
            for i,sec in enumerate(mooring.ss.lineList):
                mooring.dd['connectors'].append({'CdA':0,'m':0,'v':0})
                mooring.dd['sections'].append({'type':mooring.ss.lineList[i].type,
                                               'L':mooring.ss.lineList[i].L})
            mooring.dd['connectors'].append({'CdA':0,'m':0,'v':0})
            mooring.adjuster = mooringAdjuster   # set the designer/adjuster function
        
        
        # ----- Platforms ----- 
        for i in range(self.nt):
            r = [self.turb_coords[i][0],self.turb_coords[i][1],0]
            self.platformList[i] = Platform(id=i, r=r, heading=0, mooring_headings=[0,120,240],rFair=ss.rad_fair,zFair=ss.z_fair)
            self.platformList[i].entity = 'FOWT'
            
            for j in range(3):
                self.platformList[i].attach(self.mooringList[i*3+j], end='b')
        
        
        # ---- Anchors ----
        self.anchorList = {}
        if 'anchor_settings' in kwargs:
            anchor_settings = True
        else:
            anchor_settings = False
        # set up anchor design dictionary 
        ad = {'design':{}, 'cost':{}} 
        if anchor_settings and 'anchor_design' in kwargs['anchor_settings']:
            anchor_design_initial = kwargs['anchor_settings']['anchor_design']
            ad['type'] = kwargs['anchor_settings']['anchor_type']
        else:
            print('No anchor type given, defaulting to suction bucket anchor.')
            anchor_design_initial = {'D':3.0,'L':16.5,'zlug':10}
            ad['type'] = 'suction_pile'
        ad['design'] = anchor_design_initial # INPUT or not???
        for i, moor in enumerate(self.mooringList.values()):
            if self.soil_x is not None: # get soil conditions at anchor location if soil info available
                name, props = self.getSoilAtLocation(moor.rA[0], moor.rA[1])

            # create anchor object
            anch = Anchor(dd=ad,aNum=i,id=moor.id)
            anch.soilProps = {name:props}
            self.anchorList[anch.id] = anch
            # attach to mooring line
            moor.attachTo(anch,end='a')
            if 'mass' in ad:
                anch.mass = ad['mass']               
            elif anchor_settings and 'mass' in kwargs['anchor_settings']:
                anch.mass = kwargs['anchor_settings']['mass']
            
        
        # --- develop anchor types ---
        self.anchorTypes = {}
        self.anchorMasses = {}
        self.anchorCosts = {}
        # pull out mean depth
        meandepth = np.mean(-self.grid_depth)
        pf = self.platformList[0]
        # artificially set platform at 0,0
        pf.setPosition([0,0],project=self) # put in a random place and reposition moorings
        # create ms for this platform
        msPF = pf.mooringSystem()
        # set depth artificially to mean depth
        msPF.depth = -meandepth
        # set mooring object depth artificially for now
        for moor in pf.getMoorings().values():
            moor.dd['zAnchor'] = meandepth
            moor.z_anch = meandepth
            moor.ss.depth = -meandepth
            moor.rad_fair = 58
            moor.z_fair = -14
        # call set position function again to use adjuster function on all moorings
        pf.setPosition([0,0])
        
        # get anchors connected to this platform
        anchors = pf.getAnchors()
        # choose one (all should be same)
        anch = anchors[0]

        # keep zlug constant?
        if anchor_settings and 'fix_zlug' in kwargs['anchor_settings']:
            fix_zlug=kwargs['anchor_settings']['fix_zlug']
        else:
            fix_zlug=False
        # set minimum allowable FS
        if anchor_settings and 'FS_min' in kwargs['anchor_settings']:
            minfs = kwargs['anchor_settings']['FS_min']
        else:
            minfs = {'Ha':2,'Va':2}
        # set FSdiff_max if provided
        if anchor_settings and 'FSdiff_max' in kwargs['anchor_settings']:
            FSdiff_max = kwargs['anchor_settings']['FSdiff_max']
        else:
            FSdiff_max = None
        # create anchor for each soil type
        for name, soil in self.soilProps.items():
            if anchor_settings and 'anchor_resize' in kwargs['anchor_settings'] and kwargs['anchor_settings']['anchor_resize']:
                # get anchor forces from array watch circle
                pf.getWatchCircle(ms = msPF)
                # get loads dictionary but get rid of any Ha Va loads that might already be there
                anch.loads = {'Hm':anch.loads['Hm'],'Vm':anch.loads['Vm'],
                              'thetam':anch.loads['thetam'], 'mudline_load_type':'max'}
                # update soil type for anchor
                anch.soilProps = {name:soil}
                geom = [val for val in anch.dd['design'].values()]
                geomKeys = [key for key in anch.dd['design'].keys()]
                anch.getSize(geom,geomKeys, FSdiff_max=FSdiff_max,
                             fix_zlug=fix_zlug, minfs=minfs)
                
            self.anchorTypes[name] = deepcopy(anch.dd) if anch.dd else {}
            self.anchorMasses[name] = deepcopy(anch.mass) if anch.mass else 0
            try:
                self.anchorCosts[name] = deepcopy(anch.getCost())
            except:
                self.anchorCosts[name] = 0
                       
                
                
        self.ms_na = 3  # Number of anchors per turbine. For now ONLY 3 point mooring system.
        #self.ms_anchor_depth = np.zeros((self.nt*self.ms_na))  # depths of anchors
        self.anchor_coords= np.zeros((self.nt*self.ms_na,2))  # anchor x-y coordinate list [m]
        self.ms_bufferzones_pos = np.zeros((self.nt,), dtype=object)    # Buffer zones for moorign system
        self.ms_bufferzones_rout = np.zeros((self.nt,), dtype=object)  
        self.ms_bufferzones_rout_points = np.zeros((self.nt,), dtype=object)


        # ----- Initialize the FLORIS interface fi -----
        self.use_FLORIS = getFromDict(kwargs,'use_FLORIS', default = False)
        if self.use_FLORIS:  # If using FLORIS, initialize it
            print("initializing FLORIS")
            # How to do this more elegant?
            dirname = '' #'./_input/'
            #flName = 'gch_floating.yaml'
            if self.parallel:
                from floris import ParFlorisModel
                self.floris_file = getFromDict(kwargs, 'floris_file', dtype = str, default = '')
                self.flow = ParFlorisModel(self.floris_file)
                
            else:
                self.floris_file = getFromDict(kwargs, 'floris_file', dtype = str, default = '')
                self.flow = FlorisModel(self.floris_file) #FlorisInterface
            
            # FLORIS inputs x y positions in m
            self.flow.set(layout_x=self.turb_coords[:,0],
                                   layout_y=self.turb_coords[:,1],
                                   wind_data = self.wind_rose
                                   )
            #run floris simulation
            # self.flow.run()
            
            # # SAVE INITIAL AEP
            # self.aep0 = self.flow.get_farm_AEP()

            # ----- Wind Turbine Data ----- 
            # https://nrel.github.io/floris/turbine_interaction.html
            # self.ti = TurbineInterface.from_internal_library("iea_15MW.yaml")

            if self.display > 0:
                self.plotWakes(wind_spd = 10, wind_dir = 270, ti = 0.06)

        else:  # if not using FLORIS, indicate it with a None
            self.flow = None

        print("updating layout")
        if np.size(Xu) != 0:
            self.updateLayoutUG(Xu)
        elif np.size(Xdb) != 0:
            self.db_ext_spacing = getFromDict(kwargs, 'db_ext_spacing', default = [0, 1, 0, 1])
            self.updateLayoutDB(Xdb)
        else:
            self.updateLayoutOPT(X)
        

      


    def generateGridPoints(self, Xu, trans_mode, boundary_index=-1):
        ''' Generate uniform grid points and save resulting coordinates into vector X.
            This transforms the uniform grid (UG) design variables into the design variables of 
            the free layout optimization.
            
            trans_mode = 'x': Shear transformation in x direction only
            trans_mode = 'xy': Shear transformation in x and y direction            
        '''
        grid_spacing_x = Xu[0]
        grid_spacing_y = Xu[1]
        grid_trans_x = Xu[2]
        grid_trans_y = Xu[3]
        grid_rotang = Xu[4]  
        grid_skew = Xu[5]  
        
        if boundary_index >= 0:
            boundary = self.sub_boundary_sh[boundary_index]
            bound_centroid_y = self.sub_boundary_centroid_y[boundary_index]
            bound_centroid_x = self.sub_boundary_centroid_x[boundary_index]

        else:
            boundary = self.boundary_sh
            bound_centroid_y = self.boundary_centroid_y
            bound_centroid_x = self.boundary_centroid_x
        
        if self.rotation_mode:
            if len(Xu) != 7:
                raise ValueError('If rotation mode is True, Xu[6] is turbine rotation')
            self.turb_rot = np.radians(Xu[6])
        
        # Check if self.grid_spacing_x/y is equal to 0, if so, set it to 1000 m
        if grid_spacing_x == 0:
            grid_spacing_x = self.turb_minrad*0.5
        if grid_spacing_y == 0:
            grid_spacing_y = self.turb_minrad*0.5
            
        # Shear transformation
        # Calculate trigonometric values
        cos_theta = np.cos(np.radians(grid_rotang))
        sin_theta = np.sin(np.radians(grid_rotang))
        tan_phi = np.tan(np.radians(grid_skew))
    
        # Transmoration matrix, considering shear transformatio and rotation
        # Default: shear direction in x direction only
        # xy: shear direction in x and direction
        if trans_mode == 'xy':
            # Compute combined x and y shear
            transformation_matrix = np.array([[cos_theta-sin_theta*tan_phi, -sin_theta + tan_phi * cos_theta],
                                               [sin_theta+cos_theta*tan_phi, sin_theta*tan_phi+cos_theta]])
        else:
            # default transformation: x shear only
            transformation_matrix = np.array([[cos_theta, -sin_theta + tan_phi * cos_theta],
                                               [sin_theta, sin_theta*tan_phi+cos_theta]])
        
        # Generate points in the local coordinate system
        points = []
        
        # Lease area shape: Get min and max xy coordinates and calculate width
        min_x, min_y, max_x, max_y = boundary.bounds # self.boundary_sh.bounds
        xwidth = abs(max_x-min_x)
        ywidth = abs(max_y-min_y)
        
        
        # LOCAL COORDINATE SYSTEM WITH (0,0) LEASE AREA CENTROID
        # Therefore, +/- self.boundary_centroid_y/x cover the entire area
        # Loop through y values within the boundary_centroid_y range with grid_spacing_y increments
        column_count = 0
        rotations = []
        grid_position =[]
        for y in np.arange(-bound_centroid_y-ywidth, bound_centroid_y+ywidth, grid_spacing_y):
            column_count += 1
            row_count = 0
            # Loop through x values within the boundary_centroid_x range with grid_spacing_x increments
            for x in np.arange(-bound_centroid_x-xwidth, bound_centroid_x+xwidth, grid_spacing_x):
                
                row_count += 1
                # Apply transformation matrix to x, y coordinates
                local_x, local_y = np.dot(transformation_matrix, [x, y])
                # Add grid translation offsets to local coordinates
                local_x += grid_trans_x
                local_y += grid_trans_y
                # Create a Point object representing the transformed coordinates
                # Transform back into global coordinate system with by adding centroid to local coordinates
                point = Point(local_x + bound_centroid_x, local_y + bound_centroid_y)
                points.append(point)
                
                if self.alternate_rows:
                    rotations.append(self.turb_rot + np.radians(180 * (column_count % 2)))
                #store column, row for each turbine
                grid_position.append([column_count, row_count])
        
         
        # remove points that are not in boundaries
        bound_lines = boundary.boundary # get boundary lines for shapely analysis
        out_lines = [bound_lines]
        # keep only points inside bounds
        points_ib = [pt for pt in points if (boundary.contains(pt))]
        if self.alternate_rows:
            self.turb_rot = [rotations[ind] for ind in range(0, len(points)) if boundary.contains(points[ind])]
        self.grid_positions = [grid_position[ind] for ind in range(0, len(points)) if boundary.contains(points[ind])]
        
        points_ibe = points_ib
        # remove points in exclusion zones
        if self.exclusion_polygons_sh:
            for ie in range(len(self.exclusion)):
                points_ibe = [pt for pt in points_ibe if not self.exclusion_polygons_sh[ie].contains(pt)]
                out_lines.append(self.exclusion_polygons_sh[ie].boundary) # get boundary lines for exclusion zones
        
        return(points_ibe)

  
    def pareGridPoints(self,points_ibe):
        '''
        Function to pare number of grid points down to desired amount, place oss 
        at closest grid points (if substations allowed to move) and return 
        array of x, y(, rotation) values. Sorts points by distance from all borders
        (lease boundary, inner boundaries, exclusion zones) and keeps the nt points
        furthest from all boundaries

        Parameters
        ----------
        points_ibe : list
            List of shapely point objects that are inside all boundaries and outside all exclusion zones

        Returns
        -------
        X : np.ndarray
            1D array of concatenated x, y(, rotation) for each turbine

        '''
        # determine number of points to keep (usually # turbines + # substations)
        if self.static_substations:
            # in this case, keep substations where they are
            nt = self.nt 
        else:
            nt = self.nt + self.noss
        
        # create list of boundary lines from outside boundary, exclusion zones, and inner boundaries
        out_lines = [self.boundary_sh.boundary]
        if len(self.sub_boundary_sh) > 0:
            for sub in self.sub_boundary_sh:
                out_lines.append(sub.boundary)
        for ie in range(len(self.exclusion)):
            out_lines.append(self.exclusion_polygons_sh[ie].boundary)
        
        lines = MultiLineString(out_lines)
        point_dists = [pt.distance(lines) for pt in points_ibe] # get min dist between bounds and each point
        points_ibe = np.array(points_ibe)
        # get indices of sorting by descending minimum distance
        points_sorted_idx = [int(ind) for ind in np.flip(np.argsort(point_dists,kind='stable'))]
        furthest_points = list([points_ibe[i] for i in range (0, len(points_ibe)) if i in points_sorted_idx[:nt]]) # pull out the points that are furthest from bounds
        self.grid_positions = list(self.grid_positions[i] for i in range (0, len(points_ibe)) if i in points_sorted_idx[:nt])
        if self.alternate_rows:
            furthest_rotations = list(self.turb_rot[i] for i in range (0, len(points_ibe)) if i in points_sorted_idx[:nt]) 

        
        # add points outside lease area if more points are needed
        min_x, min_y, max_x, max_y = self.boundary_sh.bounds
        if len(points_sorted_idx)< nt:
            # determine remaining number of turbines to add
            leftover = nt-len(points_sorted_idx)
            # choose point outside bounds for leftovers
            leftover_loc = Point(min_x-1,min_y-1)
            furthest_points.extend([leftover_loc]*leftover)
            if self.alternate_rows:
                furthest_rotations.extend([0]*leftover)
            
        # put substation(s) in place closest to oss_coords if substations can move
        if not self.static_substations:
            for oo in range(self.noss):
                # make a multipoint fro
                turb_multipoint = sh.MultiPoint(furthest_points)
                oss_point_start = Point(self.oss_coords_initial[oo])
                # find point closest to initial oss coord & set as new oss position
                oss_point = nearest_points(turb_multipoint,oss_point_start)[0]
                # remove turbine from new oss position (extra turbines have been placed already)
                if oss_point in furthest_points:
                    if self.alternate_rows:
                        del furthest_rotations[furthest_points.index(oss_point)]
                        
                    index = furthest_points.index(oss_point)
                    furthest_points.remove(oss_point)
                    self.grid_positions.remove(self.grid_positions[index])
                    self.oss_coords[oo] = [oss_point.x, oss_point.y]
                else:
                    print('Could not find nearby point for oss, setting oss to initial coords')
                    self.oss_coords[oo] = self.oss_coords_initial[oo]
                    
        # save points furthest from bounds into turb_coords
        x_coords = np.array([point.x for point in furthest_points])#/1000 
        y_coords = np.array([point.y for point in furthest_points])#/1000   
        for i,coord in enumerate(self.turb_coords):
            coord[0] = x_coords[i]
            coord[1] = y_coords[i]
        
        #update grid_positions row and column coordinates based on minimum
        self.grid_positions = np.array(self.grid_positions)
        self.grid_positions[:,0] = self.grid_positions[:,0] - min(self.grid_positions[:,0])
        self.grid_positions[:,1] = self.grid_positions[:,1] - min(self.grid_positions[:,1])
        
        # Return Design Vector X with x,y coordinates, same as used for the free layout optimization.
        # Coordinates in (km)
        # This completes the interface
        if self.rotation_mode:
            
            if self.alternate_rows:
                self.turb_rot = furthest_rotations
                X = np.concatenate((self.turb_coords[:,0], self.turb_coords[:,1], self.turb_rot))
            else:
                X = np.concatenate((self.turb_coords[:,0], self.turb_coords[:,1], self.turb_rot*np.ones((nt))))
        else:
            X = np.concatenate((self.turb_coords[:,0], self.turb_coords[:,1]))

        return X
    

    def updateLayout(self, X, level=0, refresh=False):
        '''Update the layout based on the specified design vector, X. This
        will adjust the turbine positions stored in the Layout object as 
        well as those in the FLORIS and any other sub-objects.
        
        Parameters
        ----------
        X
            Design vector.
        level
            Analysis level to use. Simplest is 0.
        refresh : bool
            If true, forces a re-analysis, even if this design vector is old.
        '''
        if len(X)==0: # if any empty design vector is passed (useful for checking constraints quickly)
            if refresh:
                X = np.array(self.Xlast)
            else:
                return
        if np.array_equal(X, self.Xlast) and not refresh:
        #if all(X == self.Xlast) and not refresh:  # if X is same as last time
            #breakpoint()
            pass  # just continue, skip the update steps
        
       
        elif any(np.isnan(X)):
            raise ValueError("NaN value found in design vector")

        else:  # Update things iff the design vector is valid and has changed
            if self.display > 1:
                print("Updated design")
               
            self.iter += 1  # update internal iteration counter
            
            # Parse the design vector and store updated positions internally
            if self.rotation_mode:
                x_pos, y_pos, rot_rad = X[:self.nt], X[self.nt:2*self.nt], X[2*self.nt:]
                #self.turb_rot = np.radians(rot_deg)
                self.turb_rot = rot_rad
            else:
                x_pos, y_pos = X[:len(X)//2], X[len(X)//2:]
                #self.turb_rot = self.turb_rot_const
            
            self.turb_coords[:,0]= x_pos
            self.turb_coords[:,1]= y_pos
             
            # Update things for each turbine
            #breakpoint()
            #print(self.nt, len(self.turb_depth), X)
            #print(self.turb_coords)
          
            # Update Paltform class
            for i in range(self.nt):
                self.platformList[i].setPosition(self.turb_coords[i], heading=self.turb_rot[i], degrees=False, project = self)  
                # switch anchor type
                anchs = self.platformList[i].getAnchors() 
                for anch in anchs.values():
                    name, props = self.getSoilAtLocation(anch.r[0],anch.r[1])
                    atype = self.anchorTypes[name]
                    anch.dd.update(atype)
                    anch.mass = self.anchorMasses[name]
                    anch.cost['materials'] = self.anchorCosts[name]
                    anch.soilProps = {name:props}
                       
                # Get depth at turbine postion          
                self.turb_depth[i] = -self.getDepthAtLocation(
                          self.turb_coords[i,0], self.turb_coords[i,1])
            # update substation platform location(s)
            for oo in range(self.noss):
                self.platformList[self.nt+oo].setPosition(self.oss_coords[oo],
                                                          heading=self.turb_rot[i],
                                                          degrees=False,project=self)
                
                
            # Anchor locations - to be repalced when integration is further advanced
            #for j in range(3):
            #    im = i*3 + j  # index in mooringList
            #    self.ms_anchor_depth[im] = self.mooringList[im].z_anch#self.mooringList[im].rA[2]  OLD, not needed anymore
            #    self.anchor_coords[im,:] = self.mooringList[im].rA[:2]
                
            '''
            # Calculate anchor position based on headings
            
            theta = self.turb_rot[i] #  turbine heading
            #R = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
        
            headings = np.radians([60,180,300])
            
            for j in range(len(headings)):
            
                im = i*3 + j  # index in mooringList
                
                # heading of the mooring line
                heading_i = headings[j] + theta
                
                # adjust the whole Mooring
                #self.mooringList[im].reposition(self.turb_coords[i,:], 
                #                      heading=heading_i, project=self, level=level)
                #self.mooringList[im].reposition(r_center=self.turb_coords[i,:], 
                #                      heading=heading_i, project=self, level=level)
              '''  
                # get the anchor location from the mooring
                #self.anchor_coords[im,:] = self.mooringList[im].rA[:2]
                #self.ms_anchor_depth[im] = self.mooringList[im].rA[2] 
                #self.mooringList[0].z_anch
                #self.mooringList[0].rA
                #self.mooringList[0].rA
                    
                    
                    
            # ----- evaluate constraints -----
            
            # ----- Calculate buffer zone shape around mooring lines and anchors. -----    
            # ISO 19901-7: 100 m safety zone to other offshore assets, therefore 50 m per mooring line is recommended

            # SAFE BUFFERZONES IN PLATFORM OBJECT?
            
                
            
      
            # Create LineString geometries and buffer them
            for i in range(self.nt):
                # Buffer group for turbine positioning
                buffer_group_pos = []
                # Buffer group for cable routing
                buffer_group_rout = []
                
                for j in range(self.ms_na):
                    im = 3*i + j  # global index of mooring/anchor
                                  
                    moor_bf_start = get_point_along_line(self.turb_coords[i,:], self.mooringList[im].rA[:2], self.turb_minrad)
                    # Buffer zone mooring line
                    #line = LineString([self.turb_coords[i,:], self.mooringList[im].rA[:2]])
                    line = LineString([moor_bf_start, self.mooringList[im].rA[:2]])
                    mooringline_buffer = line.buffer(self.moor_minrad)
                    
                    # Buffer zone anchor
                    # Create a point at coordinates (x, y)
                    point = Point(self.mooringList[im].rA[:2])
                    # Create a buffer around the anchor with a radius of X
                    anchor_buffer = point.buffer(self.anchor_minrad)
                    
                    # Buffer zone turbine
                    # Create a point at coordinates (x, y)
                    point = Point(self.turb_coords[i,:],)
                    # Create a buffer around the anchor with a radius of X
                    turb_buffer = point.buffer(self.turb_minrad)
                    
                    # Buffer group for turbine positioning
                    buffer_group_pos.append(mooringline_buffer)
                    buffer_group_pos.append(anchor_buffer)
                    buffer_group_pos.append(turb_buffer)
                    
                    # Buffer group for cable routing
                    buffer_group_rout.append(mooringline_buffer)
                    buffer_group_rout.append(anchor_buffer)
                               
                # Combine the buffered lines connected to the same turbine into one polygon
                polygon = unary_union(buffer_group_pos)  # Combine buffers for each turbine
                if isinstance(polygon, MultiLineString):
                    # Convert MultiLineString to Polygon
                    polygon = Polygon(polygon)
                self.ms_bufferzones_pos[i] = polygon
                
                polygon = unary_union(buffer_group_rout)  # Combine buffers for each turbine
                if isinstance(polygon, MultiLineString):
                    # Convert MultiLineString to Polygon
                    polygon = Polygon(polygon)
                self.ms_bufferzones_rout[i] = polygon
                
                #envelopes['buffer_zones']['shape']
 
 

            # ----- Overlap between mooring zones -----
            # Create an empty 2D array to store the areas of intersection
            intersection_areas = np.zeros((self.nt, self.nt))
            # Calculate and fill the array with the areas of intersection
            for i in range(self.nt):
                for j in range(i + 1, self.nt):
                    polygon1 = self.ms_bufferzones_pos[i]
                    polygon2 = self.ms_bufferzones_pos[j]
                    # Calculate intersection
                    intersection = polygon1.intersection(polygon2)
                    # Fill the array with the area of intersection
                    intersection_areas[i, j] = intersection.area*(-1)
                    intersection_areas[j, i] = intersection.area*(-1)
            
            self.con_moor_moor = getLower(intersection_areas) # get lower diagonal
            
            
            # ----- Overlap between mooring zones and boundary -----        
            # Calculate areas of the parts of polygons outside the boundary
            self.con_moor_boundary = np.zeros(self.nt) 
            # Iterate over polygons and fill the array with areas
            for i, polygon in enumerate(self.ms_bufferzones_pos):
                if isinstance(polygon, (Polygon, MultiPolygon)) and polygon.intersects( self.boundary_sh):
                    # Calculate the intersection with the boundary polygon
                    intersection = polygon.intersection( self.boundary_sh)
                    # Calculate the area of the parts outside the boundary
                    outside_area = polygon.difference(intersection).area
                    # Fill the array with the area
                    self.con_moor_boundary[i] = -outside_area
            

            # ----- Between exclusion zones and turbines ----- 
            self.con_moor_ez_m2 = np.zeros((self.nt*len(self.exclusion)))
            r = 0
            for ie in range(len(self.exclusion)):
                #exclusion_polygon = sh.Polygon(self.exclusion[ie])                
                # Iterate over polygons and fill the array with areas
                for i, polygon in enumerate(self.ms_bufferzones_pos):
                    if isinstance(polygon, (Polygon, MultiPolygon)) and polygon.intersects(self.exclusion_polygons_sh[ie]):
                        # Calculate the intersection with the boundary polygon
                        intersection = polygon.intersection(self.exclusion_polygons_sh[ie])
                        # Calculate the area of the parts inside exclusion areas
                        inside_area = polygon.difference(intersection).area
                        # Fill the array with the area
                        self.con_moor_ez_m2[r] = -inside_area
                        r += 1


            # ----- Margin between turbines -----

            # Distance matrix between turbines
            distances = cdist(self.turb_coords, self.turb_coords)
            
            dists = distances [np.tril_indices_from( distances , k=-1)]  # get lower diagonal
            
            # Reduce by buffer radius (for each turbine) then store
            self.con_turb_turb = dists - 2*self.turb_minrad
            
            # ----- Margin between turbines and OSS -----
            # Distance matrix between turbines
            #distances = cdist(self.turb_coords,self.oss_coords)
            #print(self.oss_coords)
            #print(self.turb_coords)
            #dists = distance.cdist(self.oss_coords, self.turb_coords, 'euclidean')
            #dists = distances [np.tril_indices_from( distances , k=-1)]  # get lower diagonal
            dists = []
            for oo in range(self.noss):
                dists.extend(np.linalg.norm(self.turb_coords - self.oss_coords[oo], axis=1))
            # Reduce by buffer radius (for each turbine) then store
            self.con_turb_oss = np.array(dists) - self.oss_minrad

            # ----- margin between turbines and lease area boundary -----
            r = 0
            self.con_turb_ez_m = np.zeros((self.nt*len(self.exclusion)))
            self.con_oss_boundary = np.zeros(self.noss)
            self.con_oss_ez_m = np.zeros((self.noss*len(self.exclusion)))
            coords = np.zeros((self.nt+self.noss,2))
            coords[:self.nt] = self.turb_coords
            coords[self.nt:] = self.oss_coords
            isturb=True
            for i in range(self.nt+self.noss):
                if i>=self.nt:
                    isturb = False
                # Create a Shapely Point for the given xy of turbine or oss
                p_turb = Point(coords[i,0], coords[i,1])
                
                #breakpoint()
                # Find the nearest point on the shape to the given point
                p_bound = nearest_points(self.boundary_sh.exterior, p_turb)[0]
                
                # Calculate the Euclidean distance between point and nearest point on boundary
                distance_within = p_turb.distance(p_bound)
                
                # If point is outside boundary, give the distance a negative sign
                if not p_turb.within(self.boundary_sh):
                    distance_within = -abs(distance_within)

                # Reduce by buffer radius, then add to constraint list
                if isturb:
                    self.con_turb_boundary[i] = distance_within - self.turb_minrad
                else:
                    self.con_oss_boundary[i-self.nt] = distance_within - self.oss_minrad
                
                for ie in range(len(self.exclusion)):
                    p_exclusion = nearest_points(self.exclusion_polygons_sh[ie].exterior, p_turb)[0]
                    dist_outside = p_turb.distance(p_exclusion)
                    # if turbine is inside exclusion zone, give distance - sign
                    if p_turb.within(self.exclusion_polygons_sh[ie]):
                        dist_outside = -abs(dist_outside)
                        
                    if isturb:                    
                        self.con_turb_ez_m[r] = dist_outside
                    else:
                        self.con_oss_ez_m[r-self.nt*len(self.exclusion)] = dist_outside
                    r += 1
                
                # could handle exclusion zones in this same loop
            # # ----- margin between turbines and exclusion zones -----
            # # Optimize: creat point once, together with above
            # if len(self.exclusion) > 0:
            #     r = 0
            #     for ie in range(len(self.exclusion)):
            #         #exclusion_polygon = sh.Polygon(self.exclusion[ie])   
            #         #breakpoint()
                    
            #         for i in range(self.nt):
            #             # Create a Shapely Point for the given xy
            #             #point = Point(x_pos[i], y_pos[i])
            #             point = Point(self.turb_coords[i,0], self.turb_coords[i,1])
            #             # Find the nearest point on the shape to the given point
            #             nearest_point = nearest_points(self.exclusion_polygons_sh[ie].exterior, point)[0]
            #             # Calculate the Euclidean distance between WTG anchor radius and nearest point on shape
            #             # Reduce distance by radius (distance has to be equal or greater than anchor radius)
            #             self.turb_dist_tez1_m[r] = point.distance(
            #                 nearest_point) - self.turb_minrad
            #             # Calculate the Euclidean distance between WTG center and shape
            #             self.turb_dist_tez2_m[r] = point.distance(nearest_point)
            #             # Check if turbine is outside the boundary
            #             # Ensure if point is outside shape, distance is always negative
            #             if point.within(self.exclusion_polygons_sh[ie]):
            #                 self.turb_dist_tez1_m[r] = -abs(self.turb_dist_tez1_m[r])
            #                 # Weight the contraints so that the turbines stay within the specifified area
            #                 self.turb_dist_tez2_m[r] = -abs(self.turb_dist_tez2_m[r])
            #             r =+1

                
            
                     
            
            
            # ----- Concatenate constraints vector -----
            
            # Note: exclusions are temporarily skipped, but can be added back in to the below
            
            #!! QUESTION MB: Should this be considered at all as a constraint? I think it is more important that
            #               anchor buffer zones do not exceed the lease boundaries, but not a wind turbine spacing parameter.
            
            # distances
            constraint_vals_m = np.concatenate([self.con_turb_turb, self.con_turb_boundary, 
                                                self.con_turb_oss, self.con_turb_ez_m,
                                                self.con_oss_boundary, self.con_turb_ez_m]) 
            constraint_vals_km = constraint_vals_m/1000
     
            # areas
            constraint_vals_m2 = np.concatenate([self.con_moor_moor, self.con_moor_boundary, self.con_moor_ez_m2])            
            constraint_vals_km2 = constraint_vals_m2/(1000**2)
            # Combine constraint values (scaling to be around 1)
            self.con_vals = 10*np.concatenate([constraint_vals_km, constraint_vals_km2])


            # Sum of Constraint values
            negative_values = [val for val in self.con_vals if val < 0]
            
            
            if not negative_values:
                self.con_sum = 0
                # ----- Cable Layout - ONLY FOR FEASIBLE LAYOUT  
                if self.cable_mode:

                    self.iac_dic,_,_ = getCableLayout(self.turb_coords, self.oss_coords, self.iac_typical_conductor, 
                                                      self.iac_type, self.turb_rating_MW, turb_cluster_id=[], 
                                                      n_cluster_sub=self.n_cluster, n_tcmax=self.n_tcmax, plot=False, oss_rerouting=1)

                    # Save cables in cable objects
                    self.addCablesConnections(self.iac_dic,cableType_def=self.iac_type)    
                 
            else:
                self.con_sum = sum(negative_values) # sum of all values below zero
            

            if self.optimizer == 'PSO':            
                # PSO constraints only
                # Constraints above zero 0: satisfied (often it is g < 0 for satisfied constraints for a PSO)
                # Solution: Sum of negative constraint values, because it has to be one value only 
                self.con_vals = self.con_sum
            
              
            # Penalty factor: (1+abs(self.con_vals)) or 1 
            if self.obj_penalty == 1: # penalty ON
                f_pentalty = (1+abs(self.con_sum)) 
            else: # penalty OFF
                f_pentalty = 1
            

            # ----- evaluate objective function -----
            # compute the objective function value
            # objective function includes a constraint term, leading to a penalty when constraints are not satisfied
            #       (1+abs(self.con_vals))
            # objective funciton 
            if not negative_values or self.infeasible_obj_update or not self.obj_value:
                if self.mode == 'LCOE':  # minimize LCOE (this LCOE version focuses on mooring and cable costs/AEP)
                    self.getLCOE()
                    #self.constraintFuns_penalty(X)                
                    self.obj_value = self.lcoe*1e5*f_pentalty #(1+abs(self.con_vals)) #+ self.cost_penalty / self.aep#self.getLCOE() #+ self.constraintFuns_penalty(X)/self.aep
                elif self.mode == 'LCOE2': # minimize LCOE (this LCOE version includes opex estimates and platform/turbine cost estimates)
                    self.getLCOE2()
                    self.obj_value = self.lcoe*f_pentalty 
                elif self.mode == 'AEP':  # maximize AEP
                    self.getAEP(display = self.display) 
                    self.obj_value = -self.aep/1e12/f_pentalty #-self.getAEP() #+ self.constraintFuns_penalty(X) # minus, because algorithm minimizes the objective function
                elif self.mode == 'CAPEX':  # maximize AEP
                    self.getCost()
                    #self.constraintFuns_penalty(X)
                    self.obj_value = (self.cost_total/1e7)*f_pentalty #+ self.cost_penalty#+ abs(self.con_vals)#self.constraintFuns_penalty#(X)/1e7 #self.getCAPEX() #+ self.constraintFuns_penalty(X)
    
                else:
                    raise Exception(
                        "The layout 'mode' must be either LCOE, AEP or CAPEX.")

            ''' 
            # ----- write to log -----
            # only log if the design has significantly changed
            if np.linalg.norm(X - self.Xlast) > 100:  # <<< threshold should be customized
                # log the iteration number, design variables, objective, and constraints
                self.log['x'].append(list(X))
                self.log['f'].append(list([self.obj_value]))
                # Check if self.con_vals is an integer - Different optimizer require different constraints
                if isinstance(self.con_vals, int):
                    # Convert self.con_vals to a list before appending to self.log['g']
                    self.log['g'].append([self.con_vals])
                else:
                    # If self.con_vals is already iterable, directly append it to self.log['g']
                    self.log['g'].append(list(self.con_vals))
            '''
        
        self.Xlast = np.array(X)   # record the current design variables
       

    def updateLayoutUG(self, Xu, level=0, refresh=False):
        '''Interface from uniform grid design variables to turbine coordinates.'''
        
        X_points = []
        # create grid points
        if len(self.sub_boundary_sh) > 0:
            # determine # of grid variables per sub boundary
            nXu = 7 if self.rotation_mode else 6

            # create grid points for each sub grid
            for ind in range(len(self.sub_boundary_sh)):
                # pull out relevant design variables
                Xus = Xu[nXu*ind:nXu*(ind+1)]
                # convert km to m for first 4 variables
                Xum = np.hstack([[x*1000 for x in Xus[0:4]], Xus[4:]])
                # generate grid points
                X_points.extend(self.generateGridPoints(Xum,trans_mode='x',boundary_index=ind))
        else:
            # create grid points for entire grid
            Xum = np.hstack([[x*1000 for x in Xu[0:4]], Xu[4:]])  # convert first 4 entries from km to m
            # generate grid points
            X_points.extend(self.generateGridPoints(Xum,trans_mode='x'))
            
        # pare down grid points to those furthest from boundaries & optionally add substation(s) in grid
        X = self.pareGridPoints(X_points)
        
        self.updateLayout(X, level, refresh)  # update each turbine's position

    #def updateLayoutOPTUG(self, Xu):
     #    '''Interface from uniform grid design variables to turbine coordinates.'''              
    #     X = self.generateGridPoints(Xu) 
    #     X2 = np.array(X)   # make a copy of the design vector
    #     X2[:2*self.nt] = X[:2*self.nt]#*1000  # convert coordinates from km to m
    #     self.updateLayout(X2)


    def updateLayoutDB(self, Xdb, level=0, refresh=False):
        '''Interface for Dogger Bank style layouts.'''
        
        ### Xdb[0] and Xdb[1] are exterior spacings. db_ext_spacing allows the user to set what spacing each side uses (in order of coordinates)
        interior =self.boundary_sh.buffer(-self.mooringList[0].rad_anch - self.anchor_minrad) ### this buffer should ensure anchor stays within boundary --- need to check
        coords = list(interior.exterior.coords)
        
        from shapely.geometry import LineString
        
        #iterate through boundaries
        points = []
        for i in range(0, len(coords) - 1):
            
            # connect exterior coordinates in order
            line = LineString([coords[i], coords[i+1]])

            # db_ext_spacing input allows the user to set which boundaries use which outer spacing
            # determine number of turbines that will fit
            num = math.floor(line.length/Xdb[self.db_ext_spacing[i]])
            
            #interpolate along the side for the num turbines
            if i == 0:
                points.extend([line.interpolate(i/num , normalized = True) for i in range(num)])
            else:
                
                #after the first side, start turbines at +spacing so there isn't overlap at the corner
                points.extend([line.interpolate(i/num , normalized = True) for i in range(1, num)])
        
        
        xs = [point.coords[0][0] for point in points]
        ys = [point.coords[0][1] for point in points]
        
        
        
        #fill the interio using generateGridPoints
        interiorinterior = interior.buffer(-self.mooringList[0].rad_anch - self.anchor_minrad) ### again this buffer needs to be checked
        
        #store original exterior boundary and numturbines
        boundary_sh_int = self.boundary_sh_int
        nt = self.nt
        
        interior_nt = nt - len(xs)
        if interior_nt < 0:
            interior_nt = 0
        
        self.boundary_sh_int = interiorinterior
        self.nt = interior_nt
        
        if self.nt > 0:
            
            X = self.generateGridPoints(Xdb[2:],trans_mode='x') 
            #combined exterior and interior turbines into X vector
            Xall = list(X[:self.nt]) + xs + list(X[self.nt:]) + ys
            x_coords =  list(X[:self.nt]) + xs
            y_coords = list(X[self.nt:]) + ys
        else:
            print('Exterior coords filled the required number of turbines')
            
            xs = xs[:nt]
            ys = ys[:nt]
        
            Xall =  xs + ys
            x_coords =   xs
            y_coords = ys
        
        #revert boundary and nt
        self.boundary_sh_int = boundary_sh_int  
        self.nt = nt
        self.turb_coords = np.zeros((self.nt,2))
        
        
        #create buffers for exterior points (generateGridPoints did this for interior already)
        for i in range(0, len(xs)):
            point = Point(xs[i], ys[i])
            self.platformList[interior_nt + i].setPosition([point.x,point.y], heading=None, degrees=False, project = self)  
            atts = [x['obj'] for x in self.platformList[interior_nt + i].attachments.values()]
            mList = [x for x in atts if type(x)==Mooring]

            # switch anchor type
            anchs = self.platformList[i].getAnchors() 
            for anch in anchs.values():
                name, props = self.getSoilAtLocation(anch.r[0],anch.r[1])
                atype = self.anchorTypes[name]
                anch.dd.update(atype)
                anch.mass = self.anchorMasses[name]
                anch.cost['materials'] = self.anchorCosts[name]
                anch.soilProps = {name:props}

            # Get depth at turbine postion          
            self.turb_depth[interior_nt + i] = -self.getDepthAtLocation(
                      point.x, point.y)
            buffer_group_pos = []    
            
            
            for j in range(self.ms_na):
                # im = 3*len(points) + j  # global index of mooring/anchor        
                moor_bf_start = get_point_along_line([point.x, point.y], mList[j].rA[:2],self.turb_minrad)
                # Buffer zone mooring line
                #line = LineString([self.turb_coords[i,:], self.mooringList[im].rA[:2]])
                line = LineString([moor_bf_start, mList[j].rA[:2]])
                mooringline_buffer = line.buffer(self.moor_minrad)
                
                # Buffer zone anchor
                # Create a point at coordinates (x, y)
                point1 = Point(mList[j].rA[:2])
                # Create a buffer around the anchor with a radius of X
                anchor_buffer = point1.buffer(self.anchor_minrad)
                
                # Buffer zone turbine
                # Create a buffer around the anchor with a radius of X
                turb_buffer = point.buffer(self.turb_minrad)
                
                # Buffer group for turbine positioning
                buffer_group_pos.append(mooringline_buffer)
                buffer_group_pos.append(anchor_buffer)
                buffer_group_pos.append(turb_buffer)
            polygon = unary_union(buffer_group_pos)  # Combine buffers for each turbine
            
            if self.boundary_sh_int.contains(polygon):
                # If the point is within the shape, append it to the list of bufferzones
                
                self.ms_bufferzones_pos[interior_nt + i] = polygon
                
        
        
        if len(x_coords) < self.nt:
            for i in range(len(points)):
                self.turb_coords[i,0] = x_coords[i]
                self.turb_coords[i,1] = y_coords[i]
        else:    
            self.turb_coords[:,0] = x_coords
            self.turb_coords[:,1] = y_coords
        
        self.updateLayout(Xall, level, refresh)
        
        

    def updateLayoutOPT(self, X):
        '''Wrapper for updateLayout that uses km instead of m.'''
        X2 = np.array(X)   # make a copy of the design vector
        X2[:2*self.nt] = X[:2*self.nt]*1000  # convert coordinates from km to m
        self.updateLayout(X2)
    


    # ----- OBJECTIVE FUNCTION -----
    def objectiveFunUG(self, Xu):
        '''The general objective function. Will behave differently depending 
        on settings. Only input is the design variable vector, Xu.'''
        # print('Xu in objective function: ',Xu)
        # X = self.generateGridPoints(Xu,trans_mode='x') 

        # update the layout with the specified design vector
        # self.updateLayoutUG(X)
        #Xum = np.hstack([[x*1000 for x in Xu[0:4]], Xu[4:]])  # convert first 4 entries from km to m  
        self.updateLayoutUG(Xu)
        #self.updateLayoutOPTUG(X)
        return self.obj_value
    
    
    def objectiveFunDB(self, Xdb):
        '''The general objective function. Will behave differently depending 
        on settings. Only input is the design variable vector, Xu.'''

        # update the layout with the specified design vector
        # self.updateLayoutUG(X)

        self.updateLayoutDB(Xdb)
        
        #self.updateLayoutOPTUG(X)

        return self.obj_value
    
    
    def objectiveFun(self, X):
        '''The general objective function. Will behave differently depending 
        on settings. Only input is the design variable vector, X.'''

        # update the layout with the specified design vector
        self.updateLayoutOPT(X)

        return self.obj_value


    # ----- ANCHORS ----- 


                
    # ----- AEP / FLORIS ----- 
    def getAEP(self, display = 0):
        '''Compute AEP using FLORIS, based on whatever data and turbine
        positions are already stored in the Layout object. 
        (updateLayout should have been called before this method.'''
        
        # FLORIS inputs positions in m
        self.flow.set(layout_x=self.turb_coords[:,0],
                               layout_y=self.turb_coords[:,1] )

        #run floris simulation
        self.flow.set(wind_data = self.wind_rose)
        self.flow.run()
        
        #frequencies must be in list 
        self.aep = self.flow.get_farm_AEP()
        
        if display > 0:
            self.plotWakes(wind_spd = 10, wind_dir = 270, ti = 0.06)
        #return self.aep
    """
    # ----- CAPEX ----- 
    def getCAPEXMooring(self):
        '''Compute CAPEX of mooring systems. Currently test function only'''
        self.capex_mooring = sum(abs(self.ms_anchor_depth))*1000
        #capex_mooring = sum(self.turb_depth**2)
        #return self.capex_mooring
   
    # ----- TOTAL CAPEX function
    def getCAPEX(self):
        '''Compute TOTAL CAPEX, adding sub-methods together. Currently test function only'''
        self.getCAPEXMooring()
        self.capex_total = self.capex_mooring
    """
    # ----- TOTAL CAPEX function
    def getCost(self):
        '''Compute TOTAL CAPEX, adding sub-methods together. Currently test function only'''
        
        CapEx_mooring = 0
        for mooring in self.mooringList.values():
            CapEx_mooring += mooring.getCost()

        CapEx_anchors = 0
        for anchor in self.anchorList.values():
            CapEx_anchors += anchor.getCost()
        
        CapEx_cables = 0
        for cable in self.cableList.values():
            CapEx_cables  += cable.getCost()
            
            
            
        # Include cable costs for feasible layouts
        #if self.con_sum == 0:
            #self.cost_cable = sum(self.iac_cost)
        #    self.cost_total = CapEx+self.cost_cable
        #else:
        self.cost_total = CapEx_mooring + CapEx_cables + CapEx_anchors
        return self.cost_total     
    
    
    # ----- LCOE ----- 
    def getLCOE(self):
        '''Compute LCOE = CAPEX / AEP. Currently test function and based on CAPEX only.'''
        self.getAEP(display = self.display)
        self.getCost()
        self.lcoe = self.cost_total/self.aep#self.getCOST()/(self.getAEP()/ 1.0e6) # [$ / MWh]self.getAEP() 
        #return self.cost
        
    # ----- LCOE ----- 
    def getLCOE2(self):
        '''updated LCOE function using capex, opex, and fcr assumptions from previous projects'''
 
        farm_capacity = self.turb_rating_MW * self.nt * 1000 # kW
        
        capex = 3748.8 * farm_capacity # $ does NOT include moorings/cables. 
        #from DeepFarm LCOE report GW scale individual wind farm, substracted mooring system and array system costs
        opex = 62.51 *farm_capacity # $ annually.  from DeepFarm LCOE report
        fcr = 5.82/100 # fixed charge rate %. from DeepFarm LCOE report
        
        self.getAEP()
        self.getCost()
        self.lcoe = ((self.cost_total+capex)*fcr+ opex)/self.aep*1e6 # [$ / MWh]

        #return self.cost
   
    # ----- PENALTY FUNCTION -----    
    def constraintFuns_penalty(self, X):
        '''Penalty function to better guide the optimization. Only input is the design variable vector, X.'''  
        self.getCost()          
        self.constraintFuns(X) 
       
        #con_vals =  self.con_vals#self.constraintFuns(X)  
        # Get the indices of negative values
        #negative_indices = np.where(con_vals < 0)[0]
        #return self.getCAPEX()*0.1*abs(np.sum(con_vals[negative_indices]))#*1e3*self.nt**2 
        #self.cost_penalty = self.cost_total*0.5*abs(np.sum(con_vals[negative_indices]))#*1e3*self.n
        
        self.cost_penalty =  self.con_vals
        
        #return self.cost_penalty



    # ----- CONSTRAINTS FUNCTION -----
    # --------------------------------
    def constraintFunsUG(self, Xu):
        '''The general constraints function. Will behave differently depending 
        on settings. Only input is the design variable vector, X.'''
        
        #X = self.generateGridPoints(Xu) 
        #Xum = np.hstack([[x*1000 for x in Xu[0:4]], Xu[4:]])  # convert first 4 entries from km to m 
        # print(Xu)
        # if any([x>2500 for x in Xu]):
        #     breakpoint()
        # update the layout with the specified design vector
        self.updateLayoutUG(Xu)
        #self.updateLayoutOPTUG(Xu)
        return self.con_vals

    def constraintFunsDB(self, Xdb):
        '''The general constraints function. Will behave differently depending 
        on settings. Only input is the design variable vector, X.'''

        # update the layout with the specified design vector
        self.updateLayoutDB(Xdb)
        
        return self.con_vals
        
    
    def constraintFuns(self, X):
        '''The general constraints function. Will behave differently depending 
        on settings. Only input is the design variable vector, X.'''

        # update the layout with the specified design vector
        self.updateLayoutOPT(X)
        
        return self.con_vals
        

    def calcDerivatives(self):
        '''Compute the derivatives about the current state of the layout,
        for use with optimizers that accept a Jacobian function.
        This is explicitly designed for when variables are x, y, and h.
        >>> PLACEHOLDER <<<
        '''
        
        nDOF = 3*self.nt
        '''
        # Perturb each DOF in turn and compute AEP results
        J_AEP = np.zeros([nt,nt])

        # Perturp each turbine and figure out the effects on cost and constraint
        for i in range(nt):
            J_CONS_i = np.zeros([ng, 3])   # fill in each row of this (or each column?)
        
        
        # then combine them into overall matrices
        
        J_cost
        
        J_constraints....
        
        # >>>> need to straighten out constraint vectors...
        '''

    def saveLOG(self, filename):
     
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            # Write header
            writer.writerow(['it','x', 'f', 'g'])
            # Write data
            for i in range(len(self.log['x'])):
                it = [i]
                x_values = self.log['x'][i]    # Design variables
                f_value = self.log['f'][i]  # Result of objective function
                g_values = self.log['g'][i] # Scalar, either -1 or 1

                writer.writerow([it, x_values, f_value, g_values])


    # ----- Plot wind farm layout -----
    def plotLayout(self, ax=None, bare=False, save=False):
        '''Plot wind farm layout.'''

        # if axes not passed in, make a new figure
        if ax == None:
            fig, ax = plt.subplots(1,1, figsize=[6,6])
        else:
            fig = ax.get_figure()
        
        # Set font sizes
        fsize_legend = 12    # Legend 
        fsize_ax_label = 12  # Ax Label 
        fsize_ax_ticks = 12  # Ax ticks
        fsize_title = 16     # Title 
        
        x0 = self.turb_coords[:,0]
        y0 = self.turb_coords[:,1]

        
        # Plot the layout, using the internally stored information.
        
        #breakpoint()
        # ----- Bathymetry / contourf

        #num_levels = 10  # Adjust this value as needed        
        X, Y = np.meshgrid(self.grid_x, self.grid_y)
        #breakpoint()
        depth_min =np.min(self.grid_depth)
        depth_min=math.floor(depth_min / 10) * 10 
        depth_max =np.max(self.grid_depth)
        depth_max=math.ceil(depth_min / 10) * 10 
        
        depth_range = depth_max- depth_min
        
        if depth_range < 100:
            steps_m = 10
        else:
            steps_m = 100
        
        num_levels = round((depth_max- depth_min)/steps_m)        
       

        if depth_min != depth_max:
            contourf = ax.contourf(X, Y, self.grid_depth, num_levels, cmap='Blues', vmin=depth_min, vmax=depth_max)
            #contourf = ax.contourf(X, Y, self.grid_depth[x_indices, y_indices], num_levels, cmap='Blues', vmin=0, vmax=1000)
            #contourf.norm.autoscale([0,1])
            
            #contourf.set_clim(0, 1000)
            
            # Add colorbar with label
            if not bare:            
                cbar = plt.colorbar(contourf, ax=ax, fraction=0.04, label='Water Depth (m)')
                # Set the font size for the colorbar label and ticks
                #cbar.ax.yaxis.label.set_fontsize(fsize_ax_label)
                #cbar.ax.tick_params(axis='y', labelsize=fsize_ax_ticks)
           
        
        # seabed
        X, Y = np.meshgrid(self.soil_x, self.soil_y)
        ax.scatter(X, Y, s=4, cmap='cividis_r', vmin=-0.5, vmax=1.5)
        
        # ----- OSS        
        for oo in self.oss_coords:
            ax.scatter(oo[0],oo[1], color='red', marker='*', label='OSS', s=100)
            circle = plt.Circle((oo[0], oo[1]), self.oss_minrad, edgecolor=[.5,0,0,.8],
                     facecolor='none', linestyle='dashed', lw=0.8)

        # (AEP: {aep / 1.0e9:.2f} GWh,\n CAPEX: M$ {cost/1.0e6:.2f},\n LCOE: {lcoe:.2f} $/MWh)'
        
        # plt.scatter(x, y, color='blue', marker='D')
        # plt.scatter(optimized_x_pos, optimized_y_pos, label=f'Optimized Positions (AEP: {optimized_aep / 1.0e9:.2f} GWh)', color='red', marker='D')

        # Anchors
        #plt.scatter(self.anchor_coords[:,0], self.anchor_coords[:,1],
        #    label='Anchor Positions', color='red', marker='.')

        # Plot mooring buffer zones
        for i, polygon in enumerate(self.ms_bufferzones_pos):
            if isinstance(polygon, MultiPolygon):
                for poly in polygon:
                    x, y = poly.exterior.xy
                    ax.plot(x, y,color='red')
            else:
                x, y = polygon.exterior.xy
                #ax.plot(x, y,color='red')
                ax.fill(x, y,color=[.6,.3,.3,.6])
        # Add a single legend entry outside the loop
        if not bare:
            legend_entry = ax.fill([], [], color=[.6,.3,.3,.6], label='Mooring Buffer Zone')
        
        # Add a legend with fontsize
        if not bare:
            ax.legend(handles=legend_entry) #, fontsize=fsize_legend)

        # ----- mooring lines
        for i in range(self.nt):
            for j in range(3):
                plt.plot([self.turb_coords[i,0], self.mooringList[3*i+j].rA[0]], 
                         [self.turb_coords[i,1], self.mooringList[3*i+j].rA[1]], 'k', lw=0.5)

               # plt.plot([self.turb_coords[i,0], self.anchor_coords[3*i+j,0]], 
                #         [self.turb_coords[i,1], self.anchor_coords[3*i+j,1]], 'k', lw=0.5)



        # ----- Minimum distance
        i = 0
        for x, y in zip(x0, y0):
            if i == 0:
                circle = plt.Circle((x, y), self.turb_minrad, edgecolor=[.5,0,0,.8],
                                facecolor='none', linestyle='dashed', label='Turbine Buffer Zone', lw=0.8)
            else:
                circle = plt.Circle((x, y), self.turb_minrad, edgecolor=[.5,0,0,.8],
                                facecolor='none', linestyle='dashed', lw=0.8)
            i =+ 1 
            ax.add_patch(circle)
            # Add a legend to the axes with fontsize
            if not bare:
                ax.legend() #fontsize=fsize_legend)
            # plt.gca().add_patch(circle)

        # ----- Lease area boundary
        #shape_polygon = sh.Polygon(self.boundary)
        x, y =  self.boundary_sh.exterior.xy
        ax.plot(x, y, label='Boundary', linestyle='dashed', color='black')
        
        # ----- Sub boundaries
        for subb in self.sub_boundary_sh:
            x,y = subb.exterior.xy
            ax.plot(x,y, label='Sub-boundary', linestyle=':', color='blue')

      
        # ----- Exclusion zones
        if len(self.exclusion) !=0:
            for ie in range(len(self.exclusion)):
                shape_polygon = self.exclusion_polygons_sh[ie]#sh.Polygon(self.exclusion[i])
                x, y = shape_polygon.exterior.xy
                ax.plot(x, y, linestyle='dashed', color='orange', label='Exclusion Zone')   
            #ax.plot([], [], linestyle='dashed', color='orange', label='Exclusion Zone')
            
        # turbine locations
        ax.scatter(x0, y0, c='black', s=12, label='Turbines')
         
        
        
        if self.cable_mode:
            # ----- Cables                
            # Create a colormap and a legend entry for each unique cable section
            # Find unique values 
            unique_cables = np.unique([x['conductor_area'] for x in self.iac_dic]) #(self.iac_dic['minimum_con'].values)
            colors = plt.cm.viridis(np.linspace(0, 1, len(unique_cables)))  # Create a colormap based on the number of unique sections
            section_to_color = {sec: col for sec, col in zip(unique_cables, colors)}
                  
     
            # ----- Cables in Cluster
            # Cable array
            iac_array = self.iac_dic
            count = 0
            # Loop over each cluster
            for ic in range(self.n_cluster*self.noss):      
                # Plot vertices
                #plt.scatter(self.cluster_arrays[ic][:, 0], self.cluster_arrays[ic][:, 1], color='red', label='Turbines')
                
                # Annotate each point with its index
                #for i, point in enumerate(self.cluster_arrays[ic]):
                    #plt.annotate(str(i), (point[0], point[1]), textcoords="offset points", xytext=(0, 10), ha='center')
                 
                # Get index of cluster
                #ind_cluster = np.where(iac_array[:, 0] == 0)[0]
                # Loop over edges / cable ids
                len_cluster = len(np.where(np.array([x['cluster_id']==ic for x in iac_array]))[0])
                for i in range(len_cluster):
                    ix = np.where((np.array([x['cluster_id']== ic for x in iac_array])) & (np.array([y['cable_id']== count for y in iac_array]) ))[0]
                    if len(ix)<1:
                        breakpoint()
                    ind = ix[0]
                    #ind = np.where((iac_array[:, 0] == ic) & (iac_array[:, 2] == i))[0][0]
                    # Plot edge
                    #edge = self.iac_edges[ic][i]    
                    start = iac_array[ind]['coordinates'][0]#self.cluster_arrays[ic][edge[0]]
                    end = iac_array[ind]['coordinates'][1]
                    # Cable selection
                    color = section_to_color[iac_array[ind]['conductor_area']]
                    ax.plot([start[0], end[0]], [start[1], end[1]], color=color, label=f'Section {int(iac_array[ind]["conductor_area"])} mm²' if int(iac_array[ind]["conductor_area"]) not in plt.gca().get_legend_handles_labels()[1] else "")
                    #plt.text((start[0] + end[0]) / 2, (start[1] + end[1]) / 2, str(i), fontsize=9, color='black')
                    # for sid in oss_ids:
                    #     if iac_array[ix]['turbineA_glob_id'] == sid or iac_array[ix]['turbineB_glob_id'] == sid:
                    #         iac_array_oss.append(iac_array[ix])
                    #         iac_oss_id.append(sid)
                    
                    count += 1
                    
                # Plot gate as a diamond marker
                #plt.scatter(self.gate_coords[ic][0], self.gate_coords[ic][1], marker='D', color='green', label='Gate')
            
            
            ## ----- Cables Gates to OSS

            # for i in range(self.n_cluster):
            #     cable_section_size = int(iac_array_oss[i]['conductor_area'])  # Assuming cable section size is in the 7th column
            #     color = section_to_color.get(cable_section_size, 'black')  # Default to black if section size not found
            #     oss_coord = self.substationList[iac_oss_id[i]].r
            #     ax.plot([iac_array_oss[i]['coordinates'][1][0], oss_coord[0]], [iac_array_oss[i]['coordinates'][1][1],oss_coord[1]], color=color, label=f'Section {cable_section_size} mm²' if cable_section_size not in plt.gca().get_legend_handles_labels()[1] else "")
                
              
        
        '''
        # NEW: TURBINE CLUSTER AND CABLES
        # Plot turbines by cluster
        for label in set(self.cluster_labels):
            cluster_turbines = [self.turb_coords[i] for i, lbl in enumerate(self.cluster_labels) if lbl == label]
            if cluster_turbines:  # Check if list is not empty
                x, y = zip(*cluster_turbines)
                ax.scatter(x, y, label=f'Cluster {label}')
        
        # Plot edges
        for i in range(len(self.cluster_edges)):
            P = self.cluster_arrays[i]
            for edge in self.cluster_edges[i]:
                i, j = edge
                plt.plot([P[i, 0], P[j, 0]], [P[i, 1], P[j, 1]], color ='black')
        
        # Plot OSS and gates
        ax.scatter(*self.oss_coords, color='red', marker='*', label='OSS')
        ax.scatter(self.gate_coords[:, 0], self.gate_coords[:, 1], color='black', marker='d', label='Gates')
        
        # Legend adjustment might be needed depending on the number of elements
        #ax.legend(loc='upper center', fancybox=True, ncol=2)
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, ncol=2)
        '''
        
        
        
        # ----- Additional plot customization 
        # Set x- and y-axis ticks fontsize
        if not bare:
            ax.set_xticks(ax.get_xticks())
            ax.set_yticks(ax.get_yticks())
            ax.set_xticklabels(ax.get_xticklabels()) #, fontsize=fsize_ax_ticks)
            ax.set_yticklabels(ax.get_yticklabels()) #, fontsize=fsize_ax_ticks)
     
            # # Define a custom formatter to divide ticks by 1000
            # def divide_by_1000(value, tick_number):
            #     return f'{value/1000:.0f}'
        
            # # Apply the custom formatter to the x and y axis ticks
            # ax.xaxis.set_major_formatter(FuncFormatter(divide_by_1000))
            # ax.yaxis.set_major_formatter(FuncFormatter(divide_by_1000))
                       
            #ax.axis("equal")
            ax.set_aspect('equal')
            
            # # Use AutoLocator for major ticks
            # ax.xaxis.set_major_locator(AutoLocator())
            # ax.yaxis.set_major_locator(AutoLocator())     
            # # Use AutoMinorLocator for minor ticks
            # ax.xaxis.set_minor_locator(AutoMinorLocator())
            # ax.yaxis.set_minor_locator(AutoMinorLocator())
            
            # ax.set_xlim([self.grid_x[0], self.grid_x[-1]])
            # ax.set_ylim([self.grid_y[0], self.grid_y[-1]])
            #ax.set_xlim([x_min_bounds-1000, x_max_bounds+1000])
            #ax.set_ylim([y_min_bounds-1000, y_max_bounds+1000])
            
            
            #plt.title('Optimized Wind Farm Layout',fontsize=fsize_title)
            plt.xlabel('x (km)') #,fontsize=fsize_ax_label)
            plt.ylabel('y (km)') #,fontsize=fsize_ax_label)
            #plt.legend(loc='upper center', bbox_to_anchor=(
            #    0.5, -0.2), fancybox=True, ncol=3)
            #plt.legend(loc='upper center', fancybox=True, ncol=2)
            handles, labels = plt.gca().get_legend_handles_labels()
            unique_labels = list(set(labels))  # Get unique labels
            unique_labels.sort()  # Sort the unique labels alphabetically
            unique_handles = [handles[labels.index(label)] for label in unique_labels]  # Get handles corresponding to unique labels
            plt.legend(unique_handles, unique_labels, loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, ncol=2)
            plt.gca().set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
            #ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, ncol=2)
            
            #breakpoint()
            # Set plot area to around the lease area
            
            
            # Calc plot bounds
            #offset = 1000
            #offset_polygon = translate(self.boundary_sh, xoff=offset, yoff=offset)
            # Get bounds
            #x_min_bounds, y_min_bounds, x_max_bounds, y_max_bounds = offset_polygon.bounds
            # Round to next 100       
            #x_min_bounds, y_min_bounds = [math.floor(v / 1000) * 1000 for v in (x_min, y_min)]
            #x_max_bounds, y_max_bounds = [math.ceil(v / 1000) * 1000 for v in (x_max, y_max)]   

            
            
            # ----- Save plot with an incremented number if it already exists
            if save:
                counter = 1
                output_filename = f'wind farm layout_{counter}.png'
                while os.path.exists(output_filename):
                    counter += 1
                    output_filename = f'wind farm layout_{counter}.png'
                
                # Increase the resolution when saving the plot
                plt.savefig(output_filename, dpi=300, bbox_inches='tight')  # Adjust the dpi as needed
            
            # also print some output
            
            if self.flow:  # if FLORIS
                print('AEP:', self.aep)
                
            self.getCost()
            print('Cost:', self.cost_total)
            
            # for mooring in self.mooringList.values():
            #     print(mooring.cost)
            
    """        
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

        # try icnraesing dpeht grid density for nicer plot
        xs = np.arange(-1000,8000,500)
        ys = np.arange(-1000,9500,500)
        #self.setGrid(xs, ys)
        zs = np.zeros([len(ys), len(xs)])
        for i in range(len(ys)):
            for j in range(len(xs)):
                zs[i,j] = self.getDepthAtLocation(xs[j], ys[i])
        X, Y = np.meshgrid(xs, ys)  # 2D mesh of seabed grid
        

        # plot the bathymetry in matplotlib using a plot_surface
        #X, Y = np.meshgrid(self.grid_x, self.grid_y)  # 2D mesh of seabed grid
        ax.plot_surface(X, Y, -zs, **args_bath)
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
        #bath = ax.plot_surface(X, Y, -self.grid_depth, **args_bath)
        #
        
        
        # also if there are rocky bits... (TEMPORARY)
        '''
        X, Y = np.meshgrid(self.soil_x, self.soil_y)
        Z = np.zeros_like(X)
        xs = self.soil_x
        ys = self.soil_y
        for i in range(len(ys)):
            for j in range(len(xs)):
                Z[i,j] = -self.getDepthAtLocation(xs[j], ys[i])
        ax.scatter(X, Y, Z+5, c=self.soil_rocky, s=6, cmap='cividis_r', vmin=-0.5, vmax=1.5, zorder=0)
        '''
        
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
                mooring.subsystem.drawLine(0, ax, shadow=False)
        
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
    """        

    def playOptimization(self):
        '''A very slow clunky way to animate the optimization'''
        fig, ax = plt.subplots(1,1)
        
        self.updateLayout(self.log['x'][0])
        self.plotLayout(ax=ax)

        def animate(i):
            ax.clear()
            self.updateLayout(self.log['x'][i])
            self.plotLayout(ax=ax, bare=True)
        
        ani = FuncAnimation(fig, animate, frames=len(self.log['x']),
                            interval=500, repeat=True)
    
        return ani
  
    
  
    
  
    
    
    
    
      

  
    
  
    
  
    
  
    
  
    
    def plotOptimization(self):
    
        if len(self.log['x']) == 0:
            print("No optimization trajectory saved (log is empty). Nothing to plot.")
            return
        
        
        
        fig, ax = plt.subplots(5,1, sharex=True, figsize=[6,8])
        fig.subplots_adjust(left=0.4)
        
        X = np.array(self.log['x'])
        Fs = np.array(self.log['f'])
        Gs = np.array(self.log['g'])
        
        
        if self.rotation_mode:
            x_pos, y_pos, rot_rad = X[:,:self.nt], X[:,self.nt:2*self.nt], X[:,2*self.nt:]
        else:
            x_pos, y_pos = X[:,:len(X)//2], X[:,len(X)//2:]
            rot_rad = np.zeros_like(x_pos)
        
        for i in range(self.nt):
            ax[0].plot(x_pos[:,i])
            ax[1].plot(y_pos[:,i])
            ax[2].plot(rot_rad[:,i])

        ax[3].plot(Fs)
        ax[3].set_ylabel("cost", rotation='horizontal')
        
        Gs_neg = Gs*(Gs < 0)
        ax[4].plot(np.sum(Gs_neg, axis=1))
        ax[4].set_ylabel("constaint violation sum", rotation='horizontal')
        '''
        for i, con in enumerate(self.constraints):
            j = i+1+len(X)
            ax[j].axhline(0, color=[0.5,0.5,0.5])
            ax[j].plot(Gs[:,i])
            ax[j].set_ylabel(f"{con['name']}({con['threshold']})", 
                           rotation='horizontal', labelpad=80)
        '''
        ax[-1].set_xlabel("iteration roughly")
        
        
        
        """
        nX = len(self.log['x'][0])
        fig, ax = plt.subplots(nX+1+1,1, sharex=True, figsize=[6,8])
        fig.subplots_adjust(left=0.4)
        Xs = np.array(self.log['x'])
        Fs = np.array(self.log['f'])
        Gs = np.array(self.log['g'])
        
        for i in range(nX):
            ax[i].plot(Xs[:,i])
            #ax[i].axhline(self.Xmin[i], color=[0.5,0.5,0.5], dashes=[1,1])
            #ax[i].axhline(self.Xmax[i], color=[0.5,0.5,0.5], dashes=[1,1])

        ax[nX].plot(Fs)
        ax[nX].set_ylabel("cost", rotation='horizontal')
        
        Glist = Gs.ravel()
        
        ax[nX+1].plot(np.sum(Glist[Glist<0]))
        ax[nX+1].set_ylabel("constaint violation sum", rotation='horizontal')
        '''
        for i, con in enumerate(self.constraints):
            j = i+1+len(X)
            ax[j].axhline(0, color=[0.5,0.5,0.5])
            ax[j].plot(Gs[:,i])
            ax[j].set_ylabel(f"{con['name']}({con['threshold']})", 
                           rotation='horizontal', labelpad=80)
        '''
        ax[-1].set_xlabel("iteration roughly")
        """
    
    def plotCost(self):
        '''Makes a bar chart of the cost breakdown.'''
    
    
    def plotWakes(self, wind_spd, wind_dir, ti):
        '''uses floris tools to plot wakes'''
        import floris.layout_visualization as layoutviz
        from floris.flow_visualization import visualize_cut_plane
        
        fmodel = self.flow
        
        # Create the plotting objects using matplotlib
        fig, ax = plt.subplots()

      
        layoutviz.plot_turbine_points(fmodel, ax=ax)
        layoutviz.plot_turbine_labels(fmodel, ax=ax)
        ax.set_title("Turbine Points and Labels")
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        
     
        
        fmodel.set(wind_speeds=[wind_spd], wind_directions=[wind_dir], turbulence_intensities=[ti])
        horizontal_plane = fmodel.calculate_horizontal_plane(
            x_resolution=200,
            y_resolution=100,
            height=90.0,
        )
        
        # Plot the flow field with rotors
        fig, ax = plt.subplots()
        visualize_cut_plane(
            horizontal_plane,
            ax=ax,
            label_contours=False,
            title="Horizontal Flow with Turbine Rotors and labels",
        )
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        
        # Plot the turbine rotors
        layoutviz.plot_turbine_rotors(fmodel, ax=ax)
        
        plt.show()
                        
# Calculate offset from the turbine to create buffer zones for cable routing
def get_point_along_line(start, end, diste):
    # Convert inputs to numpy arrays
    start = np.array(start)
    end = np.array(end)
    # Calculate the direction vector from start to end
    direction = end - start
    # Normalize the direction vector
    length = np.linalg.norm(direction)
    unit_direction = direction / length
    # Calculate the new point at the specified distance along the direction vector
    new_point = start + unit_direction * diste
    return new_point       
    
# def mooringAdjuster1(mooring, project, r, u, level=0):
#     '''Custom function to adjust a mooring, called by
#     Mooring.adjust. Fairlead point should have already
#     been adjusted.'''
    
#     ss = mooring.ss  # shorthand for the mooring's subsystem
    
#     T_target = 1e6  # target mooring line pretension [N] (hardcoded example)
#     i_line = 0  # line section to adjust (if multiple) (hardcoded example)
    
#     #>>> pit in better prpfile <<<
    
#     # Find anchor location based on desired relation
#     r_i = np.hstack([r + 58*u, -14]) # fairlead point
#     slope = 0.58  # slope from horizontal
#     u_a = np.hstack([u, -slope])  # direct vector from r_i to anchor
#     r_anch = project.seabedIntersect(r_i, u_a)  # seabed intersection
    
#     # save some stuff for the heck of it
#     mooring.z_anch = r_anch[2]
#     mooring.anch_rad = np.linalg.norm(r_anch[:2]-r)
    
#     mooring.setEndPosition(r_anch, 'a')  # set the anchor position

#     # Estimate the correct line length to start with
#     ss.lineList[0].setL(np.linalg.norm(mooring.rB - mooring.rA))
        
#     # Next we could adjust the line length/tension (if there's a subsystem)
#     if level==1:  # level 1 analysis (static solve to update node positions)
#         ss.staticSolve()
        
#     elif level==2:  # adjust pretension (hardcoded example method for now)
        
#         def eval_func(X, args):
#             '''Tension evaluation function for different line lengths'''
#             ss.lineList[i_line].L = X[0]  # set the first line section's length
#             ss.staticSolve(tol=0.0001)  # solve the equilibrium of the subsystem
#             return np.array([ss.TB]), dict(status=1), False  # return the end tension

#         # run dsolve2 solver to solve for the line length that matches the initial tension
#         X0 = [ss.lineList[i_line].L]  # start with the current section length
#         L_final, T_final, _ = dsolve2(eval_func, X0, Ytarget=[T_target], 
#                               Xmin=[1], Xmax=[1.1*np.linalg.norm(ss.rB-ss.rA)],
#                               dX_last=[1], tol=[0.1], maxIter=50, stepfac=4)
#         ss.lineList[i_line].L = L_final[0]
    
    
#     # Compute anchor size and cost
#     soilr = project.getSoilAtLocation(*r_anch[:2])
#     if 'rock' in soilr:
#         rocky = 1 
#     else:
#         rocky = 0
    
#     anchor_cost = 300e3 + rocky*200e3
#     mooring.cost['anchor'] = anchor_cost
    
    
    # getWatchCircle() method
    # getMudlineForces(, max_forces=True)
    
    
    
    

if __name__ == '__main__':

    # Wind rose
    from floris import WindRose
    wind_rose = WindRose.read_csv_long(
        'humboldt_rose.csv', wd_col="wd", ws_col="ws", freq_col="freq_val", ti_col_or_value=0.06)
    
    
    # ----- LEASE AREA BOUNDARIES -----
    WestStart = 10000
    NorthStart = 10000
    boundary_coords = np.array([
    (0, 0),
    (WestStart, 0),
    (WestStart, NorthStart),
    (0,NorthStart)
    ])
 

    
    # Make a sample Subsystem to hold the mooring design (used for initialization)
    print("Making subsystem")
    newFile = '..\scripts\input_files\GoMxOntology.yaml'
    project = Project(file=newFile,raft=0)
    project.getMoorPyArray()
    ss = deepcopy(project.ms.lineList[0])      

    # ----- Set optimization mode
    opt_mode = 'CAPEX'
    #opt_mode = 'AEP'
    #opt_mode = 'LCOE'
    # remember to set use_FLORIS accordingly when initializing Layout
      
    # set substation location
    oss_coords = np.array([0, 0])
    
    
   
    #layouttype = 'freelayout' 
    layouttype = 'uniformgridlayout'
    
    # ----- UNIFORM GRID -----
    if layouttype == 'uniformgridlayout':
        
        
        
        #[grid_spacing_x, grid_spacing_y, grid_trans_x, grid_trans_y, grid_rotang, grid_skew, optional turb_rotation]
        Xu = [1000/1000, 1000/1000, 500/1000, 150/1000, 45, 0, 0]



        # Amount of wind turbines
        nt = 20
        
        #rotation mode and turbine rotation
        rotation_mode = True
        rot_rad=np.zeros((nt))
     
                               
        # Boundaries of design variables for PSO
        boundaries_UG=np.array([[0.5, 3],[0.5, 3], [-.5, .5], [-.5, .5], [0, 180], [0,0.2],[0, 180] ])
        
        #cable routing
        cable_mode = True
    
    # ----- FREE LAYOUT -----
    elif layouttype == 'freelayout':
        
        #first iteration turbine coordinates 
        gulfofmaine_int = np.array([
            [3000, 2000],
            [2000, 2000],
            [0, 2000],
            [1000, 2000],
            [2000, 2000],
            [1000, 100],
            [2000, 100],
            [0, 4000],
            [1000, 4000],
            [2000, 4000]
        ])
        
        
        x_coords = gulfofmaine_int[:, 0]  # x coordinates
        y_coords = gulfofmaine_int[:, 1]  # y coordinates
     
        nt = len(x_coords) # Number of wind turbines
        
        #first iteration rotations
        rot_deg = np.zeros((nt))
        
        
        # ----- Bounds vectors for design variables -----
        #  FOR OPTIMIZER ONLY
        # Lease area boundaries 
        boundaries_x = np.tile([(min(boundary_coords[:,0]), max(boundary_coords[:,0]))], (nt, 1))
        boundaries_y = np.tile([(min(boundary_coords[:,1]), max(boundary_coords[:,1]))], (nt, 1))
        # Rotation in rad 0 - 360*pi/180
        boundaries_rot = np.tile([(0.001, 6.283)], (nt, 1))  
        # Combine into one array

        
            
        # ----- Set rotation mode
        # If True, rotations are considered as design variable, therefore included 
        # into same vector as x and y. Otherwise not.
        rotation_mode = True
        rot_rad = np.deg2rad(rot_deg) # Rotations need to be in rad for the optimization
        x = np.array(x_coords/1000) #km
        y = np.array(y_coords/1000) #km
        
        # Create flattened array xy for initial positions for Layout [km, rad]
        if rotation_mode:
            xy = np.concatenate((x, y, rot_rad))
            boundary_xy = np.concatenate((boundaries_x/1000, boundaries_y/1000, boundaries_rot))
        else:
            xy = np.concatenate((x, y)) 
            boundary_xy = np.concatenate((boundaries_x/1000, boundaries_y/1000))

        # cable routing 
        cable_mode = True



    # ----- Initialize LAYOUT class -----
    print("Initializing Layout")
    
    settings = {}
    settings['n_turbines'] = nt
    settings['turb_rot'] = rot_rad
    settings['rotation_mode'] = rotation_mode
    settings['cable_mode'] = cable_mode
    settings['oss_coords'] = oss_coords
    settings['boundary_coords'] = boundary_coords
    settings['bathymetry_file'] = '..\scripts\input_files\GulfOfMaine_bathymetry_100x100.txt'
    settings['soil_file'] = '..\scripts\input_files\soil_sample.txt'
    settings['floris_file']='gch_floating.yaml'
    #settings['exclusion_coords'] = exclusion_coords
    settings['use_FLORIS'] = False
    settings['mode'] = opt_mode
    settings['optimizer'] ='PSO'
    settings['obj_penalty'] = 1
    settings['parallel'] = False
    settings['n_cluster'] = 3
    
    # set up anchor dictionary
    anchor_settings = {}
    anchor_settings['anchor_design'] = {'L':20,'D':4.5,'zlug':13.3} # geometry of anchor
    anchor_settings['anchor_type'] = 'suction' # anchor type
    anchor_settings['anchor_resize'] = True # bool to resize the anchor or not
    anchor_settings['fix_zlug'] = False # bool to keep zlug the same when resizing anchor
    anchor_settings['FSdiff_max'] = {'Ha':.2,'Va':.2} # max allowed difference between FS and minimum FS
    anchor_settings['FS_min'] = {'Ha':2,'Va':2} # horizontal and vertical minimum safety factors
    
    settings['anchor_settings'] = anchor_settings
    
    
    if layouttype == 'freelayout':
        layout1 = Layout(X=xy, Xu=[], wind_rose = wind_rose, ss=ss, **settings)
    elif layouttype == 'uniformgridlayout':
        layout1 = Layout(X=[], Xu=Xu, wind_rose = wind_rose, ss=ss, **settings)
    


    '''
    # ----- Sequential Least Squares Programming (SLSQP)
    if layouttype == 'freelayout':
        res = minimize(fun=layout1.objectiveFun, x0=xy, method='SLSQP',
                    bounds = boundary_xy,
                    constraints={'type': 'ineq', 'fun': layout1.constraintFuns},
                    options={'maxiter':100, 'eps':0.02,'ftol':1e-6, 'disp': True, 'iprint': 99})
                    #options={'maxiter': 5000,'eps': 0.2, 'finite_diff_rel_step': '2-point', 'ftol': 1e-6, 'disp': True, 'iprint': 99})
    elif layouttype == 'uniformgridlayout': 
        res = minimize(fun=layout1.objectiveFunUG, x0=Xu, method='SLSQP',
                    bounds = boundaries_UG,
                    constraints={'type': 'ineq', 'fun': layout1.constraintFunsUG},
                    options={'maxiter':1, 'eps':0.02,'ftol':1e-6, 'disp': True, 'iprint': 99})
    '''
    
    '''
    # ----- Constrained Optimization BY Linear Approximation (COBYLA)
    res = minimize(fun=layout1.objectiveFun, x0=xy, method='COBYLA',
                    constraints={'type': 'ineq', 'fun': layout1.constraintFuns},
                    options={'maxiter':2000,'catol': 1e-6, 'tol': 1e-6, 'disp': True})
    '''


    '''
    # ----- Differential Evolution (DE)
    # NonlinearConstraint
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.NonlinearConstraint.html#scipy.optimize.NonlinearConstraint
    #cons_fun = NonlinearConstraint(fun, lb, ub, jac='2-point', hess=<scipy.optimize._hessian_update_strategy.BFGS object>, keep_feasible=False, finite_diff_rel_step=None, finite_diff_jac_sparsity=None)
    cons_fun = NonlinearConstraint(layout1.constraintFuns, lb = 0, ub = np.inf)
    
    # takes FOREVER and is not the best solution
    res = differential_evolution(func=layout1.objectiveFun, bounds=boundary_xy, args=(), strategy='best1bin', 
                                 maxiter=1000, popsize=25, tol=0.01, mutation=(0.5, 1.0), 
                                 recombination=0.7, seed=None, callback=None, disp=True, polish=True, 
                                 init='latinhypercube', atol=0, updating='immediate', workers=1, 
                                 constraints=cons_fun, x0=xy, integrality=None, vectorized=False)
    '''


    # ----- Particle Swarm Optimization  
    
    # Other PSO (PSO with Scipy interface, but not that elaborated?)
    # https://github.com/jerrytheo/psopy?tab=readme-ov-file

    
    # Pyswarm (NOT pyswarms)
    # https://pythonhosted.org/pyswarm/
    
    if layouttype == 'freelayout':
        res, fopt = pso(layout1.objectiveFun, lb=boundary_xy[:,0], ub=boundary_xy[:,1], f_ieqcons=layout1.constraintFuns,  
                        swarmsize=20, omega=0.72984, phip=0.6, phig=0.8, maxiter=20, minstep=1e-8, minfunc=1e-8, debug=True)
    elif layouttype == 'uniformgridlayout':
        res, fopt = pso(layout1.objectiveFunUG, lb=boundaries_UG[:,0], ub=boundaries_UG[:,1], f_ieqcons=layout1.constraintFunsUG,  
                        swarmsize=20, omega=0.72984, phip=0.6, phig=0.8, maxiter=20, minstep=1e-8, minfunc=1e-8, debug=True)
   

    if layouttype == 'freelayout':
        layout1.updateLayoutOPT(res)  # make sure it is using the optimized layout => ONLY NEEDED WHEN INPUT WAS in km
        layout1.updateLayout(X=[], level=2, refresh=True)  # do a higher-fidelity update
        layout1.plotLayout(save=True)   
        
    elif layouttype == 'uniformgridlayout':  
        
        #optimized_xy_m = [1400, 1400, 500, 1000, 45, 0]
        
        layout1.updateLayoutUG(Xu=res, level=2, refresh=True)  # do a higher-fidelity update
        layout1.plotLayout(save=True) 


    
    plt.show()
    
   
    
    

########################## END
# ARCHIVE
