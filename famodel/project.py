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
from .platform.platform import Platform
from .mooring.anchor import Anchor

class Project():
    '''
    The overall object that defines a floating array for analysis and design 
    purposes. Its main model component is the RAFT model but it also includes
    more general information such as seabed info, metocean info, and...
    
    Ideally, this class can function even if RAFT and MoorPy are not loaded,
    at least for some high-level site processing functions.
    
    '''
    
    def __init__(self, lon=0, lat=0, file=None, depth=202):
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
        self.array = None  # RAFT model for dynamics analysis
        self.flow = None  # FLORIS interface instance for wake analysis
        
        # Lists describing the array, divided by structure type
        self.turbineList = []
        self.platformList = []
        self.mooringList = []  # A list of Mooring objects
        self.anchorList  = []
        self.cables = None  # CableSystem
        
        # Dictionaries of component/product properties used in the array
        self.turbineTypes = None
        self.lineTypes = None
        self.anchorTypes = None
        self.cableTypes = None
        
        
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
        self.soil_x     = None
        
        
        
        
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
        if not isinstance(d,dict):#if the input is not a dictionary, it is a yaml
            self.load(d)#load yaml into dictionary
        #d = 
        
        # ===== load FAM-specific model parts =====
        
        # cable types
        
        # dynamic cable basic properties (details are later via MoorPy)
        
        # ----- table of cables -----
        if 'array_cables' in d:
        
            cableInfo = [dict(zip( d['array_cables']['keys'], row))
                         for row in d['array_cables']['data']]
            
            # for ci in cableInfo:
            #     ...
                
            #     self.cables.addCable(...)
        
        # ----- cables one-by-one -----
        if 'cables' in d:
        
            for ci in d['cables']:
                ...
                
                self.cables.addCable(...)
        
        # ----- array mooring -----
        arrayMooring = {}
        arrayAnchor = {}
        if 'array_mooring' in d:
            # for mooring lines: save a list of dictionaries from each row in the data section
            arrayMooring = [dict(zip(d['array_mooring']['line_keys'], row)) for row in d['array_mooring']['line_data']]
            # for anchors: save a list of dictionaries from each row in the data section
            arrayAnchor = [dict(zip(d['array_mooring']['anchor_keys'], row)) for row in d['array_mooring']['anchor_data']]
        # ----- mooring systems ------
        mSystems = {}
        if 'mooring_systems' in d:
            for k, v in d['mooring_systems'].items():
                # set up mooring systems dictionary
                mSystems[k] = v
        
        # # load in shared mooring
        # if 'array_mooring' in d:
            
        
        # ----- mooring line section types ----- 
        self.lineTypes = {}
        
        if 'mooring_line_types' in d:
            # check if table format was used at all
            if 'keys' and 'data' in d['mooring_line_types']: # table-based
                dt = d['mooring_line_types'] # save location for code clarity
                # save a list of dictionaries from each row in the data section
                ms_info = [dict(zip(dt['keys'], row)) for row in dt['data']]
                # save the list into lineTypes dictionary and rename the index as the linetype name
                for k in range(0,len(ms_info)):
                    self.lineTypes[ms_info[k]['name']] = ms_info[k]
            # read in line types from list format as well(will overwrite any repeats from table)
            for k, v in d['mooring_line_types'].items():
                # set up line types dictionary
                self.lineTypes[k] = v
        
        # ----- mooring connectors -----
        connectorTypes = {}
        
        if 'mooring_connector_types' in d:
            for k, v in d['mooring_connector_types'].items():
                connectorTypes[k] = v
        
        # ----- anchor types -----
        self.anchorTypes = {}
        
        if 'anchor_types' in d:
            for k, v in d['anchor_types'].items():
                self.anchorTypes[k] = v
        
        # ----- mooring line configurations -----
        lineConfigs = {}
        
        if 'mooring_line_configs' in d:
            for k, v in d['mooring_line_configs'].items():
                # set up mooring config... (k, v)
                lineConfigs[k] = v
                # check line types listed in line configs matches those in linetypes section
                if self.lineTypes: # if linetypes section is included in dictionary
                    for j in range(0,len(v['sections'])): # loop through each line config section 
                        if not v['sections'][j]['type'] in self.lineTypes: # check if they match
                            raise Exception(f"Mooring line type '{v['sections'][j]['type']}' listed in mooring_line_configs is not found in mooring_line_types")
            # check line configurations listed in mooring systems matches those in line configs list
            if mSystems: # if mooring_systems section is included in dictionary
                for j,m_s in enumerate(mSystems): # loop through each mooring system
                    for i in range(0,len(mSystems[m_s]['data'])): # loop through each line listed in the system
                        if not mSystems[m_s]['data'][i][0] in lineConfigs: # check if they match
                            raise Exception(f"Mooring line configuration '{mSystems[m_s]['data'][i][0]}' listed in mooring_systems is not found in mooring_line_configs")

                # I think we want to just store it as a dictionary
                # validate it,
                # then have it available for use when making Mooring objects and subsystems
        def getMoorings(lineconfig):
            '''

            Parameters
            ----------
            lineconfig : string
                Line configuration type

            Returns
            -------
            m_config : dict
                mooring configuration dictionary

            '''
            # set up dictionary of information on the mooring configurations
            m_config = {'sections':[],'connectors':{},'anchor':{},'rAnchor':{},'zAnchor':{},'rFair':{},'zFair':{}}#,'EndPositions':{}}
            for k in range(0,len(lineConfigs[lineconfig]['sections'])): # loop through each section in the line
                
                lc = lineConfigs[lineconfig]['sections'][k] # set location for code clarity later
                lt = self.lineTypes[lc['type']] # set location for code clarity and brevity later
                # set up sub-dictionaries that will contain info on the line type
                m_config['sections'].append({'type':{'name':str(k)+'_'+lc['type'],'d_nom':lt['d_nom'],'material':lt['material'],'d_vol':lt['d_vol'],'m':lt['m'],'EA':float(lt['EA']),'MBL':lt['MBL'],'cost':lt['cost']}})
                # need to calculate the submerged weight of the line (not currently available in ontology yaml file)
                m_config['sections'][k]['type']['w'] = (lt['m']-np.pi/4*lt['d_vol']**2*1025)*9.81
                
                # add dynamic stretching if there is any
                if 'EAd' in lt: 
                    m_config['sections'][k]['type']['EAd'] = lt['EAd']
                    m_config['sections'][k]['type']['EAd_Lm'] = lt['EAd_Lm']
                    
                # set line length
                m_config['sections'][k]['length'] = lc['length']
                # make a connector dictionary, set info
                m_config['sections'][k]['connector'] = connectorTypes[lc['connector']]           
            
            # set general information on the whole line (not just a section/line type)
            # set to general depth first (will adjust to depth at anchor location after repositioning finds new anchor location)
            m_config['zAnchor'] = -self.depth 
            m_config['rAnchor'] = lineConfigs[lineconfig]['anchoring_radius']
            m_config['zFair'] = lineConfigs[lineconfig]['fairlead_depth']
            m_config['rFair'] = lineConfigs[lineconfig]['fairlead_radius']
   
            return(m_config)
        
        def getAnchors(lineAnch, mc=None,aNum=0):
            '''Create anchor design dictionary based on a given anchor type

            Parameters
            ----------
            lineAnch : string
                anchor type, for reference in the 'anchor_types' dictionary
            mc : mooring class instance, optional
                mooring class that the anchor is a part of (used only if not a shared anchor). The default is None.
            aNum : int, optional
                anchor ID in the anchor data table (used only if it is a shared anchor). The default is 0.

            Returns
            -------
            ad : dict
                anchor design dictionary

            '''
            ad = {'design':{}, 'cost':{}} 
            ad['design'] = self.anchorTypes[lineAnch]
            if self.soil_x: # get soil conditions at anchor location if soil info available
                if mc:
                    ad['soil_type'] = self.getSoilAtLocation(mc.rA[0], mc.rA[1])
                else:
                    ad['soil_type'] = self.getSoilAtLocation(arrayAnchor[aNum-1]['x'],arrayAnchor[aNum-1]['y'])
            ad['type'] = lineAnch[0:-2]
            ad['name'] = lineAnch
            
            return(ad)
        
        # ----- set up dictionary for each individual mooring line, create anchor, mooring, and platform classes ----
                
        # check that all necessary sections of design dictionary exist to create non-shared lines
        if self.lineTypes and lineConfigs and mSystems:
                    
            for i in range(0, len(d['array']['data'])): # loop through each platform in array
                
                # create platform instance (even if it only has shared moorings / anchors)
                self.platformList.append(Platform(r=[d['array']['data'][i][3],d['array']['data'][i][4]],heading=d['array']['data'][i][5]))
                # remove pre-set headings (need to append to this list so list should start off empty)
                self.platformList[-1].mooring_headings = []
                if d['array']['data'][i][2] > 0: #if no shared mooring on this platform
                    m_s = 'ms'+str(d['array']['data'][i][2]) # get mooring system number
                    
                    # get mooring headings (need this for platform class)
                    headings = []
                    for ii in range(0,len(mSystems[m_s]['data'])):
                        headings.append(np.radians(mSystems[m_s]['data'][ii][1])) 
                    
                    # add mooring headings to platform class instance
                    self.platformList[i].mooring_headings = headings
                    
                    # get the mooring line information 
                    for j in range(0,len(mSystems[m_s]['data'])): # loop through each line in the mooring system
                        # get the configuration for that line in the mooring system
                        lineconfig = mSystems[m_s]['data'][j][0] 
                   
                        # create mooring dictionary
                        m_config = getMoorings(lineconfig)
                        
                        # create mooring class instance as part of mooring list in the project class instance
                        mc = (Mooring(dd=m_config, rA=[m_config['rAnchor'],0,m_config['zAnchor']], rB=[m_config['rFair'],0,m_config['zFair']], rad_anch=m_config['rAnchor'], rad_fair=m_config['rFair'], z_anch=m_config['zAnchor'], z_fair=m_config['zFair']))
                        # adjust end positions based on platform location and mooring heading
                        mc.reposition(r_center=self.platformList[i].r, heading=headings[j], project=self)
                        # adjust anchor z location and rA based on location of anchor
                        zAnew, nAngle = self.getDepthAtLocation(mc.rA[0], mc.rA[1], return_n=True)
                        mc.rA[2] = -zAnew
                        mc.dd['zAnchor'] = -zAnew
                        mc.z_anch = -zAnew
                        
                        # set anchor info
                        lineAnch = mSystems[m_s]['data'][j][2] # get the anchor type for the line
                        ad = getAnchors(lineAnch, mc=mc) # call method to create anchor dictionary
                        ad['angle'] = nAngle
                        # add anchor class instance to anchorList in project class
                        self.anchorList.append(Anchor(dd=ad, r=mc.rA))
                        # add mooring class instance to list in anchor class
                        self.anchorList[-1].mooringList.append(mc)
                        # add mooring class instance to mooringlist in project class
                        self.mooringList.append(mc)                        
                        # # record index in project mooring list to mooring object
                        # self.mooringList[-1].mNum = len(self.mooringList)-1
                        # add mooring class instance to mooring list in the platform class instance
                        self.platformList[i].mooringList.append(mc)
                        # add 0 to boolean list (platform not connected to line end A)
                        self.platformList[i].endA.append(0)
                        # add anchor class instance to anchor list in platform class instance
                        self.platformList[i].anchorList.append(self.anchorList[-1])
        
        # ----- set up dictionary for each shared mooring line or shared anchor, create mooring and anchor classes ----
    
        # # create any shared anchor objects
        # if arrayAnchor:
        #     for i in range(0,len(arrayAnchor)):
        aNum = []
        
        # create any shared mooring lines / lines connected to shared anchors
        if arrayMooring:
            # get mooring line info for all lines 
            for j in range(0, len(arrayMooring)): # run through each line            
                # get configuration for that line 
                lineconfig = arrayMooring[j]['MooringConfigID']                       
                # create mooring dictionary for that line
                m_config = getMoorings(lineconfig)
                
                PFNum = [] # platform ID(s) connected to the mooring line
                
                # determine if end A is an anchor or a platform
                if arrayMooring[j]['end A'][0:4] == 'FOWT': # shared mooring line (no anchor)
                    # get ID of platforms connected to line
                    PFNum.append(int(arrayMooring[j]['end B'][-1])-1)
                    PFNum.append(int(arrayMooring[j]['end A'][-1])-1)
                    # set locations of row for platform connected to end A and end B for simpler code
                    rowB = d['array']['data'][PFNum[0]]
                    rowA = d['array']['data'][PFNum[1]]
                    # calculate fairlead locations (can't use reposition method because both ends need separate repositioning)
                    Aloc = [rowA[3]+np.cos(np.radians(arrayMooring[j]['headingA']))*m_config['rFair'], rowA[4]+np.sin(np.radians(arrayMooring[j]['headingA']))*m_config['rFair'], m_config['zFair']]
                    Bloc = [rowB[3]+np.cos(np.radians(arrayMooring[j]['headingB']))*m_config['rFair'], rowB[4]+np.sin(np.radians(arrayMooring[j]['headingB']))*m_config['rFair'], m_config['zFair']]
                    # create mooring class instance
                    mc = (Mooring(dd=m_config, rA=Aloc, rB=Bloc, rad_fair=m_config['rFair'], z_fair=m_config['zFair'], rad_anch=m_config['rAnchor'], z_anch=m_config['zAnchor']))
                    # mc.sPF = PFNum
                else: # end A is an anchor
                    # get ID of platform connected to line
                    PFNum.append(int(arrayMooring[j]['end B'][-1])-1)
                    # create mooring class instance
                    mc = (Mooring(dd=m_config, rA=[m_config['rAnchor'],0,m_config['zAnchor']], rB=[m_config['rFair'],0,m_config['zFair']], rad_anch=m_config['rAnchor'], rad_fair=m_config['rFair'], z_anch=m_config['zAnchor'], z_fair=m_config['zFair']))
                    # adjust end positions based on platform location and mooring heading
                    mc.reposition(r_center=self.platformList[PFNum[0]].r, heading=arrayMooring[j]['headingB'], project=self, degrees=True)

                    # check if anchor instance already exists
                    if any(tt == int(arrayMooring[j]['end A'][-1]) for tt in aNum): # anchor exists
                        # find anchor class instance
                        for k in range(0,len(self.anchorList)):
                            if self.anchorList[k].aNum:
                                if self.anchorList[k].aNum == int(arrayMooring[j]['end A'][-1]):
                                    sharedAnch = k
                        # add mooring object to list in anchor class
                        self.anchorList[sharedAnch].mooringList.append(mc)
                        # add anchor object to list in platform class
                        self.platformList[PFNum[0]].anchorList.append(self.anchorList[sharedAnch])
                    else:
                        # set line anchor type and get dictionary of anchor information
                        lineAnch = arrayAnchor[j]['type']
                        aNum.append(int(arrayMooring[j]['end A'][-1])) # get anchor number
                        ad = getAnchors(lineAnch,aNum=aNum[-1]) # call method to create dictionary
                        # get exact anchor x,y location from arrayAnchor 
                        aloc = [arrayAnchor[aNum[-1]-1]['x'],arrayAnchor[aNum[-1]-1]['y']] 
                        # adjust anchor location and rA based on location of anchor
                        zAnew, nAngle = self.getDepthAtLocation(aloc[0], aloc[1], return_n=True)
                        mc.rA = [aloc[0],aloc[1],-zAnew]
                        mc.dd['zAnchor'] = -zAnew
                        mc.z_anch = -zAnew
                        # create anchor object
                        self.anchorList.append(Anchor(dd=ad, r=[aloc[0],aloc[1],-zAnew], aNum=aNum[-1]))
                        # add mooring object to anchor mooring list
                        self.anchorList[-1].mooringList.append(mc)
                        # add anchor object to platform anchor list
                        self.platformList[PFNum[0]].anchorList.append(self.anchorList[-1])
                        # mc.anchor = self.anchorList[aNum]
                # add mooring object to project mooring list
                self.mooringList.append(mc)
                # record index in project mooring list to mooring object
                # self.mooringList[-1].mNum = len(self.mooringList)-1
                # add mooring object to platform mooring list    
                self.platformList[PFNum[0]].mooringList.append(mc)
                self.platformList[PFNum[0]].endA.append(0)
                # add heading
                self.platformList[PFNum[0]].mooring_headings.append(np.radians(arrayMooring[j]['headingB']))
                if len(PFNum)>1: # if shared line
                    # record shared line on the other platform as well
                    self.platformList[PFNum[1]].mooringList.append(mc)
                    self.platformList[PFNum[1]].mooring_headings.append(np.radians(arrayMooring[j]['headingA'])) # add heading
                    self.platformList[PFNum[1]].endA.append(1)
        # ===== load RAFT model parts =====
        

    # ----- Site conditions processing functions -----

    def loadSite(self, site):
        '''Load site information from a dictionary or YAML file
        (specified by input). This should be the site portion of
        the floating wind array ontology.'''
        # standard function to load dict if input is yaml
        
        # load general information
        self.depth = getFromDict(site['general'], 'water_depth', default=100)
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
       
    def getMoorPyArray(self,plt=0):
        '''

        Parameters
        ----------
        plt : boolean, optional
            Controls whether to create a plot of the MoorPy array. 1=create plot, 0=no plot The default is 0.

        Returns
        -------
        ms : class instance
            MoorPy system for the whole array based on the subsystems in the mooringList

        '''
        ms = mp.System(depth=self.depth)
        
        
        for i in range(0,len(self.platformList)): # make all the bodies up front
            PF = self.platformList[i]
            # add a moorpy body at the correct location
            r6 = [PF.r[0],PF.r[1],0,0,0,0]
            ms.addBody(0,r6,m=1.784e7,v=20206,rM=100,AWP=1011)
        
        # create anchor points and all mooring lines connected to the anchors (since all connected to anchors, can't be a shared mooring)
        for i in range(0,len(self.anchorList)):
            ssloc = []
            for j in range(0,len(self.anchorList[i].mooringList)):
                # create subsystem
                self.anchorList[i].mooringList[j].createSubsystem()
                # set location of subsystem for simpler coding
                ssloc.append(self.anchorList[i].mooringList[j].subsystem)
                # add subsystem as a line in moorpy system
                ms.lineList.append(ssloc[j])
                ssloc[j].number = len(ms.lineList)
                # add fairlead point and attach the line to it
                ms.addPoint(1,ssloc[j].rB)
                ms.pointList[-1].attachLine(ssloc[j].number,1)
                # find associated platform and attach body to point (since not a shared line, should only be one platform with this mooring object)
                for k in range(0,len(self.platformList)):
                    if self.anchorList[i].mooringList[j] in self.platformList[k].mooringList:
                        PF = self.platformList[k]
                        PFNum = k
                # attach rB point to platform
                ms.bodyList[PFNum].attachPoint(len(ms.pointList),[ssloc[j].rB[0]-PF.r[0],ssloc[j].rB[1]-PF.r[1],ssloc[j].rB[2]])#attach to fairlead
            # create anchor point
            self.anchorList[i].makeMoorPyAnchor(ms)
            # attach line(s) to anchor point
            for j in range(0,len(ssloc)):
                ms.pointList[-1].attachLine(ssloc[j].number,0)
                
        # now create and attach any shared lines
        for i in range(0,len(self.mooringList)): # loop through all lines
            if self.mooringList[i].subsystem in ms.lineList: # check if line has already been put in ms
                pass
            else: # new shared line
                # create subsystem for shared line
                self.mooringList[i].createSubsystem(case=1)
                # set location of subsystem for simpler coding
                ssloc = self.mooringList[i].subsystem
                # add subsystem as a line in moorpy system
                ms.lineList.append(ssloc)
                ssloc.number = len(ms.lineList)               
                
                # find associated platforms
                PF = []
                PFNum = []
                idx = []
                for k in range(0,len(self.platformList)):
                    if self.mooringList[i] in self.platformList[k].mooringList:
                        PF.append(self.platformList[k])
                        PFNum.append(k)
                        # find index of mooring object in platform mooring list
                        idx.append(PF[-1].mooringList.index(self.mooringList[i]))
                        
                # add fairlead point A and attach the line to it
                ms.addPoint(1,ssloc.rA)
                ms.pointList[-1].attachLine(ssloc.number,0)
                if PF[0].endA[idx[0]] == 1: # end A connected to PF[0], end B connected to PF[1]
                    # connect line end A to the body
                    ms.bodyList[PFNum[0]].attachPoint(len(ms.pointList),[ssloc.rA[0]-PF[0].r[0],ssloc.rA[1]-PF[0].r[1],ssloc.rA[2]])
                    # add fairlead point 2 and attach the line to it
                    ms.addPoint(1,ssloc.rB)
                    ms.pointList[-1].attachLine(ssloc.number,1)
                    # connect line end B to the body
                    ms.bodyList[PFNum[1]].attachPoint(len(ms.pointList),[ssloc.rB[0]-PF[1].r[0],ssloc.rB[1]-PF[1].r[1],ssloc.rB[2]])
                else: # end A connected to PF[1], end B connected to PF[0]
                    # connect line end A to the body
                    ms.bodyList[PFNum[1]].attachPoint(len(ms.pointList),[ssloc.rA[0]-PF[1].r[0],ssloc.rA[1]-PF[1].r[1],ssloc.rA[2]])
                    # add fairlead point 2 and attach the line to it
                    ms.addPoint(1,ssloc.rB)
                    ms.pointList[-1].attachLine(ssloc.number,1)
                    # connect line end B to the body
                    ms.bodyList[PFNum[0]].attachPoint(len(ms.pointList),[ssloc.rB[0]-PF[0].r[0],ssloc.rB[1]-PF[0].r[1],ssloc.rB[2]])

        # initialize, solve equilibrium, and plot the system       
        ms.initialize()
        ms.solveEquilibrium()
        if plt:
            settings = {}
            # settings["linelabels"] = True
            # settings["pointlabels"] = True
            ms.plot( **settings)
            
        return(ms)                            

                
                
                
            
        
        # for i in range(0,len(self.platformList)): # loop through for each body/platform
        #     PF = self.platformList[i]
        #     # get number of mooring lines for this platform
        #     nLines = PF.n_mooring
        #     # loop through each line on the body
        #     for j in range(0,nLines):
                
        #         # value = self.platformList[i].mooringList[j]
        #         # pfs = [x for x in self.platformList[0:-1].mooringList if value in x]
                
        #         # check if there's an anchor class associated with the line
        #         if self.platformList[i].mooringList[j].anchor: # not a shared line
        #             # record line index from project level
        #             mCheck.append(self.platformList[i].mooringList[j].mNum)
        #             # get anchor number from project level
        #             aNum = self.platformList[i].mooringList[j].anchor.aNum
                    
        #             # create subsystem
        #             self.mooringList[mCheck[-1]].createSubsystem()
        #             # set the location of subsystem to make next few lines shorter
        #             ssloc = self.mooringList[mCheck[-1]].subsystem
        #             # add subsystem to mp system line list
        #             ms.lineList.append(ssloc)
        #             ssloc.number = len(mCheck)
                    
        #             if any(a==aNum for a in aCheck): # check if anchor has already been made
        #                 # shared anchor, anchor has already been created in moorpy array
        #                 # find point list index for the shared anchor
        #                 pass
        #             else:
        #                 # record anchor index 
        #                 aCheck.append(aNum)
        #                 # create anchor moorpy object for anchor point
        #                 ms = self.mooringList[mCheck[-1]].anchor.makeMoorPyAnchor(ms)
                        
                    
        #             ms.pointList[-1].attachLine(ssloc.number,0)
        #             ms.addPoint(1,ssloc.rB) # fairlead
        #             ms.pointList[-1].attachLine(ssloc.number,1)
        #             # attach body to fairlead line
        #             # Note: must use an rB location that is centered on 0, otherwise the rB location will be incorrect after initialization
        #             ms.bodyList[i].attachPoint(len(ms.pointList),[ssloc.rB[0]-PF.r[0],ssloc.rB[1]-PF.r[1],ssloc.rB[2]])#attach to fairlead
        #         else: # this is a shared mooring line then
        #             mNum = self.platformList[i].mooringList[j].mNum
        #             if any(m==mNum for m in mCheck):
        #                 # shared line has already been created
        #                 pass
        #             else:
        #                 # record mooring index
        #                 mCheck.append(mNum)
        #                 # get platform number attached to point A
        #                 PFA = self.mooringList[mCheck[-1]].sPF
        #                 # create subsystem
        #                 self.mooringList[mCheck[-1]].createSubsystem(case=1)
        #                 # add subsystem to mp system line list
        #                 ms.lineList.append(ssloc)
        #                 ssloc.number = len(mCheck)
        #                 # add fairlead 1 point (instead of anchor point)
        #                 ms.addPoint(1,ssloc.rA)
        #                 ms.pointList[-1].attachLine(ssloc.number,0)
        #                 ms.bodyList[PFA].attachPoint(len(ms.pointList),[ssloc.rA[0]-self.platformList[PFA].r[0],ssloc.rA[1]-self.platformList[PFA].r[1],ssloc.rA[2]])
        #                 ms.addPoint(1,ssloc.rB) # fairlead 2 point
        #                 ms.pointList[-1].attachLine(ssloc.number,1)
        #                 # attach body to fairlead line
        #                 # Note: must use an rB location that is centered on 0, otherwise the rB location will be incorrect after initialization
        #                 ms.bodyList[i].attachPoint(len(ms.pointList),[ssloc.rB[0]-PF.r[0],ssloc.rB[1]-PF.r[1],ssloc.rB[2]])#attach to fairlead
                
                
        #         # increase counter
        #         counter = counter + 1

        

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
