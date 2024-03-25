"""Project class for FAModel, containing information and key methods for
the site information and design information that make up a project."""

import os
import numpy as np
import matplotlib.pyplot as plt

import moorpy as mp
from moorpy.helpers import set_axes_equal
from moorpy import helpers
import yaml
from copy import deepcopy
import raft

#from shapely.geometry import Point, Polygon, LineString

from .anchors.anchor_capacity import anchorCapacity
from .seabed import seabed_tools as sbt
from .mooring.mooring import Mooring
from .platform.platform import Platform
from .mooring.anchor import Anchor
from .mooring.connector import Connector

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
        self.mooringListPristine = [] # A list of Mooring objects in initial condition (no marine growth, corrosion, etc)
        self.anchorList  = []
        self.connectorList = []
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
        self.soil_y     = None
        
        # MoorPy system associated with the project
        self.ms  = None
        
        
        
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
            if 'line_data' in d['array_mooring']:
                if d['array_mooring']['line_data']:                    
                    arrayMooring = [dict(zip(d['array_mooring']['line_keys'], row)) for row in d['array_mooring']['line_data']]
            # for anchors: save a list of dictionaries from each row in the data section
            if 'anchor_data' in d['array_mooring']:
                if d['array_mooring']['anchor_data']:
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
                # set up mooring config
                lineConfigs[k] = v
                # check line types listed in line configs matches those in linetypes section
                if self.lineTypes: # if linetypes section is included in dictionary
                    for j in range(0,len(v['sections'])): # loop through each line config section
                        if 'type' in v['sections'][j]: # check if it is a connector or line config
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
            c_config : dict
                connector configuration dictionary

            '''
            # set up dictionary of information on the mooring configurations
            m_config = {'sections':[],'anchor':{},'rAnchor':{},'zAnchor':{},'rFair':{},'zFair':{}}#,'EndPositions':{}}
            # set up connector dictionary
            c_config = []
                        
            lineLast = 1    # boolean whether item with index k-1 is a line. Set to 1 for first run through of for loop
            ct = 0   # counter for number of line types
            for k in range(0,len(lineConfigs[lineconfig]['sections'])): # loop through each section in the line
            
                lc = lineConfigs[lineconfig]['sections'][k] # set location for code clarity later
                # determine if it's a line type or a connector listed
                if 'type' in lc: 
                    # this is a line
                    if lineLast: # previous item in list was a line (or this is the first item in a list)
                        # no connector was specified for before this line - add an empty connector
                        c_config.append(None)                        
                    # set line information                                                
                    lt = self.lineTypes[lc['type']] # set location for code clarity and brevity later
                    # set up sub-dictionaries that will contain info on the line type
                    m_config['sections'].append({'type':{'name':str(ct)+'_'+lc['type'],'d_nom':lt['d_nom'],'material':lt['material'],'d_vol':lt['d_vol'],'m':lt['m'],'EA':float(lt['EA'])}})
                    # need to calculate the submerged weight of the line (not currently available in ontology yaml file)
                    m_config['sections'][ct]['type']['w'] = (lt['m']-np.pi/4*lt['d_vol']**2*1025)*9.81
                    # add cost if given
                    if 'cost' in lt:
                        m_config['sections'][ct]['type']['cost'] = lt['cost']
                    # add MBL if given
                    if 'MBL' in lt:
                        m_config['sections'][ct]['type']['MBL'] = lt['MBL']
                    # add dynamic stretching if there is any
                    if 'EAd' in lt: 
                        m_config['sections'][ct]['type']['EAd'] = lt['EAd']
                        m_config['sections'][ct]['type']['EAd_Lm'] = lt['EAd_Lm']
                        
                    # set line length
                    m_config['sections'][ct]['length'] = lc['length']
                    # update counter for line types 
                    ct = ct + 1
                    # update line last boolean
                    lineLast = 1
                    
                elif 'connectorType' in lc:
                    # this is a connector
                    if lineLast == 0:
                        # last item in list was a connector
                        raise Exception(f"Two connectors were specified in a row for line configuration '{lineconfig}', please remove one of the connectors")
                    else:
                        # last item in list was a line
                        c_config.append(connectorTypes[lc['connectorType']]) # add connector to list
                        c_config[-1]['type'] = lc['connectorType']
                        # update lineLast boolean
                        lineLast = 0
                else:
                    # not a connector or a line
                    raise Exception(f"Please make sure that all section entries for line configuration '{lineconfig}' are either line sections (which must have a 'type' key) or connectors (which must have a 'connectorType' key")

            # check if line is a shared symmetrical configuration
            if 'symmetric' in lineConfigs[lineconfig] and lineConfigs[lineconfig]['symmetric']:
                if not lineLast: # check if last item in line config list was a connector
                    for ii in range(0,ct):
                        # set mooring configuration 
                        m_config['sections'].append(m_config['sections'][-1-2*ii])
                        # set connector (since it's mirrored, connector B becomes connector A)
                        c_config.append(c_config[-2-2*ii])
                else: # double the length of the end line
                    m_config['sections'][-1]['length'] = m_config['sections'][-1]['length']*2
                    # set connector B for line same as previous listed connector
                    c_config.append(c_config[-1])
                    for ii in range(0,ct-1): # go through every line config except the last (since it was doubled already)
                        # set mooring configuration
                        m_config['sections'].append(m_config['sections'][-2-2*ii])
                        # set connector
                        c_config.append(c_config[-3-2*ii])
            else: # if not a symmetric line, check if last item was a line (if so need to add another empty connector)
                if lineLast:
                    # add an empty connector object
                    c_config.append(None)
            # set general information on the whole line (not just a section/line type)
            # set to general depth first (will adjust to depth at anchor location after repositioning finds new anchor location)
            m_config['zAnchor'] = -self.depth 
            m_config['rAnchor'] = lineConfigs[lineconfig]['anchoring_radius']
            m_config['zFair'] = lineConfigs[lineconfig]['fairlead_depth']
            m_config['rFair'] = lineConfigs[lineconfig]['fairlead_radius']
   
            return(m_config, c_config)
        
        def getConnectors(c_config):
            '''

            Parameters
            ----------
            c_config : dict
                Dictionary of connector configurations for a mooring line.

            Returns
            -------
            None.

            '''
            # make connector objects for all sections of a mooring line configuration in order
            for i in range(0,len(c_config)):
                # check if connector is a none-type
                if c_config[i] == None:
                    # create empty connector object
                    self.connectorList.append(Connector())
                else:
                    # create connector object
                    self.connectorList.append(Connector(dd=c_config[i]))
                # add connector to mooring object list
                self.mooringList[-1].connectorList.append(self.connectorList[-1])
        
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
            ad['type'] = lineAnch[0:-1]
            ad['name'] = lineAnch
            
            return(ad)
        
        
        
        # ----- set up dictionary for each individual mooring line, create anchor, mooring, and platform classes ----
                
        # check that all necessary sections of design dictionary exist to create non-shared lines
        if self.lineTypes and lineConfigs and mSystems:
                    
            for i in range(0, len(d['array']['data'])): # loop through each platform in array
                
                # create platform instance (even if it only has shared moorings / anchors)
                self.platformList.append(Platform(r=[d['array']['data'][i][4],d['array']['data'][i][5]],heading=d['array']['data'][i][6]))
                # remove pre-set headings (need to append to this list so list should start off empty)
                self.platformList[-1].mooring_headings = []
                if not d['array']['data'][i][3] == 0: #if no shared mooring on this platform
                    m_s = d['array']['data'][i][3] # get mooring system number
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
                   
                        # create mooring and connector dictionary
                        m_config, c_config = getMoorings(lineconfig)
                        
                        # create mooring class instance as part of mooring list in the project class instance
                        mc = (Mooring(dd=m_config, rA=[m_config['rAnchor'],0,m_config['zAnchor']], rB=[m_config['rFair'],0,m_config['zFair']], rad_anch=m_config['rAnchor'], rad_fair=m_config['rFair'], z_anch=m_config['zAnchor'], z_fair=m_config['zFair']))
                        # adjust end positions based on platform location and mooring and platform headings
                        mc.reposition(r_center=self.platformList[i].r, heading=headings[j]+self.platformList[i].phi, project=self)
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
                        # create connector dictionaries and objects for line 
                        getConnectors(c_config)                        
                        # add mooring class instance to mooring list in the platform class instance
                        self.platformList[i].mooringList.append(mc)
                        # add 0 to boolean list (platform not connected to line end A)
                        self.platformList[i].endA.append(0)
                        # add anchor class instance to anchor list in platform class instance
                        self.platformList[i].anchorList.append(self.anchorList[-1])
                        
        
        # ----- set up dictionary for each shared mooring line or shared anchor, create mooring and anchor classes ----
    
        aNum = []
        
        # create any shared mooring lines / lines connected to shared anchors
        if arrayMooring:
            # get mooring line info for all lines 
            for j in range(0, len(arrayMooring)): # run through each line            
                # get configuration for that line 
                lineconfig = arrayMooring[j]['MooringConfigID']                       
                # create mooring and connector dictionary for that line
                m_config, c_config = getMoorings(lineconfig)
                
                PFNum = [] # platform ID(s) connected to the mooring line
                
                # Error check for putting an anchor (or something else) at end B
                if not arrayMooring[j]['end B'][0:4].upper() == 'FOWT':
                    raise Exception(f"input for end B of line_data table row '{j}' in array_mooring must be 'FOWT #'. Any anchors should be listed as end A.")
                # determine if end A is an anchor or a platform
                if arrayMooring[j]['end A'][0:4].upper() == 'FOWT': # shared mooring line (no anchor)
                    # get ID of platforms connected to line
                    PFNum.append(int(arrayMooring[j]['end B'][-1])-1)
                    PFNum.append(int(arrayMooring[j]['end A'][-1])-1)
                    # set locations of row for platform connected to end A and end B for simpler code
                    rowB = d['array']['data'][PFNum[0]]
                    rowA = d['array']['data'][PFNum[1]]
                    # get headings (mooring heading combined with platform heading)
                    headingA = np.radians(arrayMooring[j]['headingA']) + self.platformList[PFNum[1]].phi - np.pi
                    headingB = np.radians(arrayMooring[j]['headingB']) + self.platformList[PFNum[0]].phi - np.pi
                    # calculate fairlead locations (can't use reposition method because both ends need separate repositioning)
                    Aloc = [rowA[4]+np.cos(headingA)*m_config['rFair'], rowA[5]+np.sin(headingA)*m_config['rFair'], m_config['zFair']]
                    Bloc = [rowB[4]+np.cos(headingB)*m_config['rFair'], rowB[5]+np.sin(headingB)*m_config['rFair'], m_config['zFair']]
                    # create mooring class instance
                    mc = (Mooring(dd=m_config, rA=Aloc, rB=Bloc, rad_fair=m_config['rFair'], z_fair=m_config['zFair'], rad_anch=m_config['rAnchor'], z_anch=m_config['zAnchor']))
                    mc.shared = 1
                elif arrayMooring[j]['end A'][0:4].upper() == 'ANCH': # end A is an anchor
                    # check if anchor number exists in anchor_data table
                    if int(arrayMooring[j]['end A'][-1]) > len(arrayAnchor):
                        raise Exception(f"anchor number {arrayMooring[j]['end A'][-1]} listed in array_mooring line_data table row {j} does not exist in anchor_data table")
                    # get ID of platform connected to line
                    PFNum.append(int(arrayMooring[j]['end B'][-1])-1)
                    # create mooring class instance
                    mc = (Mooring(dd=m_config, rA=[m_config['rAnchor'],0,m_config['zAnchor']], rB=[m_config['rFair'],0,m_config['zFair']], rad_anch=m_config['rAnchor'], rad_fair=m_config['rFair'], z_anch=m_config['zAnchor'], z_fair=m_config['zFair']))
                    # adjust end positions based on platform location and mooring and platform heading
                    mc.reposition(r_center=self.platformList[PFNum[0]].r, heading=np.radians(arrayMooring[j]['headingB'])+self.platformList[PFNum[0]].phi, project=self)

                    # check if anchor instance already exists
                    if any(tt == int(arrayMooring[j]['end A'][-1]) for tt in aNum): # anchor exists
                        # find anchor class instance
                        for k in range(0,len(self.anchorList)):
                            if self.anchorList[k].aNum and self.anchorList[k].aNum == int(arrayMooring[j]['end A'][-1]):
                                # add mooring object to list in anchor class
                                self.anchorList[k].mooringList.append(mc)
                                # add anchor object to list in platform class
                                self.platformList[PFNum[0]].anchorList.append(self.anchorList[k])
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
                else: # error in input
                    raise Exception(f"end A input in array_mooring line_data table line '{j}' must be either 'Anch #' or 'FOWT #', where  # refers to the row number in the anchor_data or array data table respectively.")

                # add mooring object to project mooring list
                self.mooringList.append(mc)
                
                # create connector dictionaries and objects for line 
                getConnectors(c_config)

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
        
        # create a deepcopy of the mooring list to preserve original in case marine growth, corrosion, or other changes made
        self.mooringListPristine = deepcopy(self.mooringList)    
        
        # ===== load RAFT model parts =====
        # load info into RAFT dictionary and create RAFT model
        RAFTDict = {}
        # load turbine dictionary into RAFT dictionary
        if 'turbine' in d and d['turbine']:            
            # check that there is only one turbine
            if 'turbines' in d and d['turbines']:
                raise Exception("Cannot read in items for both 'turbines' and 'turbine' keywords. Use either 'turbine' keyword for one turbine or 'turbines' keyword for a list of turbines.")
            elif type(d['turbine']) is list and len(d['turbine'])>1:
                raise Exception("'turbine' section keyword must be changed to 'turbines' if multiple turbines are listed")
            else:
                RAFTDict['turbine'] = d['turbine']
        # load list of turbine dictionaries into RAFT dictionary
        elif 'turbines' in d and d['turbines']:
            RAFTDict['turbines'] = d['turbines']
        
        # load platform dictionary into RAFT dictionary if only one platform      
        if 'platform' in d and d['platform']:
            # check that there is only one platform
            if 'platforms' in d and d['platforms']:
                raise Exception("Cannot read in items for both 'platforms' and 'platform' keywords. Use either 'platform' keyword for one platform or 'platforms' keyword for a list of platforms.")
            elif type(d['platform']) is list and len(d['platform'])>1:
                raise Exception("'platform' section keyword must be changed to 'platforms' if multiple platforms are listed")
            else:
                RAFTDict['platform'] = d['platform']
        # load list of platform dictionaries into RAFT dictionary
        elif 'platforms' in d and d['platforms']:
            RAFTDict['platforms'] = d['platforms']
        
        # load global RAFT settings into RAFT dictionary
        if 'RAFT_Settings' in d['site'] and d['site']['RAFT_Settings']:
            RAFTDict['settings'] = d['site']['RAFT_Settings']
        # load RAFT cases into RAFT dictionary
        if 'RAFT_cases' in d['site'] and d['site']['RAFT_cases']:
            RAFTDict['cases'] = d['site']['RAFT_cases']
        
        # load array information into RAFT dictionary
        RAFTDict['array'] = deepcopy(d['array']) # need to change items so make a deepcopy
        RAFTDict['array']['keys'].pop(0) # remove key for ID because this doesn't exist in RAFT array table
        for i in range(0,len(d['array']['data'])):
            RAFTDict['array']['data'][i][3] = 0 # make mooringID = 0 (mooring data will come from MoorPy)
            RAFTDict['array']['data'][i].pop(0) # remove ID column because this doesn't exist in RAFT array data table
        # load general site info to RAFT dictionary
        RAFTDict['site'] = {'water_depth':self.depth,'rho_water':self.rho_water,'rho_air':self.rho_air,'mu_air':self.mu_air}
        RAFTDict['site']['shearExp'] = getFromDict(d['site']['general'],'shearExp',default=0.12)
        
        # create a name for the raft model
        RAFTDict['name'] = 'Project_Array'
        RAFTDict['type'] = 'input file for RAFT'
        self.array = raft.Model(RAFTDict)
        # create moorpy array if it doesn't exist
        if not self.ms:
            self.getMoorPyArray()
        # assign moorpy array to RAFT object
        self.array.ms = self.ms
        
        # connect RAFT fowt to the correct moorpy body
        for i in range(0,len(self.ms.bodyList)):
            self.array.fowtList[i].body = self.ms.bodyList[i]
        

        
        


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
            if 'file' in site['bathymetry'] and site['bathymetry']['file']: # make sure there was a file provided even if the key is there
                self.loadBathymetry(site['bathymetry']['file'])
            elif 'x_y_z' in site['bathymetry'] and site['bathymetry']['x_y_z']:
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
            else:
                # assume a flat bathymetry
                self.grid_depth  = np.array([[self.depth]])

        # Load project boundary, if provided
        if 'boundaries' in site:
            if 'file' in site['boundaries'] and site['boundaries']['file']:  # load boundary data from file if filename provided
                    self.loadBoundary(site['boundaries']['file'])
            elif 'x_y' in site['boundaries'] and site['boundaries']['x_y']:  # process list of boundary x,y vertices
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
        xs = np.arange(-2000,2000,500)
        ys = np.arange(-2000,2000,500)
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
        
        
        # # also if there are rocky bits... (TEMPORARY)
        # X, Y = np.meshgrid(self.soil_x, self.soil_y)
        # ax.scatter(X, Y, c=self.soil_rocky, s=4, cmap='cividis_r', vmin=-0.5, vmax=1.5)
        
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
       
    def getMoorPyArray(self,plt=0, mgDict=None):
        '''Creates an array in moorpy from the mooring, anchor, connector, and platform objects in the array.

        Parameters
        ----------
        plt : boolean, optional
            Controls whether to create a plot of the MoorPy array. 1=create plot, 0=no plot The default is 0.

        mgDict : dictionary, optional
            Dictionary of marine growth information for the location
    
        Returns
        -------
        ms : class instance
            MoorPy system for the whole array based on the subsystems in the mooringList

        '''
        def createArray(pristineLines=0):
            '''Creates the moorpy array
            Parameters
            ----------
            pristineLines : boolean, optional
                Describes whether the created subsystems are pristine and should be added the the mooringListPristine objects
                Default is 0 (do not add subsystems to pristine mooring objects)

            Returns
            -------
            None.

            '''
            # create MoorPy system
            self.ms = mp.System(depth=self.depth)       
            
            for i in range(0,len(self.platformList)): # make all the bodies up front
                PF = self.platformList[i]
                # add a moorpy body at the correct location
                r6 = [PF.r[0],PF.r[1],0,0,0,0]
                self.ms.addBody(-1,r6,m=19911423.956678286,rCG=np.array([ 1.49820657e-15,  1.49820657e-15, -2.54122031e+00]),v=19480.104108645974,rM=np.array([2.24104273e-15, 1.49402849e-15, 1.19971829e+01]),AWP=446.69520543229874)
            
            # create anchor points and all mooring lines connected to the anchors (since all connected to anchors, can't be a shared mooring)
            for i in range(0,len(self.anchorList)):
                ssloc = []
                for j in range(0,len(self.anchorList[i].mooringList)):
                    # create subsystem
                    self.anchorList[i].mooringList[j].createSubsystem()
                    # set location of subsystem for simpler coding
                    ssloc.append(self.anchorList[i].mooringList[j].subsystem)
                    self.ms.lineList.append(ssloc[j])
                    ssloc[j].number = len(self.ms.lineList)
                    # add fairlead point
                    self.ms.addPoint(1,ssloc[j].rB)
                    # add connector info for fairlead point
                    self.ms.pointList[-1].m = self.ms.lineList[-1].pointList[-1].m 
                    self.ms.pointList[-1].v = self.ms.lineList[-1].pointList[-1].v
                    self.ms.pointList[-1].CdA = self.ms.lineList[-1].pointList[-1].CdA
                    # attach the line to point
                    self.ms.pointList[-1].attachLine(ssloc[j].number,1)
                    # find associated platform and attach body to point (since not a shared line, should only be one platform with this mooring object)
                    for k in range(0,len(self.platformList)):
                        if self.anchorList[i].mooringList[j] in self.platformList[k].mooringList:
                            PF = self.platformList[k]
                            PFNum = k
                    # attach rB point to platform (need to subtract out location of platform from point for subsystem integration to work correctly)
                    self.ms.bodyList[PFNum].attachPoint(len(self.ms.pointList),[ssloc[j].rB[0]-PF.r[0],ssloc[j].rB[1]-PF.r[1],ssloc[j].rB[2]])#attach to fairlead
                # create anchor point
                self.anchorList[i].makeMoorPyAnchor(self.ms)
                # attach line(s) to anchor point
                for j in range(0,len(ssloc)):
                    self.ms.pointList[-1].attachLine(ssloc[j].number,0)
            
            check = np.ones((len(self.mooringList),1))
            # now create and attach any shared lines
            for i in range(0,len(self.mooringList)): # loop through all lines
                for j in range(0,len(self.anchorList)):
                    if self.mooringList[i] in self.anchorList[j].mooringList: # check if line has already been put in ms
                        check[i] = 0     
                if check[i] == 1: # mooring object not in any anchor lists
                    # new shared line
                    # create subsystem for shared line
                    self.mooringList[i].createSubsystem(case=1) # we doubled all symmetric lines so any shared lines should be case 1
                    # set location of subsystem for simpler coding
                    ssloc = self.mooringList[i].subsystem
                    # add subsystem as a line in moorpy system
                    self.ms.lineList.append(ssloc)
                    ssloc.number = len(self.ms.lineList)               
                    
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
                    self.ms.addPoint(1,ssloc.rA)
                    self.ms.pointList[-1].attachLine(ssloc.number,0)
                    if PF[0].endA[idx[0]] == 1: # end A connected to PF[0], end B connected to PF[1]
                        # connect line end A to the body
                        self.ms.bodyList[PFNum[0]].attachPoint(len(self.ms.pointList),[ssloc.rA[0]-PF[0].r[0],ssloc.rA[1]-PF[0].r[1],ssloc.rA[2]])
                        # add fairlead point 2
                        self.ms.addPoint(1,ssloc.rB)
                        # add connector info for fairlead point
                        self.ms.pointList[-1].m = self.ms.lineList[-1].pointList[-1].m 
                        self.ms.pointList[-1].v = self.ms.lineList[-1].pointList[-1].v
                        self.ms.pointList[-1].CdA = self.ms.lineList[-1].pointList[-1].CdA
                        # attach the line to point
                        self.ms.pointList[-1].attachLine(ssloc.number,1)
                        # connect line end B to the body
                        self.ms.bodyList[PFNum[1]].attachPoint(len(self.ms.pointList),[ssloc.rB[0]-PF[1].r[0],ssloc.rB[1]-PF[1].r[1],ssloc.rB[2]])
                    else: # end A connected to PF[1], end B connected to PF[0]
                        # connect line end A to the body
                        self.ms.bodyList[PFNum[1]].attachPoint(len(self.ms.pointList),[ssloc.rA[0]-PF[1].r[0],ssloc.rA[1]-PF[1].r[1],ssloc.rA[2]])
                        # add fairlead point 2 and attach the line to it
                        self.ms.addPoint(1,ssloc.rB)
                        self.ms.pointList[-1].attachLine(ssloc.number,1)
                        # connect line end B to the body
                        self.ms.bodyList[PFNum[0]].attachPoint(len(self.ms.pointList),[ssloc.rB[0]-PF[0].r[0],ssloc.rB[1]-PF[0].r[1],ssloc.rB[2]])

            # initialize, solve equilibrium, and plot the system 
            self.ms.initialize()
            self.ms.solveEquilibrium(DOFtype='coupled')
            
            if pristineLines:
                for i in range(0,len(self.mooringList)):
                    # add subsystems to pristine mooring objects
                    self.mooringListPristine[i].subsystem = deepcopy(self.mooringList[i].subsystem)
                    
        createArray()
        
        # Plot array if requested
        if plt:
            settings = {}
            settings["linelabels"] = True
            settings["pointlabels"] = True                          
            self.ms.plot( **settings)
            
        # return the mooring system   
        return(self.ms)

    def getFLORISArray(self, config, turblist, windSpeeds, thrustForces):
        '''
        Sets up FLORIS interface and stores thrust/windspeed curve

        Parameters
        ----------
        config : str
            Filename for FLORIS wake input file
        turblist : list of str
            List of FLORIS turbine input files (if 1, use for all turbines)
        windSpeeds : list of float
            List of wind speeds for thrust curve
        thrustForces : list of float
            List of turbine thrust forces at each wind speed in windSpeeds

        Returns
        -------
        None.

        '''
        from floris.tools import FlorisInterface
        
        # Setup FLORIS interface using base yaml file
        self.fi = FlorisInterface(config)
        
        self.fi.reinitialize(layout_x=[PF.r[0] for PF in self.platformList], layout_y=[PF.r[1] for PF in self.platformList])
        
        #right now, point to FLORIS turbine yaml. eventually should connect to ontology
        self.fi.reinitialize(turbine_type= turblist)       
        
        #store ws and thrust data... eventually store for each turbine
        self.windSpeeds = windSpeeds
        self.thrustForces = thrustForces
        
    def getFLORISMPequilibrium(self, ws, wd, ti, cutin, hubht, plotting = True):
        '''
        Function to find array equilibrium with FLORIS wake losses and MoorPy platform offsets

        Parameters
        ----------
        ws : float
            Wind speed (m/s)
        wd : float
            Wind direction (heading direction in deg where due West is 0 and due North is 90)
        ti : float
            Turbulence intenstiy (input to floris)
        cutin : float
            Cut in wind speed
        hubht : float
            Hub height
        plotting : bool
            True plots wakes and mooring systems. The default is True.

        Returns
        -------
        winds : array of float
            Initial and updated wind speed at each turbine
        xpositions : array of float
            Initial and updated x position at each turbine
        ypositions : array of float
            Initial and updated y position at each turbine
        turbine_powers : array of float
            Final power at each turbine

        '''
        # Wind directional convention ??? 
        
        from scipy import interpolate
        import time
        
        if ws < min(self.windSpeeds) or ws > max(self.windSpeeds):
            return ValueError("Wind speed outside of stored range")
        
        #FLORIS inputs the wind direction as direction wind is coming from (where the -X axis is 0)
        self.fi.reinitialize(wind_directions = [-wd+270], wind_speeds = [ws], turbulence_intensity= ti)
        
        # wind speed thrust curve for interpolation
        f = interpolate.interp1d(self.windSpeeds, self.thrustForces)
        
        #initialize list of wind speeds (1 wind speed for each turbine)
        nturbs = len(self.platformList)
        ws = [ws] * nturbs
        
        winds = []
        xpositions = []
        ypositions = []
        
        #iterate twice through updating wind speeds/platform positions
        for n in range(0, 2):
            
            # interpolate thrust force from speed/thrust curve
            for i in range(0, nturbs):
                if ws[i] < cutin:
                    T = 0
                else:
                    T = f(ws[i])
                
                # apply thrust force/moments (split into x and y components)
                self.ms.bodyList[i].f6Ext = np.array([T*np.cos(np.radians(wd)), T*np.sin(np.radians(wd)), 0, T*np.cos(np.radians(wd))*hubht, T*np.sin(np.radians(wd))*hubht, 0])       # apply an external force on the body 
                
            
            #solve statics to find updated turbine positions
            self.ms.initialize()
            self.ms.solveEquilibrium(DOFtype='both')

            #update floris turbine positions and calculate wake losses
            self.fi.reinitialize(layout_x=[body.r6[0] for body in self.ms.bodyList], layout_y=[body.r6[1] for body in self.ms.bodyList])
            self.fi.calculate_wake()
    
          
           
            #update wind speed list for RAFT
            ws = list(self.fi.turbine_average_velocities[0][0])
            winds.append(ws)
            xpositions.append([body.r6[0] for body in self.ms.bodyList])
            ypositions.append([body.r6[1] for body in self.ms.bodyList])


        #return FLORIS turbine powers (in order of turbine list)
        if min(self.fi.turbine_effective_velocities[0][0]) > cutin:
            turbine_powers = self.fi.get_turbine_powers()[0]

           
        else:
            turbine_powers = np.zeros((1,4))
        
        
        if plotting:
            
            # plot wakes
            import floris.tools.visualization as wakeviz
            horizontal_plane = self.fi.calculate_horizontal_plane(
                x_resolution=200,
                y_resolution=100,
                height=90.0,
                #yaw_angles=yaw_angles, 
            )
    
            y_plane = self.fi.calculate_y_plane(
                x_resolution=200,
                z_resolution=100,
                crossstream_dist=0.0,
                #yaw_angles=yaw_angles,
            )
            cross_plane = self.fi.calculate_cross_plane(
                y_resolution=100,
                z_resolution=100,
                downstream_dist=630.0,
                #yaw_angles=yaw_angles,
            )
    
            # Create the plots
            fig, ax = plt.subplots(1, 1, figsize=(10, 8))
            #ax_list = ax_list.flatten()
            wakeviz.visualize_cut_plane(horizontal_plane, ax=ax)
    
            cmap = plt.cm.get_cmap('viridis_r')
            self.ms.plot2d(ax = ax, Yuvec = [0,1,0])
     
            #return turbines to neutral positions **** only done if plotting - this reduces runtime for AEP calculation
            for i in range(0, nturbs):
    
                self.ms.bodyList[i].f6Ext = np.array([0, 0, 0, 0, 0, 0])       # apply an external force on the body 
                
            #solve statics to find updated turbine positions
            self.ms.initialize()
            self.ms.solveEquilibrium(DOFtype='both')
            self.ms.plot2d(ax = ax, Yuvec = [0,1,0], color = 'darkblue')
        
        
        winds = np.array(winds)
        xpositions = np.array(xpositions)
        ypositions = np.array(ypositions)
        return(winds,xpositions, ypositions, turbine_powers)  

    def calcoffsetAEP(self, windrose, ti, cutin, hubht):
        '''
        Function to calculate AEP in FLORIS with moorpy platform offsets

        Parameters
        ----------
        windrose : str
            Wind rose filename
        ti : float
            Turbulence intensity (input to FLORIS for all cases)
        cutin : float
            Tut in wind speed
        hubht : float
            Hub height

        Returns
        -------
        aeps : list of floats
            AEP for each turbine in farm

        '''

        
        # remove windrose entries below cut in or with zero frequency
        wr = np.loadtxt(windrose, delimiter=',', skiprows = 1)
        wr = wr[wr[:,0] >= cutin]
        wr = wr[wr[:,2] > 0]

        # iterate through windrose directions/speeds
        for ind in range(0, len(wr)):
            
            ws = wr[ind, 0]
            wd = wr[ind, 1]
            fq = wr[ind, 2]
            
            # solve equilibrium with moorpy and floris at given wind direction/speed
            winds,xpositions, ypositions, turbine_powers = self.getFLORISMPequilibrium(ws, wd, ti, cutin, hubht, plotting = False)
            if ind == 0:
                powers = turbine_powers
            else:
                powers = np.vstack([powers,turbine_powers])

        aeps = np.matmul(wr[:,2], powers) * 365 * 24
        return(aeps)



    # def getMarineGrowth(self,changeDepth,rho_mg=1325):
    #     ################DOES NOT CURRENTLY WORK for shared lines!!!!#################################
    #     '''

    #     Parameters
    #     ----------
    #     changeDepthSt : list
    #         [m] z elevation(s) to change from one thickness to the next (from sea floor to surface, negative below MWL)
    #     th : list
    #         [m] Marine growth thicknesses (from sea floor to surface)
    #     (optional) rho_mg : float
    #         [kg/m^3] density of marine growth (default 1325 in accordance with DNV standards)

    #     Returns
    #     -------
    #     ptCheck : list
    #         point numbers of split points in a subsystem
    #     ssCheck : list
    #         subsystem numbers associated with split points in ptCheck
    #     dCheck : list
    #         change depth values associated with the split point numbers in ptCheck
    #     '''
    #     ssCheck = []
    #     ptCheck = []
    #     dCheck = []
               
    #     for i in range(0,len(self.mooringListPristine)): # go through each line/subsystem 
            
    #         # Set up list variables for new lines
    #         LineLengths = [] # lengths of new lines
    #         LineMats = [] # line material list for new lines
    #         LType = [] # line type for new lines
            
    #         LThick = [] # thickness for each line segment
    #         bb = 0 # counter
    #         ch = [] # stores counter to see if line does not cross any cutoff depths
    #         splitNum = [] # number of the split line segment
    #         nodeD = [] # node closest to the split depth
    #         ln_raw = [] # length from point A to nodeD
    #         connList = [] # list of connectors (need to add empty connectors for between split segments)
    #         xChange = np.zeros((len(changeDepth['depth'])-1,1)) # difference between actual x-location of split and the closest node
    #         yChange = np.zeros((len(changeDepth['depth'])-1,1)) # difference between actual y-location of split and the closest node
            
    #         ch.append(bb)
    #         # set first connector
    #         connList.append(self.mooringListPristine[i].connectorList[0])
    #         for j in range(0,len(self.mooringListPristine[i].dd['sections'])): # go through each line type in the subsystem
    #             # add line material, type to list
    #             LineMats.append(self.mooringListPristine[i].dd['sections'][j]['type']['material'])
    #             LType.append(self.mooringListPristine[i].dd['sections'][j]['type']['name'])
    #             # set location for ease of coding
    #             ss = self.mooringListPristine[i].subsystem.lineList[j]
                
    #             lenseg = ss.L/ss.nNodes # length of line between nodes
    #             # get thicknesses and depths (may need to reverse order if end A is above end B)
    #             th = deepcopy(changeDepth['th'])
    #             dpt = deepcopy(changeDepth['depth'])
    #             if ss.rA[2]>ss.rB[2]: # case where rA is above rB (some shared lines could be)
    #                 low = ss.rB[2]
    #                 high = ss.rA[2]
    #                 th.reverse() # need to reverse order of changedepth
    #                 dpt.reverse()
    #                 flip = 1 # denotes rA is above rB
    #             else: # regular case
    #                 low = ss.rA[2]
    #                 high = ss.rB[2]
    #                 flip = 0
                
    #             # look up what thickness this line section starts at
    #             rAth = 0 # exit while loop when we found thickness at rA
    #             count1 = 0 # counter
    #             while rAth==0 and count1 <= len(dpt):
    #                if flip:
    #                    if high > dpt[count1]:
    #                        LThick.append(th[count1-1])
    #                        rAth = 1 
    #                else:
    #                    if low < dpt[count1]:
    #                        LThick.append(th[count1])
    #                        rAth = 1 
    #                count1 = count1 + 1

                
    #             # determine if this line section will be split
    #             for k in range(0,len(changeDepth['depth']) - 1):
    #                 if flip: # adjust depth choice
    #                     kk = k + 1 
    #                 else:
    #                     kk = k

    #                 if low < dpt[kk] and high >= dpt[kk]: # line is split by this depth
    #                     bb = bb + 1
    #                     ssCheck.append(j)
    #                     if flip:
    #                         LThick.append(th[kk])
    #                     else:
    #                         LThick.append(th[kk+1])

    #                     LType.append(ss.type['name'])
    #                     LineMats.append(ss.type['material'])
    #                     splitNum.append(j) # record the line number
    #                     # add an empty connector object to list for split location
    #                     connList.append(Connector())
    #                     # add point number to the list to check for correct depth
    #                     ptCheck.append(len(connList)-1)
    #                     if flip:
    #                         dCheck.append(len(dpt)-k-1)
    #                     else:
    #                         dCheck.append(k)
    #                     old_line = ss.getLineCoords(Time=0) # get the coordinates of the line
    #                     #find length of each new section by finding node at changeDepth
    #                     for ii in range(0, ss.nNodes-1): # go through each node in the line
    #                         if flip:
    #                             if old_line[2][ii+1]<=dpt[kk] and old_line[2][ii]>dpt[kk]:
    #                                 nodeD.append(ii) # find the node closest to the changeDepth (the node right below this depth)
    #                                 # interpolate to find x & y coordinates at chosen depth (since node might not be exactly at the right depth)
    #                                 xChange = float(np.interp(dpt[kk], old_line[2][:], old_line[0][:]))
    #                                 yChange = float(np.interp(dpt[kk], old_line[2][:], old_line[1][:]))
    #                         else:
    #                             if old_line[2][ii]<=dpt[k] and old_line[2][ii+1]>dpt[k]:
    #                                 nodeD.append(ii) # find the node closest to the changeDepth (the node right below this depth)
    #                                 # interpolate to find x & y coordinates at chosen depth (since node might not be exactly at the right depth)
    #                                 xChange = float(np.interp(dpt[k], old_line[2][:], old_line[0][:]))
    #                                 yChange = float(np.interp(dpt[k], old_line[2][:], old_line[1][:]))
                        
    #                     ln_raw.append(lenseg*nodeD[-1] + np.sqrt((xChange-old_line[0][nodeD[-1]])**2 + (yChange-old_line[1][nodeD[-1]])**2 + (dpt[kk]-old_line[2][nodeD[-1]])**2))
                        
    #                     if len(splitNum)>1 and splitNum[-1]==splitNum[-2]: # line has multiple cuts (upper cut sections have to count the length only from previous nodeD)
    #                         LineLengths.append(float(ln_raw[-1]-ln_raw[-2]))
                            
    #                     else: # first split
    #                         LineLengths.append(float(ln_raw[-1]))
                    
    #                 #print('SplitNum: ',len(splitNum),' changeDepth-1: ',len(changeDepth['depth'])-1)                            
                    
    #             if splitNum and splitNum[-1]==j: # you're on the last split - need to get the top length and number of nodes
    #                 LineLengths.append(float(ss.L-ln_raw[-1]))
                
    #             ch.append(bb)
    #             if ch[-1]==ch[-2]: # not a split line (leave length and number of nodes alone)
    #                 LineLengths.append(ss.L)
    #             # add connector at end of section to list
    #             connList.append(self.mooringListPristine[i].connectorList[j+1])
            
            
    #         # Set up list variables for pristine line info
    #         EA = []
    #         m = []
    #         d_ve_old = []
    #         cd = []
    #         cdAx = []
                                                
    #         # create arrays
    #         d_nom_old = np.zeros((len(LType),1))        
    #         ve_nom_adjust = np.zeros((len(LType),1))
    #         mu_mg = np.zeros((len(LType),1))
        
    #         nd = [] # list of dictionaries for new design dictionary sections part
            
    #         for j,ltyp in enumerate(LType):
    #             # add in information for each line type without marine growth
    #             st =  self.mooringListPristine[i].subsystem.lineTypes
    #             EA.append(st[ltyp]['EA'])
    #             m.append(st[ltyp]['m'])
    #             d_ve_old.append(st[ltyp]['d_vol'])
                
                
    #             # get ratio between ve and nom diameter from MoorProps yaml
    #             opt = helpers.loadLineProps(None)
    #             ve_nom_adjust[j] = opt[LineMats[j]]['dvol_dnom']
    #             # get cd and cdAx
    #             if not 'Cd' in st[ltyp]:
    #                 cd.append(opt[LineMats[j]]['Cd'])
    #             else:
    #                 cd.append(st[LType[j]]['Cd'])
    #             if not 'CdAx' in st[ltyp]:
    #                 cdAx.append(opt[LineMats[j]]['CdAx'])
    #             else:
    #                 cdAx.append(st[LType[j]]['CdAx'])
    #             if LineMats[j] == 'chain' or LineMats[j] == 'chain_studlink':
    #                 mu_mg[j] = 2
    #             else:
    #                 mu_mg[j] = 1
                
    #             # re-form dictionaries with marine growth values
    #             nd.append({'type':{}, 'length':{}}) # new design dictionary
    #             ndt = nd[j]['type']
                
    #             # calculate nominal diameter
    #             d_nom_old[j] = d_ve_old[j]/ve_nom_adjust[j] # m
                
    #             # calculate new line diameter that includes marine growth
    #             ndt['d_nom'] = float(d_nom_old[j]+2*LThick[j]) #m
                
    #             # calculate the new mass per meter including marine growth
    #             growthMass = np.pi/4*(ndt['d_nom']**2-d_nom_old[j]**2)*rho_mg*mu_mg[j] # marine growth mass
    #             ndt['m'] =  float(growthMass + m[j]) # kg/m (total mass)
                
    #             # calculate the submerged weight per meter including marine growth
    #             ndt['w'] = float(growthMass*(1-1025/rho_mg)*9.81 + (m[j]-np.pi/4*d_ve_old[j]**2*1025)*9.81) # N/m
                
    #             # calculate new volume-equivalent diameter (cannot use regular chain/polyester conversion because of added marine growth)
    #             ndt['d_vol'] = np.sqrt(4*((ndt['m']*9.81-ndt['w'])/1025/9.81)/np.pi)
                
    #             # calculate new increased drag coefficient from marine growth
    #             # convert cd to cd for nominal diameter, then multiply by inverse of new ve_nom_adjust (ratio of d_nom with mg to d_ve with mg) to return to cd for volume equivalent diameter
    #             ndt['Cd'] = float(cd[j]*ve_nom_adjust[j]*(ndt['d_nom']/ndt['d_vol']))
    #             ndt['CdAx'] = float(cdAx[j]*ve_nom_adjust[j]*(ndt['d_nom']/ndt['d_vol']))
                
    #             ndt['material'] = LineMats[j]
    #             ndt['name'] = str(j)
    #             ndt['MBL'] = self.mooringListPristine[i].subsystem.lineTypes[ltyp]['MBL']
    #             ndt['cost'] = self.mooringListPristine[i].subsystem.lineTypes[ltyp]['cost']
    #             ndt['EA'] = EA[j]
    #             if 'EAd' in self.mooringListPristine[i].subsystem.lineTypes[ltyp]:
    #                 ndt['EAd'] = self.mooringListPristine[i].subsystem.lineTypes[ltyp]['EAd']
                                  
    #             nd[j]['length'] = LineLengths[j]
          
    #         self.mooringList[i].dd['sections'] = nd
    #         self.mooringList[i].connectorList = connList
        
        # return(ptCheck,ssCheck,dCheck)
            
    def getMarineGrowth(self,mgDict_start,lines='all',tol=2):
        '''Calls the addMarineGrowth mooring object method for the chosen mooring objects
           and applies the specified marine growth thicknesses at the specified depth ranges
           for the specified marine growth densities.

        Parameters
        ----------
        mgDict_start : dictionary
            Provides marine growth thicknesses and the associated depth ranges
            'rho' entry is optional. If no 'rho' entry is created in the dictionary, addMarineGrowth defaults to 1325 kg/m^3
            {
                th : list with 3 values in each entry - thickness, range lower z-cutoff, range higher z-cutoff [m]
                    *In order from sea floor to surface*
                    example, if depth is 200 m: - [0.00,-200,-100]
                                                - [0.05,-100,-80]
                                                - [0.10,-80,-40]
                                                - [0.20,-40,0]
                rho : list of densities for each thickness, or one density for all thicknesses, [kg/m^3] (optional - default is 1325 kg/m^3)
                }
        lines : list, optional
            List of indices from self.mooringList to add marine growth to. Default is the string
            'all', which triggers the addition of marine growth to all mooring lines.
                
        tol : float, optional [m]
            Tolerance for marine growth cutoff depth values. Default is 2 m.
        Returns
        -------
        None.

        '''
        
        # get indices of lines to add marine growth
        if lines == 'all':
            idx = []
            for i in range(0,len(self.mooringList)):
                idx.append(i)
        else:
            idx = lines
        
        for i in idx:
            ct = 0 # counter variable
            cEq = [10,10] # make cEq mean large enough to enter the while loop for the first time.
            # get a deepcopy of the mgDict_start (may need to change the depths)
            mgDict = deepcopy(mgDict_start)
            while np.absolute(sum(cEq)/len(cEq)) > tol and ct<10: # while mean difference between actual and desired change depth is larger than tolerance and count is less than 10
                cEq = [] # reset cEq
                cD,cP = self.mooringList[i].addMarineGrowth(mgDict,project=self,idx=i)
                for j in range(0,len(cP)):
                    cEq.append(mgDict_start['th'][cD[j][0]][cD[j][1]] - self.mooringList[i].subsystem.pointList[cP[j]].r[2])
                # adjust depth to change based on difference between actual and desired change depth
                if cEq:
                    mgDict['th'][0][2] = mgDict['th'][0][2] + sum(cEq)/len(cEq)
                    for j in range(1,len(mgDict['th'])):
                        for k in range(1,3):
                            if ct < 4:
                                mgDict['th'][j][k] = mgDict['th'][j][k] + sum(cEq)/len(cEq)
                            elif ct >= 4 and ct < 9:
                                # could be ping-ponging between two different things, try adding half
                                mgDict['th'][j][k] = mgDict['th'][j][k] + 0.5*sum(cEq)/len(cEq)
                    print('average difference between expected and actual change depth is: ',sum(cEq)/len(cEq))
                else: # there were no change depths in the line (could be the case for a shared line)
                    cEq = [0,0] # kick out of the while loop
                ct = ct + 1 # add to counter
                    
                if ct == 10:
                    raise Exception(f"Unable to produce marine growth at the indicated change depths within the depth tolerance provided for mooring line index {i}. Please check for errors or increase tolerance.")
                # assign the newly created subsystem into the right place in the line list
                self.ms.lineList[i] = self.mooringList[i].subsystem
                
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
