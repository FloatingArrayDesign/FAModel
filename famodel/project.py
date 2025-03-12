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
try: 
    import raft as RAFT
except:
    pass

#from shapely.geometry import Point, Polygon, LineString
from famodel.seabed import seabed_tools as sbt
from famodel.mooring.mooring import Mooring
from famodel.platform.platform import Platform
from famodel.anchors.anchor import Anchor
from famodel.mooring.connector import Connector
from famodel.substation.substation import Substation
from famodel.cables.cable import Cable
from famodel.cables.dynamic_cable import DynamicCable
from famodel.cables.static_cable import StaticCable
from famodel.cables.cable_properties import getCableProps, getBuoyProps, loadCableProps,loadBuoyProps
from famodel.cables.components import Joint
from famodel.turbine.turbine import Turbine
from famodel.famodel_base import Node

# Import select required helper functions
from famodel.helpers import (check_headings, head_adjust, getCableDD, getDynamicCables, 
                            getMoorings, getAnchors, getFromDict, cleanDataTypes, 
                            getStaticCables, getCableDesign)


class Project():
    '''
    The overall object that defines a floating array for analysis and design 
    purposes. Its main model component is the RAFT model but it also includes
    more general information such as seabed info, metocean info, and...
    
    Ideally, this class can function even if RAFT and MoorPy are not loaded,
    at least for some high-level site processing functions.
    
    '''
    
    def __init__(self, lon=0, lat=0, file=None, depth=202,raft=1):
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
        
        # Dictionaries describing the array, divided by structure type
        self.turbineList = {}
        self.platformList = {}
        self.mooringList = {}  # A dictionary of Mooring objects
        self.anchorList  = {}
        self.cableList = {}  # CableSystem
        self.substationList = {}
        self.midConnList = {} # multi-line connectors
        
        # Dictionaries of component/product properties used in the array
        self.turbineTypes = None # list of turbine designs (RAFT input file style)
        self.platformTypes = None # list of platform designs (RAFT input file style, with rFair and zFair added)
        self.lineTypes = None
        self.anchorTypes = None
        self.cableTypes = None
        
        
        # ----- site information -----
        self.lat0  = lat  # lattitude of site reference point [deg]
        self.lon0  = lon  # longitude of site reference point [deg]
        self.g = 9.81

        # Project boundary (vertical stack of x,y coordinate pairs [m])
        self.boundary = np.zeros([0,2])
        
        # Seabed grid
        self.grid_x      = np.array([0])  # coordinates of x grid lines [m]
        self.grid_y      = np.array([0])  # coordinates of y grid lines [m]
        self.grid_depth  = np.array([[depth]])  # depth at each grid point [iy, ix]
        self.depth = depth
        
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
            self.load(file,raft=raft)
    


    def load(self, project_yaml,raft=1):
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
            if not project:
                raise Exception(f'File {file} does not exist or cannot be read. Please check filename.')
        
        # look for site section
        # call load site method
        self.loadSite(project['site'])
        
        # look for design section
        # call load design method
        self.loadDesign(project,raft=raft)
        
        
    

    # ----- Design loading/processing methods -----
    
    def loadDesign(self, d, raft=1):
        '''Load design information from a dictionary or YAML file
        (specified by input). This should be the design portion of
        the floating wind array ontology.'''
        
        print('Loading design')
        # standard function to load dict if input is yaml
        if not isinstance(d,dict):#if the input is not a dictionary, it is a yaml
            self.load(d)#load yaml into dictionary
        #d = 
        
        # ===== load FAM-specific model parts =====
        
        # array table
        if 'array' in d and d['array']['data']:
            arrayInfo = [dict(zip(d['array']['keys'], row)) for row in d['array']['data']]
        elif 'uniform_array' in d and d['uniform_array']:
            # build array info dictionary from uniform array
            ua = d['uniform_array']
            # pull out information
            WestStart = ua['west_start'] 
            NorthStart = ua['north_start']
            xSpacing = ua['spacing_x']
            ySpacing = ua['spacing_y']
            topID = ua['topsideID'] 
            pfID = ua['platformID']
            moorID = ua['mooringID']
            pfhead = ua['heading_adjust']
            
            # get locations of platforms
            arrayInfo = []
            xs = WestStart + np.arange(0, ua['n_cols']) * xSpacing
            ys = NorthStart + np.arange(0, ua['n_rows']) * ySpacing
            
            xlocs,ylocs = np.meshgrid(xs,ys)

            outx = np.hstack(xlocs)
            outy = np.hstack(ylocs)

            for i in range(ua['n_rows']*ua['n_cols']):
                arrayInfo.append({'ID':'fowt'+str(i), 'topsideID':topID, 'platformID':pfID,
                                  'mooringID':moorID, 'x_location':outx[i], 'y_location':outy[i],
                                  'heading_adjust':pfhead})
            
            
        
        # cable types
        
        # dynamic cable basic properties (details are later via MoorPy)
        
        # ----- table of cables -----
        arrayCableInfo = []
        if 'array_cables' in d and d['array_cables'] and d['array_cables']['data']:
        
            arrayCableInfo = [dict(zip( d['array_cables']['keys'], row))
                         for row in d['array_cables']['data']]
            
            
            # for ci in cableInfo:
            #     ...
                
            #     self.cables.addCable(...)
        
        # ----- cables info -----
        cableInfo = {}
        if 'cables' in d and d['cables']:
            
            cableInfo = d['cables']
        
            # for ci in d['cables']:
            #     for k, v in d['cables'].items():
            #         cableInfo[k] = v

        # ----- cable configurations -----
        dyn_cable_configs = {}
        if 'dynamic_cable_configs' in d and d['dynamic_cable_configs']:
            for k, v in d['dynamic_cable_configs'].items():
                dyn_cable_configs[k] = v
                
        # ----- cable types -----
        cable_types = {}
        if 'cable_types' in d and d['cable_types']:
            for k, v in d['cable_types'].items():
                cable_types[k] = v
                
        # ----- cable appendages -----
        cable_appendages = {}
        if 'cable_appendages' in d and d['cable_appendages']:
            for k,v in d['cable_appendages'].items():
                cable_appendages[k] = v

        
        # ----- array mooring -----
        arrayMooring = {}
        arrayAnchor = {}
        if 'array_mooring' in d and d['array_mooring']:
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
        if 'mooring_systems' in d and d['mooring_systems']:
            for k, v in d['mooring_systems'].items():
                # set up mooring systems dictionary
                mSystems[k] = v
        
        # # load in shared mooring
        # if 'array_mooring' in d:
            
        
        # ----- mooring line section types ----- 
        self.lineTypes = {}
        
        if 'mooring_line_types' in d and d['mooring_line_types']:
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
        
        if 'mooring_connector_types' in d and d['mooring_connector_types']:
            for k, v in d['mooring_connector_types'].items():
                connectorTypes[k] = v
        
        # ----- anchor types -----
        self.anchorTypes = {}
        
        if 'anchor_types' in d and d['anchor_types']:
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
                    for i in range(0, len(arrayInfo)): # loop through each entry in array
                        if m_s == arrayInfo[i]['mooringID']:
                            msys = [dict(zip(d['mooring_systems'][m_s]['keys'], row)) for row in d['mooring_systems'][m_s]['data']]
                            for i in range(0,len(msys)): #len(mSystems[m_s]['data'])): # loop through each line listed in the system
                                 if not msys[i]['MooringConfigID'] in lineConfigs: # check if they match
                                    
                                    raise Exception(f"Mooring line configuration '{msys[i]['MooringConfigID']}' listed in mooring_systems is not found in mooring_line_configs")
                            
        # ----- platforms -----
        
        RAFTDict = {} # dictionary for raft platform information
        
        if 'platform' in d and d['platform']:
            # check that there is only one platform
            if 'platforms' in d and d['platforms']:
                raise Exception("Cannot read in items for both 'platforms' and 'platform' keywords. Use either 'platform' keyword for one platform or 'platforms' keyword for a list of platforms.")
            elif type(d['platform']) is list and len(d['platform'])>1:
                raise Exception("'platform' section keyword must be changed to 'platforms' if multiple platforms are listed")
            else:
                platforms = [] # dictionary of platform information
                platforms.append(d['platform'])
                RAFTDict['platform'] = d['platform']
                self.platformTypes = d['platform']
        # load list of platform dictionaries into RAFT dictionary
        elif 'platforms' in d and d['platforms']:
            platforms = [] # list of dictionaries of platform information
            platforms.extend(d['platforms'])
            self.platformTypes = d['platforms']
            RAFTDict['platforms'] = d['platforms']
            
        # ----- turbines & topsides -----
        turbines = []
        substations = []
        if 'topsides' in d and d['topsides']:
            topsides = d['topsides']
            for ts in d['topsides']:
                if 'TURBINE' in ts['type'].upper():
                    turbines.append(ts)
                elif 'SUBSTATION' in ts['type'].upper():
                    substations.append(ts)
                    
        self.turbineTypes = turbines
        # ----- set up dictionary for each individual mooring line, create anchor, mooring, and platform classes ----
        # make platforms first if they exist, as there may be no moorings called out
        if arrayInfo:
            for i in range(len(arrayInfo)):
                # create platform instance (even if it only has shared moorings / anchors), store under name of ID for that row
                if 'z_location' in arrayInfo[i]:
                    self.platformList[arrayInfo[i]['ID']] = Platform(arrayInfo[i]['ID'],
                                                                     r=[arrayInfo[i]['x_location'],
                                                                        arrayInfo[i]['y_location'],
                                                                        arrayInfo[i]['z_location']],
                                                                     heading=arrayInfo[i]['heading_adjust'])
                elif 'z_location' in platforms[arrayInfo[i]['platformID']-1]:
                    self.platformList[arrayInfo[i]['ID']] = Platform(arrayInfo[i]['ID'],
                                                                     r=[arrayInfo[i]['x_location'],
                                                                        arrayInfo[i]['y_location'],
                                                                        platforms[arrayInfo[i]['platformID']]['z_location']],
                                                                     heading=arrayInfo[i]['heading_adjust'])
                else: # assume 0 depth
                    self.platformList[arrayInfo[i]['ID']] = Platform(arrayInfo[i]['ID'],
                                                                     r=[arrayInfo[i]['x_location'],
                                                                        arrayInfo[i]['y_location'],
                                                                        0],
                                                                     heading=arrayInfo[i]['heading_adjust'])
        # check that all necessary sections of design dictionary exist
        if arrayInfo and lineConfigs:
            
            mct = 0 # counter for number of mooring lines
            # set up a list of the alphabet for assigning names purposes
            import string
            alph = list(string.ascii_lowercase)
                               
            for i in range(0, len(arrayInfo)): # loop through each platform in array
                
                           
                # add fairlead radius and fairlead depth of this platform type from platform information section
                # get index of platform from array table
                pfID = arrayInfo[i]['platformID']-1
                if type(platforms) == list:
                    self.platformList[arrayInfo[i]['ID']].rFair = platforms[pfID]['rFair']
                    self.platformList[arrayInfo[i]['ID']].zFair = platforms[pfID]['zFair']
                else:
                    self.platformList[arrayInfo[i]['ID']].rFair = platforms['rFair']
                    self.platformList[arrayInfo[i]['ID']].zFair = platforms['zFair']
                # platform ID correlating to the ID in the platformTypes list
                self.platformList[arrayInfo[i]['ID']].dd['type'] = pfID
                # remove pre-set headings (need to append to this list so list should start off empty)
                self.platformList[arrayInfo[i]['ID']].mooring_headings = []
                entity = platforms[pfID]['type']
                # create topside instance as needed
                if arrayInfo[i]['topsideID'] > 0:
                    # create topside design dictionary
                    if isinstance(topsides,list):
                        topside_dd = topsides[arrayInfo[i]['topsideID']-1]
                    else:
                        topside_dd = topsides
                    # determine what type of topside it is
                    # create topside object and attach to platform as needed
                    if entity.upper() == 'FOWT':
                        topside_name = 'T'+str(arrayInfo[i]['topsideID'])+'_'+str(i)
                        blade_diameter = topside_dd['blade']['Rtip']*2
                        self.turbineList[topside_name] = Turbine(topside_dd,topside_name,D=blade_diameter)
                        self.turbineList[topside_name].dd['type'] = arrayInfo[i]['topsideID']-1
                        self.platformList[arrayInfo[i]['ID']].attach(self.turbineList[topside_name])
                        
                    elif entity.upper() == 'SUBSTATION':
                        topside_name = 'S'+str(arrayInfo[i]['topsideID'])+'_'+str(i)
                        self.substationList[topside_name] = Substation(topside_dd, topside_name)
                        self.platformList[arrayInfo[i]['ID']].attach(self.substationList[topside_name]) 
                    else:
                        # unknown topside - just create a basic node to attach to platform
                        topside_name = 'X'+str(arrayInfo[i]['topsideID'])+'_'+str(i)
                        node = Node(topside_name)
                        self.platformList[arrayInfo[i]['ID']].attach(node)

                # assign the topside type as the associated platform's entity
                self.platformList[arrayInfo[i]['ID']].entity = entity
                    
                
                if mSystems and not arrayInfo[i]['mooringID'] == 0: #if not fully shared mooring on this platform
                    m_s = arrayInfo[i]['mooringID'] # get mooring system ID
                    # mySys_unsorted = [dict(zip(d['mooring_systems'][m_s]['keys'], row)) for row in d['mooring_systems'][m_s]['data']]
                    # # sort the mooring lines in the mooring system by heading from 0 (North)
                    # mySys = sorted(mySys_unsorted,key=lambda x:x['heading'])
                    mySys = [dict(zip(d['mooring_systems'][m_s]['keys'], row)) for row in d['mooring_systems'][m_s]['data']]
                    # get mooring headings (need this for platform class)
                    headings = []
                    for ii in range(0,len(mySys)):
                        headings.append(np.radians(mySys[ii]['heading']))
                    # for ii in range(0,len(mSystems[m_s]['data'])):
                    #     headings.append(np.radians(mSystems[m_s]['data'][ii][1])) 
                    
                    # add mooring headings to platform class instance
                    self.platformList[arrayInfo[i]['ID']].mooring_headings = headings
                    
                    # get the mooring line information 
                    for j in range(0,len(mySys)): # loop through each line in the mooring system
                    # for j in range(0,len(mSystems[m_s]['data'])): # loop through each line in the mooring system
                        # get the configuration for that line in the mooring system
                        lineconfig = mySys[j]['MooringConfigID']
                        # lineconfig = mSystems[m_s]['data'][j][0] 
                   
                        # create mooring and connector dictionary
                        m_config = getMoorings(lineconfig, lineConfigs, connectorTypes, arrayInfo[i]['ID'], self)
                        
                        
                        # create mooring class instance as part of mooring list in the project class instance
                        mc = (Mooring(dd=m_config, id=str(arrayInfo[i]['ID'])+alph[j]))

                        # set anchor info
                        lineAnch = mySys[j]['anchorType'] # get the anchor type for the line
                        ad = getAnchors(lineAnch, arrayAnchor, self, mc=mc) # call method to create anchor dictionary
                        # add anchor class instance to anchorList in project class
                        name = str(arrayInfo[i]['ID'])+alph[j]
                        self.anchorList[name] = (Anchor(dd=ad, r=mc.rA, id=name))
                        # add mooring class instance to mooringlist in project class
                        self.mooringList[name] = mc
                        # attach mooring object to anchor and platform
                        mc.attachTo(self.anchorList[name],end='A')
                        mc.attachTo(self.platformList[arrayInfo[i]['ID']],end='B')
                        
                        # adjust mooring end positions based on platform location and mooring and platform headings, adjust anchor position
                        mc.reposition(r_center=self.platformList[arrayInfo[i]['ID']].r, 
                                      heading=headings[j]+self.platformList[arrayInfo[i]['ID']].phi, 
                                      project=self)
                        # adjust anchor z location and rA based on location of anchor
                        zAnew, nAngle = self.getDepthAtLocation(mc.rA[0], 
                                                                mc.rA[1], 
                                                                return_n=True)
                        mc.rA[2] = -zAnew
                        mc.dd['zAnchor'] = -zAnew
                        #mc.z_anch = -zAnew
                
                        
                        # update counter
                        mct += 1
                        
        
        # ----- set up dictionary for each shared mooring line or shared anchor, create mooring and anchor classes ----
        
        # create any shared mooring lines / lines connected to shared anchors / lines called out in arrayMooring table
        if arrayMooring:
            # get mooring line info for all lines 
            for j in range(0, len(arrayMooring)): # run through each line            
                
                PFNum = [] # platform ID(s) connected to the mooring line
                
                # Error check for putting an anchor (or something else) at end B
                if not any(ids['ID'] == arrayMooring[j]['endB'] for ids in arrayInfo):
                    raise Exception("Input for end B must match an ID from the array table.")
                if any(ids['ID'] == arrayMooring[j]['endB'] for ids in arrayAnchor):
                    raise Exception(f"input for end B of line_data table row '{j}' in array_mooring must be an ID for a FOWT from the array table. Any anchors should be listed as end A.")
                # Make sure no anchor IDs in arrayAnchor table are the same as IDs in array table
                for k in range(0,len(arrayInfo)):
                    if any(ids['ID'] == arrayInfo[k] for ids in arrayAnchor):
                        raise Exception(f"ID for array table row {k} must be different from any ID in anchor_data table in array_mooring section")
                # determine if end A is an anchor or a platform
                if any(ids['ID'] == arrayMooring[j]['endA'] for ids in arrayInfo): # shared mooring line (no anchor)
                    # get ID of platforms connected to line
                    PFNum.append(arrayMooring[j]['endB'])
                    PFNum.append(arrayMooring[j]['endA'])
                    # find row in array table associated with these platform IDs and set locations
                    for k in range(0, len(arrayInfo)):
                        if arrayInfo[k]['ID'] == PFNum[0]:
                            rowB = arrayInfo[k]
                            Bnum = k
                        elif arrayInfo[k]['ID'] == PFNum[1]:
                            rowA = arrayInfo[k]
                    # get headings (mooring heading combined with platform heading)
                    headingB = np.radians(arrayMooring[j]['headingB']) + self.platformList[PFNum[0]].phi
                    # get configuration for the line 
                    lineconfig = arrayMooring[j]['MooringConfigID']       
                    
                    # create mooring and connector dictionary for that line
                    m_config = getMoorings(lineconfig, lineConfigs, connectorTypes, arrayInfo[Bnum]['ID'], self)
                    # get letter number for mooring line
                    ind = len(self.platformList[PFNum[1]].getMoorings())
                    # create mooring class instance
                    mc = (Mooring(dd=m_config, id=str(PFNum[1])+'-'+str(PFNum[0])))
                    mc.shared = 1

                    # add mooring object to project mooring list           
                    self.mooringList[str(PFNum[1])+'-'+str(PFNum[0])] = mc
                    # attach mooring object to platforms
                    mc.attachTo(self.platformList[PFNum[0]],end='B')
                    mc.attachTo(self.platformList[PFNum[1]],end='A')
                    
                    # reposition both ends
                    mc.reposition(r_center=[self.platformList[PFNum[1]].r,self.platformList[PFNum[0]].r],heading=headingB,project=self)

                elif any(ids['ID'] == arrayMooring[j]['endA'] for ids in arrayAnchor): # end A is an anchor
                    # get ID of platform connected to line
                    PFNum.append(arrayMooring[j]['endB'])
                    for k in range(0,len(arrayInfo)):
                        if arrayInfo[k]['ID'] == PFNum[0]:
                            Bnum = k
                    # get configuration for that line 
                    lineconfig = arrayMooring[j]['MooringConfigID']                       
                    # create mooring and connector dictionary for that line
                    m_config = getMoorings(lineconfig, lineConfigs, connectorTypes, arrayInfo[Bnum]['ID'], self)
                    # get letter number for mooring line
                    ind = len(self.platformList[PFNum[0]].getMoorings())
                    # create mooring class instance
                    mc = (Mooring(dd=m_config, id=str(PFNum[0])+alph[ind]))
                    mc.rA = [m_config['span']+self.platformList[PFNum[0]].rFair,0,m_config['zAnchor']]
                    mc.rB = [self.platformList[PFNum[0]].rFair,0,self.platformList[PFNum[0]].zFair]

                    # check if anchor instance already exists
                    #if any(tt == 'shared_'+ arrayMooring[j]['end A'] for tt in self.anchorList): # anchor name exists already in list
                    if any(tt == arrayMooring[j]['endA'] for tt in self.anchorList): # anchor name exists already in list
                        # find anchor class instance
                        for anch in self.anchorList: #range(0,len(self.anchorList)):
                            if anch == arrayMooring[j]['endA']:
                                mc.attachTo(self.anchorList[anch],end='A')
                                mc.rA[-1] = self.anchorList[anch].r[-1]
                                mc.dd['zAnchor'] = self.anchorList[anch].r[-1]
                                mc.z_anch = self.anchorList[anch].r[-1]
                    else:
                        # find location of anchor in arrayAnchor table
                        for k in range(0,len(arrayAnchor)):
                            if arrayAnchor[k]['ID'] == arrayMooring[j]['endA']:
                                aloc = [arrayAnchor[k]['x'],arrayAnchor[k]['y']] 
                                aNum = k # get anchor row number
                                # set line anchor type and get dictionary of anchor information
                                lineAnch = arrayAnchor[k]['type']
                        ad = getAnchors(lineAnch, arrayAnchor, self, mc=mc, aNum=aNum) # call method to create dictionary
                        #mc.z_anch = -zAnew
                        # create anchor object
                        zAnew, nAngle = self.getDepthAtLocation(aloc[0], aloc[1], return_n=True)
                        self.anchorList[arrayAnchor[aNum]['ID']] = Anchor(dd=ad, r=[aloc[0],aloc[1],-zAnew], aNum=aNum,id=arrayAnchor[aNum]['ID'])
                        # attach mooring object to anchor
                        mc.attachTo(self.anchorList[(arrayAnchor[aNum]['ID'])],end='A')
                        mc.rA = [aloc[0],aloc[1],-zAnew]
                        mc.dd['zAnchor'] = -zAnew
                        mc.z_anch = -zAnew
                    # add mooring object to project mooring list
                    self.mooringList[str(PFNum[0])+alph[ind]] = mc
                    # attach mooring object to platform
                    mc.attachTo(self.platformList[PFNum[0]],end='B')
                    
                    # adjust end positions based on platform location and mooring and platform heading, adjust anchor depth
                    mc.reposition(r_center=self.platformList[PFNum[0]].r, heading=np.radians(arrayMooring[j]['headingB'])+self.platformList[PFNum[0]].phi, project=self)
                    # adjust anchor location and rA based on location of anchor
                    zAnew, nAngle = self.getDepthAtLocation(aloc[0], aloc[1], return_n=True)
                    mc.rA = [aloc[0],aloc[1],-zAnew]
                    mc.dd['zAnchor'] = -zAnew

                else: # error in input
                    raise Exception(f"end A input in array_mooring line_data table line '{j}' must be either an ID from the anchor_data table (to specify an anchor) or an ID from the array table (to specify a FOWT).")
                                           
                # add heading
                self.platformList[PFNum[0]].mooring_headings.append(np.radians(arrayMooring[j]['headingB']))
                if len(PFNum)>1: # if shared line
                    self.platformList[PFNum[1]].mooring_headings.append(np.radians(arrayMooring[j]['headingA'])) # add heading
                    
                # increment counter
                mct += 1
                
        # ===== load Cables ======
        
        # load in array cables from table (no routing, assume dynamic-static-dynamic or dynamic suspended setup)
        if arrayCableInfo:
            for i,cab in enumerate(arrayCableInfo):
                A=None
                rJTubeA = None; rJTubeB = None
                # create design dictionary for subsea cable
                dd = {'cables':[],'joints':[]}
                
                # build out specified sections
                dyn_cabA = cab['DynCableA'] if not 'NONE' in cab['DynCableA'].upper() else None
                dyn_cabB = cab['DynCableB'] if not 'NONE' in cab['DynCableB'].upper() else None
                stat_cab = cab['cableType'] if not 'NONE' in cab['cableType'].upper() else None

                if dyn_cabA:
                    dyn_cab = cab['DynCableA']
                    Acondd, jAcondd = getDynamicCables(dyn_cable_configs[dyn_cab],
                                       cable_types, cable_appendages, 
                                       self.depth, rho_water=self.rho_water, g=self.g)
                    Acondd['headingA'] = np.radians(90-cab['headingA'])
                    # only add a joint if there's a cable section after this
                    if stat_cab or dyn_cabB: 
                        dd['joints'].append(jAcondd)
                    else:
                        # this is a suspended cable - add headingB
                        Acondd['headingB'] = np.radians(90-cab['headingB'])
                    
                    rJTubeA = dyn_cable_configs[dyn_cabA]['rJTube']
                    Acondd['rJTube'] = rJTubeA
                    dd['cables'].append(Acondd)
                    # get conductor area to send in for static cable
                    A = Acondd['A']
                    
                if stat_cab:    
                    # add static cable
                    dd['cables'].append(getStaticCables(stat_cab, cable_types, 
                                                        rho_water=self.rho_water, 
                                                        g=self.g, A=A))
                    
                if dyn_cabB:
                    # assume dynamic-static-dynamic configuration
                    dyn_cabB = cab['DynCableB']
                    Bcondd, jBcondd = getDynamicCables(dyn_cable_configs[dyn_cabB],
                                                      cable_types, cable_appendages,
                                                      self.depth, rho_water=self.rho_water, 
                                                      g=self.g)
                    
                    rJTubeB = dyn_cable_configs[dyn_cabB]['rJTube']
                    Bcondd['rJTube'] = rJTubeB
                    # add heading for end A to this cable
                    Bcondd['headingB'] = np.radians(90-arrayCableInfo[i]['headingB'])
                    dd['cables'].append(Bcondd)
                    # add joint (even if empty)
                    dd['joints'].append(jBcondd)
                    

                    
                # create subsea cable object
                cableID = 'cable'+str(len(self.cableList))
                self.cableList[cableID] = Cable(cableID,d=dd)
                # attach ends
                if cab['AttachA'] in self.platformList.keys():
                    # connect to platform
                    self.cableList[cableID].attachTo(self.platformList[cab['AttachA']],end='A')
                elif cab['AttachA'] in cable_appendages:
                    pass
                else:
                    raise Exception(f'AttachA {arrayCableInfo[i]["AttachA"]} for array cable {i} does not match any platforms or appendages.')
                if cab['AttachB'] in self.platformList.keys():
                    # connect to platform
                    self.cableList[cableID].attachTo(self.platformList[cab['AttachB']],end='B')
                elif cab['AttachB'] in cable_appendages:
                    pass     
                else:
                    raise Exception(f'AttachB {arrayCableInfo[i]["AttachB"]} for array cable {i} does not match any platforms or appendages.')
                  
                # reposition the cable
                self.cableList[cableID].reposition(project=self, rad_fair=[rJTubeA,rJTubeB]) 
                        
        # create any cables from cables section (this is a more descriptive cable format that may have routing etc)           
        if cableInfo:
            
            for cab in cableInfo:
                
                rJTubeA = None; rJTubeB = None
                
                # create design dictionary for subsea cable
                dd = {'cables':[],'joints':[]}

                # pull out cable sections (some may be 'NONE')
                dyn_cabA = cab['endA']['dynamicID'] if not 'NONE' in cab['endA']['dynamicID'].upper() else None
                dyn_cabB = cab['endB']['dynamicID'] if not 'NONE' in cab['endB']['dynamicID'].upper() else None
                stat_cab = cab['type'] if not 'NONE' in cab['type'].upper() else None                       
                
                # load in end A cable section type
                if dyn_cabA:
                    Acondd, jAcondd = getDynamicCables(dyn_cable_configs[dyn_cabA],
                                                      cable_types, cable_appendages,
                                                      self.depth, rho_water=self.rho_water, 
                                                      g=self.g)
                    # only add a joint if there's a cable section after this
                    if stat_cab or dyn_cabB: 
                        dd['joints'].append(jAcondd)
                    else:
                        # this is a suspended cable - add headingB
                        Acondd['headingB'] = np.radians(90-cab['endB']['heading'])
                       
                    # add headingA
                    Acondd['headingA'] = np.radians(90-cab['endA']['heading'])
                    rJTubeA = dyn_cable_configs[dyn_cabA]['rJTube']
                    Acondd['rJTube'] = rJTubeA
                    # append to cables list
                    dd['cables'].append(Acondd)
                    
                    A = dyn_cable_configs[dyn_cabA]['A']
                    
                    
                # load in static cable design incl. routing and burial
                if stat_cab:
                    statcondd = getStaticCables(stat_cab, cable_types, 
                                                rho_water=self.rho_water, 
                                                g=self.g, A=A)
                    # load in any routing / burial
                    if 'routing_x_y_r' in cab and cab['routing_x_y_r']:
                        statcondd['routing'] = cab['routing_x_y_r']
                        
                    if 'burial' in cab and cab['burial']:
                        statcondd['burial'] = cab['burial']
                        
                    dd['cables'].append(statcondd)
                    
                # load in end B cable section type
                if dyn_cabB:
                    Bcondd, jBcondd = getDynamicCables(dyn_cable_configs[dyn_cabB],
                                                      cable_types, cable_appendages,
                                                      self.depth, rho_water=self.rho_water,
                                                      g=self.g)
                    # add headingB
                    Bcondd['headingB'] = np.radians(90-cab['endB']['heading'])
                    
                    rJTubeB = dyn_cable_configs[dyn_cabB]['rJTube']
                    Bcondd['rJTube'] = rJTubeB
                    # append to cables list
                    dd['cables'].append(Bcondd)
                    # append to joints list
                    dd['joints'].append(jBcondd)
                    
                    
                # cable ID
                cableID = cab['name'] + str(len(self.cableList))
                # create cable object
                self.cableList[cableID] = Cable(cableID,d=dd)
                
                # attach end A
                if cab['endA']['attachID'] in self.platformList.keys():
                    # connect to platform
                    self.cableList[cableID].attachTo(self.platformList[cab['endA']['attachID']],end='A')
                elif cab['endA']['attachID'] in cable_appendages:
                    pass
                else:
                    raise Exception(f"AttachA {cab['endA']['attachID']} for cable {cab['name']} does not match any platforms or appendages.") 
                # attach end B  
                if cab['endB']['attachID'] in self.platformList.keys():
                    # connect to platform
                    self.cableList[cableID].attachTo(self.platformList[cab['endB']['attachID']],end='B')
                elif cab['endB']['attachID'] in cable_appendages:
                    pass
                else:
                    raise Exception(f"AttachB {cab['endB']['attachID']} for cable {cab['name']} does not match any platforms or appendages.")
                
                # reposition the cable
                self.cableList[cableID].reposition(project=self, rad_fair=[rJTubeA,rJTubeB])       
        
        
        # ===== load RAFT model parts =====
        # load info into RAFT dictionary and create RAFT model
        if raft:
            # load turbine dictionary into RAFT dictionary
            nt = len(self.turbineTypes)
            if nt==1:            
                RAFTDict['turbine'] = self.turbineTypes[0]
            # load list of turbine dictionaries into RAFT dictionary
            else:
                RAFTDict['turbines'] = self.turbineTypes
            
            # load global RAFT settings into RAFT dictionary
            if 'RAFT_settings' in d['site'] and d['site']['RAFT_settings']:
                RAFTDict['settings'] = d['site']['RAFT_settings']
            # load RAFT cases into RAFT dictionary
            if 'RAFT_cases' in d['site'] and d['site']['RAFT_cases']:
                RAFTDict['cases'] = d['site']['RAFT_cases']
            
            # load array information into RAFT dictionary
            RAFTDict['array'] = {}
            RAFTDict['array']['keys'] = deepcopy(list(arrayInfo[0].keys())) # need to change items so make a deepcopy
            RAFTDict['array']['data'] = deepcopy([list(x.values()) for x in arrayInfo])
            # load general site info to RAFT dictionary
            RAFTDict['site'] = {'water_depth':self.depth,'rho_water':self.rho_water,'rho_air':self.rho_air,'mu_air':self.mu_air}
            RAFTDict['site']['shearExp'] = getFromDict(d['site']['general'],'shearExp',default=0.12)
            
            # create a name for the raft model
            RAFTDict['name'] = 'Project_Array'
            RAFTDict['type'] = 'input file for RAFT'
    
            # create RAFT model if necessary components exist
            if 'platforms' in RAFTDict or 'platform' in RAFTDict:

                    self.getRAFT(RAFTDict,pristine=1)
                    

        
        


    # ----- Site conditions processing functions -----

    def loadSite(self, site):
        '''Load site information from a dictionary or YAML file
        (specified by input). This should be the site portion of
        the floating wind array ontology.'''
        # standard function to load dict if input is yaml
        
        # load general information
        self.depth = getFromDict(site['general'], 'water_depth', default=self.depth)
        self.rho_water = getFromDict(site['general'], 'rho_water', default=1025.0)
        self.rho_air = getFromDict(site['general'], 'rho_air', default=1.225)
        self.mu_air = getFromDict(site['general'], 'mu_air', default=1.81e-5)
        
        # load bathymetry information, if provided
        if 'bathymetry' in site and site['bathymetry']:
            if 'file' in site['bathymetry'] and site['bathymetry']['file']: # make sure there was a file provided even if the key is there
                self.loadBathymetry(site['bathymetry']['file'])
            elif 'x' in site['bathymetry'] and 'y' in site['bathymetry']:
                self.grid_x = np.array(site['bathymetry']['x'])
                self.grid_y = np.array(site['bathymetry']['y'])
                self.grid_depth = np.array(site['bathymetry']['depths'])
                # breakpoint()
                # self.setGrid(xs,ys)
            else:
                # assume a flat bathymetry
                self.grid_depth  = np.array([[self.depth]])
                self.setGrid([0],[0])
                
                
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
                        self.boundary[i,0] = float(xy[i][0])
                        self.boundary[i,1] = float(xy[i][1])
        
        if 'seabed' in site and site['seabed']:
            # if there's a file listed in the seabed dictionary
            if 'file' in site['seabed'] and site['seabed']['file']:
                # without reading the file to tell whether it has soil property information listed, check to see if soil property information is given
                if 'soil_types' in site['seabed'] and site['seabed']['soil_types']:     # if the yaml does have soil property information
                    self.loadSoil(filename=site['seabed']['file'], yaml=site['seabed'])
                else:       # if the yaml doesn't have soil property information, read in just the filename to get all the information out of that
                    self.loadSoil(filename=site['seabed']['file'])
            # if there's no file listed in the seabed dictionary, load in just the yaml information (assuming the ['x', 'y', and 'type_array'] information are there)
            else:
                self.loadSoil(yaml=site['seabed'])

        # and set the project boundary/grid based on the loaded information
        # TBD, may not be necessary in the short term. self.setGrid(xs, ys)
        
        # load metocean portions
        
        # load marine growth info
        if 'marine_growth' in site and site['marine_growth']['data']:
            mg = site['marine_growth']['data']
            mgDict = {'th':[]}
            # adjust to match dictionary requirements for project.getMarineGrowth
            for i in range(0,len(mg)):
                if len(mg[i])>3:
                    # pull out density
                    if not 'density' in mgDict:
                        mgDict['density'] = [mg[i][3]]
                    else:
                        mgDict['density'].append(mg[i][3])
                mgDict['th'].append(mg[i][:3])
            if 'buoys' in site['marine_growth']:
                mgDict['buoy_th'] = site['marine_growth']['buoys']
                
            self.marine_growth = mgDict
                    


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
        depths = np.zeros([len(ys), len(xs)])  # note: indices are iy, ix
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
                    z[i,j] = self.getDepthAtLocation(x[j], y[i])
            
        return z



    # METHODS TO USE WITH ANCHOR TOOLS

    def loadSoil(self, filename=None, yaml=None):
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
        xs = None
        ys = None
        soil_names = None
        if filename is not None:    # if the filename option was selected, then that means there is at least a grid in the file, and maybe soil type information
            if filename[-3:]=='shp':
                raise ValueError("Geography-related operations not directly supported in Project class")
            
            elif filename[-3:]=='txt' or filename[-3:]=='dat':

                # load in the grid portion of the soil input file
                xs, ys, soil_names = sbt.readBathymetryFile(filename, dtype=str)  # read MoorDyn-style file

                soilProps = sbt.getSoilTypes(filename)     # load in the soil property information (if there is any)
            
            # regardless of whether there is soil type information in the file, if there is soil information in the yaml, read that in
            if yaml:
                soilProps = yaml['soil_types']     # if there is a yaml file as input, load in the soil props that way (overwrites the other one)


        elif filename is None:  # if the filename option was not selected
            if yaml:            # and if there was a yaml option selected, simply read in that yaml information
                xs = yaml['x']
                ys = yaml['y']
                soil_names = yaml['type_array']
                soilProps = yaml['soil_types']
            else:               # but if there was no yaml option selected (and no file option selected) -> set default values
                print('Warning: No soil grid nor soil properties were selected, but this function ran -> use preprogrammed default values')
                xs = [0]
                ys = [0]
                soil_names = ['mud']
                soilProps = dict(mud={'Su0':[2.39], 'k':[1.41], 'gamma':[10], 'depth':[0]},
                                rock={'UCS':[5], 'Em':[7], 'depth':[0]})
        
        else:
            raise ValueError("Something is wrong")
        
        '''
        # check that correct soil properties are being provided for the different soil types
        for soil in soilProps:
            if 'rock' in soil or 'hard' in soil:
                if not 'UCS' in soilProps[soil] or not 'Em' in soilProps[soil]:
                    raise ValueError('Rock soil type requires UCS and Em values')
            elif 'sand' in soil:
                if not 'phi' in soilProps[soil] or not 'gamma' in soilProps[soil]:
                    raise ValueError('Sand soil type requires phi and gamma values')
            elif 'clay' in soil:
                if not 'Su0' in soilProps[soil] or not 'k' in soilProps[soil]:
                    raise ValueError('Clay soil type requires Su0 and k values')
            elif 'mud' in soil or 'mud_soft':
                if not 'Su0' in soilProps[soil] or not 'k' in soilProps[soil]:
                    raise ValueError('Mud soil type requires Su0 and k values')
            else:
                raise ValueError(f'Soil type {soil} not recognized. Soil type key must contain one of the following keywords: rock, sand, clay, mud')
        '''
        
        # make sure the soilProps dictionary has all the required information (should be updated later with exact properties based on anchor capacity functions)
        # setting each soil type dictionary with all the values, just in case they need them for whatever reason - here are the default values
        # the default types (and values) are set if there is no other information provided
        for key,props in soilProps.items():
            props['Su0']   = getFromDict(props, 'Su0'  , shape=-1, dtype=list, default=[2.39], index=None)
            props['k']     = getFromDict(props, 'k'    , shape=-1, dtype=list, default=[1.41], index=None)
            props['alpha'] = getFromDict(props, 'alpha', shape=-1, dtype=list, default=[0.7] , index=None)
            props['gamma'] = getFromDict(props, 'gamma', shape=-1, dtype=list, default=[4.7] , index=None)
            props['phi']   = getFromDict(props, 'phi'  , shape=-1, dtype=list, default=[0.0] , index=None)
            props['UCS']   = getFromDict(props, 'UCS'  , shape=-1, dtype=list, default=[7.0] , index=None)
            props['Em']    = getFromDict(props, 'Em'   , shape=-1, dtype=list, default=[50.0], index=None)
        
        
        self.soilProps = soilProps


        
        
        if xs is not None:
            self.soil_x = np.array(xs)
            self.soil_y = np.array(ys)
            self.soil_names = np.array(soil_names)
        
        
        # update soil info for anchor if needed
        if self.anchorList:
            for anch in self.anchorList.values():
                anch.dd['soil_type'], anch.dd['soil_properties'] = self.getSoilAtLocation(anch.r[0],anch.r[1])
        
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
        
        # NEW: finds the specific soil grid point that the xy point is closest to and assigns it that soil type
        ix = np.argmin([abs(x-soil_x) for soil_x in self.soil_x])
        iy = np.argmin([abs(y-soil_y) for soil_y in self.soil_y])

        soil_name = self.soil_names[iy, ix]

        soil_info = self.soilProps[soil_name]

        return soil_name, soil_info
        
        '''
        # SIMPLE HACK FOR NOW        
        rocky, _,_,_,_ = sbt.interpFromGrid(x, y, self.soil_x, self.soil_y, self.soil_rocky)
        
        return rocky
        '''
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

    # # ----- Anchor capacity calculation functions -----


    # def calcAnchorCapacity(self, anchor):
    #     '''Compute holding capacity of a given anchor based on the soil
    #     info at its position. The anchor object's anchor properties and
    #     location will be used to determine the holding capacity, which
    #     will be saved to the anchor object.
        
    #     Parameters
    #     ----------
    #     anchor : MoorPy Anchor object (derived from Point)
    #         The anchor object in question.
    #     '''

    #     # interpolate soil properties/class based on anchor position
    #     anchor.soilProps = self.getSoilAtLocation(anchor.r[0], anchor.r[1])
        
    #     # fill in generic anchor properties if anchor info not provided
    #     if not type(anchor.anchorProps) == dict:
    #         anchor.anchorProps = dict(type='suction', diameter=6, length=12)
        
    #     # apply anchor capacity model
    #     capacity, info = anchorCapacity(anchorProps, soilProps)
        
    #     # save all information to the anchor (attributes of the Point)
    #     anchor.soilProps = soilProps
    #     anchor.anchorCapacity = capacity
    #     anchor.anchorInfo = info
        
    #     # also return it
    #     return capacity

    
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
        
        self.setBoundary(Xs, Ys)


    def setBoundary(self, Xs, Ys):
        '''Set the boundaries of the project based on x-y polygon vertices.'''
        
        # check compatibility with project grid size
        
        # save as project boundaries
        self.boundary = np.vstack([[Xs[i],Ys[i]] for i in range(len(Xs))])
        # self.boundary = np.vstack([Xs, Ys])
        
        # if the boundary doesn't repeat the first vertex at the end, add it
        if not all(self.boundary[0,:] == self.boundary[-1,:]):
            self.boundary = np.vstack([self.boundary, self.boundary[0,:]])
        
        # figure out masking to exclude grid data outside the project boundary
    
    
    def trimGrids(self, buffer=100):
        '''Trims bathymetry and soil grid information that is outside the
        project boundaries, for faster execution and plotting.'''

        # boundary extents
        xmin = np.min(self.boundary[:,0]) - buffer
        xmax = np.max(self.boundary[:,0]) + buffer
        ymin = np.min(self.boundary[:,1]) - buffer
        ymax = np.max(self.boundary[:,1]) + buffer
        
        # figure out indices to trim at if needed
        i_x1 = np.max((np.argmax(self.grid_x > xmin) - 1, 0))  # start x index
        i_y1 = np.max((np.argmax(self.grid_y > ymin) - 1, 0))  # start y index
        if xmax < np.max(self.grid_x):           
            i_x2 = np.max((np.argmin(self.grid_x < xmax) + 1, 0))  # end x index+1
        else:
            i_x2 = len(self.grid_x)
        if ymax < np.max(self.grid_y):           
            i_y2 = np.max((np.argmin(self.grid_y < ymax) + 1, 0))  # end y index+1
        else:
            i_y2 = len(self.grid_y)
        
        # trim things
        self.grid_x     = self.grid_x    [i_x1:i_x2]
        self.grid_y     = self.grid_y    [i_y1:i_y2]
        self.grid_depth = self.grid_depth[i_y1:i_y2, i_x1:i_x2]
    
    
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
                    self.grid_depth[i,j], _ = sbt.getDepthAtLocation(x, y, Xs, Ys, Zs)
                    
        # update anchor depths and the mooring ends connected to these anchors
        if self.anchorList:
            for anch in self.anchorList.values():
                anch.r[2] = -self.getDepthAtLocation(anch.r[0],anch.r[1])
                # now update any connected moorings
                for att in anch.attachments.values():
                    if isinstance(att['obj'], Mooring):
                        moor = att['obj']
                        if att['end'] == 'a' or att['end']== 'A' or att['end']== 0:
                            moor.rA[2] = -self.getDepthAtLocation(moor.rA[0],moor.rA[1])
                        else:
                            moor.rB[2] = -self.getDepthAtLocation(moor.rB[0],moor.rB[1])
        # update any joint depths
        if self.cableList:               
            for cab in self.cableList.values():
                for sub in cab.subcomponents():
                    if isinstance(sub,Joint):
                        sub.r[2] = -self.getDepthAtLocation(sub.r[0],sub.r[1])
                
        
        
        # also save in RAFT, in its MoorPy System(s)

    def cableDesignInterpolation(self,depth,cables):
        '''Interpolates between dynamic cable designs for different depths to produce 
        a design for the given depth
        '''
        # grab list of values for all cables
        cabdesign = {}
        cabdesign['span'] = [x.dd['span'] for x in cables]
        depths = [-x.z_anch for x in cables]
        cabdesign['n_buoys'] = [x.dd['buoyancy_sections']['N_modules'] for x in cables]
        cabdesign['spacings'] = [x.dd['buoyancy_sections']['spacing'] for x in cables]
        cabdesign['L_mids'] = [x.dd['buoyancy_sections']['L_mid'] for x in cables]
        cabdesign['L'] = [x.dd['L'] for x in cables]

        # sort and interp all lists by increasing depths
        sorted_indices = np.argsort(depths)
        depths_sorted = [depths[i] for i in sorted_indices]
        newdd = deepcopy(cable[0].dd)
        newdd['span'] = np.interp(depth,depths_sorted,[cabdesign['span'][i] for i in sorted_indices])
        newdd['buoyancy_sections']['N_modules'] = np.interp(depth, depths_sorted,
                                                            [cabdesign['n_buoys'][i] for i in sorted_indices])
        newdd['buoyancy_sections']['spacing'] = np.interp(depth,depths_sorted, 
                                                          [cabdesign['spacings'][i] for i in sorted_indices])
        newdd['buoyancy_sections']['L_mids'] = np.interp(depth,depths_sorted,
                                                         [cabdesign['L_mids'][i] for i in sorted_indices])
        newdd['L'] = np.interp(depth,depths,[cabdesign['L'][i] for i in sorted_indices])

        
        return(newdd)



    
    def addCablesConnections(self,connDict,cableType_def='dynamic_cable_66',oss=False,
                             substation_r=[None],ss_id=200,id_method='location',
                             keep_old_cables=False, connect_ss=True, 
                             cableConfig=None, configType=0,heading_buffer=30):
        '''Adds cables and connects them to existing platforms/substations based on info in connDict
        Designed to work with cable optimization output designed by Michael Biglu

        Parameters
        ----------
        connDict : dict
            Connection dictionary that describes the cables to create and their connections
        cableType_def : str, optional
            Cable family name corresponding to a name in the CableProps_default file. Default is 'dynamic_cable_66'
        oss : bool, optional
            Controls whether to create an offshore substation object. Default is False.
        substation_r : list, optional
            x-y location of substation object to be created. Default is None
        id_method : str, optional
            Method to identify which platforms/substations to connect to each cable. Options are 'location' or 'id'.
            'location' option connects cable to platforms/substation that has the same location as the start and end coordinates listed 
            provided for that cable in the connDict. 'id' option connects cable to platforms/substation that has the same id 
            as the global ids listed for that cable in the connDict. Default is 'location'.
        keep_old_cables : bool, optional
            Controls whether to keep any existing cables in the project object or disconnect and remove all of them.
            Default is False.
        connect_ss : bool, optional
            Controls whether to add the cable connections to the substation
        cableConfig : dict, optional
            Dictionary listing details on cable configurations to apply to the cable objects. 
        configType : int, optional
            0 = default to dynamic-static-dynamic cables, 1 = default to suspended cable systems
        heading_buffer : float, optional
            Minimum buffer between moorings and cables (degrees). Default is 30

        Returns
        -------
        None.

        '''
        
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # create substation object with id 
        if oss:
            if None in substation_r:
                dd = {'r':[5000,1000]}
            else:
                dd = {'r':substation_r}
            self.substationList[ss_id] = Substation(dd,id=ss_id)
            self.substationList[ss_id].rFair = 5 ##### TEMPORARY #####
            self.substationList[ss_id].zFair = -14
            self.substationList[ss_id].phi = 0
            
        # create reference cables for resizing
        ref_cables = []
        
        # detach and delete existing cable list unless specified to keep old cables
        if keep_old_cables:
            lcab = len(self.cableList)
        else:
            if self.cableList:
                for j,cab in enumerate(self.cableList.values()):
                    cab.detachFrom('a')
                    cab.detachFrom('b')
            self.cableList = {} 
            lcab = 0
        #### - - - -  Find cable attachments  - - - - #### 
        for i in range(0,len(connDict)): # go through each cable
            if connDict[i]['cable_id']>=100 and not connect_ss:
                # this is a substation, and we don't want to connect to it
                attachCable = 0
            else:
                attachCable = 1
                    
            if attachCable:        
                
                # attach to platforms/substations
                for pf in self.platformList.values():
                    if id_method == 'location':
                        # find platform associated with ends
                        if np.allclose(pf.r[:2],connDict[i]['coordinates'][0],atol=.01): 
                            attA = pf
                        elif np.allclose(pf.r[:2],connDict[i]['coordinates'][-1],atol=.01):
                            attB = pf
                    elif id_method == 'id':
                        # find platform associated with global id
                        if connDict[i]['turbineA_glob_id'] == pf.id:
                            attA = pf
                            # update platform location
                            pf.r[:2] = connDict[i]['coordinates'][0]
                        elif connDict[i]['turbineB_glob_id'] == pf.id:
                            attB = pf
                            # update platform location
                            pf.r[:2] = connDict[i]['coordinates'][-1]
    
                # get heading of cable from attached object coordinates 
                headingA = np.radians(90) - np.arctan2((connDict[i]['coordinates'][-1][0]-connDict[i]['coordinates'][0][0]),
                                                       (connDict[i]['coordinates'][-1][1]-connDict[i]['coordinates'][0][1]))
                headingB = np.radians(90) - np.arctan2((connDict[i]['coordinates'][0][0]-connDict[i]['coordinates'][-1][0]),
                                                       (connDict[i]['coordinates'][0][1]-connDict[i]['coordinates'][-1][1]))
    
                # figure out approx. depth at location
                initial_depths = []
                if cableConfig:
                    avgspan = np.mean([x['span'] for x in cableConfig['cableTypes'].values() if 'span' in x])
                    # depth at avg span (this is just a rough estimate!)
                    endLocA = connDict[i]['coordinates'][0] + avgspan*np.array([np.cos(headingA),np.sin(headingA)])
                    endLocB = connDict[i]['coordinates'][-1] + avgspan*np.array([np.cos(headingB),np.sin(headingB)])
                    # get depth at these locs
                    initial_depths.append(self.getDepthAtLocation(*endLocA))
                    initial_depths.append(self.getDepthAtLocation(*endLocB))
                else:
                    initial_depths = self.depth
                
                # choose cable object
                selected_cable, dd = getCableDesign(connDict[i], cableType_def, 
                                                    cableConfig, configType, 
                                                    depth=np.mean(initial_depths))
    
                # create cable object
                cab = Cable(cableType_def+str(i+lcab),d=dd)
                self.cableList[cab.id] = cab
    
                # update upstream turbines property
                cab.upstream_turb_count = getFromDict(connDict[i],'upstream_turb_count',default=0)
                
                # attach cable
                cab.attachTo(attA,end='a')
                cab.attachTo(attB,end='b')
                
                if cableConfig:
                    if 'head_offset' in selected_cable:
                        headingA += np.radians(selected_cable['head_offset'])
                        headingB -= np.radians(selected_cable['head_offset'])
                    
                    # adjust heading if too close to moorings
                    rad_buff = np.radians(heading_buffer)
                    dc0s = cab.subcomponents[0].span
                    moors = attA.getMoorings() 
                    msp = list(moors.values())[0].span + attA.rFair
                    # consider mooring headings from both ends if close enough
                    pfsp = np.sqrt((attA.r[0]-attB.r[0])**2+(attA.r[1]-attB.r[1])**2)
                    if pfsp-2*attA.rFair < msp+dc0s:
                        headingA = head_adjust([attA,attB],headingA,rad_buff=rad_buff)
                        headingB = head_adjust([attB,attA],headingB,rad_buff=rad_buff,endA_dir=-1)
                    else:
                        headingA = head_adjust([attA],headingA,rad_buff=rad_buff)
                        headingB = head_adjust([attB],headingB,rad_buff=rad_buff)
                        
                heads = [headingA,headingB]
                # reposition cable
                cab.reposition(project=self,headings=heads,rad_fair=[5,5])
                
                coords = []
                if cableConfig:
                    ref_cables = None
                    # add routing for static cable to continue along adjusted heading for total of 500m (inluding dynamic cable span) & adjust dynamic cable depths as needed
                    if len(cab.subcomponents)>1:
                        inds = [0,-1]
                        for ii,ind in enumerate(inds):
                            # adjust for depth as needed
                            dc = cab.subcomponents[ind]
                            if ref_cables:
                                dc.dd = self.cableDesignInterpolation(dc.z_anch,ref_cables) #***** left off here 1/28/25 need to check z_anch sign here and figure out how to send cables options *****
                            # add static routing
                            if cab.subcomponents[ind].span < 500:
                                spandiff = 500 - cab.subcomponents[ind].span
                                ind_of_stat = 2-4*ii # 2 for end A, -2 for end B -- relative loc of static cable compared to dynamic cable
                                stat_cable = cab.subcomponents[ind+ind_of_stat]
                                # get new coordinate routing point
                                stat_cable_end = stat_cable.rA if ind==0 else stat_cable.rB
                                coord = [stat_cable_end[0] + np.cos(heads[ii])*spandiff,
                                            stat_cable_end[1] + np.sin(heads[ii])*spandiff]
                                # append it to static cable object coordinates
                                coords.append(coord)
                
                # update lengths & spans of any static cables as needed
                cts = np.where([isinstance(a,StaticCable) for a in cab.subcomponents])[0]
                for cs in cts:
                    # update routing
                    cab.subcomponents[cs].updateRouting(coords) # also updates static and general cable lengths

                              
    
    def updatePositions(self):
        '''Temporary quick-fix to update Platform object positions based on
        RAFT or MoorPy body positions.'''
        
        for platform in self.platformList.values():
            platform.r[0] = platform.body.r6[0]
            platform.r[1] = platform.body.r6[1]
        
    
    def plot2d(self, ax=None, plot_seabed=True,plot_bathymetry=True, plot_boundary=True, bare=False, axis_equal=True,save=False,**kwargs):
        '''Plot aspects of the Project object in matplotlib in 3D.
        
        TODO - harmonize a lot of the seabed stuff with MoorPy System.plot...
        
        Parameters
        ----------
        ...
        bare : bool
            If True, supress display of extra labeling like the colorbar.
        '''
     
        # Handle extra keyword arguments or use default values
        figsize = kwargs.get('figsize', (8,8))  # the dimensions of the figure to be plotted
        edgecolor = kwargs.get('env_color',[.5,0,0,.8])
        color = kwargs.get('fenv_color',[.6,.3,.3,.6])
        alpha = kwargs.get('alpha',0.5)
        
        
        # if axes not passed in, make a new figure
        if ax == None:
            fig, ax = plt.subplots(1,1, figsize=figsize)
        else:
            fig = ax.get_figure()
        
        
        # Bathymetry 
        if plot_bathymetry:
            if len(self.grid_x) > 1 and len(self.grid_y) > 1:
                num_levels = 100  # Adjust this value as needed
                X, Y = np.meshgrid(self.grid_x, self.grid_y)
                
                contourf = ax.contourf(X, Y, self.grid_depth, num_levels, cmap='Blues', vmin=np.min(self.grid_depth), vmax=np.max(self.grid_depth))
                #contourf.norm.autoscale([0,1])
                #contourf.set_clim(0, 1000)
            
                if not bare:  # Add colorbar with label
                    cbar = plt.colorbar(contourf, ax=ax, fraction=0.04, label='Water Depth (m)')
                    
        if plot_boundary:
            ax.plot(self.boundary[:,0], self.boundary[:,1], 'b-.',label='Lease Boundary')
            
        
        # Seabed ground/soil type (to update)
        #X, Y = np.meshgrid(self.soil_x, self.soil_y)
        #ax.scatter(X, Y, c=self.soil_rocky, s=4, cmap='cividis_r', vmin=-0.5, vmax=1.5)
        # or if we have a grid of soil types, something like
        # pcolormesh([X, Y,] C, **kwargs)  wjere C is [x, y, 3 rgba]
        
        
        # Plot any object envelopes
        from shapely import Point
        for platform in self.platformList.values():
            
            for name, env in platform.envelopes.items():                 
                ax.fill(env['x'], env['y'], edgecolor=edgecolor, facecolor='none', linestyle='dashed', lw=0.8, label='Platform envelope')
        
        for mooring in self.mooringList.values():
            for name, env in mooring.envelopes.items():
                #if 'shape' in env:  # if there's a shapely object
                #    pass  # do nothing for now...
                #elif 'x' in env and 'y' in env:  # otherwise just use coordinates
                ax.fill(env['x'], env['y'], color=color,label='Mooring envelope',alpha=alpha)
        
        
        # Plot moorings one way or another (eventually might want to give Mooring a plot method)
        for mooring in self.mooringList.values():
        
            if mooring.ss:  # plot with Subsystem if available
                mooring.ss.drawLine2d(0, ax, color="k", endpoints=False, 
                                      Xuvec=[1,0,0], Yuvec=[0,1,0],label='Mooring Line')        
            else: # simple line plot
                ax.plot([mooring.rA[0], mooring.rB[0]], 
                        [mooring.rA[1], mooring.rB[1]], 'k', lw=0.5)
        
        # Plot cables one way or another (eventually might want to give Mooring a plot method)
        if self.cableList:
            maxcableSize = max([cab.dd['cables'][0].dd['A'] for cab in self.cableList.values()])
        for cable in self.cableList.values():
            # get cable color
            import matplotlib.cm as cm
            cmap = cm.get_cmap('plasma_r')
            cableSize = int(cable.dd['cables'][0].dd['A'])
            Ccable = cmap(cableSize/(1.1*maxcableSize))
            # # simple line plot for now
            # ax.plot([cable.subcomponents[0].rA[0], cable.subcomponents[-1].rB[0]], 
            #         [cable.subcomponents[0].rA[1], cable.subcomponents[-1].rB[1]],'--',color = Ccable, lw=1,label='Cable '+str(cableSize)+' mm$^{2}$')
            
            # add in routing if it exists
            for sub in cable.subcomponents:
                if isinstance(sub,StaticCable):
                    if sub.coordinates:
                        # has routing  - first plot rA to sub.coordinate[0] connection
                        ax.plot([sub.rA[0],sub.coordinates[0][0]],
                                [sub.rA[1],sub.coordinates[0][1]],':',color = Ccable,
                                lw=1,label='Buried Cable '+str(cableSize)+' mm$^{2}$')
                        # now plot route
                        if len(sub.coordinates) > 1:
                            for i in range(1,len(sub.coordinates)):
                                ax.plot([sub.coordinates[i-1][0],sub.coordinates[i][0]],
                                        [sub.coordinates[i-1][1],sub.coordinates[i][1]],
                                        ':',color=Ccable,lw=1,label='Buried Cable '+str(cableSize)+' mm$^{2}$')
                        # finally plot sub.coordinates[-1] to rB connection
                        ax.plot([sub.coordinates[-1][0],sub.rB[0]],
                                [sub.coordinates[-1][1],sub.rB[1]],':',color=Ccable,
                                lw=1,label='Buried Cable '+str(cableSize)+' mm$^{2}$')
                    else:
                        # if not routing just do simple line plot
                        ax.plot([sub.rA[0],sub.rB[0]], 
                                [sub.rA[1], sub.rB[1]],':',color = Ccable, lw=1,
                                label='Buried Cable '+str(cableSize)+' mm$^{2}$')
                elif isinstance(sub,DynamicCable):
                        ax.plot([sub.rA[0],sub.rB[0]], 
                                [sub.rA[1], sub.rB[1]],'--',color = Ccable, lw=1,
                                label='Cable '+str(cableSize)+' mm$^{2}$')
            
                # ax.plot([cable.subcomponents[0].rA[0], cable.subcomponents[-1].rB[0]], 
                #         [cable.subcomponents[0].rA[1], cable.subcomponents[0].rB[1]], 'r--', lw=0.5)
        
        # Plot platform one way or another (might want to give Platform a plot method)
        for platform in self.platformList.values():
            entity = platform.entity
            if 'FOWT' in entity.upper():
                plotstring = 'ko'
            elif 'SUBSTATION' in entity.upper():
                plotstring = 'go'
            elif 'WEC' in entity.upper():
                plotstring = 'ro'
                
            ax.plot(platform.r[0], platform.r[1], plotstring ,label=entity)

            
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        if axis_equal:
            ax.axis('equal')
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))  # Removing duplicate labels
        ax.legend(by_label.values(), by_label.keys(),loc='upper center',bbox_to_anchor=(0.5, -0.1), fancybox=True, ncol=4)
        if save:
            plt.savefig('2dfarm.png', dpi=600, bbox_inches='tight')  # Adjust the dpi as needed
            
            # TODO - add ability to plot from RAFT FOWT
        return(ax)   
        
        

    def plot3d(self, ax=None, figsize=(10,8), fowt=False, save=False,
               draw_boundary=True, boundary_on_bath=True, args_bath={}, 
               draw_axes=True, draw_bathymetry=True, draw_soil=False,
               colorbar=True, boundary_only=False):
        '''Plot aspects of the Project object in matplotlib in 3D.
        
        TODO - harmonize a lot of the seabed stuff with MoorPy System.plot...
        
        Parameters
        ----------
        ...
        
        Returns
        -------
        ax : figure axes
        '''
        
        # color map for soil plotting
        import matplotlib.cm as cm
        # from matplotlib.colors import Normalize
        
            
        

        # if axes not passed in, make a new figure
        if ax == None:    
            fig = plt.figure(figsize=figsize)
            ax = plt.axes(projection='3d')
        else:
            fig = ax.get_figure()

        if draw_bathymetry:
            if draw_soil:
                raise ValueError('The bathymetry grid and soil grid cannot yet be plotted at the same time')
            # # try increasing grid density
            if len(self.grid_x)<=1:
                pass
            else:
                if not args_bath:
                    cmap = cm.gist_earth
                    # set vmax based on bathymetry, if > 0 set max = 0
                    vmax = max([max(x) for x in -self.grid_depth])
                    if vmax > 0:
                        vmax = 0
                    args_bath = {'cmap': cmap, 'vmin':min([min(x) for x in -self.grid_depth]),
                                 'vmax': vmax}
                
                # >>> below code sections seem to do unnecessary interpolation and regridding... >>>
                if boundary_only:   # if you only want to plot the bathymetry that's underneath the boundary, rather than the whole file
                    boundary = np.vstack([self.boundary, self.boundary[0,:]])
                    xs = np.linspace(min(boundary[:,0]),max(boundary[:,0]),len(boundary[:,0]))
                    ys = np.linspace(min(boundary[:,1]),max(boundary[:,1]),len(boundary[:,1]))

                    self.setGrid(xs, ys)
            
                    # plot the bathymetry in matplotlib using a plot_surface
                    X, Y = np.meshgrid(xs, ys)  # 2D mesh of seabed grid

                    plot_depths = np.zeros([len(ys), len(xs)])
                    for i in range(len(ys)):
                        for j in range(len(xs)):
                            plot_depths[i,j], nvec = self.getDepthAtLocation(xs[j], ys[i])

                else:
                    xs = np.linspace(min(self.grid_x),max(self.grid_x),len(self.grid_x))
                    ys = np.linspace(min(self.grid_y),max(self.grid_y),len(self.grid_y))

                    self.setGrid(xs, ys)
            
                    # plot the bathymetry in matplotlib using a plot_surface
                    X, Y = np.meshgrid(self.grid_x, self.grid_y)  # 2D mesh of seabed grid
                    plot_depths = -self.grid_depth
                '''
                # interpolate soil rockyness factor onto this grid
                xs = self.grid_x
                ys = self.grid_y
                rocky = np.zeros([len(ys), len(xs)])
                for i in range(len(ys)):
                    for j in range(len(xs)):
                        rocky[i,j], _,_,_,_ = sbt.interpFromGrid(xs[j], ys[i], 
                                self.soil_x, self.soil_y, self.soil_rocky)
                                
                # or if we have a grid of soil types, something like
                ax.plot_surface(X, Y, h, rstride=1, cstride=1, facecolors = soil grid converted to colors <<<,
                            linewidth=0, antialiased=False)
                                
                                
                # apply colormap
                rc = cmap(norm(rocky))
                bath = ax.plot_surface(X, Y, -self.grid_depth, facecolors=rc, **args_bath)
                '''
                #################
                # from matplotlib import cm
                # args_bath = {'color':'#C8A2C8'}
                ####################
                bath = ax.plot_surface(X, Y, plot_depths, rstride=1, cstride=1, **args_bath)
                if colorbar:
                    fig.colorbar(bath,ax=ax,shrink=0.5,aspect=10)
        
        if draw_soil:
            if draw_bathymetry:
                raise ValueError('The bathymetry grid and soil grid cannot yet be plotted at the same time')
            X, Y = np.meshgrid(self.soil_x, self.soil_y)
            Z = np.ones([len(X), len(X[0])])*-10000
            colors = np.empty_like(self.soil_names, dtype='<U16')
            for i in range(len(self.soil_names)):
                for j in range(len(self.soil_names[0])):
                    if self.soil_names[i,j]=='mud':
                        colors[i,j] = 'g'
                    elif self.soil_names[i,j]=='hard' or self.soil_names[i,j]=='rock':
                        colors[i,j] = 'r'
                    else:
                        colors[i,j] = 'w'
            xplot = np.hstack(X)
            yplot = np.hstack(Y)
            zplot = np.hstack(Z)
            cplot = np.hstack(colors)
            ax.scatter(xplot, yplot, zplot, c=cplot, marker='o', s=50)
        
        # # also if there are rocky bits... (TEMPORARY)
        # X, Y = np.meshgrid(self.soil_x, self.soil_y)
        # ax.scatter(X, Y, c=self.soil_rocky, s=4, cmap='cividis_r', vmin=-0.5, vmax=1.5)
        
        # plot the project boundary
        if draw_boundary:
            ax.plot(self.boundary[:,0], self.boundary[:,1], 
                    np.zeros(self.boundary.shape[0]), 
                    'b--', zorder=100, lw=0.5, alpha=0.5)
            
        # plot the projection of the boundary on the seabed, if desired
        if boundary_on_bath:
            boundary_z = -self.projectAlongSeabed(self.boundary[:,0], self.boundary[:,1])
            ax.plot(self.boundary[:,0], self.boundary[:,1], boundary_z, 
                    'k--', zorder=10, lw=0.5, alpha=0.7)

        lw=1 #0.5
        # find max cable size as applicable
        if self.cableList:
            maxA = max([a.subcomponents[0].dd['A'] for a in self.cableList.values()])
        for cable in self.cableList.values():
            # get cable color
            import matplotlib.cm as cm
            cmap = cm.get_cmap('plasma_r')
            cableSize = cable.dd['cables'][0].dd['A']
            Ccable = cmap(cableSize/maxA)
            
            for j,sub in enumerate(cable.subcomponents):
                if isinstance(sub,DynamicCable):
                    if sub.ss:
                        for ii in range(len(sub.ss.lineList)):
                            sub.ss.lineList[ii].color=Ccable
                            sub.ss.lineList[ii].lw=lw
                        sub.ss.drawLine(0,ax,color=Ccable)
                        
                elif isinstance(sub,StaticCable):
                    # add static cable routing if it exists
                    if hasattr(sub,'x'):
                        # burial = sub.burial
                        # bur = []
                        # if burial and 'station' in burial:
                        #     # replace any NA with 0
                        #     for i,b in enumerate(burial['station']):
                        #         if b.upper() == 'NONE':
                        #            bur[i] = 0 
                        #     if not burial:
                        #         burial.append(0)
                        burial = []
                        burial.append(0) #### UPDATE LATER!!!
                        # get joint locations
                        jointA = cable.subcomponents[j-1]['r']
                        jointB = cable.subcomponents[j+1]['r']
                        try:
                            if len(sub.x)==1:
                                soil_z= [self.getDepthAtLocation(sub.x[0],sub.y[0])]
                            else:
                                soil_z = self.projectAlongSeabed(sub.x,sub.y)
                        except:
                            soil_z = [self.depth]
                        # plot connections from joints to first and last routing point
                        ax.plot([jointA[0],sub.x[0]],[jointA[1],sub.y[0]],[-soil_z[0],-soil_z[0]-burial[0]],
                                ':',color=Ccable,zorder=5,lw=lw)
                        ax.plot([jointB[0],sub.x[-1]],[jointB[1],sub.y[-1]],[-soil_z[-1],-soil_z[-1]-burial[-1]],
                                ':',color=Ccable,zorder=5,lw=lw)
                        
                        if len(sub.x)> 1:
                            # plot in 3d along soil_z
                            ax.plot(sub.x,sub.y,-soil_z-burial,':',color=Ccable,zorder=5,lw=lw)
                    else:
                        # no routing - just plot a straight line
                        ax.plot([sub.rA[0],sub.rB[0]],[sub.rA[1],sub.rB[1]],[sub.rA[2],sub.rB[2]],':',color=Ccable,zorder=5,lw=lw)
        
        # plot the Moorings
        ct = 0
        for mooring in self.mooringList.values():
            #mooring.subsystem.plot(ax = ax, draw_seabed=False)
            if mooring.ss:
                for line in mooring.ss.lineList:
                    if 'chain' in line.type['material']:
                        line.color = 'k'
                    elif 'polyester' in line.type['material']:
                        line.color = [.3,.5,.5]
                    else:
                        line.color = [0.5,0.5,0.5]
                    line.lw = lw
                mooring.ss.drawLine(0, ax, color='self')
                        
        
        # plot the FOWTs using a RAFT FOWT if one is passed in (TEMPORARY)
        if fowt:
            for pf in self.array.fowtList:
                pf.plot(ax,zorder=20)
            # for i in range(self.nt):
            #     xy = self.turb_coords[i,:]
            #     fowt.setPosition([xy[0], xy[1], 0,0,0,0])
            #     fowt.plot(ax, zorder=20)
        
        # Show full depth range
        if len(self.grid_depth)>1:
            ax.set_zlim([-np.max(self.grid_depth), 0])
        else:
            ax.set_zlim([-self.depth,0])

        set_axes_equal(ax)
        if not draw_axes:
            ax.axis('off')
        
        ax.view_init(20, -130)
        
        # # Remove x-axis tick labels
        # ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        # ax.set_xticks([])
        # ax.set_yticks([])
        # # # # Remove y-axis tick labels
        # # # ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
        # ax.set_zticks([])
        ax.set_axis_off()
        # ax.tick_params(axis='z', labelsize=0)
        # ax.dist -= 3
        fig.tight_layout()
        
        # ----- Save plot with an incremented number if it already exists
        if save:
            counter = 1
            output_filename = f'wind farm 3d_{counter}.png'
            while os.path.exists(output_filename):
                counter += 1
                output_filename = f'wind farm 3d_{counter}.png'
            
            # Increase the resolution when saving the plot
            plt.savefig(output_filename, dpi=600, bbox_inches='tight')  # Adjust the dpi as needed
            
        return(ax)
        
       
    def getMoorPyArray(self, plt=0, pristineLines=True, cables=True):
        '''Creates an array in moorpy from the mooring, anchor, connector, and platform objects in the array.

        Parameters
        ----------
        plt : number
            Controls whether to create a plot of the MoorPy array. 2=create plot with labels, 1=create plot, 0=no plot The default is 0.

        pristineLines : boolean, optional
            A boolean describing if the mooringList objects have been altered (0/False) or not (1/True) from their pristine, initial condition.
    
        cables : boolean, optional
            Controls whether to include cables in the moorpy array
        
        Returns
        -------
        ms : class instance
            MoorPy system for the whole array based on the subsystems in the mooringList

        '''
        print('Creating MoorPy Array')
        # create MoorPy system
        if len(self.grid_x)>1: # add bathymetry if available
            bath = {'x':self.grid_x,'y':self.grid_y,'depth':self.grid_depth}
            self.ms = mp.System(depth=self.depth,bathymetry=bath)
        else:
            self.ms = mp.System(depth=self.depth)
        # reset all anchor objects listed in self.anchorList to None in case there was a previous ms
        for i in self.anchorList:
            self.anchorList[i].mpAnchor = None
        for i in self.platformList:
            self.platformList[i].body = None
        
        wflag = 0 # warning flag has not yet been printed (prevent multiple printings of same hydrostatics warning)
        for i,body in enumerate(self.platformList.values()): # make all the bodies up front - i is index in dictionary, body is key (name of platform)
            # add a moorpy body at the correct location
            r6 = [body.r[0],body.r[1],body.r[2],0,0,0]
            # use model_details portion of dictionary to create moorpy body if given
            if 'hydrostatics' in body.dd and body.dd['hydrostatics']:
                bodyInfo = body.dd['hydrostatics']
                self.ms.addBody(-1,r6,m=bodyInfo['m'],v=bodyInfo['v'],rCG=np.array(bodyInfo['rCG']),rM=np.array(bodyInfo['rM']),AWP=bodyInfo['AWP'])
            elif not 'hydrostatics' in body.dd and wflag == 0: # default to UMaine VolturnUS-S design hydrostatics info
                print('No hydrostatics information given, so default body hydrostatics from UMaine VolturnUS-S will be used.')
                wflag = 1
                self.ms.addBody(-1,r6,m=19911423.956678286,rCG=np.array([ 1.49820657e-15,  1.49820657e-15, -2.54122031e+00]),v=19480.104108645974,rM=np.array([2.24104273e-15, 1.49402849e-15, 1.19971829e+01]),AWP=446.69520543229874)
            else:
                self.ms.addBody(-1,r6,m=19911423.956678286,rCG=np.array([ 1.49820657e-15,  1.49820657e-15, -2.54122031e+00]),v=19480.104108645974,rM=np.array([2.24104273e-15, 1.49402849e-15, 1.19971829e+01]),AWP=446.69520543229874)
            body.body = self.ms.bodyList[-1]
        # create anchor points and all mooring lines connected to the anchors (since all connected to anchors, can't be a shared mooring)
        for i in self.anchorList: # i is key (name) of anchor
            ssloc = []
            for j in self.anchorList[i].attachments: # j is key (name) of mooring object in anchor i
                # create subsystem
                if pristineLines:
                    self.anchorList[i].attachments[j]['obj'].createSubsystem(pristine=1, project=self, mooringSys=self.ms)
                    # set location of subsystem for simpler coding
                    ssloc.append(self.anchorList[i].attachments[j]['obj'].ss)
                else:
                    self.anchorList[i].attachments[j]['obj'].createSubsystem(project=self, mooringSys=self.ms)
                    # set location of subsystem for simpler coding
                    ssloc.append(self.anchorList[i].attachments[j]['obj'].ss_mod)
                self.ms.lineList.append(ssloc[-1])
                ssloc[-1].number = len(self.ms.lineList)
                # create anchor point if it doesn't already exist
                if self.anchorList[i].mpAnchor:
                    # get point number of anchor
                    num = self.anchorList[i].mpAnchor.number
                    # attach line to anchor point
                    self.ms.pointList[num-1].attachLine(ssloc[-1].number,0)
                else:
                    self.anchorList[i].makeMoorPyAnchor(self.ms)
                    # attach line to anchor point
                    self.ms.pointList[-1].attachLine(ssloc[-1].number,0)
                
                # find associated platform and attach body to point (since not a shared line, should only be one platform with this mooring object)
                for ii,k in enumerate(self.platformList): # ii is index in dictionary, k is key (name) of platform
                    if j in self.platformList[k].attachments: # j is key (name) of mooring object in anchor i checking if that same mooring object name is attached to platform k
                        PF = self.platformList[k] # platform object associated with mooring line j and anchor i
                        body = PF.body
                # attach rB point to platform 
                # add fairlead point
                self.ms.addPoint(1,ssloc[-1].rB)
                # add connector info for fairlead point
                self.ms.pointList[-1].m = self.ms.lineList[-1].pointList[-1].m 
                self.ms.pointList[-1].v = self.ms.lineList[-1].pointList[-1].v
                self.ms.pointList[-1].CdA = self.ms.lineList[-1].pointList[-1].CdA
                # attach the line to point
                self.ms.pointList[-1].attachLine(ssloc[-1].number,1)
                body.attachPoint(len(self.ms.pointList),[ssloc[-1].rB[0]-PF.r[0],ssloc[-1].rB[1]-PF.r[1],ssloc[-1].rB[2]]) # attach to fairlead (need to subtract out location of platform from point for subsystem integration to work correctly)

        
        check = np.ones((len(self.mooringList),1))
        # now create and attach any shared lines or hybrid lines attached to buoys
        for ii,i in enumerate(self.mooringList): # loop through all lines - ii is index of mooring object in dictionary, i is key (name) of mooring object
            for j in self.anchorList: # j is key (name) of anchor object
                if i in self.anchorList[j].attachments: # check if line has already been put in ms
                    check[ii] = 0     
            if check[ii] == 1: # mooring object not in any anchor lists
                # new shared line
                # create subsystem for shared line
                if hasattr(self.mooringList[i],'shared'):
                    self.mooringList[i].createSubsystem(case=self.mooringList[i].shared,pristine=pristineLines, project=self)
                else:
                    self.mooringList[i].createSubsystem(case=1,pristine=pristineLines, project=self) # we doubled all symmetric lines so any shared lines should be case 1
                # set location of subsystem for simpler coding
                if pristineLines:
                    ssloc = self.mooringList[i].ss
                else:
                    ssloc = self.mooringList[i].ss_mod
                # add subsystem as a line in moorpy system
                self.ms.lineList.append(ssloc)
                ssloc.number = len(self.ms.lineList)               
                
                # find associated platforms/ buoys
                att = self.mooringList[i].attached_to
                
                # connect line ends to the body/buoy
                ends = [ssloc.rA,ssloc.rB]
                for ki in range(0,2):
                    if isinstance(att[ki],Platform):
                        if att[ki]:
                            # add fairlead point and attach the line to it
                            self.ms.addPoint(1,ends[ki])
                            self.ms.pointList[-1].attachLine(ssloc.number,ki)
                            att[ki].body.attachPoint(len(self.ms.pointList),[ends[ki][0]-att[ki].r[0],ends[ki][1]-att[ki].r[1],ends[ki][2]])
                        else:
                            # this end is unattached
                            pass

                        
        # add in cables if desired
        if cables:

            for i in self.cableList:           
                # determine if suspended cable or not - having a static cable as a subcomponent means this is not a suspended cable
                for j,comp in enumerate(self.cableList[i].subcomponents):
                    # # check for all attachments (may be failurs enacted preventing all attachments) 
                    # don't make a subsystem for a joint - make a point                   
                    if isinstance(comp,Joint):
                        if not comp.mpConn:
                            comp.makeMoorPyConnector(self.ms)
                    elif isinstance(comp,StaticCable):
                        # don't make a subsystem for a static cable (yet...)
                        pass
                    elif isinstance(comp,DynamicCable):
                        # create subsystem for dynamic cable
                        comp.createSubsystem(pristine=pristineLines)
                        if pristineLines:                           
                            ssloc = comp.ss
                        else:
                            ssloc = comp.ss_mod
                        ssloc.number = len(self.ms.lineList)+1 
                        # add subsystem to line list
                        self.ms.lineList.append(ssloc)
                        
                        # attach each end to correct bodies
                        attach = comp.attached_to #                                 
                        
                        if j==0 and attach[0]: # only attach if cable is attached to something
                            self.ms.addPoint(1,ssloc.rA)
                            self.ms.pointList[-1].attachLine(ssloc.number,0)
                            body = comp.attached_to[0].body
                            body.attachPoint(len(self.ms.pointList),[ssloc.rA[0]-body.r6[0],ssloc.rA[1]-body.r6[1],ssloc.rA[2]])
                        elif isinstance(attach[0],Joint):
                            # connect to joint at end A
                            if not comp.attached_to[0].mpConn:
                                comp.attached_to[0].makeMoorPyConnector(self.ms)
                            joint = comp.attached_to[0].mpConn
                            joint.attachLine(ssloc.number,0)
                            
                        if j==len(self.cableList[i].subcomponents)-1 and attach[1]: # last subcomponent could be first subcomponent
                            self.ms.addPoint(1,ssloc.rB)
                            self.ms.pointList[-1].attachLine(ssloc.number,1)
                            body = comp.attached_to[-1].body
                            body.attachPoint(len(self.ms.pointList),[ssloc.rB[0]-body.r6[0],ssloc.rB[1]-body.r6[1],ssloc.rB[2]])
                        elif isinstance(attach[1],Joint):
                            # connect to joint at end B
                            if not comp.attached_to[-1].mpConn:
                                comp.attached_to[-1].makeMoorPyConnector(self.ms)
                            joint = comp.attached_to[-1].mpConn
                            joint.attachLine(ssloc.number,1)                                   

        # initialize, solve equilibrium, and plot the system 
        self.ms.initialize()
        self.ms.solveEquilibrium()       
        
        # Plot array if requested
        if plt>0:
            if plt==1:
                self.ms.plot()
            else:
                settings = {}
                settings["linelabels"] = True
                settings["pointlabels"] = True                          
                self.ms.plot( **settings)
            
        # return the mooring system   
        return(self.ms)

    def getFLORISArray(self, config, turblist, windSpeeds = [], thrustForces = []):
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
        import floris
        from floris import FlorisModel

        
        # Setup FLORIS interface using base yaml file
        self.flow = FlorisModel(config)


        x = []
        y = []
        for pf in self.platformList.values():
                if pf.entity == 'FOWT':
                    x.append(pf.r[0])
                    y.append(pf.r[1])
        self.flow.set(layout_x=x, layout_y=y)

        #right now, point to FLORIS turbine yaml. eventually should connect to ontology
        self.flow.set(turbine_type= turblist)       
        self.flow.reset_operation() # Remove any previously applied yaw angles
       # Since we are changing the turbine type, make a matching change to the reference wind height
        self.flow.assign_hub_height_to_ref_height()
        
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
        self.flow.set(wind_directions = [-wd+270], wind_speeds = [ws], turbulence_intensities= [ti])
        
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
            self.flow.set(layout_x=[body.r6[0] for body in self.ms.bodyList], layout_y=[body.r6[1] for body in self.ms.bodyList])
            self.flow.run()
    
          
           
            #update wind speed list for RAFT
            ws = list(self.flow.turbine_average_velocities[0])
            winds.append(ws)
            xpositions.append([body.r6[0] for body in self.ms.bodyList])
            ypositions.append([body.r6[1] for body in self.ms.bodyList])


        #return FLORIS turbine powers (in order of turbine list)
        if min(self.flow.turbine_average_velocities[0]) > cutin:
            turbine_powers = self.flow.get_turbine_powers()[0]

           
        else:
            turbine_powers = np.zeros((1,4))
        
        
        if plotting:
            
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
            
         
            
            #fmodel.set(wind_speeds=[wind_spd], wind_directions=[wind_dir], turbulence_intensities=[ti])
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
            
       
            #self.ms.plot2d(ax = ax, Yuvec = [0,1,0])
     
            #return turbines to neutral positions **** only done if plotting - this reduces runtime for AEP calculation
            for i in range(0, nturbs):
    
                self.ms.bodyList[i].f6Ext = np.array([0, 0, 0, 0, 0, 0])       # apply an external force on the body 
                
            #solve statics to find updated turbine positions
            # self.ms.initialize()
            # self.ms.solveEquilibrium(DOFtype='both')
            # self.ms.plot2d(ax = ax, Yuvec = [0,1,0], color = 'darkblue')
        
        
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

    
    def getRAFT(self,RAFTDict,pristine=1):
        '''Create a RAFT object and store in the project class
        
        Parameters
        ---------
        RAFTDict : dictionary
            Provides information needed to create RAFT dictionary. Need settings, cases, turbines, and platforms sub-dictionaries.
            See RAFT documentation for requirements for each sub-dictionary
        '''
        print('Creating RAFT object')
        # create RAFT model if necessary components exist
        if 'platforms' in RAFTDict or 'platform' in RAFTDict:
            # set up a dictionary with keys as the table names for each row (ease of use later)
            RAFTable = [dict(zip(RAFTDict['array']['keys'], row)) for row in RAFTDict['array']['data']]

            if 'ID' in RAFTDict['array']['keys']:
                IDindex = np.where(np.array(RAFTDict['array']['keys'])=='ID')[0][0]
                RAFTDict['array']['keys'].pop(IDindex) # remove key for ID because this doesn't exist in RAFT array table
            if 'topsideID' in RAFTDict['array']['keys']:
                ts_loc = RAFTDict['array']['keys'].index('topsideID')
                RAFTDict['array']['keys'][ts_loc] = 'turbineID'
            mooringIDindex = RAFTDict['array']['keys'].index('mooringID')
            headindex = RAFTDict['array']['keys'].index('heading_adjust')
            tsIDindex = RAFTDict['array']['keys'].index('turbineID')
            nonturbID = []
            for i in range(0,len(RAFTDict['array']['data'])):
                entity = self.platformList[RAFTDict['array']['data'][i][IDindex]].entity
                RAFTDict['array']['data'][i].pop(IDindex) # remove ID column because this doesn't exist in RAFT array data table
                if entity.upper() != 'FOWT':
                    nonturbID.append(deepcopy(RAFTDict['array']['data'][i][tsIDindex]))
                    RAFTDict['array']['data'][i][tsIDindex] = 0
                RAFTDict['array']['data'][i][mooringIDindex] = 0 # make mooringID = 0 (mooring data will come from MoorPy)
                RAFTDict['array']['data'][i][headindex] = - RAFTDict['array']['data'][i][headindex] # convert heading to cartesian from compass
            # adjust topside ids that are turbines as necessary
            for turb in RAFTDict['array']['data']:
                for ti in nonturbID:
                    if turb[tsIDindex] > ti:
                        turb[tsIDindex] -= 1

                
            # create empty mooring dictionary
            RAFTDict['mooring'] = {}
            #RAFTDict['mooring']['currentMod']
            # create raft model
            from raft.raft_model import Model
            self.array = Model(RAFTDict)
            # create dictionary of dictionaries of body hydrostatics for MoorPy bodies
            bodyInfo = {}
            for i,body in enumerate(self.array.fowtList):
                # set position (required before you can calcStatics)
                body.setPosition([RAFTable[i]['x_location'],RAFTable[i]['y_location'],0,0,0,0])
                # get body hydrostatics info for MoorPy bodies
                body.calcStatics()
                # populate dictionary of body info to send to moorpy
                if 'ID' in RAFTable[i]:
                    self.platformList[RAFTable[i]['ID']].dd['hydrostatics'] = {'m':body.m,'rCG':body.rCG,'v':body.V,'rM':body.rM,'AWP':body.AWP}
            # create moorpy array if it doesn't exist
            if not self.ms:
                if self.cableList:
                    self.getMoorPyArray(pristineLines=pristine)
                else:
                    self.getMoorPyArray(pristineLines=pristine)
            # assign moorpy array to RAFT object
            self.array.ms = self.ms
            # connect RAFT fowt to the correct moorpy body
            for i in range(0,len(self.platformList)): # do not include substations (these are made last)
                self.array.fowtList[i].body = self.ms.bodyList[i]
        else:
            raise Exception('Platform(s) must be specified in YAML file')
            
    def getMarineGrowth(self,mgDict_start=None,lines='all',tol=2,display=False):
        '''Calls the addMarineGrowth mooring object method for the chosen mooring objects
           and applies the specified marine growth thicknesses at the specified depth ranges
           for the specified marine growth densities.

        Parameters
        ----------
        mgDict_start : dictionary, optional
            Provides marine growth thicknesses and the associated depth ranges
            Default is None, which triggers the use of the marine growth dictionary saved in the project class
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
            List of keys from self.mooringList or self.cableList to add marine growth to. For cables, each entry must be a list
            of length 2 with the cable id in the first position and the subcomponent number of the dynamic cable in the 2nd position.
            Default is the string 'all', which triggers the addition of marine growth to all mooring lines.
                Ex: lines = ['fowt1a',['suspended_cable10',0]] to add mg just to mooring line 'fowt1a' 
                and the first subcomponent of cable 'suspended_cable10'
                
        tol : float, optional [m]
            Tolerance for marine growth cutoff depth values. Default is 2 m.
        display : bool, optional
            Whether or not to print out difference between expected and actual depth of change 
            for marine growth thicknesses
        
        Returns
        -------
        None.

        '''
        # make sure moorpy system exists already
        if not self.ms:
            if self.cableList:
                self.getMoorPyArray(cables=1)
            else:
                self.getMoorPyArray()
        # get marine growth dictionary
        if not mgDict_start:
            mgDict_start = self.marine_growth
        
        # get indices of lines to add marine growth
        if lines == 'all':
            idx = []
            for i in self.mooringList:
                idx.append(i)
            for i in self.cableList:
                for j,cab in enumerate(self.cableList[i].subcomponents):
                    if isinstance(cab,DynamicCable):
                        idx.append([i,j])
        else:
            idx = lines
        
        for ii,i in enumerate(idx):
            ct = 0 # counter variable
            cEq = [10,10] # make cEq mean large enough to enter the while loop for the first time.
            # get a deepcopy of the mgDict_start (may need to change the depths)
            mgDict = deepcopy(mgDict_start)
            while np.absolute(sum(cEq)/len(cEq)) > tol and ct<10: # while mean difference between actual and desired change depth is larger than tolerance and count is less than 10
                cEq = [] # reset cEq
                if isinstance(i,list):
                    # this is a cable
                    cD,cP = self.cableList[i[0]].subcomponents[i[1]].addMarineGrowth(mgDict)
                    for j in range(0,len(cP)):
                        cEq.append(mgDict_start['th'][cD[j][0]][cD[j][1]] - self.cableList[i[0]].subcomponents[i[1]].ss_mod.pointList[cP[j]].r[2])
                else:
                    cD,cP = self.mooringList[i].addMarineGrowth(mgDict,project=self,idx=i)
                    for j in range(0,len(cP)):
                        cEq.append(mgDict_start['th'][cD[j][0]][cD[j][1]] - self.mooringList[i].ss_mod.pointList[cP[j]].r[2])
                # adjust depth to change based on difference between actual and desired change depth
                if cEq:
                    mcEq = sum(cEq)/len(cEq)
                    mgDict['th'][0][2] = mgDict['th'][0][2] + mcEq
                    for j in range(1,len(mgDict['th'])):
                        for k in range(1,3):
                            if ct < 4 and abs(mcEq)<12:
                                mgDict['th'][j][k] = mgDict['th'][j][k] + mcEq
                            elif (ct >= 4 and ct < 9) or abs(mcEq)>=12:
                                # could be ping-ponging between two different things, try adding half
                                mgDict['th'][j][k] = mgDict['th'][j][k] + 0.5*mcEq
                    if display:
                        print('average difference between expected and actual change depth is: ',mcEq)
                else: # there were no change depths in the line (could be the case for a shared line)
                    cEq = [0,0] # kick out of the while loop
                ct = ct + 1 # add to counter
                    
                if ct == 10:
                    raise Exception(f"Unable to produce marine growth at the indicated change depths within the depth tolerance provided for mooring line index {i}. Please check for errors or increase tolerance.")
                # assign the newly created subsystem into the right place in the line list
                if isinstance(i,list):
                    self.ms.lineList[ii] = self.cableList[i[0]].subcomponents[i[1]].ss_mod
                else:
                    self.ms.lineList[ii] = self.mooringList[i].ss_mod
                
    
    def getCorrosion(self,corr_th=10,lines='all'):
        '''
        Function to reduce MBL of specified lines based on corrosion thickness

        Parameters
        ----------
        corr_th : float, optional
            Thickness of corrosion in mm. The default is 10.
        lines : list or string, optional
            List of line ids to add corrosion to. The default is 'all', which adds to every line.

        Returns
        -------
        None.

        '''
        if lines == 'all':
            idx = []
            for i in self.mooringList:
                idx.append(i)
        elif isinstance(lines,list):
            idx = lines
        
        for ii,i in enumerate(idx):
            self.mooringList[i].addCorrosion(corrosion_mm=corr_th)
    
    
    def updateUniformArray(self,n_rows,pf_rows,spacingXY,grid_rotang=0,grid_skew_x=0,grid_skew_y=0,grid_trans_x=0,grid_trans_y=0,phis=[0,0],center=[0,0]):
        '''
        Function to update the array spacing, angles, platform locations and headings, etc for a uniform array
        
        Parameters
        ----------
        n_rows : int
            Number of rows
        pf_rows : int or list
            Number of platforms in each row. If an int, the number of platforms is constant (essentially equivalent to number of columns). 
            If a list, each entry represents the number of platforms for an individual row
        spacingXY : list
            2 entry list that provides the spacing in the x and y direction respectively [m]. If an integer is provided 
            that is less than 20, it is assumed the value given is in multiples of the turbine diameter instead of [m].
            If the x and y entries are lists themselves, this specifies the spacing between each row/col individually
        grid_rotang : float
            Rotation of the grid [deg] from 0 North, optional. Default is 0
        grid_skew_x : float
            Angle of the parallelogram between adjacent rows [deg], optional. Default is 0
        grid_skew_y : float
            Angle of the parallelogram between adjacent columns [deg], optional. Default is 0
        grid_trans_x : float
            x offset [m] from Western-most boundary point (if no boundary in project, x offset from [0,0] origin) for all turbine positions.
        grid_trans_y : float
            y offset [m] from Southern-most boundary point (if no boundary in project, y offset from [0,0] origin) for all turbine positions.
        phis : list
            Platform headings for rows. The length of the list dictates the repeat pattern 
            (a length of 2 indicates that odd rows will have phis[0] heading and even rows will have phis[1] heading)
            If each list entry is also a list, this sublist indicates the repeat pattern of phi within the row from west to east
        center : list
            x and y coordinate of array layout center, around which any rotations/skews would occur. Optional, default is [0,0]
        
        Returns
        -------
        None.

        '''
        r = np.zeros([len(self.platformList),3])   # xy location of each platform
        r_no_off = np.zeros([len(self.platformList),3]) # xy location of each platform not considering offsets
        phi = np.zeros(len(self.platformList)) # rotation heading of each platform
        D = 0
        for i,turb in enumerate(self.turbineList.values()):
            D += turb.D 
        D = D/i/2
        # get a list of platform names for use later
        pfIDs = [None]*len(self.platformList)
        
        for i,pf in enumerate(self.platformList.values()):
            if pf.rc:
                # convert rc to number (as if array were hstacked)
                numh = pf.rc[0]*n_rows + pf.rc[1]
                pfIDs[numh] = pf.id              
            else:    
                # if not row-col indices, assume dictionary order is in rows starting from north west corner
                pfIDs[i] = pf.id

        # calculate grid rotation and skew (from fadesign layout_helpers create_rotated_layout function)
        # Shear transformation in X
        # Calculate trigonometric values
        cos_theta = np.cos(np.radians(-grid_rotang))
        sin_theta = np.sin(np.radians(-grid_rotang))
        tan_phi_x = np.tan(np.radians(grid_skew_x))
        tan_phi_y = np.tan(np.radians(grid_skew_y))

        # Compute combined rotation and skew transformation matrix
        transformation_matrix = np.array([[cos_theta - sin_theta*tan_phi_y, cos_theta*tan_phi_x - sin_theta],
                                          [sin_theta + cos_theta*tan_phi_y, sin_theta*tan_phi_x + cos_theta]])
        
        west_bound = None
        north_bound = None
        # determine if a boundary exists for the project
        if self.boundary.size > 0:
            # find western-most boundary point
            west_bound = self.boundary[0][0]
            north_bound = self.boundary[0][1]
            for xy in self.boundary:
                if xy[0] < west_bound:
                    west_bound = xy[0]
                if xy[1] > north_bound:
                    north_bound = xy[1]
        
        # set up grid
        ct = 0 # counter variable
        p = 0  # counter variable for phi rows
        xx = 0 # store previous row location
        yy = 0 # store previous column location
        for i in range(0,n_rows):
            # find number of platforms in the row
            if isinstance(pf_rows,list):
                num = pf_rows[i]
            else:
                num = pf_rows     
            
            for j in range(0,num):
                # find spacing from previous rows + cols
                # if first row or first col, spacing = 0
                if j == 0:
                    spacing_x = 0 
                else:
                    if isinstance(spacingXY[0],list):
                        spacing_x = spacingXY[0][j-1]
                    else:
                        spacing_x = spacingXY[0]
                    if spacing_x < 20:
                        spacing_x = spacing_x*D
                if i == 0:
                    spacing_y = 0 
                else:
                    if isinstance(spacingXY[1],list):
                        spacing_y = spacingXY[1][i-1]
                    else:
                        spacing_y = spacingXY[1]                    
                    if spacing_y < 20:
                        spacing_y = spacing_y*D
                
                # assign location of platform in the row (temporary, does not include rotation or skew)
                r_no_off[ct][0] = spacing_x + xx
                r_no_off[ct][1] = -spacing_y + yy
                
                xx = r_no_off[ct][0]
                                               
                if west_bound:
                    # set offsets from north and west boundaries
                    xoff = west_bound + grid_trans_x
                    yoff = north_bound - grid_trans_y
                else:
                    # use grid_trans_x and grid_trans_y as offsets
                    xoff = grid_trans_x
                    yoff = grid_trans_y
                    
                r[ct][0] = r_no_off[ct][0] + xoff + center[0]
                r[ct][1] = r_no_off[ct][1] + yoff + center[1]
                
                # include grid rotation around center
                x,y = np.dot(transformation_matrix,[r[ct][0]-center[0],r[ct][1]-center[1]])
                
                # get platform heading
                if p == len(phis):
                    p = 0

                if isinstance(phis[p],list):
                    phi[ct] = phis[p][j]
                else:
                    phi[ct] = phis[p]
                    
                # call platform function to set location, phi, moorings, cables, and anchors
                self.platformList[pfIDs[ct]].setPosition([x,y],heading=phi[ct],degrees=True)
                self.platformList[pfIDs[ct]].rc = [i,j]
                # # update cable lengths as needed (assumes the attachments are correct)
                # cabs = self.platformList[pfIDs[ct]].getCables()
                # for cab in cabs.values():
                #     # determine what the cable is connected to
                #     rr = [0,0]
                #     for k,att in enumerate(cab.attached_to):
                #         if isinstance(att,Platform):
                #             rr[k] = att.r
                #     newSpan = np.sqrt((rr[0][0]-rr[1][0])**2 + (rr[0][1]-rr[1][1])**2)
                #     cab.updateSpan(newSpan)
                #     # update cable end locations
                #     cab.reposition()
                # update anchor locations
                anchs = self.platformList[pfIDs[ct]].getAnchors()
                for anch in anchs.values():
                    for att in anch.attachments:
                        if isinstance(anch.attachments[att]['obj'],Mooring):
                            anch.r = anch.attachments[att]['obj'].rA
                            # update depth of anchor as needed
                            aDepth = self.getDepthAtLocation(anch.r[0],anch.r[1])
                            anch.attachments[att]['obj'].rA[2] = -aDepth                                                                       
                
                # increase counter
                ct += 1
           
            xx = 0
            yy = r_no_off[ct-1][1]
            # increase phi counter
            p += 1
                
        # update moorpy
        self.getMoorPyArray(plt=1,cables=1)
        
    def duplicate(self,pf, r=None,heading=None):
        '''
        Function to duplicate a platform object and all
        of its associated moorings and anchors (NOT CABLES)

        Parameters
        ----------
        pf : Platform object
            Platform object to duplicate
        r : list, optional
            Location of new platform object. Default is None,
            in which case platform is not moved

        Returns
        -------
        pf2 : Platform object
            Duplicated object

        '''
        from copy import deepcopy
        # get name for new platform based on length of platformList
        if self.platformList:
            lp = len(self.platformList)
        else:
            lp = 0
        newid = 'fowt'+str(lp)
        import string
        alph = list(string.ascii_lowercase)
        
        # copy platform object and its attachments and disconnect from attachments
        pf2 = deepcopy(pf)
        pf2.id = newid
        self.platformList[newid] = pf2
        count = 0 
        
        for att in pf.attachments.values():
            if isinstance(att['obj'],Mooring):
                if att['end'] == 'a':
                    endB = 0 
                else:
                    endB = 1
                # grab all info from mooring object
                md = deepcopy(att['obj'].dd)
                mhead = att['obj'].heading
                # detach mooring object from platform
                pf2.detach(att['obj'],end=endB)
                # create new mooring object
                newm = Mooring(dd=md,id=newid+alph[count])
                self.mooringList[newm.id] = newm
                newm.heading = mhead
                # attach to platform
                pf2.attach(newm,end=endB)
                # grab info from anchor object and create new one
                ad = deepcopy(att['obj'].attached_to[1-endB].dd)
                newa = Anchor(dd=ad,id=newid+alph[count])
                self.anchorList[newa.id] = newa
                # attach anchor to mooring
                newm.attachTo(newa,end=1-endB)
                newm.reposition(r_center=r,project=self)
                zAnew, nAngle = self.getDepthAtLocation(newm.rA[0], newm.rA[1], return_n=True)
                newm.rA[2] = -zAnew
                newm.dd['zAnchor'] = -zAnew
                newa.r = newm.rA
                
                count += 1
                
            elif isinstance(att['obj'],Turbine):
                pf2.detach(att['obj'])
                turb = deepcopy(att['obj'])
                turb.id = newid+'turb'
                self.turbineList[turb.id] = turb
                pf2.attach(turb)
                
            else:
                # could be cable, just detach for now
                pf2.detach(att['obj'],att['end'])
        
        # reposition platform as needed
        pf2.setPosition(r,heading=heading,project=self)
 
        
        
        # delete body object from pf2
        pf2.body = None
        
        return(pf2)
            
        
    def addPlatformMS(self,ms,r=[0,0]):
        '''
        Create a platform object, along with associated mooring and anchor objects
        from a moorpy system
        Currently only works for regular (non-shared) moorings.

        Parameters
        ----------
        ms : moorpy system, optional
            Moorpy system representing the platform, its moorings, and its anchors. The default is None.
        
        Returns
        -------
        new platform object

        '''
        # create platform, moorings, and anchors from ms
        ix = len(self.platformList)
        # check there is just one body
        if len(ms.bodyList) > 1:
            raise Exception('This function only works with a 1 body system')
            
        # switch to subsystems if lineList doesn't already have them
        if not isinstance(ms.lineList[0],mp.subsystem.Subsystem):
            from moorpy.helpers import lines2ss
            lines2ss(ms)
            
        # get lines attached to platform and headings
        
        mhead = []
        mList = []
        endB = []
        count = 0
        pfid = 'fowt'+str(ix)
        import string
        alph = list(string.ascii_lowercase)
        for point in ms.bodyList[0].attachedP:
            for j,line in enumerate(ms.pointList[point-1].attached):
                md = {'sections':[],'connectors':[]} # start set up of mooring design dictionary
                rA = ms.lineList[line-1].rA
                rB = ms.lineList[line-1].rB
                pfloc = ms.bodyList[0].r6
                if ms.pointList[point-1].attachedEndB[j]:
                    vals = rA[0:2]-rB[0:2]
                    zFair = rB[2]
                    rFair = np.hypot(rB[0]-pfloc[0],rB[1]-pfloc[1])
                    endB.append(1)
                else:
                    vals = rB[0:2]-rA[0:2]
                    zFair = rA[2]
                    rFair = np.hypot(rA[0]-pfloc[0],rA[1]-pfloc[1])
                    endB.append(0)
                    
                # pull out mooring line info
                md['rad_fair'] = rFair
                md['z_fair'] = zFair
                md['span'] = np.hypot(vals[0],vals[1])
                if not endB[-1]:
                    md['zAnchor'] = -self.getDepthAtLocation(rA[0],rA[1])
                else:
                    md['zAnchor'] = -self.getDepthAtLocation(rB[0],rB[1])
                    
                # # add section and connector info
                # md['sections'].append({'type':line.type})
                # md['sections'][-1]['L'] = line.L
                # md['connectors'].append({'m':point.m,'v':point.v,'Ca':point.Ca,'CdA':point.CdA})
                  
                # anline = True
                # for pt in ms.pointList:
                #     if line in pt.attached and pt != point:
                #         n_att = len(pt.attached)
                #         nextloc = np.where([x!=line for x in pt.attached])[0][0]
                #         if n_att == 1:
                #             # this is the anchor point
                #             ad = {'design':{}}
                #             ad['design']['m'] = pt.m
                #             ad['design']['v'] = pt.v
                #             ad['design']['CdA'] = pt.CdA
                #             ad['design']['Ca'] = pt.Ca
                #             if 'anchor_type' in pt.entity:
                #                 ad['type'] = pt.entity['anchor_type']
                #             self.anchorList[mList[-1].id] = Anchor(dd=ad,r=pt.r,id=mList[-1].id)
                #             self.anchorList[mList[-1].id].attach(mList[-1],end=1-endB[-1])
                #             # reposition mooring and anchor
                #             mList[-1].reposition(r_center=r)
                #             zAnew = self.getDepthAtLocation(mList[-1].rA[0], 
                #                                             mList[-1].rA[1])
                #             mList[-1].rA[2] = -zAnew
                #             mList[-1].dd['zAnchor'] = -zAnew
                #             self.anchorList[mList[-1].id].r = mList[-1].rA
                #             anline = False
                #         else:
                #             # add section and connector info
                #             md['sections'].append({'type':sline.type})
                #             md['sections'][-1]['L'] = sline.L
                #             spt = ms.lineList[line-1].pointList[k]
                #             md['connectors'].append({'m':spt.m,'v':spt.v,'Ca':spt.Ca,'CdA':spt.CdA})
                            
                            
                
                for k,sline in enumerate(ms.lineList[line-1].lineList):
                    # add section and connector info
                    md['sections'].append({'type':sline.type})
                    md['sections'][-1]['L'] = sline.L
                    spt = ms.lineList[line-1].pointList[k]
                    md['connectors'].append({'m':spt.m,'v':spt.v,'Ca':spt.Ca,'CdA':spt.CdA})
                spt = ms.lineList[line-1].pointList[k+1]
                md['connectors'].append({'m':spt.m,'v':spt.v,'Ca':spt.Ca,'CdA':spt.CdA})
                mhead.append(90 - np.degrees(np.arctan2(vals[1],vals[0])))
                mList.append(Mooring(dd=md,id=pfid+alph[count]))
                mList[-1].heading = mhead[-1]
                self.mooringList[mList[-1].id] = mList[-1]
                
                # pull out anchor info
                for pt in ms.pointList:
                    if line in pt.attached:
                        loc = np.where([x==line for x in pt.attached])[0][0]
                        if pt.attachedEndB[loc]:
                            ad = {'design':{}}
                            ad['design']['m'] = pt.m
                            ad['design']['v'] = pt.v
                            ad['design']['CdA'] = pt.CdA
                            ad['design']['Ca'] = pt.Ca
                            if 'anchor_type' in pt.entity:
                                ad['type'] = pt.entity['anchor_type']
                            self.anchorList[mList[-1].id] = Anchor(dd=ad,r=pt.r,id=mList[-1].id)
                            self.anchorList[mList[-1].id].attach(mList[-1],end=1-endB[-1])
                            # reposition mooring and anchor
                            mList[-1].reposition(r_center=r,project=self)
                            zAnew = self.getDepthAtLocation(mList[-1].rA[0], 
                                                            mList[-1].rA[1])
                            mList[-1].rA[2] = -zAnew
                            mList[-1].dd['zAnchor'] = -zAnew
                            self.anchorList[mList[-1].id].r = mList[-1].rA
                        
                count += 1
                            
                
        # add platform at ms.body location       
        self.platformList[pfid] = Platform(pfid, r=ms.bodyList[0].r6[0:2],
                                                     mooring_headings=mhead,
                                                     rFair=rFair, zFair=zFair)
        # attach moorings
        for i,moor in enumerate(mList):
            self.platformList[pfid].attach(moor,end=endB[i])
            
        # update location
        self.platformList[pfid].setPosition(r=r)
            
        return(self.platformList[pfid])
            
    def addPlatformConfig(self,configDict,r=[0,0]):
        '''
        Add a platform, anchors, and mooring lines based on a configuration dictionary
        
        Parameters
        ----------
        configDict : dict
            Dictionary describing design of platform, moorings, and anchors
            The dictionary layout is as follows:
                platform:
                    rFair: 
                    zFair:
                    moor_headings:
                    platform_heading:
                mooring_config:
                    span:
                    segment list: # in order from end A to end B
                        - <segment name>:
                            material: # polyester,chain,etc
                            length:
                            d_nom: # can list material properties below this, or it will be imported from MoorProps
                                    # based on the material type and nominal diameter
                        - <c
                        - <segment name>:
                            material: 
                            length:
                            d_nom:
                anchor_config:
                    geometry:
                    type:
        r : list, optional
            x,y location of platform. The default is [0,0]
                        
        Returns
        -------
        pf : new platform object
                    
        '''
        
        if self.platformList:
            lp = len(self.platformList)
        else:
            lp = 0
            
        pfid = 'fowt'+str(lp)
        import string
        alph = list(string.ascii_lowercase)
            
        # pull out platform info
        pfinfo = configDict['platform']
        
        # create platform object
        self.platformList[pfid] = Platform(pfid, r=r,
                                           mooring_headings=pfinfo['mooring_headings'],
                                           rFair=pfinfo['rFair'], zFair=pfinfo['zFair'],
                                           phi=pfinfo['platform_heading'])
        
        # pull out mooring info
        minfo = configDict['mooring_config']
        # create mooring objects
        for i in range(len(pfinfo['mooring_headings'])):
            head = pfinfo['mooring_headings'][i]+pfinfo['platform_heading']
            md = {'span':minfo['span'],'sections':[],'connectors':[]}
      
    def arrayWatchCircle(self,plot=False, ang_spacing=45, RNAheight=150,
                         shapes=True,thrust=1.95e6,SFs=True,moor_envelopes=True,
                         DAF=1):
        '''
        Method to get watch circles on all platforms at once

        Parameters
        ----------
        plot : bool, optional
            Controls whether to plot watch circles at the end. The default is False.
        ang_spacing : int/float, optional
            Spacing between angles to check the platform offsets at. The default is 45.
        RNAheight : int/float, optional
            Height of RNA from mean sea level. The default is 150.
        shapes : bool, optional
            Controls whether to create shapely objects. The default is True.
        thrust : float, optional
            Thrust force on turbine. The default is 1.95e6 (thrust on IEA 15 MW reference
                                                            turbine at rated wind speed.
        SFs : bool, optional
            Controls whether to output safety factors of moorings and cables,
            and calculate anchor loads. The default is True.
        moor_envelopes : bool, optional
            Controls whether to create and store motion envelopes for moorings
        DAF : float, optional
            Dynamic amplification factor. Multiplies max forces by the DAF to account for 
            under-estimation of loads when using static method

        Returns
        -------
        x : np.array
            matrix of platform x-locations for watch circles
        y : np.array
            matrix of platform y-locations for watch circles
        maxVals : dict
            dictionary of safety factors for mooring line tensions for each turbine

        '''
            
        # get angles to iterate over
        angs = np.arange(0,360+ang_spacing,ang_spacing)
        n_angs = len(angs)
        
        # lists to save info in       
        minSag = [None]*len(self.cableList)
        minCurvSF = [None]*len(self.cableList)
        CminTenSF = [None]*len(self.cableList)
        minTenSF = [None]*len(self.mooringList)
        F = [None]*len(self.anchorList) 
        x = np.zeros((len(self.platformList),n_angs))
        y = np.zeros((len(self.platformList),n_angs))
        
        if not self.ms:
            self.getMoorPyArray()
             
        # apply thrust force to platforms at specified angle intervals
        for i,ang in enumerate(angs):
            print('Analyzing platform offsets at angle ',ang)
            fx = thrust*np.cos(np.radians(ang))
            fy = thrust*np.sin(np.radians(ang))
            
            # add thrust force and moment to the body
            for pf in self.platformList.values():
                if pf.entity == 'FOWT':
                    pf.body.f6Ext = np.array([fx, fy, 0, fy*RNAheight, fx*RNAheight, 0])       # apply an external force on the body [N]  
                else:
                    pass                     
            # solve equilibrium 
            self.ms.solveEquilibrium3(DOFtype='both')
        
            # save info if requested
            if SFs:
                # get loads on anchors (may be shared)
                for j,anch in enumerate(self.anchorList.values()):
                    atts = [att['obj'] for att in anch.attachments.values()]
                    F1 = [None]*len(atts)
                    for jj,moor in enumerate(atts):
                        if isinstance(moor.attached_to[0],Anchor):
                            # anchor attached to end A
                            F1[jj] = moor.ss.fA*DAF
                        else:
                            F1[jj] = moor.ss.fB*DAF
                    # add up all tensions on anchor in each direction (x,y,z)
                    F2 = [sum([a[0] for a in F1]),sum([a[1] for a in F1]),sum([a[2] for a in F1])]
                    H = np.hypot(F2[0],F2[1]) # horizontal force
                    T = np.sqrt(F2[0]**2+F2[1]**2+F2[2]**2) # total tension force
                    if not F[j] or T>np.sqrt(F[j][0]**2+F[j][1]**2+F[j][2]**2):
                        F[j] = F2 # max load on anchor                         
                        # save anchor load information
                        anch.loads['Hm'] = H
                        anch.loads['Vm'] = F[j][2]
                        anch.loads['thetam'] = np.degrees(np.arctan(anch.loads['Vm']/anch.loads['Hm'])) #[deg]
                        anch.loads['mudline_load_type'] = 'max'
                            
                # get tensions on mooring line
                for j, moor in enumerate(self.mooringList.values()):
                    MBLA = float(moor.ss.lineList[0].type['MBL'])
                    MBLB = float(moor.ss.lineList[-1].type['MBL'])
                    # print(MBLA,MBLB,moor.ss.TA,moor.ss.TB,MBLA/moor.ss.TA,MBLB/moor.ss.TB,abs(MBLA/moor.ss.TA),abs(MBLB/moor.ss.TB))
                    MTSF = min([abs(MBLA/(moor.ss.TA*DAF)),abs(MBLB/(moor.ss.TB*DAF))])
                    # atenMax[j], btenMax[j] = moor.updateTensions()
                    if not minTenSF[j] or minTenSF[j]>MTSF:
                        minTenSF[j] = deepcopy(MTSF)
                        moor.loads['TAmax'] = moor.ss.TA*DAF
                        moor.loads['TBmax'] = moor.ss.TB*DAF
                        moor.loads['info'] = f'determined from arrayWatchCircle() with DAF of {DAF}'
                        moor.safety_factors['tension'] = minTenSF[j]
                        
                                
                # get tensions and curvature on cables
                for j,cab in enumerate(self.cableList.values()):
                    dcs = [a for a in cab.subcomponents if isinstance(a,DynamicCable)] # dynamic cables in this cable 
                    ndc = len(dcs) # number of dynamic cable objects in this single cable object
                    CminTenSF[j] = [None]*ndc
                    minCurvSF[j] = [None]*ndc
                    minSag[j] = [None]*ndc
                    if dcs[0].ss:
                        for jj,dc in enumerate(dcs):               
                            MBLA = dc.ss.lineList[0].type['MBL']
                            MBLB = dc.ss.lineList[-1].type['MBL']
                            CMTSF = min([abs(MBLA/dc.ss.TA),abs(MBLB/dc.ss.TB)])
                            if not CminTenSF[j][jj] or CminTenSF[j][jj]>CMTSF:
                                CminTenSF[j][jj] = deepcopy(CMTSF)
                                dc.loads['TAmax'] = dc.ss.TA*DAF
                                dc.loads['TBmax'] = dc.ss.TB*DAF
                                dc.loads['info'] = f'determined from arrayWatchCircle() with DAF of {DAF}'
                                dc.safety_factors['tension'] = CminTenSF[j][jj]
                            # CatenMax[j], CbtenMax[j] = cab.updateTensions()
                            dc.ss.calcCurvature()
                            mCSF = dc.ss.getMinCurvSF()
                            if not minCurvSF[j][jj] or minCurvSF[j][jj]>mCSF:
                                minCurvSF[j][jj] = mCSF
                                dc.safety_factors['curvature'] = minCurvSF[j][jj]
                        # # determine number of buoyancy sections
                        # nb = len(dc.dd['buoyancy_sections'])
                        # m_s = []
                        # for k in range(0,nb):
                        #     m_s.append(dc.ss.getSag(2*k))
                        # mS = min(m_s)
                        # if not minSag[j][jj] or minSag[j][jj]<mS:
                        #     minSag[j][jj] = deepcopy(mS)
                        #     dc.
                        
                
            # save location of each platform for envelopes
            for k, pf in enumerate(self.platformList.values()):                       
                x[k,i] = pf.body.r6[0]     
                y[k,i] = pf.body.r6[1]
                
                    
        for k, body in enumerate(self.platformList.values()):
            # save motion envelope in the correct platform 
            body.envelopes['mean'] = dict(x=x[k,:], y=y[k,:])
        
            if shapes:  # want to *optionally* make a shapely polygon
                from shapely import Polygon
                body.envelopes['mean']['shape'] = Polygon(list(zip(x[k,:],y[k,:])))
            
        for body in self.ms.bodyList:
            body.f6Ext = np.array([0, 0, 0, 0, 0, 0])
        self.ms.solveEquilibrium3(DOFtype='both')  
        
        if moor_envelopes:
            for moor in self.mooringList.values():
                moor.getEnvelope()
        
        if SFs:
            maxVals = {'minTenSF':minTenSF,'minTenSF_cable':CminTenSF,'minCurvSF':minCurvSF,'maxF':F}# np.vstack((minTenSF,CminTenSF,minCurvSF,minSag))    
            return(x,y,maxVals)     
        else:
            return(x,y)
        
        
        
        
    def getArrayCost(self):
        '''
        Function to sum all available costs for the array components and produce a 
        spreadsheet with itemized costs for each component.

        Returns
        -------
        total cost

        '''
        total_cost = 0
        anch_costs = {}
        pf_costs = {}
        mooring_costs = {}
        cable_costs = {}
        turbine_costs = {}
        substation_costs = {}
        # anchors
        anch_cost = 0
        for anch in self.anchorList.values():
            anch_cost += anch.getCost()
            
            anch_costs[anch.id] = anch.cost
        # maxrows_anch = max([len(anch.cost) for anch in self.anchorList.values()])
        #headings_anch = set([anch.cost.keys() for anch in self.anchorList.values()])
        # platforms
        pf_cost = 0
        for pf in self.platformList.values():
            if pf.cost:
                pf_cost += sum(pf.cost.values())
            pf_costs[pf.id] = pf.cost
        # maxrows_pf = max([len(pf.cost) for pf in self.platformList.values()])
        # moorings
        moor_cost = 0
        for moor in self.mooringList.values():
            moor_cost += moor.getCost()
            mooring_costs[moor.id] = moor.cost
        # maxrows_moor = max([len(moor.cost) for moor in self.mooringList.values()])
        # cables
        cab_cost = 0
        for cab in self.cableList.values():
            cab_cost += cab.getCost()
            cable_costs[cab.id] = cab.cost
        # maxrows_cab = max([len(cab.cost) for cab in self.cableList.values()])
        # turbine
        turb_cost = 0
        for turb in self.turbineList.values():
            if turb.cost:
                turb_cost += sum(turb.cost.values())
            turbine_costs[turb.id] = turb.cost
        # maxrows_turb = max([len(turb.cost) for turb in self.turbineList.values()])
        # substation
        oss_cost = 0
        for oss in self.substationList.values():
            if oss.cost:
                oss_cost += sum(oss.cost.values())
            substation_costs[oss.id] = oss.cost
            
        # add up total cost
        total_cost = anch_cost + pf_cost + moor_cost + cab_cost + turb_cost + oss_cost
        # maxrows_oss = max([len(oss.cost) for oss in self.substationList.values()])
            
        # now write the costs to a spreadsheet
        
        # let's do a new sheet for each component type and one overview sheet
        # for now let's write out a dictionary
        cost_dict = {'total cost':total_cost,'cable cost':cab_cost,
                     'anchor cost':anch_cost, 'pf cost':pf_cost,
                     'moor cost':moor_cost,'turbine cost':turb_cost,
                     'substation cost':oss_cost}
        
        return(cost_dict)
    
    def unload(self,file='project.yaml'):
        '''
        Function to unload information to a yaml file
        
        Parameters
        ----------
        file : str
            File name for output yaml file

        Returns
        -------
        None.

        '''
        print(f'Unloading project to yaml file {file}')
                   
        # build out array table
        arrayKeys = ['ID','topsideID','platformID','mooringID','x_location','y_location','z_location','heading_adjust']
        arrayData = [] #np.zeros((len(arrayKeys),len(self.platformList)))
        pf_type = []
        
        # build out platform info
        topList = []        
        for i,pf in enumerate(self.platformList.values()):
            ts_loc = 0
            # determine any connected topsides
            for att in pf.attachments.values():
                if not type(att['obj']).__name__ in ['Mooring','Cable','Fairlead']:
                    dd = att['obj'].dd
                    if isinstance(att['obj'],Turbine):
                        entity = 'Turbine'
                    elif isinstance(att['obj'],Substation):
                        entity = 'Substation'
                    else:
                        entity = att['obj'].entity
                    dd['type'] = entity
                    if att['obj'].dd in topList:
                        ts_loc = topList.index(att['obj'].dd) + 1
                    else:
                        topList.append(att['obj'].dd)
                        ts_loc = len(topList)
                

            arrayData.append([pf.id, ts_loc, pf.dd['type']+1, 0, float(pf.r[0]), 
                              float(pf.r[1]), float(pf.r[2]), float(np.degrees(pf.phi))])
            pf_type.append(pf.dd['type'])
            
        # get set of just platform types used in this project
        pf_type = set(pf_type)

        # figure out keys
        if len(pf_type) > 1:
            pfkey = 'platforms'
            pfTypes = [self.platformTypes[x] for x in pf_type]
        else:
            pfkey = 'platform'
            pfTypes = self.platformTypes[list(pf_type)[0]]
           
        # build out site info
        site = {}
        if hasattr(self,'soilProps'):                       
            if len(self.soil_x)>1:
                site['seabed'] = {'x':[float(x) for x in self.soil_x],'y':[float(x) for x in self.soil_y],'type_array':self.soil_names.tolist(),
                                  'soil_types':self.soilProps} #[[[float(v[0])] for v in x.values()] for x in self.soilProps.values()]}
                    
        if len(self.grid_x)>1:
            site['bathymetry'] = {'x':[float(x) for x in self.grid_x],'y':[float(y) for y in self.grid_y],'depths':[[float(y) for y in x] for x in self.grid_depth]}
        if len(self.boundary)>0:
            site['boundaries'] = {'x_y':[[float(y) for y in x] for x in self.boundary]}
        site['general'] = {'water_depth':float(self.depth),'rho_air':float(self.rho_air),
                           'rho_water':float(self.rho_water),'mu_air':float(self.mu_air)}
        
         
            
        # build out array mooring and array anchor section
        arrayMoor = []
        allconfigs = []
        arrayAnch = []
        anchConfigs = {}
        anchKeys = ['ID','type','x','y','embedment']
        lineKeys = ['MooringConfigID','endA','endB','headingA','headingB','lengthAdjust']
        
        for moor in self.mooringList.values():
            newcon = True
            newanch = True
            # get connected objects
            endA = moor.attached_to[0]
            endB = moor.attached_to[1]
            # get heading(s)
            if not moor.shared:  #  and type(endA) != Connector:
                headA = 'NA'
                # add anchor
                arrayAnch.append([endA.id, endA.dd['name'], float(endA.r[0]), float(endA.r[1]),0])
                if anchConfigs:
                    if any([endA.dd['name']==k for k in anchConfigs]):
                        newanch = False
                        current_anch = endA.dd['name']
                if newanch:
                    anchConfigs[endA.dd['name']] = dict(endA.dd['design'])
            # elif type(endA)==Connector:
            #     # get connector info & store like anchors
            #     pass
            
            else:
                # shared line - get end A heading
                # relloc = np.array(endA.r) - np.array(moor.rA)
                # fairleadA = np.where(endA.fairleads[ii] == relloc[ii] for ii in range(3))
                ang = np.pi/2 - np.arctan2(moor.rA[1]-endA.r[1],moor.rB[0]-endA.r[0])
                headA = float(np.degrees(ang - endA.phi))
                
            # get end B heading
            angB = np.pi/2 - np.arctan2(moor.rB[1]-endB.r[1],moor.rB[0]-endB.r[0])
            headB = np.degrees(angB - endB.phi)

            # if type(endB)==Connector:
            #     # get connector info & store like anchors
            #     pass
            # else:
            #     relloc = np.array(endB.r) - np.array(moor.rB)
            #     fairleadB = np.where(endB.fairlead[ii] == relloc[ii] for ii in range(3))
            # get mooring configuration
            config = {'span':float(moor.dd['span']),'sections':moor.dd['sections'],'connectors':moor.dd['connectors']}
            
            # check if an existing mooring configuration matches the current configuration
            if allconfigs:
                pc = np.where([config['span']==x['span'] for x in allconfigs] and [len(y['sections'])==len(config['sections']) for y in allconfigs])[0]
                for j in pc:
                    if all([allconfigs[j]['sections'][k]==config['sections'][k] for k in range(len(config['sections']))]):
                        if all([allconfigs[j]['connectors'][k]==config['connectors'][k] for k in range(len(config['connectors']))]):
                            current_config = str(j)
                            newcon = False
            if newcon:
                allconfigs.append(config)
                current_config = str(len(allconfigs) - 1)

            arrayMoor.append([current_config,endA.id, endB.id, headA,headB,int(0)])
            # arrayMoor.append([current_config,endA.id, endB.id, fairleadA,fairleadB,int(0)])

        # set up mooring configs, connector and section types dictionaries
        connTypes = {}  
        secTypes = {}
        mooringConfigs = {}
        sIdx = 0
        sUnique = []
        for j,conf in enumerate(allconfigs):
            sections = []
            # iterate through sections
            for i in range(len(conf['sections'])):
                # add connector if it isn't empty
                if not conf['connectors'][i]['m'] == 0 or not conf['connectors'][i]['CdA'] == 0 or not conf['connectors'][i]['v'] == 0:
                    # this is not an empty connector
                    if not 'type' in conf['connectors'][i]:
                        # make a new connector type
                        connTypes[str(int(len(connTypes)))] = dict(conf['connectors'][i])
                        ctn = str(int(len(connTypes)-1)) # connector type name
                    else:
                        ctn = str(conf['connectors'][i]['type'])
                        connTypes[ctn] = dict(conf['connectors'][i])
                            
                    sections.append({'connectorType':ctn})
                # add section info
                stm = conf['sections'][i]['type']['material'] # section type material
                stw = conf['sections'][i]['type']['w']        # section type weight

                sKey = (stm, stw)
                if sKey not in sUnique:
                    sUnique.append(sKey)
                    conf['sections'][i]['type']['name'] = sIdx
                    stn = conf['sections'][i]['type']['name'] # section type name
                    secTypes[stn] = dict(conf['sections'][i]['type'])
                    #secTypes[stn] = cleanDataTypes(secTypes[stn])
                    sIdx += 1
                    
                stn = sUnique.index(sKey)
                sections.append({'type':stn,'length':float(conf['sections'][i]['L'])})

            # add last connector if needed
            if not conf['connectors'][i+1]['m'] == 0 or not conf['connectors'][i+1]['CdA'] == 0 or not conf['connectors'][i+1]['v'] == 0:
                # this is not an empty connector
                if not 'type' in conf['connectors'][i+1]:
                    # make a new connector type
                    #conf['connectors'][i+1] = cleanDataTypes(conf['connectors'][i+1])
                    connTypes[str(len(connTypes))] = conf['connectors'][i+1]
                    ctn = str(int(len(connTypes)-1))
                else:
                    ctn = conf['connectors'][i+1]['type']
                    connTypes[ctn] = dict(conf['connectors'][i+1])    
                sections.append({'connectorType':ctn})
            # put mooring config dictionary together
            mooringConfigs[str(j)] = {'name':str(j),'span':float(conf['span']),'sections':sections}

                    
        # cables setup
        cables = []
        cableTypes = {}
        cableConfigs = {}
        cUnique = []
        bUnique = []
        appendageTypes = {}
        jUnique = []
        jIdx = 0
        for jj,cab in enumerate(self.cableList.values()):
            endA = cab.attached_to[0]
            endB = cab.attached_to[1]
            angA = np.pi/2 - np.arctan2(cab.subcomponents[0].rB[1]-cab.rA[1],cab.subcomponents[0].rB[0]-cab.rA[0])
            headA = float(np.degrees(angA - endA.phi))
            angB = np.pi/2 - np.arctan2(cab.subcomponents[-1].rA[1]-cab.rB[1],cab.subcomponents[-1].rA[0]-cab.rB[0])
            headB = float(np.degrees(angB - endB.phi))
            coords = []
            statcab = 'None'
            dynCabs = [None,None]
            burial = None
            
            for kk,sub in enumerate(cab.subcomponents):
                currentConfig = {}
                    
                if isinstance(sub,StaticCable):
                    # pull out cable config and compare it to existing cableConfigs
                    ctw = sub.dd['cable_type']['w']
                    ctA = sub.dd['cable_type']['A'] 
                    cKey = (ctw,ctA)
                    ctf = False
                    # check if made with getCableProps (then we can skip writing out cable type info)
                    if 'notes' in sub.dd['cable_type']:
                        if 'made with getCableProps' in sub.dd['cable_type']['notes']:
                            ctk = 'cableFamily'
                            ctn = 'static_cable_'+str(int(sub.voltage))
                            ctf = True
                    
                    # create current cable config dictionary
                    if not ctf:
                        if not cKey in cUnique:
                            cUnique.append(cKey)                        
                            ctn = 'stat_cab_'+str(len(cUnique)-1)
                            cableTypes[ctn] = sub.dd['cable_type']
                        else:
                            cIdx = cUnique.index(cKey)
                            ctn = 'stat_cab_'+str(cIdx)
                    

                    
                    # check for routing coordinates
                    if hasattr(sub,'coordinates'):
                        coords.extend(sub.coordinates)
                        
                    if hasattr(sub,'burial'):
                        burial = sub.burial
                        
                    statcab = ctn

                        
                elif isinstance(sub,DynamicCable):
                    # pull out cable config and compare it to existing cableConfigs
                    ct = sub.dd['type'] # static or dynamic
                    ctw = sub.dd['cable_type']['w']
                    ctA = sub.dd['A'] 
                    cKey = (ctw,ctA)
                    ctf = False; ctk = 'cable_type'
                    # check if made with getCableProps (then we can skip writing out cable type info)
                    if 'notes' in sub.dd['cable_type']:
                        if 'made with getCableProps' in sub.dd['cable_type']['notes']:
                            ctn = ct+'_cable_'+str(int(sub.voltage))
                            ctf = True
                    # check if cable type has already been written
                    if not ctf:
                        if not cKey in cUnique:
                            cUnique.append(cKey)                        
                            ctn = 'dyn_cab_'+str(len(cUnique)-1)
                            cableTypes[ctn] = sub.dd['cable_type']
                        else:
                            cIdx = cUnique.index(cKey)
                            ctn = 'dyn_cab_'+str(cIdx)
                    # collect buoyancy sections info if applicable
                    bs = []
                    if 'buoyancy_sections' in sub.dd:
                        for b in sub.dd['buoyancy_sections']:
                            btw = b['module_props']['w']; btv = b['module_props']['volume']
                            if not (btw,btv) in bUnique:
                                bUnique.append((btw,btv))
                                btn = 'buoy_'+str(len(bUnique)-1)
                                appendageTypes[btn] = b['module_props']
                                appendageTypes[btn]['type'] = 'buoy'
                            else:
                                bid = bUnique.index((btw,btv))
                                btn = 'buoy_'+str(bid)
                            bs.append({'L_mid':b['L_mid'],'N_modules':b['N_modules'],
                                      'spacing':b['spacing'],'V':b['module_props']['volume'],
                                      'type':btn})
                    if 'appendages' in sub.dd:
                        for app in sub.dd['appendages']:
                            pass # UPDATE TO PULL OUT APPENDAGE INFO AND STORE
                            
                    # grab joint info
                    if kk == 0 and len(cab.subcomponents)>1:
                        sc = cab.subcomponents[kk+1] 
                    else:
                        sc = cab.subcomponents[kk-1]
                    jsub = isinstance(sc,Joint)
                    if jsub:
                        # grab joint info and add
                        if 'm' in sc or 'v' in sc and (sc['m']!=0 or sc['v']!=0):
                            jKey = (getFromDict(sc,'m',default=0),getFromDict(sc,'v',default=0))
                            if not jKey in jUnique:
                                jUnique.append(jKey)
                                jtn = 'joint_'+str(len(jUnique)-1)
                                appendageTypes[jtn] = dict(deepcopy(sc))
                                appendageTypes[jtn]['type'] = 'joint'
                                if 'r' in sc:
                                    appendageTypes[jtn].pop('r')
                            else:
                                jtd = deepcopy(sc)
                                if 'r' in sc:
                                    jtd.pop('r')
                                jIdx = jUnique.index(jKey)
                                jtn = 'joint_'+str(jIdx)
                            bs.append({'type':jtn})
                    # create current cable config dictionary
                    currentConfig = {ctk:ctn,'A':ctA,'rJTube':sub.dd['rJTube'],
                                     'span':sub.dd['span'],'length':sub.L,
                                     'voltage':sub.dd['cable_type']['voltage'],'sections':bs}
                    # check if current cable config already exists in cable configs dictionary
                    if currentConfig in cableConfigs.values():
                        ccn = [key for key,val in cableConfigs.items() if val==currentConfig][0] # get cable config key
                    else:
                        # create new cable config entry in dictionary
                        ccn = 'dynamic_'+str(len(cableConfigs))
                        cableConfigs[ccn] = currentConfig
                        
                    if kk>0:
                        didx = 1
                    else:
                        didx = 0 
                    dynCabs[didx] = ccn
                    # currentCable.append({'type':ccn})
            if not coords:
                route = None
                
            cid = 'array_cable'+str(len(cables))
            endAdict = {'attachID':endA.id,
                        'heading':headA,
                        'dynamicID':dynCabs[0] if dynCabs[0] else 'None'}
            endBdict = {'attachID':endB.id,
                        'heading':headB,
                        'dynamicID':dynCabs[1] if dynCabs[1] else 'None'}
            
            cables.append({'name':cid,'endA':endAdict,'endB':endBdict,'type':statcab})
            
            if route:
                cables[-1]['routing_x_y_r'] = coords
            if burial:
                cables[-1]['burial'] = burial
                
            
            
        
        
            
         
        # create master output dictionary for yaml
        output = {'site':site, 'array':{'keys':arrayKeys,'data':arrayData}, 
                  pfkey:pfTypes, 
                  'topsides': topList, 
                  'array_mooring':{'anchor_keys':anchKeys, 'anchor_data':arrayAnch,
                                   'line_keys':lineKeys, 'line_data':arrayMoor},
                  'mooring_line_configs':mooringConfigs,
                  'mooring_line_types':secTypes, 
                  'mooring_connector_types':connTypes,
                  'anchor_types':anchConfigs,
                  'cables':cables,'dynamic_cable_configs':cableConfigs,'cable_types':cableTypes, 
                  'cable_appendages':appendageTypes}

        output = cleanDataTypes(output)
        import ruamel.yaml
        yaml = ruamel.yaml.YAML()
        # write out to file
        with open(file,'w') as f:    
            yaml.dump(output,f)
        
    def extractFarmInfo(self, cmax=5, fmax=10/6, Cmeander=1.9):
        '''
        Function to extract farm-level information required to create FAST.Farm case simulations. [Under developement]:

        Parameters
        ----------
        cmax : float, optional
            maximum blade chord (m)
        fmax: maximum excitation frequency (Hz)
        Cmeander: Meandering constant (-)
        
        Returns
        -------
        wts : dict
            General farm-level information needed for FAST.Farm from project class
        yaw_init : list 
            initial yaw offset values (for not it's set as just the platform orientation adjusted for rotational convention variation between FAM and FF)
        '''      


        # ----------- Extract Wind Farm Data
        wts = {}
        i = 0
        yaw_init = np.zeros((1, len(self.platformList.items())))
        for _, pf in self.platformList.items():
            x, y, z   = pf.body.r6[0], pf.body.r6[1], pf.body.r6[2]
            phi_deg       = np.degrees(pf.phi)  # float((90 - np.degrees(pf.phi)) % 360)  # Converting FAD's rotational convention (0deg N, +ve CW) into FF's rotational convention (0deg E, +ve CCW)
            phi_deg       = (phi_deg + 180) % 360 - 180  # Shift range to -180 to 180
            for att in pf.attachments.values():
                if isinstance(att['obj'],Turbine):
                    D    = 240   # att['obj'].D         (assuming 15MW)
                    zhub = att['obj'].dd['hHub']
                
            wts[i] = {
                'x': x, 'y': y, 'z': z, 'phi_deg': phi_deg, 'D': D, 'zhub': zhub, 
                'cmax': cmax, 'fmax': fmax, 'Cmeander': Cmeander
                }
            yaw_init[0, i] = -phi_deg
            i += 1

        # store farm-level wind turbine information
        self.wts = wts

        return wts, yaw_init  
    
    def FFarmCompatibleMDOutput(self, filename, unrotateTurbines=True, renameBody=True, removeBody=True, MDoptionsDict={}):
        '''
        Function to create FFarm-compatible MoorDyn input file:

        Parameters
        ----------
        filename : str
            Name of the MoorDyn output file (.dat)
        unrotateTurbines: bool, optional
            A flag to unrotate turbine (body) objects when passing it to MoorPy unload function [FFarm takes fairlead points in the local-unrotated reference frame]
        renameBody: bool, optional
            A flag to rename `Body` objects in the output MD file into `Turbine` to be compatible with FFarm. 
        removeBody: boo., optional
            A flag to remove 'Body' objects in the Bodies list in the output MD file to be compatible with FFarm.
        MDoptionsDict: dictionary, optional
            MoorDyn Options. If not given, default options are considered.            
        '''          
        from moorpy.helpers import subsystem2Line    
        
        # convert SS to lines
        ms_temp = deepcopy(self.ms)  # copy to avoid affecting self.ms
        lineCount = len(ms_temp.lineList)
        for _ in range(lineCount):
            subsystem2Line(ms_temp, 0)
        
        # Unrotate turbines if needed
        if unrotateTurbines:
            if self.wts:
                phiV = [wt['phi_deg'] for wt in self.wts.values()]  # to unrotate the platforms when unloading MoorDyn
            else:
                raise ValueError("wts is empty. Please run project.extractFarmInfo first before extracting MoorDyn")
        else:
            phiV = None
        
        ms_temp.unload(fileName=filename, phiV=phiV, MDoptionsDict=MDoptionsDict)
        
        # rename Body to Turbine if needed
        if renameBody:
            # Rename Body to Turbine:
            with open(filename, 'r') as f:
                filedata = f.read()

                filedata = filedata.replace('Body', 'Turbine')
                with open(filename, 'w') as f:
                    f.write(filedata)

                f.close()       
        
        if removeBody:
            with open(filename, 'r') as f:
                lines = f.readlines()

            newLines = []
            skipCount = 0

            for i, line in enumerate(lines):
                if '---' in line and ('BODIES' in line.upper() or 'BODY LIST' in line.upper() or 'BODY PROPERTIES' in line.upper()):
                    newLines.append(line)
                    newLines.append(next(iter(lines[i+1:])))  # Append 2 lines
                    newLines.append(next(iter(lines[i+2:]))) 

                    skipCount = 2 + len(self.platformList)  # Skip the number of platforms and the already appended lines above
                    continue

                if skipCount > 0:
                    skipCount -= 1
                else:
                    newLines.append(line)

            with open(filename, 'w') as f:
                f.writelines(newLines)
           

    def updateFailureProbability(self):
        '''
        Function to populate (or update) failure probability dictionaries in each object 
        based on failure probability calculations developed by Emma Slack
        
        To be filled in...

        Returns
        -------
        None.

        '''
                

    


'''
Other future items:
Cost calc functions
System Reliability/failure analysis functions
Full scenario visualization functions
Load case setup and constraint eval?
'''

if __name__ == '__main__':

    project = Project()
    project.loadSoil(filename='tests/soil_sample.txt')
    # create project class instance from yaml file
    #project = Project(file='OntologySample600m.yaml')
    project = Project(file='../tests/simple_farm.yaml')
    
    project.getMoorPyArray(cables=1,plt=1,pristineLines=1)

    # make envelopes and watch circles
    # (Mooring.getEnvelope will call Platform.getWatchCircle when needed)
    for moor in project.mooringList.values():
        moor.getEnvelope()

    project.plot2d(plot_boundary=False)  # this should also plot the watch circles/envelopes!
    
    plt.show()