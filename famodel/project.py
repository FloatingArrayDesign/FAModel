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
from famodel.anchors.anchor_capacity import anchorCapacity
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
        
        # Dictionaries of component/product properties used in the array
        self.turbineTypes = None
        self.lineTypes = None
        self.anchorTypes = None
        self.cableTypes = None
        
        
        # ----- site information -----
        self.lat0  = lat  # lattitude of site reference point [deg]
        self.lon0  = lon  # longitude of site reference point [deg]
        self.g = 9.81

        # Project boundary (list of x,y coordinate pairs [m])
        self.boundary = np.zeros([0,2])
        
        # Seabed grid
        self.grid_x      = np.array([2200])
        self.grid_y      = np.array([200])
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
        if 'array' in d:
            arrayInfo = [dict(zip(d['array']['keys'], row)) for row in d['array']['data']]
        
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
        
            for ci in d['cables']:
                for k, v in d['cables'].items():
                    cableInfo[k] = v

        # ----- cable configurations -----
        cable_configs = {}
        if 'cable_configs' in d and d['cable_configs']:
            for k, v in d['cable_configs'].items():
                cable_configs[k] = v
                
        # # ----- cable types -----
        # cable_types = {}
        # if 'cable_types' in d and d['cable_types']:
        #     for k, v in d['cable_types']:
        #         cable_types[k] = v

        # ----- substation -----
        if 'substation' in d and d['substation']:
            # create substation design dictionary and object for each substation
            for k in d['substation'].keys():
                subID = k
                r = [d['substation'][k]['x_location'],d['substation'][k]['y_location']]
                dd = {'r':r}
                self.substationList[subID] = Substation(dd, subID)
                self.substationList[subID].rFair = d['substation'][k]['rFair']
                self.substationList[subID].zFair = d['substation'][k]['zFair']
                self.substationList[subID].phi = d['substation'][k]['heading']

        
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
                platforms = {} # dictionary of platform information
                platforms = d['platform']
                RAFTDict['platform'] = d['platform']
        # load list of platform dictionaries into RAFT dictionary
        elif 'platforms' in d and d['platforms']:
            platforms = [] # list of dictionaries of platform information
            platforms = d['platforms']
            RAFTDict['platforms'] = d['platforms']
            
        # ----- turbines -----
        if 'turbine' in d and d['turbine']:
            # check that there is only one turbine
            if 'turbines' in d and d['turbines']:
                raise Exception("Cannot read in items for both 'turbines' and 'turbine' keywords. Use either 'turbine' keyword for one platform or 'turbines' keyword for a list of platforms.")
            elif type(d['turbine']) is list and len(d['turbine'])>1:
                raise Exception("'turbine' section keyword must be changed to 'turbines' if multiple turbines are listed")
            else:
                turbines = {} # dictionary of turbine information
                turbines = d['turbine']
                RAFTDict['turbine'] = d['turbine']
        # load list of turbine dictionaries into RAFT dictionary
        elif 'turbines' in d and d['turbines']:
            turbines = [] # list of dictionaries of turbine information
            turbines = d['turbines']
            RAFTDict['turbines'] = d['turbines']

                # I think we want to just store it as a dictionary
                # validate it,
                # then have it available for use when making Mooring objects and subsystems
        
        def MooringProps(mCon,checkType=1):
            '''
            Parameters
            ----------
            mType : dictionary
                Dictionary of mooring details from the mooring_line_configs key
                Includes type (reference to name in mooring_line_types or cable_props yaml)
            checkType : boolean
                Controls whether or not to first look for the cable type in the project yaml dictionary before
                attempting to get the cable properties from cable props yaml.

            Returns
            -------
            dd : design dictionary

            '''
            if 'type' in mCon and mCon['type'] in d['mooring_line_types']:
                dd = d['mooring_line_types'][mCon['type']]
                dd['name'] = mCon['type']
                if 'd_vol' in dd:
                    d_vol = dd['d_vol']
                # else:
                #     d_vol = dd['d']
                dd['w'] = (dd['m']-np.pi/4*d_vol**2*self.rho_water)*self.g
                if 'mooringFamily' in mCon:
                    raise Exception('type and moorFamily listed in yaml - use type to reference a mooring type in the mooring_line_types section of the yaml and mooringFamily to obtain mooring properties from MoorProps_default.yaml')
            elif 'mooringFamily' in mCon:
                from moorpy.helpers import loadLineProps, getLineProps
                if not 'd_nom' in mCon:
                    raise Exception('To use MoorProps yaml, you must specify a nominal diameter in mm for the mooring line family')
                lineprops = loadLineProps(None)
                mProps = getLineProps(mCon['d_nom']*1000,mCon['mooringFamily'],lineProps=lineprops)
                dd = mProps
                dd['name'] = mCon['mooringFamily']
                dd['d_nom'] = mProps['input_d']
            elif 'type' in mCon and not mCon['type'] in d['mooring_line_types']:
                raise Exception(f'Type {mCon["type"]} provided in mooring_line_config {mCon} is not found in mooring_line_types section. Check for errors.')

            return(dd)
        
        def getMoorings(lineconfig,i):
            '''

            Parameters
            ----------
            lineconfig : string
                Line configuration type
            i : int
                Index in array table (essentially which platform)

            Returns
            -------
            m_config : dict
                mooring configuration dictionary
            c_config : dict
                connector configuration dictionary

            '''
            # set up dictionary of information on the mooring configurations
            m_config = {'sections':[],'anchor':{},'span':{},'zAnchor':{}}#,'EndPositions':{}}
            # set up connector dictionary
            c_config = []
                        
            lineLast = 1    # boolean whether item with index k-1 is a line. Set to 1 for first run through of for loop
            ct = 0   # counter for number of line types
            for k in range(0,len(lineConfigs[lineconfig]['sections'])): # loop through each section in the line
            
                lc = lineConfigs[lineconfig]['sections'][k] # set location for code clarity later
                # determine if it's a line type or a connector listed
                if 'type' in lc or 'mooringFamily' in lc: 
                    # this is a line
                    if lineLast: # previous item in list was a line (or this is the first item in a list)
                        # no connector was specified for before this line - add an empty connector
                        c_config.append({})                        
                    # set line information
                    lt = MooringProps(lc)                                             
                    # lt = self.lineTypes[lc['type']] # set location for code clarity and brevity later
                    # set up sub-dictionaries that will contain info on the line type
                    m_config['sections'].append({'type':lt})# {'name':str(ct)+'_'+lc['type'],'d_nom':lt['d_nom'],'material':lt['material'],'d_vol':lt['d_vol'],'m':lt['m'],'EA':float(lt['EA'])}})
                    m_config['sections'][ct]['type']['name'] = str(ct)+'_'+lt['name']
                    # make EA a float not a string
                    m_config['sections'][ct]['type']['EA'] = float(lt['EA'])
                    # # need to calculate the submerged weight of the line (not currently available in ontology yaml file)
                    # m_config['sections'][ct]['type']['w'] = (lt['m']-np.pi/4*lt['d_vol']**2*self.rho_water)*self.g
                    # add cost if given
                    # if 'cost' in lt:
                    #     m_config['sections'][ct]['type']['cost'] = lt['cost']
                    # # add MBL if given
                    # if 'MBL' in lt:
                    #     m_config['sections'][ct]['type']['MBL'] = lt['MBL']
                    # # add dynamic stretching if there is any
                    # if 'EAd' in lt: 
                    #     m_config['sections'][ct]['type']['EAd'] = lt['EAd']
                    #     m_config['sections'][ct]['type']['EAd_Lm'] = lt['EAd_Lm']
                        
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
                    c_config.append({})
            # set general information on the whole line (not just a section/line type)
            # set to general depth first (will adjust to depth at anchor location after repositioning finds new anchor location)
            m_config['zAnchor'] = -self.depth 
            m_config['span'] = lineConfigs[lineconfig]['span']
            m_config['name'] = lineconfig
            # add fairlead radius and depth to dictionary
            m_config['rad_fair'] = self.platformList[arrayInfo[i]['ID']].rFair
            m_config['z_fair'] = self.platformList[arrayInfo[i]['ID']].zFair
            
            m_config['connectors'] = c_config  # add connectors section to the mooring dict
            
            return(m_config) #, c_config)
        
        def getConnectors(c_config,mName):
            '''

            Parameters
            ----------
            c_config : dict
                Dictionary of connector configurations for a mooring line.
            mName : tuple
                Key name in the mooringList dictionary for the associated mooring object

            Returns
            -------
            None.

            '''
            
            # make connector objects for all sections of a mooring line configuration in order
            for i in range(0,len(c_config)):
                # check if connector is a none-type
                if c_config[i] == None:                   
                    # create empty connector object
                    self.mooringList[mName].dd['connectors'].append(Connector())
                else:
                    # create connector object with c_config entries
                    self.mooringList[mName].dd['connectors'].append(Connector(**c_config[i]))
        
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
            if self.soil_x is not None: # get soil conditions at anchor location if soil info available
                if mc:
                    ad['soil_type'], ad['soil_properties'] = self.getSoilAtLocation(mc.rA[0], mc.rA[1])
                else:
                    ad['soil_type'], ad['soil_properties'] = self.getSoilAtLocation(arrayAnchor[aNum-1]['x'],arrayAnchor[aNum-1]['y'])
            ad['type'] = self.anchorTypes[lineAnch]['type']
            ad['name'] = lineAnch
            
            return(ad)
        
        
        
        # ----- set up dictionary for each individual mooring line, create anchor, mooring, and platform classes ----
                
        # check that all necessary sections of design dictionary exist to create non-shared lines
        if self.lineTypes and lineConfigs and mSystems:
            mct = 0 # counter for number of mooring lines
            # set up a list of the alphabet for assigning names purposes
            import string
            alph = list(string.ascii_lowercase)
                               
            for i in range(0, len(arrayInfo)): # loop through each platform in array
                
                # create platform instance (even if it only has shared moorings / anchors), store under name of ID for that row
                self.platformList[arrayInfo[i]['ID']] = Platform(arrayInfo[i]['ID'],r=[arrayInfo[i]['x_location'],arrayInfo[i]['y_location']],heading=arrayInfo[i]['heading_adjust'])
                # add fairlead radius and fairlead depth of this platform type from platform information section
                if type(platforms) == list:
                    # get index of platform from array table
                    pfID = arrayInfo[i]['platformID']-1
                    self.platformList[arrayInfo[i]['ID']].rFair = platforms[pfID]['rFair']
                    self.platformList[arrayInfo[i]['ID']].zFair = platforms[pfID]['zFair']
                else:
                    self.platformList[arrayInfo[i]['ID']].rFair = platforms['rFair']
                    self.platformList[arrayInfo[i]['ID']].zFair = platforms['zFair']
                # remove pre-set headings (need to append to this list so list should start off empty)
                self.platformList[arrayInfo[i]['ID']].mooring_headings = []
                
                # create turbine instance
                turb_name = str(arrayInfo[i]['turbineID'])+str(i)
                if isinstance(turbines,list):
                    turb_dd = turbines[arrayInfo[i]['turbineID']-1]
                else:
                    turb_dd = turbines
                self.turbineList[turb_name] = Turbine(turb_dd,turb_name,D=turb_dd['blade']['Rtip']*2)
                # attach turbine to platform
                self.platformList[arrayInfo[i]['ID']].attach(self.turbineList[turb_name])
                
                if not arrayInfo[i]['mooringID'] == 0: #if not fully shared mooring on this platform
                    m_s = arrayInfo[i]['mooringID'] # get mooring system ID
                    mySys_unsorted = [dict(zip(d['mooring_systems'][m_s]['keys'], row)) for row in d['mooring_systems'][m_s]['data']]
                    # sort the mooring lines in the mooring system by heading from 0 (North)
                    mySys = sorted(mySys_unsorted,key=lambda x:x['heading'])
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
                        m_config = getMoorings(lineconfig,i)
                        
                        
                        # create mooring class instance as part of mooring list in the project class instance
                        mc = (Mooring(dd=m_config, id=str(arrayInfo[i]['ID'])+alph[j]))
                        # mc.rA = [m_config['span']+m_config['rad_fair'],0,m_config['zAnchor']]
                        # mc.rB = [m_config['rad_fair'],0,m_config['z_fair']]
                        # adjust end positions based on platform location and mooring and platform headings
                        mc.reposition(r_center=self.platformList[arrayInfo[i]['ID']].r, heading=headings[j]+self.platformList[arrayInfo[i]['ID']].phi, project=self)
                        # adjust anchor z location and rA based on location of anchor
                        zAnew, nAngle = self.getDepthAtLocation(mc.rA[0], mc.rA[1], return_n=True)
                        mc.rA[2] = -zAnew
                        mc.dd['zAnchor'] = -zAnew
                        mc.z_anch = -zAnew
                        
                        # set anchor info
                        lineAnch = mySys[j]['anchorType'] # get the anchor type for the line
                        ad = getAnchors(lineAnch, mc=mc) # call method to create anchor dictionary
                        ad['angle'] = nAngle
                        
                        # add anchor class instance to anchorList in project class
                        self.anchorList[str(arrayInfo[i]['ID'])+alph[j]] = (Anchor(dd=ad, r=mc.rA, id=str(arrayInfo[i]['ID'])+alph[j]))
                        # add mooring class instance to mooringlist in project class
                        self.mooringList[str(arrayInfo[i]['ID'])+alph[j]] = mc
                        # attach mooring object to anchor and platform
                        mc.attachTo(self.anchorList[str(arrayInfo[i]['ID'])+alph[j]],end='A')
                        mc.attachTo(self.platformList[arrayInfo[i]['ID']],end='B')
                        
                        # update counter
                        mct += 1
                        
        
        # ----- set up dictionary for each shared mooring line or shared anchor, create mooring and anchor classes ----
    
        aNum = []
        
        # create any shared mooring lines / lines connected to shared anchors
        if arrayMooring:
            # get mooring line info for all lines 
            for j in range(0, len(arrayMooring)): # run through each line            
                
                PFNum = [] # platform ID(s) connected to the mooring line
                
                # Error check for putting an anchor (or something else) at end B
                if not any(ids['ID'] == arrayMooring[j]['end B'] for ids in arrayInfo):
                    raise Exception("Input for end B must match an ID from the array table.")
                if any(ids['ID'] == arrayMooring[j]['end B'] for ids in arrayAnchor):
                    raise Exception(f"input for end B of line_data table row '{j}' in array_mooring must be an ID for a FOWT from the array table. Any anchors should be listed as end A.")
                # Make sure no anchor IDs in arrayAnchor table are the same as IDs in array table
                for k in range(0,len(arrayInfo)):
                    if any(ids['ID'] == arrayInfo[k] for ids in arrayAnchor):
                        raise Exception(f"ID for array table row {k} must be different from any ID in anchor_data table in array_mooring section")
                # determine if end A is an anchor or a platform
                if any(ids['ID'] == arrayMooring[j]['end A'] for ids in arrayInfo): # shared mooring line (no anchor)
                    # get ID of platforms connected to line
                    PFNum.append(arrayMooring[j]['end B'])
                    PFNum.append(arrayMooring[j]['end A'])
                    # find row in array table associated with these platform IDs and set locations
                    for k in range(0, len(arrayInfo)):
                        if arrayInfo[k]['ID'] == PFNum[0]:
                            rowB = arrayInfo[k]
                            Bnum = k
                        elif arrayInfo[k]['ID'] == PFNum[1]:
                            rowA = arrayInfo[k]
                    # get headings (mooring heading combined with platform heading)
                    headingA = np.radians(90-arrayMooring[j]['headingA']) - self.platformList[PFNum[1]].phi
                    headingB = np.radians(90-arrayMooring[j]['headingB']) - self.platformList[PFNum[0]].phi
                    # print('headingA: ',headingA,'listed headingA',arrayMooring[j]['headingA'],'phi: ',self.platformList[PFNum[1]].phi)
                    # calculate fairlead locations (can't use reposition method because both ends need separate repositioning)
                    Aloc = [rowA['x_location']+np.cos(headingA)*self.platformList[PFNum[1]].rFair, rowA['y_location']+np.sin(headingA)*self.platformList[PFNum[1]].rFair, self.platformList[PFNum[1]].zFair]
                    Bloc = [rowB['x_location']+np.cos(headingB)*self.platformList[PFNum[0]].rFair, rowB['y_location']+np.sin(headingB)*self.platformList[PFNum[0]].rFair, self.platformList[PFNum[0]].zFair]
                    # get configuration for the line 
                    lineconfig = arrayMooring[j]['MooringConfigID']                       
                    # create mooring and connector dictionary for that line
                    m_config = getMoorings(lineconfig,Bnum)
                    # get letter number for mooring line
                    ind = len(self.platformList[PFNum[1]].getMoorings())
                    # create mooring class instance
                    mc = (Mooring(dd=m_config, id=str(PFNum[1])+'-'+str(PFNum[0])))
                    mc.shared = 1
                    mc.rA = Aloc
                    mc.rB = Bloc
                    # add mooring object to project mooring list
                    
                    
                    self.mooringList[str(PFNum[1])+'-'+str(PFNum[0])] = mc
                    # attach mooring object to platforms
                    mc.attachTo(self.platformList[PFNum[0]],end='B')
                    mc.attachTo(self.platformList[PFNum[1]],end='A')

                elif any(ids['ID'] == arrayMooring[j]['end A'] for ids in arrayAnchor): # end A is an anchor
                    # get ID of platform connected to line
                    PFNum.append(arrayMooring[j]['end B'])
                    for k in range(0,len(arrayInfo)):
                        if arrayInfo[k]['ID'] == PFNum[0]:
                            Bnum = k
                    # get configuration for that line 
                    lineconfig = arrayMooring[j]['MooringConfigID']                       
                    # create mooring and connector dictionary for that line
                    m_config = getMoorings(lineconfig,Bnum)
                    # get letter number for mooring line
                    ind = len(self.platformList[PFNum[0]].getMoorings())
                    # create mooring class instance
                    mc = (Mooring(dd=m_config, id=str(PFNum[0])+alph[ind]))
                    mc.rA = [m_config['span']+self.platformList[PFNum[0]].rFair,0,m_config['zAnchor']]
                    mc.rB = [self.platformList[PFNum[0]].rFair,0,self.platformList[PFNum[0]].zFair]
                    # adjust end positions based on platform location and mooring and platform heading
                    mc.reposition(r_center=self.platformList[PFNum[0]].r, heading=np.radians(arrayMooring[j]['headingB'])+self.platformList[PFNum[0]].phi, project=self)

                    # check if anchor instance already exists
                    if any(tt == 'shared_'+ arrayMooring[j]['end A'] for tt in self.anchorList): # anchor name exists already in list
                        # find anchor class instance
                        for anch in self.anchorList:#range(0,len(self.anchorList)):
                            if anch == 'shared_'+arrayMooring[j]['end A']:
                                mc.attachTo(self.anchorList[anch],end='A')

                    else:
                                               
                        # find location of anchor in arrayAnchor table
                        for k in range(0,len(arrayAnchor)):
                            if arrayAnchor[k]['ID'] == arrayMooring[j]['end A']:
                                aloc = [arrayAnchor[k]['x'],arrayAnchor[k]['y']] 
                                aNum.append(k) # get anchor row number
                                # set line anchor type and get dictionary of anchor information
                                lineAnch = arrayAnchor[k]['type']
                        ad = getAnchors(lineAnch,aNum=aNum[-1]) # call method to create dictionary
                        # adjust anchor location and rA based on location of anchor
                        zAnew, nAngle = self.getDepthAtLocation(aloc[0], aloc[1], return_n=True)
                        mc.rA = [aloc[0],aloc[1],-zAnew]
                        mc.dd['zAnchor'] = -zAnew
                        mc.z_anch = -zAnew
                        # create anchor object
                        self.anchorList['shared_'+arrayAnchor[k]['ID']] = Anchor(dd=ad, r=[aloc[0],aloc[1],-zAnew], aNum=aNum[-1],id='shared_'+arrayAnchor[k]['ID'])
                        # attach mooring object to anchor
                        mc.attachTo(self.anchorList[('shared_'+arrayAnchor[k]['ID'])],end='A')
                  
                    # add mooring object to project mooring list
                    self.mooringList[str(PFNum[0])+alph[ind]] = mc
                    # attach mooring object to platform
                    mc.attachTo(self.platformList[PFNum[0]],end='B')

                else: # error in input
                    raise Exception(f"end A input in array_mooring line_data table line '{j}' must be either an ID from the anchor_data table (to specify an anchor) or an ID from the array table (to specify a FOWT).")
                                           
                # add heading
                self.platformList[PFNum[0]].mooring_headings.append(np.radians(arrayMooring[j]['headingB']))
                if len(PFNum)>1: # if shared line
                    self.platformList[PFNum[1]].mooring_headings.append(np.radians(arrayMooring[j]['headingA'])) # add heading
                    
                # increment counter
                mct += 1
        
  
        
        # ===== load Cables ======
        def CableProps(cabType,checkType=1):
            '''
            Parameters
            ----------
            cabType : dictionary
                Dictionary of cable details from the cable_configs typeID
                Includes type (reference to name in cable_types or cable_props yaml)
            checkType : boolean
                Controls whether or not to first look for the cable type in the project yaml dictionary before
                attempting to get the cable properties from cable props yaml.

            Returns
            -------
            dd : design dictionary

            '''
            if 'typeID' in cabType and cabType['typeID'] in d['cable_types']:
                dd = d['cable_types'][cabType['typeID']]
                dd['name'] = cabType['typeID']
                if 'd_vol' in dd:
                    d_vol = dd['d_vol']
                else:
                    d_vol = dd['d']
                dd['w'] = (dd['m']-np.pi/4*d_vol**2*self.rho_water)*self.g
                if 'cableFamily' in cabType:
                    raise Exception('typeID and cableFamily listed in yaml - use typeID to reference a cable type in the cable_type section of the yaml and cableFamily to obtain cable properties from CableProps_default.yaml')
            elif 'cableFamily' in cabType:
                if not 'A' in cabType:
                    raise Exception('To use CableProps yaml, you must specify an area A for the cable family')
                cp = loadCableProps(None)
                cabProps = getCableProps(cabType['A'],cabType['cableFamily'],cableProps=cp)
                # fix units
                cabProps['power'] = cabProps['power']*1e6
                dd = cabProps
                dd['name'] = cabType['cableFamily']
            elif 'typeID' in cabType and not cabType['typeID'] in d['cable_types']:
                raise Exception(f'TypeID {cabType["typeID"]} provided in cable_config {cabType} is not found in cable_types section. Check for errors.')

            return(dd)
        
        def BuoyProps(buoyType):
            '''

            Parameters
            ----------
            buoyType : dict
                Dictionary of buoy details from the cable_configs sections list
            Returns
            -------
            dd : design dictionary

            '''
            dd = {}
            # add midpoint along length to add buoys to
            dd['L_mid'] = buoyType['L_mid']
            
            # figure out where buoy properties are located, add it in
            if buoyType['type'] in d['cable_appendages']:
                dd['module_props'] = d['cable_appendages'][buoyType['type']]
            else:
                bp = loadBuoyProps(None)
                buoyProps = getBuoyProps(buoyType['V'],buoyType['type'],buoyProps=bp)
                dd['module_props'] = buoyProps
            
            # add number of modules and spacing
            dd['N_modules'] = buoyType['N_modules']
            dd['spacing'] = buoyType['spacing']
            
            return(dd)
        
        def getCables(cabSection):
            cCondd = {}
            cC = cable_configs[cabSection['type']]
            
            if cabSection['type'] in d['cable_configs']:
                cCondd['span'] = cC['span']
                cCondd['length'] = cC['length']
                cCondd['A'] = getFromDict(cC,'A',default=0)
                cCondd['conductorSize'] = getFromDict(cC,'conductorSize',default=111)
                cCondd['type'] = cC['type']
                cCondd['powerRating'] = getFromDict(cC,'powerRating',default=0)
                if 'zJTube' in cC:
                    cCondd['zJTube'] = cC['zJTube']
                cCondd['voltage'] = getFromDict(cC,'voltage',default=66)
                
                # get cable properties for cable type (should only be one section - may change later)
                cCondd['cable_type'] = CableProps(cC)
                
                # get buoy properties
                if 'sections' in cC:
                    cCondd['buoyancy_sections'] = []
                    for i in range(0,len(cC['sections'])):
                        cCondd['buoyancy_sections'].append(BuoyProps(cC['sections'][i]))
                
                # # check for routing / burial info (generally for static cable)
                # if 'routing_x_y_r' in cC:
                #     cCondd['routing_xyr'] = cC['routing_x_y_r']
                # if 'burial' in cC:
                #     cCondd['burial'] = cC['burial']
                    
                # add depth
                cCondd['z_anch'] = -self.depth
                
            return(cCondd)
        
        
        # load in array cables
        if arrayCableInfo:
            cabLast = 0
            routing = []
            for i,cab in enumerate(arrayCableInfo):
                # create design dictionary for subsea cable
                dd = {'cables':[],'joints':[]}
                # get sections of the subsea cable 
                cable = arrayCableInfo[i]['CableID']
                # check if routing exists for this cable
                if arrayCableInfo[i]['route'].upper() == 'NA' or arrayCableInfo[i]['route'].upper() == 'NONE':
                    route = None
                else:
                    route = d['route_cables'][arrayCableInfo[i]['route']] #arrayCableInfo[i]['route']
                if cable in cableInfo:
                    if len(cableInfo[cable]['sections'])>1:
                        for j in range(0,len(cableInfo[cable]['sections'])):
                            cabSection = cableInfo[cable]['sections'][j]
                            # check for routing in this cable section
                            if route:
                                if route[0]['cable_configID'] in cabSection.values():
                                    # add routing
                                    routing = []
                                    # process any joints
                                    for f in range(0,len(route[0]['routing'])):
                                        if isinstance(route[0]['routing'][f],str):
                                            pass
                                        else:
                                            routing.append(route[0]['routing'][f])
                                else:
                                    routing = None
                                            
                            if cabLast: # last item was a cable, next should be a connector
                                if 'type' in cabSection:
                                    # no joint connecting 2 cables - add an empty joint to list
                                    dd['joints'].append({}) 
                                    # now get the sections of the cable configuration and put in dictionary
                                    cCondd = getCables(cabSection)
                                    
                                    if j == 0:
                                        # add heading for end A to this cable
                                        cCondd['headingA'] = np.radians(90-arrayCableInfo[i]['headingA'])
                                    elif j == len(cableInfo[cable]['sections']-1):
                                        # add heading for end B to this cable
                                        cCondd['headingB'] = np.radians(90-arrayCableInfo[i]['headingB'])
                                    dd['cables'].append(cCondd)
                                    if routing:
                                        cCondd['routing'] = routing
                                elif 'connectorType' in cabSection:
                                    dd['joints'].append(d['cable_joints'][cabSection['connectorType']])
                                    dd['joints'][-1]['type'] = cabSection['connectorType']

                                    cabLast = 0
                                else:
                                    # unsupported input
                                    raise Exception('Invalid section type keyword. Must be either type or connectorType')
                            else:
                                # last item was a connector
                                if 'type' in cabSection:
                                    cCondd = getCables(cabSection)
                                    if j == 0:
                                        # add heading for end A to this cable
                                        cCondd['headingA'] = np.radians(90-arrayCableInfo[i]['headingA'])
                                    elif j == len(cableInfo[cable]['sections'])-1:
                                        # add heading for end B to this cable
                                        cCondd['headingB'] = np.radians(90-arrayCableInfo[i]['headingB'])
                                    cabLast = 1
                                    dd['cables'].append(cCondd)
                                    if routing:
                                        cCondd['routing'] = routing
                                elif 'connectorType' in cabSection:
                                    raise Exception('Cannot have two connectors in a row')
                                else:
                                    # unsupported input
                                    raise Exception('Invalid section type keyword. Must be either type or connectorType')
                            
                    else:
                        # just a simple one line cable (no joints)
                        cabSection = cableInfo[cable]['sections'][0]
                        cCondd = getCables(cabSection)
                        cCondd['headingA'] = np.radians(90-arrayCableInfo[i]['headingA'])
                        cCondd['headingB'] = np.radians(90-arrayCableInfo[i]['headingB'])
                        if route:
                            if cabSection in route:
                                # add routing (no joints to process)
                                cCondd['routing'] = route[cabSection]
                        dd['cables'].append(cCondd)
                        
                else:
                    raise Exception(f'Cable configuration {cable} not found')

                            
                
                # create subsea cable object
                self.cableList[cable+str(i)] = Cable(cable+str(i),d=dd)
                
                        
                # connect cable to platform/substation
                if 'substation' in d and arrayCableInfo[i]['AttachA'] in d['substation']:
                    for j in range(0,len(arrayInfo)):
                        if arrayCableInfo[i]['AttachA'] == arrayInfo[j]['ID']:
                            raise Exception('Substation name must be different from platform ID')
                    self.cableList[cable+str(i)].attachTo(self.substationList[arrayCableInfo[i]['AttachA']],end='A')
                for j in range(0,len(arrayInfo)):
                    if arrayCableInfo[i]['AttachA'] == arrayInfo[j]['ID']:
                        # connect to platform
                        self.cableList[cable+str(i)].attachTo(self.platformList[arrayInfo[j]['ID']],end='A')
                    
                if 'substation' in d and arrayCableInfo[i]['AttachB'] in d['substation']:
                    for j in range(0,len(arrayInfo)):
                        if arrayCableInfo[i]['AttachB'] == arrayInfo[j]['ID']:
                            raise Exception('Substation name must be different from platform ID')
                    self.cableList[cable+str(i)].attachTo(self.substationList[arrayCableInfo[i]['AttachB']],end='B')
                for j in range(0,len(arrayInfo)):
                    if arrayCableInfo[i]['AttachB'] == arrayInfo[j]['ID']:
                        # connect to platform
                        self.cableList[cable+str(i)].attachTo(self.platformList[arrayInfo[j]['ID']],end='B')                            
                
                self.cableList[cable+str(i)].reposition()          
                
                # set joint positions
                for j,comp in enumerate(self.cableList[cable+str(i)].subcomponents):
                    if isinstance(comp,Joint):
                        if not 'r' in comp or comp['r'] is None:
                            jLocX,jLocY = self.cableList[cable+str(i)].estJointLoc(j)
                            depth = self.getDepthAtLocation(jLocX,jLocY)
                            comp['r']= [jLocX,jLocY,-depth]
                        # set rB of previous line and rA of next line to joint location
                        self.cableList[cable+str(i)].subcomponents[j-1].rB = comp['r']
                        self.cableList[cable+str(i)].subcomponents[j+1].rA = comp['r']
        
        
        # ===== load RAFT model parts =====
        # load info into RAFT dictionary and create RAFT model
        if raft:
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
            
            # load global RAFT settings into RAFT dictionary
            if 'RAFT_settings' in d['site'] and d['site']['RAFT_settings']:
                RAFTDict['settings'] = d['site']['RAFT_settings']
            # load RAFT cases into RAFT dictionary
            if 'RAFT_cases' in d['site'] and d['site']['RAFT_cases']:
                RAFTDict['cases'] = d['site']['RAFT_cases']
            
            # load array information into RAFT dictionary
            RAFTDict['array'] = deepcopy(d['array']) # need to change items so make a deepcopy
            # load general site info to RAFT dictionary
            RAFTDict['site'] = {'water_depth':self.depth,'rho_water':self.rho_water,'rho_air':self.rho_air,'mu_air':self.mu_air}
            RAFTDict['site']['shearExp'] = getFromDict(d['site']['general'],'shearExp',default=0.12)
            
            # create a name for the raft model
            RAFTDict['name'] = 'Project_Array'
            RAFTDict['type'] = 'input file for RAFT'
    
            # create RAFT model if necessary components exist
            if 'platforms' in RAFTDict or 'platform' in RAFTDict:
                if 'turbine' in RAFTDict or 'turbines' in RAFTDict:
                    self.getRAFT(RAFTDict,pristine=1)

        
        


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
        if 'bathymetry' in site and site['bathymetry']:
            if 'file' in site['bathymetry'] and site['bathymetry']['file']: # make sure there was a file provided even if the key is there
                self.loadBathymetry(site['bathymetry']['file'])
            elif 'x' in site['bathymetry'] and 'y' in site['bathymetry']:
                xs = np.array(site['bathymetry']['x'])
                ys = np.array(site['bathymetry']['y'])
                self.grid_depth = np.array(site['bathymetry']['depths'])
                self.setGrid(xs,ys)
            else:
                # assume a flat bathymetry
                self.grid_depth  = np.array([[self.depth]])
                
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
            if 'file' in site['seabed'] and site['seabed']['file']:
                self.loadSoil(file=site['seabed']['file'])
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
        if filename is not None:
            if filename[-3:]=='shp':
                raise ValueError("Geography-related operations not directly supported in Project class")
            
            elif filename[-3:]=='txt' or filename[-3:]=='dat':

                # load in the grid portion of the soil input file
                xs, ys, soil_names = sbt.readBathymetryFile(filename, dtype=str)  # read MoorDyn-style file

                self.soilProps = sbt.getSoilTypes(filename)

        elif filename is None:
            if yaml:
                xs = yaml['x']
                ys = yaml['y']
                soil_names = yaml['type_array']
                soilProps = yaml['soil_types']          
                
                # check that correct soil properties are being provided for the different soil types
                for soil in soilProps:
                    if 'rock' in soil:
                        if not 'UCS' in soilProps[soil] or not 'Em' in soilProps[soil]:
                            raise ValueError('Rock soil type requires UCS and Em values')
                    elif 'sand' in soil:
                        if not 'phi' in soilProps[soil] or not 'gamma' in soilProps[soil]:
                            raise ValueError('Sand soil type requires phi and gamma values')
                    elif 'clay' in soil:
                        if not 'Su0' in soilProps[soil] or not 'k' in soilProps[soil]:
                            raise ValueError('Clay soil type requires Su0 and k values')
                    elif 'mud' in soil:
                        if not 'Su0' in soilProps[soil] or not 'k' in soilProps[soil]:
                            raise ValueError('Mud soil type requires Su0 and k values')
                    else:
                        raise ValueError(f'Soil type {soil} not recognized. Soil type key must contain one of the following keywords: rock, sand, clay, mud')

                

                # make sure the soilProps dictionary has all the required information (should be updated later with exact properties based on anchor capacity functions)
                # soilProps = yaml['soil_types']
                for key,props in soilProps.items():
                    props['gamma'] = getFromDict(props, 'gamma', shape=-1, dtype=list, default=[4.7] , index=None)
                    props['Su0']   = getFromDict(props, 'Su0'  , shape=-1, dtype=list, default=[2.39], index=None)
                    props['k']     = getFromDict(props, 'k'    , shape=-1, dtype=list, default=[1.41], index=None)
                    props['alpha'] = getFromDict(props, 'alpha', shape=-1, dtype=list, default=[0.7] , index=None)
                    props['phi']   = getFromDict(props, 'phi'  , shape=-1, dtype=list, default=[0.0] , index=None)
                    props['UCS']   = getFromDict(props, 'UCS'  , shape=-1, dtype=list, default=[7.0] , index=None)
                    props['Em']    = getFromDict(props, 'Em'   , shape=-1, dtype=list, default=[50.0], index=None)
                
                self.soilProps = soilProps

        else:
            raise ValueError("Something is wrong")
        
        if xs is not None:
            self.soil_x = np.array(xs)
            self.soil_y = np.array(ys)
            self.soil_names = np.array(soil_names)
        
        
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
    
    def addCablesConnections(self,connDict,cableType='dynamic_cable_66',oss=False,substation_r=None,id_method='location',
                             keep_old_cables=False):
        '''Adds cables and connects them to existing platforms/substations based on info in connDict
        Designed to work with cable optimization output designed by Michael Biglu

        Parameters
        ----------
        connDict : dict
            Connection dictionary that describes the cables to create and their connections
        cableType : str, optional
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
            

        Returns
        -------
        None.

        '''
        
        
        # create substation object with id 
        if oss:
            if not substation_r:
                dd = {'r':[5000,1000]}
            else:
                print(substation_r)
                dd = {'r':substation_r}
            self.substationList[200] = Substation(dd,id=200)
            self.substationList[200].rFair = 58 ##### TEMPORARY #####
        
        # detach and delete existing cable list
        if keep_old_cables:
            lcab = len(self.cableList)
        else:
            if self.cableList:
                for j,cab in enumerate(self.cableList.values()):
                    cab.detachFrom('a')
                    cab.detachFrom('b')
            self.cableList = {} 
            lcab = 0
            
        # go through each index in the list and create a cable, connect to platforms
        for i in range(0,len(connDict)): # go through each cable
            # collect design dictionary info on cable
            dd = {}
            dd['cables'] = [{}]
            cd = dd['cables'][0]
            cd['span'] = connDict[i]['2Dlength']
            cd['length'] = connDict[i]['2Dlength']
            cd['A'] = connDict[i]['A_min_con']
            cd['type'] = 'dynamic'
            
            # add routing if necessary
            if len(connDict[i]['coordinates'])>2:
                cd['routing'] = []
                for coord in connDict[i]['coordinates'][1:-1]:
                    cd['routing'].append(coord)
            
            cp = loadCableProps(None)
            cabProps = getCableProps(connDict[i]['A_min_con'],cableType,cableProps=cp)
            # fix units
            cabProps['power'] = cabProps['power']*1e6
            cd['cable_type'] = cabProps
            dd['name'] = cableType

            
            if connDict[i]['cable_id']>=100 and not oss:
                # this is a substation
                pass
            else:        
                # create cable object
                self.cableList[cableType+str(i+lcab)] = Cable(cableType+str(i+lcab),d=dd)
                
                cab = self.cableList[cableType+str(i+lcab)]
                                
                # update upstream turbines property
                cab.upstream_turb_count = connDict['upstream_turb_count']
                
                # attach to platforms/substations
                for pf in self.platformList.values():
                    if id_method == 'location':
                        # find platform associated with ends
                        if np.allclose(pf.r,connDict[i]['coordinates'][0],atol=.01): 
                            attA = pf
                        elif np.allclose(pf.r,connDict[i]['coordinates'][1],atol=.01):
                            attB = pf
                    elif id_method == 'id':
                        # find platform associated with global id
                        if connDict[i]['turbineA_glob_id'] == pf.id:
                            attA = pf
                            # update platform location
                            pf.r = connDict[i]['coordinates'][0]
                        elif connDict[i]['turbineB_glob_id'] == pf.id:
                            attB = pf
                            # update platform location
                            pf.r = connDict[i]['coordinates'][1]
                        
                for substation in self.substationList.values():
                    if id_method == 'location':
                        # find substation associated with ends
                        if np.allclose(substation.r,connDict[i]['coordinates'][0],atol=.01):
                            attA = substation
                        elif np.allclose(substation.r,connDict[i]['coordinates'][1],atol=.01):
                            attB = substation
                    elif id_method == 'id':
                        # find substation associated with global id
                        if connDict[i]['turbineA_glob_id'] == substation.id:
                            attA = substation
                            # update oss location
                            substation.r = connDict[i]['coordinates'][0]
                        elif connDict[i]['turbineB_glob_id'] == substation.id:
                            attB = substation
                            # update oss location
                            substation.r = connDict[i]['coordinates']
                
                # attach cable
                cab.attachTo(attA,end='A')
                cab.attachTo(attB,end='B')
    
                # get heading of cable from attached object coordinates
                headingA = np.radians(90) - np.arctan2((attB.r[0]-attA.r[0]),(attB.r[1]-attA.r[1]))
                headingB = np.radians(90) - np.arctan2((attA.r[0]-attB.r[0]),(attA.r[1]-attB.r[1]))
                # print('headings: ',headingA,headingB)
    
                # reposition cable
                cab.reposition(headings=[headingA,headingB])                      
                
    
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
        figsize = kwargs.get('figsize', (6,4))  # the dimensions of the figure to be plotted
        
        
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
                
                contourf = ax.contourf(X, Y, self.grid_depth, num_levels, cmap='Blues', vmin=0, vmax=1000)
                #contourf.norm.autoscale([0,1])
                #contourf.set_clim(0, 1000)
            
                if not bare:  # Add colorbar with label
                    cbar = plt.colorbar(contourf, ax=ax, fraction=0.04, label='Water Depth (m)')
                    
        if plot_boundary:
            boundary = np.vstack([self.boundary, self.boundary[0,:]])
            ax.plot(boundary[:,0], boundary[:,1], 'b-.',label='Lease Boundary')
            
        
        # Seabed ground/soil type (to update)
        #X, Y = np.meshgrid(self.soil_x, self.soil_y)
        #ax.scatter(X, Y, c=self.soil_rocky, s=4, cmap='cividis_r', vmin=-0.5, vmax=1.5)
        # or if we have a grid of soil types, something like
        # pcolormesh([X, Y,] C, **kwargs)  wjere C is [x, y, 3 rgba]
        
        
        # Plot any object envelopes
        for platform in self.platformList.values():
            for name, env in platform.envelopes.items():
                ax.fill(env['x'], env['y'], edgecolor=[.5,0,0,.8], facecolor='none', linestyle='dashed', lw=0.8)
        
        for mooring in self.mooringList.values():
            for name, env in mooring.envelopes.items():
                #if 'shape' in env:  # if there's a shapely object
                #    pass  # do nothing for now...
                #elif 'x' in env and 'y' in env:  # otherwise just use coordinates
                ax.fill(env['x'], env['y'], color=[.6,.3,.3,.6])
        
        
        # Plot moorings one way or another (eventually might want to give Mooring a plot method)
        for mooring in self.mooringList.values():
        
            if mooring.ss:  # plot with Subsystem if available
                mooring.ss.drawLine2d(0, ax, color="k", endpoints=False, 
                                      Xuvec=[1,0,0], Yuvec=[0,1,0],label='Mooring Line')        
            else: # simple line plot
                ax.plot([mooring.rA[0], mooring.rB[0]], 
                        [mooring.rA[1], mooring.rB[1]], 'k', lw=0.5)
        
        # Plot cables one way or another (eventually might want to give Mooring a plot method)
        for cable in self.cableList.values():
        
            # simple line plot for now
            ax.plot([cable.subcomponents[0].rA[0], cable.subcomponents[0].rB[0]], 
                    [cable.subcomponents[0].rA[1], cable.subcomponents[0].rB[1]], 'r--', lw=0.5,label='Cable')
            # add in routing if it exists
            for sub in cable.subcomponents:
                if isinstance(sub,StaticCable):
                    if sub.coordinates:
                        # has routing  - first plot rA to sub.coordinate[0] connection
                        ax.plot([sub.rA[0],sub.coordinates[0][0]],
                                [sub.rA[1],sub.coordinates[0][1]],'r:',lw=0.6,label='Buried Cable')
                        # now plot route
                        for i in range(1,len(sub.coordinates)):
                            ax.plot([sub.coordinates[i-1][0],sub.coordinates[i][0]],
                                    [sub.coordinates[i-1][1],sub.coordinates[i][1]],'r:',lw=0.6,label='Buried Cable')
                        # finally plot sub.coordinates[-1] to rB connection
                        ax.plot([sub.coordinates[-1][0],sub.rB[0]],
                                [sub.coordinates[-1][1],sub.rB[1]],'r:',lw=0.6,label='Buried Cable')
            ax.plot([cable.subcomponents[-1].rA[0], cable.subcomponents[-1].rB[0]], 
                    [cable.subcomponents[-1].rA[1], cable.subcomponents[-1].rB[1]], 'r--', lw=0.5)
            
                # ax.plot([cable.subcomponents[0].rA[0], cable.subcomponents[-1].rB[0]], 
                #         [cable.subcomponents[0].rA[1], cable.subcomponents[0].rB[1]], 'r--', lw=0.5)
        
        # Plot platform one way or another (might want to give Platform a plot method)
        for platform in self.platformList.values():
            
            ax.plot(platform.r[0], platform.r[1], 'k*',label='Platform')
            
        for substation in self.substationList.values():
            ax.plot(substation.r[0],substation.r[1],'go',label='Substation')
            
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        if axis_equal:
            ax.axis('equal')
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))  # Removing duplicate labels
        ax.legend(by_label.values(), by_label.keys(),loc='upper center',bbox_to_anchor=(0.5, -0.1), fancybox=True, ncol=5)
        if save:
            plt.savefig('2dfarm.png', dpi=300, bbox_inches='tight')  # Adjust the dpi as needed
            
            # TODO - add ability to plot from RAFT FOWT
            
        
        

    def plot3d(self, ax=None, figsize=(10,8), fowt=False, save=False,
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
        xs = np.arange(min(self.grid_x),max(self.grid_x),50)
        ys = np.arange(min(self.grid_y),max(self.grid_y),50)
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
                           
        # or if we have a grid of soil types, something like
        ax.plot_surface(X, Y, h, rstride=1, cstride=1, facecolors = soil grid converted to colors <<<,
                       linewidth=0, antialiased=False)
                           
                           
        # apply colormap
        rc = cmap(norm(rocky))
        bath = ax.plot_surface(X, Y, -self.grid_depth, facecolors=rc, **args_bath)
        '''
        #################
        # from matplotlib import cm
        # args_bath = {'cmap':cm.GnBu_r}
        ####################
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
        ct = 0
        for mooring in self.mooringList.values():
            #mooring.subsystem.plot(ax = ax, draw_seabed=False)
            if mooring.ss:
                mooring.ss.drawLine(0,ax)
                
        for cable in self.cableList.values():
            for j,sub in enumerate(cable.subcomponents):
                if isinstance(sub,DynamicCable):
                    if sub.ss:
                        sub.ss.drawLine(0,ax)
                        
                elif isinstance(sub,StaticCable):
                    # add static cable routing if it exists
                    if sub.x:
                        burial = sub.burial
                        if 'NA' in burial:
                            # replace any NA with 0
                            for i,b in enumerate(burial):
                                if b == 'NA':
                                   burial[i] = 0 
                        # first plot from joint to start of cable route
                        jointA = cable.subcomponents[j-1].r
                        jointB = cable.subcomponents[j+1].r
                        soil_z = self.projectAlongSeabed(sub.x,sub.y)
                        ax.plot([jointA[0],sub.x[0]],[jointA[1],sub.y[0]],[-soil_z[0],-soil_z[0]-burial[0]],'k:',zorder=5,lw=1,alpha=0.7)
                        ax.plot([jointB[0],sub.x[-1]],[jointB[1],sub.y[-1]],[-soil_z[-1],-soil_z[-1]-burial[-1]],'k:',zorder=5,lw=1,alpha=0.7)
                        
                        
                        # plot in 3d along soil_z
                        ax.plot(sub.x,sub.y,-soil_z-burial,'k:',zorder=5,lw=1,alpha=0.7)
                                
            
                        
        
        # plot the FOWTs using a RAFT FOWT if one is passed in (TEMPORARY)
        if fowt:
            for pf in self.array.fowtlist:
                pf.plot(ax,zorder=20)
            # for i in range(self.nt):
            #     xy = self.turb_coords[i,:]
            #     fowt.setPosition([xy[0], xy[1], 0,0,0,0])
            #     fowt.plot(ax, zorder=20)
        
        # Show full depth range
        ax.set_zlim([-np.max(self.grid_depth), 0])

        set_axes_equal(ax)
        if not draw_axes:
            ax.axis('off')
        
        ax.view_init(20, -130)
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
            plt.savefig(output_filename, dpi=300, bbox_inches='tight')  # Adjust the dpi as needed
        
       
    def getMoorPyArray(self,bodyInfo=None,plt=0, pristineLines=True,cables=0):
        '''Creates an array in moorpy from the mooring, anchor, connector, and platform objects in the array.

        Parameters
        ----------
        bodyInfo : list of dictionaries, optional
            List of dictionaries (one list entry per body) that has information on hydrostatics for each body
        plt : boolean, optional
            Controls whether to create a plot of the MoorPy array. 1=create plot, 0=no plot The default is 0.

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
        self.ms = mp.System(depth=self.depth)
        # reset all anchor objects listed in self.anchorList to None in case there was a previous ms
        for i in self.anchorList:
            self.anchorList[i].mpAnchor = None
        for i in self.platformList:
            self.platformList[i].body = None
        
        wflag = 0 # warning flag has not yet been printed (prevent multiple printings of same hydrostatics warning)
        for i,body in enumerate(self.platformList): # make all the bodies up front - i is index in dictionary, body is key (name of platform)
            PF = self.platformList[body]
            # add a moorpy body at the correct location
            r6 = [PF.r[0],PF.r[1],0,0,0,0]
            # use bodyInfo dictionary to create moorpy body if given
            if bodyInfo:
                self.ms.addBody(-1,r6,m=bodyInfo[body]['m'],v=bodyInfo[body]['v'],rCG=np.array(bodyInfo[body]['rCG']),rM=np.array(bodyInfo[body]['rM']),AWP=bodyInfo[body]['AWP'])
            elif not bodyInfo and wflag == 0: # default to UMaine VolturnUS-S design hydrostatics info
                print('No hydrostatics information given, so default body hydrostatics from UMaine VolturnUS-S will be used.')
                wflag = 1
                self.ms.addBody(-1,r6,m=19911423.956678286,rCG=np.array([ 1.49820657e-15,  1.49820657e-15, -2.54122031e+00]),v=19480.104108645974,rM=np.array([2.24104273e-15, 1.49402849e-15, 1.19971829e+01]),AWP=446.69520543229874)
            else:
                self.ms.addBody(-1,r6,m=19911423.956678286,rCG=np.array([ 1.49820657e-15,  1.49820657e-15, -2.54122031e+00]),v=19480.104108645974,rM=np.array([2.24104273e-15, 1.49402849e-15, 1.19971829e+01]),AWP=446.69520543229874)
            PF.body = self.ms.bodyList[-1]
        # create anchor points and all mooring lines connected to the anchors (since all connected to anchors, can't be a shared mooring)
        for i in self.anchorList: # i is key (name) of anchor
            ssloc = []
            for j in self.anchorList[i].attachments: # j is key (name) of mooring object in anchor i
                # create subsystem
                if pristineLines:
                    self.anchorList[i].attachments[j]['obj'].createSubsystem(pristine=1)
                    # set location of subsystem for simpler coding
                    ssloc.append(self.anchorList[i].attachments[j]['obj'].ss)
                else:
                    self.anchorList[i].attachments[j]['obj'].createSubsystem()
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
                # add fairlead point
                self.ms.addPoint(1,ssloc[-1].rB)
                # add connector info for fairlead point
                self.ms.pointList[-1].m = self.ms.lineList[-1].pointList[-1].m 
                self.ms.pointList[-1].v = self.ms.lineList[-1].pointList[-1].v
                self.ms.pointList[-1].CdA = self.ms.lineList[-1].pointList[-1].CdA
                # attach the line to point
                self.ms.pointList[-1].attachLine(ssloc[-1].number,1)
                # find associated platform and attach body to point (since not a shared line, should only be one platform with this mooring object)
                for ii,k in enumerate(self.platformList): # ii is index in dictionary, k is key (name) of platform
                    if j in self.platformList[k].attachments: # j is key (name) of mooring object in anchor i checking if that same mooring object name is attached to platform k
                        PF = self.platformList[k] # platform object associated with mooring line j and anchor i
                        PFNum = ii # platform index
                # attach rB point to platform (need to subtract out location of platform from point for subsystem integration to work correctly)
                self.ms.bodyList[PFNum].attachPoint(len(self.ms.pointList),[ssloc[-1].rB[0]-PF.r[0],ssloc[-1].rB[1]-PF.r[1],ssloc[-1].rB[2]]) # attach to fairlead

        
        check = np.ones((len(self.mooringList),1))
        # now create and attach any shared lines
        for ii,i in enumerate(self.mooringList): # loop through all lines - ii is index of mooring object in dictionary, i is key (name) of mooring object
            for j in self.anchorList: # j is key (name) of anchor object
                if i in self.anchorList[j].attachments: # check if line has already been put in ms
                    check[ii] = 0     
            if check[ii] == 1: # mooring object not in any anchor lists
                # new shared line
                # create subsystem for shared line
                self.mooringList[i].createSubsystem(case=1,pristine=pristineLines) # we doubled all symmetric lines so any shared lines should be case 1
                # set location of subsystem for simpler coding
                if pristineLines:
                    ssloc = self.mooringList[i].ss
                else:
                    ssloc = self.mooringList[i].ss_mod
                # add subsystem as a line in moorpy system
                self.ms.lineList.append(ssloc)
                ssloc.number = len(self.ms.lineList)               
                
                # find associated platforms
                PF = [None,None]
                PFNum = [None,None]
                idx = []
                for ki,k in enumerate(self.platformList): # ki is index in dictionary, k is key (name) of platform
                    if self.platformList[k].isAttached(self.mooringList[i],end='a'):  #i in self.platformList[k].attachments:
                        PF[0] = self.platformList[k] # platform object
                        PFNum[0] = ki # platform index                    
                        # find key of mooring object in platform mooring list
                        idx.append(i)
                    elif self.platformList[k].isAttached(self.mooringList[i],end='b'):
                        PF[1] = self.platformList[k]
                        PFNum[1] = ki # platform index
                        # find key of mooring object in platform mooring list
                        idx.append(i)                          
                # add fairlead point A and attach the line to it
                self.ms.addPoint(1,ssloc.rA)
                self.ms.pointList[-1].attachLine(ssloc.number,0)
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
        
        # add in cables if desired
        if cables:
            # create a body for any substations cables are connected to 
            if self.substationList:
                for i,sub in enumerate(self.substationList.values()):
                    r6 = [sub.r[0],sub.r[1],0,0,0,0]
                    self.ms.addBody(1,r6,m=19911423.956678286,rCG=np.array([ 1.49820657e-15,  1.49820657e-15, -2.54122031e+00]),v=19480.104108645974,rM=np.array([2.24104273e-15, 1.49402849e-15, 1.19971829e+01]),AWP=446.69520543229874)
                    sub.body = self.ms.bodyList[-1]
            for i in self.cableList:
                # determine if suspended cable or not - having a static cable as a subcomponent means this is not a suspended cable
                for j,comp in enumerate(self.cableList[i].subcomponents):
                    # don't make a subsystem for a joint - make a point                   
                    if isinstance(comp,Joint):
                        if not comp.mpConn:
                            comp.makeMoorPyConnector(self.ms)
                    elif isinstance(comp,StaticCable):
                        # don't make a subsystem for a static cable (yet...)
                        pass
                    else:
                        # create subsystem for dynamic cable
                        comp.createSubsystem(pristine=pristineLines)
                        if pristineLines:                           
                            ssloc = comp.ss
                        else:
                            ssloc = comp.ss_mod
                        ssloc.number = len(self.ms.lineList)+1 
                        # add subsystem to line list
                        self.ms.lineList.append(ssloc)

                        
                        if j==0: # attach each end to correct bodies
                            self.ms.addPoint(1,ssloc.rA)
                            self.ms.pointList[-1].attachLine(ssloc.number,0)
                            body = comp.attached_to[0].body
                            body.attachPoint(len(self.ms.pointList),[ssloc.rA[0]-body.r6[0],ssloc.rA[1]-body.r6[1],ssloc.rA[2]])
                        else:
                            # connect to joint at end A
                            if not comp.attached_to[0].mpConn:
                                comp.attached_to[0].makeMoorPyConnector(self.ms)
                            joint = comp.attached_to[0].mpConn
                            joint.attachLine(ssloc.number,0)
                            
                        if j==len(self.cableList[i].subcomponents)-1: # last subcomponent could be first subcomponent
                            self.ms.addPoint(1,ssloc.rB)
                            self.ms.pointList[-1].attachLine(ssloc.number,1)
                            body = comp.attached_to[-1].body
                            body.attachPoint(len(self.ms.pointList),[ssloc.rB[0]-body.r6[0],ssloc.rB[1]-body.r6[1],ssloc.rB[2]])
                        else:
                            # connect to joint at end B
                            if not comp.attached_to[-1].mpConn:
                                comp.attached_to[-1].makeMoorPyConnector(self.ms)
                            joint = comp.attached_to[-1].mpConn
                            joint.attachLine(ssloc.number,1)                                   

        # initialize, solve equilibrium, and plot the system 
        self.ms.initialize()
        self.ms.solveEquilibrium(DOFtype='coupled')
        
        
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
            # find index for ID and mooringID
            for i in range(0,len(RAFTDict['array']['keys'])):
                if RAFTDict['array']['keys'][i] == 'ID':
                    IDindex = i
                elif RAFTDict['array']['keys'][i] =='mooringID':
                    mooringIDindex = i

            RAFTDict['array']['keys'].pop(IDindex) # remove key for ID because this doesn't exist in RAFT array table
            for i in range(0,len(RAFTDict['array']['data'])):
                RAFTDict['array']['data'][i][mooringIDindex] = 0 # make mooringID = 0 (mooring data will come from MoorPy)
                RAFTDict['array']['data'][i].pop(IDindex) # remove ID column because this doesn't exist in RAFT array data table

            # create raft model
            self.array = RAFT.Model(RAFTDict)
            # create dictionary of dictionaries of body hydrostatics for MoorPy bodies
            bodyInfo = {}
            for i,body in enumerate(self.array.fowtList):
                # set position (required before you can calcStatics)
                body.setPosition([RAFTable[i]['x_location'],RAFTable[i]['y_location'],0,0,0,0])
                # get body hydrostatics info for MoorPy bodies
                body.calcStatics()
                # populate dictionary of body info to send to moorpy
                bodyInfo[RAFTable[i]['ID']] = {'m':body.m,'rCG':body.rCG,'v':body.V,'rM':body.rM,'AWP':body.AWP}
            # create moorpy array if it doesn't exist
            if not self.ms:
                if self.cableList:
                    self.getMoorPyArray(bodyInfo,pristineLines=pristine,cables=1)
                else:
                    self.getMoorPyArray(bodyInfo,pristineLines=pristine)
            # assign moorpy array to RAFT object
            self.array.ms = self.ms
            # connect RAFT fowt to the correct moorpy body
            for i in range(0,len(self.platformList)): # do not include substations (these are made last)
                self.array.fowtList[i].body = self.ms.bodyList[i]
        else:
            raise Exception('Platform(s) must be specified in YAML file')
            
    def getMarineGrowth(self,mgDict_start=None,lines='all',tol=2):
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
        r = np.zeros([len(self.platformList),2])   # xy location of each platform
        r_no_off = np.zeros([len(self.platformList),2]) # xy location of each platform not considering offsets
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
            #print(north_bound)
        
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
                #breakpoint()
                self.platformList[pfIDs[ct]].setPosition([x,y],heading=phi[ct],degrees=True)
                self.platformList[pfIDs[ct]].rc = [i,j]
                # update cable lengths as needed (assumes the attachments are correct)
                cabs = self.platformList[pfIDs[ct]].getCables()
                for cab in cabs.values():
                    # determine what the cable is connected to
                    rr = [0,0]
                    for k,att in enumerate(cab.attached_to):
                        if isinstance(att,Platform):
                            rr[k] = att.r
                    newSpan = np.sqrt((rr[0][0]-rr[1][0])**2 + (rr[0][1]-rr[1][1])**2)
                    cab.updateSpan(newSpan)
                    # update cable end locations
                    cab.reposition()
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
        
        
    
    
    def updateFailureProbability(self):
        '''
        Function to populate (or update) failure probability dictionaries in each object 
        based on failure probability calculations developed by Emma Slack
        
        To be filled in...

        Returns
        -------
        None.

        '''
                
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

if __name__ == '__main__':

    
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