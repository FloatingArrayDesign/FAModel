
import numpy as np
import time
import yaml
import os
import re
from copy import deepcopy
from famodel.cables.cable_properties import getCableProps, getBuoyProps, loadCableProps,loadBuoyProps
import ruamel.yaml
import moorpy as mp
from moorpy.helpers import loadPointProps, getPointProps
import shapely as sh


def cart2pol(x, y):
    rho = np.hypot(x, y)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def m2nm(data):
    ''' Convert meters to nautical miles'''
    if isinstance(data,list):
        data = np.array(data)
    data = data*0.000539957
    return(data)


def printMat(mat):
    '''Prints a matrix to a format that is specified

    Parameters
    ----------
    mat : array
        Any matrix that is to be printed.

    Returns
    -------
    None.

    '''
    for i in range(mat.shape[0]):
        print( "\t".join(["{:+8.3e}"]*mat.shape[1]).format( *mat[i,:] ))
        
def printVec(vec):
    '''Prints a vector to a format that is specified

    Parameters
    ----------
    vec : array
        Any vector that is to be printed.

    Returns
    -------
    None.

    '''
    print( "\t".join(["{:+9.4e}"]*len(vec)).format( *vec ))



def unitVector(r):
    '''Returns the unit vector along the direction of input vector r.'''

    L = np.linalg.norm(r)

    return r/L


def loadYAML(filename):
    '''
    Loads a YAML file, allowing !include <filename> to include another yaml 
    in the file.

    Parameters
    ----------
    filename : str
        Filename (including path if needed) of main yaml file to load

    Returns
    -------
    info : dict
        Dictionary loaded from yaml

    '''
    
    with open(filename) as file:
        loader = yaml.FullLoader 
        loader.add_constructor('!include',yamlInclude)
        project = yaml.load(file, Loader=loader)
        if not project:
            raise Exception(f'File {file} does not exist or cannot be read. Please check filename.')
    
    return(project)

def yamlInclude(loader, node):
    '''
    Custom constructor that allows !include tag to include another yaml in
    the main yaml

    Parameters
    ----------
    loader : YAML loader object
    node : YAML node
        YAML node for include

    Returns
    -------
    None.

    '''
    # pull out f
    file_to_include = loader.construct_scalar(node)
    # pull out absolute path of file
    # if os.path.isabs(file_to_include):
    #     included_yaml = file_to_include
    # else:
    #     dir = 
    included_yaml = os.path.abspath(file_to_include)
    try:
        with open(included_yaml) as file:
            return(yaml.load(file,Loader=loader.__class__))
    except FileNotFoundError:
        raise FileNotFoundError(f"Included file {included_yaml} not found")
    except Exception as err:
        raise Exception(f"Error ocurred while loading included file {included_yaml}: {err}")
    

# ----- Cable routing support functions -----

def adjustCable(cc,project,na=None,nb=None,routeAdjustLength=500,rad_fair=None):
    cx = cc.subcomponents[2].x
    cy = cc.subcomponents[2].y
    if na==None:
        hA = np.radians(90) - np.arctan2((cc.subcomponents[0].rB[0]-cc.subcomponents[0].rA[0]),(cc.subcomponents[0].rB[1]-cc.subcomponents[0].rA[1]))
    else:
        hA = np.radians(na) #headingA
        cx[0] = [cc.attached_to[0].r[0]+500*np.cos(np.radians(na))]
        cy[0] = [cc.attached_to[0].r[1]+500*np.sin(np.radians(na))]
    if nb==None:
        hB = np.radians(90) - np.arctan2((cc.subcomponents[-1].rA[0]-cc.subcomponents[-1].rB[0]),(cc.subcomponents[-1].rA[1]-cc.subcomponents[-1].rB[1]))
    else:
        hB = np.radians(nb)
        cx[-1] = [cc.attached_to[1].r[0]+500*np.cos(np.radians(nb))]
        cy[-1] = [cc.attached_to[1].r[1]+500*np.sin(np.radians(nb))]
    cc.reposition(project=project,headings=[hA,hB],rad_fair=rad_fair)
    
# find cable(s) associated with specific platform(s) xy coordinates
def findCable(coords,project,atol=150):
    '''
    Find the cable(s) associated with specific platform(s) xy coordinates

    Parameters
    ----------
    coords : nested list or list
        xy coordinates of platforms(s) connected to this cable
    project : project object
        Associated project object
    atol : float, optional
        Absolute tolerance when looking for associated end locations. The default is 150.

    Returns
    -------
    corrcab : List or Cable object
        Returns a list of cable objects if only one platform coordinate provided
        Returns one cable object is 2 platform coordinates provided

    '''
    
    if isinstance(coords[0],list):
        allcabs = False
    else:
        allcabs = True
        corrcab = []

    for cab in project.cableList.values():
        if allcabs:
            # add any cable connected to this platform coordinate
            if np.allclose(cab.rA[0:2],coords,atol=atol) or np.allclose(cab.rB[0:2],coords,atol=atol):
                corrcab.append(cab)
        else:
            # add only the cable connected to both platform coordinates
            if np.allclose(cab.rA[0:2],coords[0],atol=atol) or np.allclose(cab.rB[0:2],coords[0],atol=atol):
                if np.allclose(cab.rA[0:2],coords[1],atol=atol) or np.allclose(cab.rB[0:2],coords[1],atol=atol):
                    corrcab = cab
    return(corrcab)


def check_headings(m_headings,c_heading,rad_buff):
    '''
    Check that cable heading is not near any mooring headings (with provided buffer zone angle)

    Parameters
    ----------
    m_headings : list
        Mooring headings [rad]
    c_heading : float
        Cable heading [rad]
    rad_buff : float
        Angular buffer zone [rad]

    Returns
    -------
    list of mooring headings that interfere with cable heading

    '''

    # convert negative headings to positive headings
    for i,mh in enumerate(m_headings):
        if mh<0:
            m_headings[i] = 2*np.pi + mh
        elif mh>2*np.pi:
            m_headings[i] = mh - 2*np.pi
    ang_diff = m_headings - c_heading
    inds_to_fix = np.where([round(abs(angd),8)<round(rad_buff,8) or round(abs(angd),8)>np.pi*2-round(rad_buff,8) for angd in ang_diff])[0]
    if len(inds_to_fix)>0:
        return([m_headings[ind] for ind in inds_to_fix])
    else:
        return([])
    
        
def head_adjust(att,heading,rad_buff=np.radians(30),endA_dir=1, adj_dir=1):
    '''
    function to adjust heading of cable based on angle buffer from mooring lines

    Parameters
    ----------
    att : list
        list of objects to attach to. 1 object if only concerned about the attached object associated with that side
    heading : float
        Cable compass heading at attachment to att in radians
    rad_buff : float
        Buffer angle in radians
    endA_dir : float, optional
        Either 1 or -1, controls sign of new heading for end B. Only altered to -1 if dynamic
        cable from end A will get close to end B moorings. Default is 1.
    adj_dir : float, optional
        Either 1 or -1, default is 1. If -1, adjusts direction heading is altered 
        to avoid mooring lines, can be used if that heading direction is more natural.
        This is a manual input to the main function adjusting cables.

    Returns
    -------
    headnew : float
        New cable heading

    '''
    if heading<0:
        headnew = np.pi*2 + heading
    else:
        headnew = heading
    attheadings = [] # complete list of mooring headings to avoid, from all platforms
    flipheads = False # whether to flip headings ( for if you are looking at mooring headings of platform on the other end)
    for at in att:
        mhs = np.radians([m.heading for m in at.getMoorings().values()])
        if flipheads:
            atmh = np.array(mhs) + np.pi
            for j,a in enumerate(atmh):
                # keep everything under 2pi angle
                if a>2*np.pi:
                    atmh[j] = a-2*np.pi
        else:
            atmh = np.array(mhs) #attached platform mooring headings array
        #attheadings.extend(atmh)
        attheadings.extend(atmh) # keep in compass heading
        flipheads = True

    interfere_h = check_headings(attheadings,headnew,rad_buff)
    # if the headings interfere, adjust them by angle buffer
    for mhead in interfere_h:
        ang_diff_dir = np.sign(headnew - mhead) if headnew != mhead else 1
        headnew = mhead - adj_dir*rad_buff*endA_dir*ang_diff_dir #headnew + np.sign(ang_diff)*(rad_buff - abs(ang_diff))*endA_dir
        interfere_hi = check_headings(attheadings,headnew,rad_buff)
        for i in interfere_hi:
            # try rotating other way
            headnew = mhead + rad_buff*endA_dir*ang_diff_dir
            # re-check offsets
            interfere_hij = check_headings(attheadings,headnew,rad_buff)
            if not interfere_hij:
                return(headnew)
            else:
                # cut buffer in half and try again
                newbuff = rad_buff/2
                headnew = mhead + newbuff*endA_dir*ang_diff_dir
                return(headnew)

    return(headnew)

def getCableDD(dd,selected_cable,cableConfig,cableType_def,connVal):
    '''
    get cable design dictionary from a cableConfig yaml. Primarily used for project.addCablesConnections()

    Parameters
    ----------
    dd : dict
        design dictionary of cable (initial information)
    selected_cable : dict
        Dictionary of information on cable configuration in cableConfig yaml format
    cableConfig : dict
        Dictionary of all cable configurations in cableConfig yaml format
    cableType_def : str
        Name of cable type referring to a cable ID key in cableProps yaml
    connVal : dict
        Dictionary of cable connection information for associated cable in layout.py iac_dic format

    Returns
    -------
    dd : dict
        Fully developed design dictionary

    '''
    
    # set up selected cable design dictionary
    if len(selected_cable['sections'])> 1:
        dd['joints'] = []
    
    # get connector and joint costs if they were given
    dd['connector_cost'] = getFromDict(selected_cable,'connector_cost',default=0)
    joint_cost = getFromDict(selected_cable,'joint_cost',default=0)
    
    for j in range(len(selected_cable['sections'])):
        dd['cables'].append(deepcopy(cableConfig['cableTypes'][selected_cable['sections'][j]]))
        cd = dd['cables'][j]
        cd['z_anch'] = -selected_cable['depth']
        # cd['cable_type'] = cableConfig['cableTypes'][selected_cable['sections'][j]] # assign info in selected cable section dict to cd
        cd['A'] = selected_cable['A']
        cd['voltage'] = cableType_def[-2:]

        
        # add joints as needed (empty for now)
        if j < len(selected_cable['sections'])-1:
            dd['joints'].append({'cost':joint_cost}) # default 0

        # add routing if necessary
        if dd['cables'][j]['type']=='static':
            cd['routing'] = []
            # if len(connDict[i]['coordinates'])>2:
            #     for coord in connDict[i]['coordinates'][1:-1]:
            #         cd['routing'].append(coord)
            cableType = 'static_cable_'+cableType_def[-2:]
        else:
            cableType = 'dynamic_cable_'+cableType_def[-2:]
            cd['rJTube'] = 5
            
        
        if not 'cable_type' in cd or not cd['cable_type']:
            cp = loadCableProps(None)
            cabProps = getCableProps(connVal['conductor_area'],cableType,cableProps=cp)
            # fix units
            cabProps['power'] = cabProps['power']*1e6
            cd['cable_type'] = cabProps

        cd['cable_type']['name'] = selected_cable['sections'][j]
        
    return(dd)

def getCableDesign(connVal, cableType_def, cableConfig, configType, depth=None):
    # go through each index in the list and create a cable, connect to platforms
    
    dd = {}
    dd['cables'] = []
    # collect design dictionary info on cable

    # create reference cables (these are not saved into the cableList, just used for reference)
    
    # find associated cable in cableConfig dict
    cableAs = []
    cableDs = []
    cable_selection = []
    for cabC in cableConfig['configs']:
        if connVal['conductor_area'] == cabC['A']:
            cableAs.append(cabC)
    if not cableAs:
        raise Exception('Cable configs provided do not match required conductor area')
    elif len(cableAs) == 1:
        cable_selection = cableAs
    else:                        
        for cabA in cableAs:                           
            # only check distance if the cable is NOT connected to substation
            if 'dist' in cabA and connVal['cable_id']<100:
                if abs(connVal['2Dlength'] - cabA['dist']) < 0.1:
                    cableDs.append(cabA)    
        
        #if there's no matching distance, assume the nonsuspended cables 
        if cableDs == []:
            for cabA in cableAs:
                if cabA['type'] == 0:
                    cableDs.append(cabA)
        
        
        for cabD in cableDs:
            if connVal['cable_id']>=100 and cabD['type']==0:
                # connected to a substation, use a dynamic-static-dynamic configuration
                cable_selection.append(cabD)
                
            elif connVal['cable_id']<100 and cabD['type']==configType:
                # not connected to substation, use default config type
                cable_selection.append(cabD)

        # if no cables are found to match, override the configType

        if cable_selection == []:
            for cabD in cableDs:
                if connVal['cable_id']<100:
                    cable_selection.append(cabD)
            
    if len(cable_selection)> 1:
        # downselect by depth
        depthdiff = np.array([x['depth']-depth for x in cable_selection])
        selected_cable = cable_selection[np.argmin(depthdiff)]
        # else:
        #     raise Exception(f"Multiple cables match selection criteria for cable {connDict[i]['cable_id']}")
    elif len(cable_selection) == 1:
        # found the correct cable
        selected_cable = cable_selection[0]

    else:
        raise Exception(f"No cable matching the selection criteria found for cable {connVal['cable_id']}")
        
    dd = getCableDD(dd,selected_cable,cableConfig,cableType_def,connVal)           
    dd['name'] = cableType_def
        
    return(selected_cable,deepcopy(dd))

def getDynamicCables(cable_config, cable_types, cable_appendages, depth, 
                     rho_water=1025, g=9.81):
    '''
    Create cable design dictionary

    Parameters
    ----------
    cable_config : dict
        Dictionary of dynamic cable configuration information in ontology yaml format
    cable_types : dict
        Dictionary of cable type information in ontology yaml format
    cable_appendages : dict
        Dictionary of cable appendage (ex: buoyancy modules) information in ontology yaml format
    depth : float
        Water depth
    rho_water : float, optional
        Water density [kg/m^3]. Default is 1025
    g : float, optional
        acceleration due to gravity [m/s^2]. Default is 9.81

    Returns
    -------
    cCondd : dict
        Dictionary of cable design
    jCondd : dict
        Dictionary of joint design

    '''
    
    cCondd = {'appendages':[]}
    jCondd = {}

    if cable_config:
        cCondd['span'] = cable_config['span']
        cCondd['L'] = cable_config['length']
        cCondd['A'] = getFromDict(cable_config,'A',default=0)
        if 'zJTube' in cable_config:
            cCondd['zJTube'] = cable_config['zJTube']
        # cCondd['voltage'] = getFromDict(cable_config,'voltage',default=66)
        
        # get cable properties for cable type (should only be one section - may change later)
        cCondd['cable_type'] = CableProps(cable_config['cable_type'], cable_types, rho_water, g, A=cCondd['A'])
        cCondd['type'] = 'dynamic'
        
        # get appendage properties (could be buoys or could be J-tubes, joints, etc)
        if 'sections' in cable_config:
            cCondd['buoyancy_sections'] = []
            for i in range(0,len(cable_config['sections'])):
                dd, appEntity = AppendageProps(cable_config['sections'][i],cable_appendages)
                if 'BUOY' in appEntity.upper():
                    cCondd['buoyancy_sections'].append(dd)
                elif 'JOINT' in appEntity.upper():
                    jCondd = dd
                else:
                    cCondd['appendages'].append(dd)
            
        # add depth
        cCondd['z_anch'] = -depth
        
    return(cCondd, jCondd)

def getStaticCables(statCabID, cable_types, routing=None, burial=None, rho_water=1025,
                    g=9.81, A=None):
    '''
    Creates design dictionary of a static cable including any routing and burial information

    Parameters
    ----------
    statCabID : str
        Name of static cable type.
    cable_types : dict
        Dictionary of cable types.
    routing : list, optional
        List of x-y-r routing coordinates (where r is radius). The default is None.
    burial : dict, optional
        Dictionary containing lists of depths and lengths along cable. The default is None.
    rho_water : float, optional
        Water density [kg/m^3]. Default is 1025
    g : float, optional
        acceleration due to gravity [m/s^2]. Default is 9.81
    A : int, optional
        conductor area in mm^2

    Returns
    -------
    dd : dict
        Design dictionary of static cable

    '''
    dd = {}
    dd['cable_type'] = CableProps(statCabID, cable_types, rho_water, g, A=A)
    dd['type'] = 'static'
    
    # check for routing / burial info
    if routing:
        dd['routing_xyr'] = routing
    if burial:
        dd['burial'] = burial
        
    return(dd)

def AppendageProps(appType,cable_appendages):
    '''
    Create appendage or buoyancy_sections portion of cable design dictionary

    Parameters
    ----------
    appType : dict
        Dictionary of appendage details from the dyn_cable_configs sections list
    Returns
    -------
    dd : buoyancy_sections portion of design dictionary (if buoy) or creates Joint design dictionary (if joint)
         or adds to appendages section of design dictionary (if other i.e. J-tube)
    entity : str
        type of appendage

    '''
    dd = {}
    if appType['type'] in cable_appendages:
        # pull out what type of appendage this is 
        entity = cable_appendages[appType['type']]['type']
        
        if 'BUOY' in entity.upper():  # this is a buoyancy module
            # add midpoint along length to add buoys to
            dd['L_mid'] = appType['L_mid']
            
            # add buoy props to dd
            dd['module_props'] = cable_appendages[appType['type']]
            # add number of modules and spacing
            dd['N_modules'] = appType['N_modules']
            dd['spacing'] = appType['spacing']
        else:
            # create dd
            dd = cable_appendages[appType['type']]
            if 'L_mid' in appType:
                dd['L_mid'] = appType['L_mid']
    elif 'V' in appType:
        
        # pull from buoy props
        bp = loadBuoyProps(None)
        buoyProps = getBuoyProps(appType['V'],appType['type'],buoyProps=bp)
        
        # assign buoyancy section info
        dd['L_mid'] = appType['L_mid']
        dd['module_props'] = buoyProps
        dd['N_modules'] = appType['N_modules']
        dd['spacing'] = appType['spacing']
        
        entity = 'buoy'
        
    else:
        raise Exception(f"appendage {appType['type']} is not found in cable_appendages dictionary. If appendage design should come from buoyProps, please provide a volume 'V'.")
        
        
    
    return(deepcopy(dd),entity)

def CableProps(cabType, cable_types, rho_water, g, checkType=1, A=None):
    '''
    Create cable_type section of design dictionary
    Parameters
    ----------
    cabType : str
        Name of cable type
    cable_types : dict
        Dictionary of cable type details in ontology yaml format
    rho_water : float
        Density of water [kg/m^3]
    g : float
        Acceleration due to gravity [m/s^2]
    checkType : boolean
        Controls whether or not to first look for the cable type in the project yaml dictionary before
        attempting to get the cable properties from cable props yaml.
    A : int
        Cable conductor area in mm^2, only needed if using cableprops yaml

    Returns
    -------
    dd : cable type section of cable design dictionary

    '''
    if cabType in cable_types:
        dd = cable_types[cabType]
        dd['name'] = cabType
        if 'd_vol' in dd:
            d_vol = dd['d_vol']
        else:
            d_vol = dd['d']
        dd['w'] = (dd['m']-np.pi/4*d_vol**2*rho_water)*g
        # if 'cableFamily' in cabType:
        #     raise Exception('typeID and cableFamily listed in yaml - use typeID to reference a cable type in the cable_type section of the yaml and cableFamily to obtain cable properties from CableProps_default.yaml')
    else:
        # cable type not listed in cable_types dictionary, default to using getCableProps
        if not A:
            raise Exception('To use CableProps yaml, you must specify an area A for the cable family')
        cp = loadCableProps(None)
        cabProps = getCableProps(A,cabType,cableProps=cp)
        # fix units
        cabProps['power'] = cabProps['power']*1e6
        dd = cabProps
        dd['name'] = cabType
        dd['voltage'] = cabProps['voltage']
    # elif 'typeID' in cabType and not cabType['typeID'] in cable_types:
    #     raise Exception(f'TypeID {cabType["typeID"]} provided in cable_config {cabType} is not found in cable_types section. Check for errors.')

    return(deepcopy(dd))

def MooringProps(mCon, lineTypes, rho_water, g, checkType=1):
    '''
    Parameters
    ----------
    mCon : dict
        Dictionary of mooring details from the mooring_line_configs key in ontology yaml
        Includes type (reference to name in mooring_line_types or cable_props yaml)
    lineTypes : dict
        Dictionary of mooring line material type information in ontology yaml format
    rho_water : float
        Water density [kg/m^3]
    g : float
        Acceleration due to gravity [m/s^2]
    checkType : boolean
        Controls whether or not to first look for the cable type in the project yaml dictionary before
        attempting to get the cable properties from cable props yaml.

    Returns
    -------
    dd : design dictionary

    '''
    if 'type' in mCon and mCon['type'] in lineTypes:
        dd = lineTypes[mCon['type']]
        dd['name'] = mCon['type']
        if 'd_vol' in dd:
            d_vol = dd['d_vol']
        # else:
        #     d_vol = dd['d']
        dd['w'] = (dd['m']-np.pi/4*d_vol**2*rho_water)*g
        dd['MBL'] = float(dd['MBL'])
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
        dd['d_nom'] = mCon['d_nom']
        dd['MBL'] = float(dd['MBL'])
    elif 'type' in mCon and not mCon['type'] in lineTypes:
        raise Exception(f'Type {mCon["type"]} provided in mooring_line_config {mCon} is not found in mooring_line_types section. Check for errors.')

    return(deepcopy(dd))

def getMoorings(lcID, lineConfigs, connectorTypes, pfID, proj):
    '''

    Parameters
    ----------
    lcID : str
        ID of line configuration in ontology yaml
    lineConfigs : dict
        Dictionary of line configuration types described in ontology yaml format
    connectorTypes : dict
        Dictionary of connector types described in ontology yaml format
    pfID : str
        Platform object ID connected to this mooring object
    proj : project class instance
        Project object associated with this mooring line

    Returns
    -------
    dd : dict
        mooring design dictionary

    '''
    # set up dictionary of information on the mooring configurations
    dd = {'span':{},'zAnchor':{}}#,'EndPositions':{}}
    # set up connector dictionary
    c_config = []
    config = [] # mooring and connector combined configuation list
                
    lineLast = 1    # boolean whether item with index k-1 is a line. Set to 1 for first run through of for loop
    ct = 0   # counter for number of line types
    nsec = len(lineConfigs[lcID]['sections']) # number of sections
    for k in range(0,nsec): # loop through each section in the line
    
        lc = lineConfigs[lcID]['sections'][k] # set location for code clarity later
        # determine if it's a line type or a connector listed
        if 'type' in lc or 'mooringFamily' in lc: 
            # this is a line
            if lineLast: # previous item in list was a line (or this is the first item in a list)
                # no connector was specified for before this line - add an empty connector
                config.append({})
                c_config.append({})                        
            # set line information
            lt = MooringProps(lc, proj.lineTypes, proj.rho_water, proj.g)                                             
            # lt = self.lineTypes[lc['type']] # set location for code clarity and brevity later
            # set up sub-dictionaries that will contain info on the line type
            config.append({'type':lt})# {'name':str(ct)+'_'+lc['type'],'d_nom':lt['d_nom'],'material':lt['material'],'d_vol':lt['d_vol'],'m':lt['m'],'EA':float(lt['EA'])}})
            config[-1]['type']['name'] = str(ct)+'_'+str(lt['name'])
            # make EA a float not a string
            config[-1]['type']['EA'] = float(lt['EA'])  
            config[-1]['type']['MBL'] = float(lt['MBL'])
            # set line length
            config[-1]['L'] = lc['length']

            # update line last boolean
            lineLast = 1
            
        elif 'connectorType' in lc:
            cID = lc['connectorType']
            # this is a connector
            if lineLast == 0:
                # last item in list was a connector
                raise Exception(f"Two connectors were specified in a row for line configuration '{lcID}', please remove one of the connectors")
            else:
                
                # last item in list was a line
                if cID in connectorTypes:
                    config.append(connectorTypes[cID]) # add connector to list
                    config[-1]['type'] = cID
                else:
                    # try pointProps
                    try:
                        props = loadPointProps(None)
                        design = {f"num_c_{cID}":1}
                        config.append(getPointProps(design, Props=props))

                    except Exception as e: 
                        raise Exception(f"Connector type {cID} not found in connector_types dictionary, and getPointProps raised the following exception:",e)
                        
                # update lineLast boolean
                lineLast = 0
        elif 'subsections' in lc:
            # TODO: LHS: ERROR CHECKING FOR ORDER OF COMPONENTS PROVIDED WITHIN SUBSECTIONS, ADD IN NEEDED CONNECTORS!!
            
            if lineLast and k != 0:
                # if this is not the first section AND last section was a line, add a empty connector first
                config.append({})
                lineLast = 0
            config.append([])
            sublineLast = [lineLast]*len(lc['subsections']) # to check if there was a connector provided before this
            for ii,sub in enumerate(lc['subsections']):
                config[-1].append([])
                for jj,subsub in enumerate(sub):
                    if 'connectorType' in subsub and sublineLast[ii]:
                        cID = subsub['connectorType']
                        if cID in connectorTypes:
                            cID = subsub['connectorType']
                            config[-1][-1].append(connectorTypes[cID])
                        else:
                            # try pointProps
                            try:
                                props = loadPointProps(None)
                                design = {f"num_c_{cID}":1}
                                config[-1][-1].append(getPointProps(design, Props=props))
                            except Exception as e: 
                                raise Exception(f"Connector type {cID} not found in connector_types dictionary, and getPointProps raised the following exception:",e)
                        sublineLast[ii] = 0
                    elif 'connectorType' in subsub and not sublineLast[ii]:
                        raise Exception('Previous section had a connector, two connectors cannot be listed in a row')
                    elif 'type' or 'mooringFamily' in subsub:
                        if sublineLast[ii]:
                            # add empty connector
                            config[-1][-1].append({})
                        lt = MooringProps(subsub,proj.lineTypes, proj.rho_water, proj.g)
                        config[-1][-1].append({'type':lt,
                                               'L': subsub['length']})
                        # make EA a float not a string
                        config[-1][-1][-1]['type']['EA'] = float(lt['EA']) 
                        config[-1][-1][-1]['type']['MBL'] = float(lt['MBL'])
                        sublineLast[ii] = 1
                    else:
                        raise Exception(f"keys in subsection line definitions must either be 'type', 'mooringFamily', or 'connectorType'")
                    # if this is the last section and the last part of the subsection in the section, it needs to end on a connector
                    # so, add a connector if last part of subsection was a line!
                    if sublineLast[ii] and k==nsec-1 and jj==len(sub)-1:
                        # end bridle needs connectors added 
                        config[-1][-1].append({})
                        sublineLast[ii] = 0
                        
            lineLast = sublineLast[-1] # TODO: LHS: think how to handle this situation for error checking...
        else:
            # not a connector or a line
            raise Exception(f"Please make sure that all section entries for line configuration '{lcID}' are either line sections (which must have a 'type' key), connectors (which must have a 'connectorType' key, or subsections")

    # check if line is a shared symmetrical configuration
    if 'symmetric' in lineConfigs[lcID] and lineConfigs[lcID]['symmetric']:
        if not lineLast: # check if last item in line config list was a connector
            for ii in range(len()):
                # set mooring configuration 
                config.append(config[-1-2*ii])
                # set connector (since it's mirrored, connector B becomes connector A)
                config.append(config[-2-2*ii])
        else: # double the length of the end line
            config[-1]['L'] =config[-1]['L']*2
            # set connector B for line same as previous listed connector
            config.append(config[-1])
            for ii in range(0,ct-1): # go through every line config except the last (since it was doubled already)
                # set mooring configuration
                config.append(config[-2-2*ii])
                # set connector
                config.append(config[-3-2*ii])
    else: # if not a symmetric line, check if last item was a line (if so need to add another empty connector)
        if lineLast:
            # add an empty connector object
            config.append({})
    # set general information on the whole line (not just a section/line type)
    # set to general depth first (will adjust to depth at anchor location after repositioning finds new anchor location)
    dd['subcomponents'] = config
    dd['zAnchor'] = -proj.depth 
    dd['span'] = lineConfigs[lcID]['span']
    dd['name'] = lcID
    # add fairlead radius and depth to dictionary
    dd['rad_fair'] = proj.platformList[pfID].rFair
    dd['z_fair'] = proj.platformList[pfID].zFair
    
    
    return(dd) #, c_config)

    

def getConnectors(c_config, mName, proj):
    '''

    Parameters
    ----------
    c_config : dict
        Dictionary of connector configurations for a mooring line.
    mName : tuple
        Key name in the mooringList dictionary for the associated mooring object
    proj : project class instance
        project object to create connectors for

    Returns
    -------
    None.

    '''
    from famodel.mooring.connector import Connector
    
    # make connector objects for all sections of a mooring line configuration in order
    for i in range(0,len(c_config)):
        # check if connector is a none-type
        if c_config[i] == None:                   
            # create empty connector object
            proj.mooringList[mName].dd['connectors'].append(Connector())
        elif c_config[i]:
            # create connector object with c_config entries
            proj.mooringList[mName].dd['connectors'].append(Connector(**c_config[i]))

def getAnchors(lineAnch, arrayAnchor, proj):
    '''Create anchor design dictionary based on a given anchor type

    Parameters
    ----------
    lineAnch : string
        anchor type, for reference in the 'anchor_types' dictionary
    arrayAnchor : list
        list of anchors listed in array_mooring anchor_data table of ontology yaml.
    proj : project class instance
        project object to develop anchors for

    Returns
    -------
    ad : dict
        anchor design dictionary

    '''
    ad = {'design':{}, 'cost':{}} 
    ad['design'] = deepcopy(proj.anchorTypes[lineAnch])
    if 'mass' in proj.anchorTypes[lineAnch]:
        mass = ad['design'].pop('mass')
    else:
        mass = 0
    ad['type'] = proj.anchorTypes[lineAnch]['type']
    ad['name'] = lineAnch
    
    return(ad, mass)

def attachFairleads(moor, end, platform, fair_ID_start=None, fair_ID=None, fair_inds=None):
    '''
    helper function for loading, attaches fairleads to mooring objects
    and runs some error checks

    Parameters
    ----------
    fair_inds : int/list
        Fairlead index/indices to attach to mooring line
    moor : Mooring class instance
        Mooring that will attach to fairlead(s)
    end : int or str
        must be in [0,a,A] for end A or [1,b,B] for end B
    platform : Platform class instance
        Platform that is associated with the fairlead
    fair_ID_start : str, optional
        start of fairlead ID, the index will be appended to this. Not needed if fair_ID provided
    fair_ID : list, optional
        fairlead ID list for each fairlead. If fair_ID_start is not provided, fair_ID must be provided


    Returns
    -------
    None.

    '''
    # convert to list if needed
    if fair_inds is not None :
        if not isinstance(fair_inds,list):
            fair_inds = list([fair_inds])
        # check lengths are the same
        if not len(moor.subcons_B)==len(fair_inds):
            raise Exception(f'Number of fairleads must equal number of parallel sections at end {end}')
    elif fair_ID is not None: 
        if not isinstance(fair_ID, list):
            fair_ID = list([fair_ID])
        # check lengths are the same
        if not len(moor.subcons_B)==len(fair_ID):
            raise Exception(f'Number of fairleads must equal number of parallel sections at end {end}')
    else:
        raise Exception('Either fairlead indices or fairlead IDs must be provided')
    # grab correct end
    end_subcons = moor.subcons_B if end in [1,'b','B'] else moor.subcons_A
    
    # put together fairlead ids as needed
    if fair_ID_start != None and fair_inds != None:
        fair_ID = []
        for i in fair_inds:
            fair_ID.append(fair_ID_start+str(i))
    # attach each fairlead to the end subcomponent
    fairs = []
    for ii,con in enumerate(end_subcons):
        fairs.append(platform.attachments[fair_ID[ii]]['obj'])
        end_subcons[ii].join(fairs[-1])
        
    return(fairs)
        
def calc_heading(pointA, pointB):
    '''calculate a compass heading from points, if pointA or pointB is a list of points,
       the average of those points will be used for that end'''
    # calculate the midpoint of the point(s) on each end first
    pointAmid = calc_midpoint(pointA) 
    pointBmid = calc_midpoint(pointB)
    dists = np.array(pointAmid) - np.array(pointBmid)
    headingB = np.pi/2 - np.arctan2(dists[1], dists[0])
    
    return(headingB)

def calc_midpoint(point):
    '''Calculates the midpoint of a list of points'''
    if isinstance(point[0],list) or isinstance(point[0],np.ndarray):
        pointx = sum([x[0] for x in point])/len(point)
        pointy = sum([x[1] for x in point])/len(point)
        # add z component if needed
        if len(point[0])==3:
            pointz = sum([x[2] for x in point])/len(point)
            return([pointx,pointy,pointz])
    else:
        pointx = point[0]
        pointy = point[1]
        # add z component if needed
        if len(point)==3:
            pointz = point[2]
            return([pointx,pointy,pointz])
        
    return([pointx,pointy])
    

def route_around_anchors(proj, anchor=True, cable=True, padding=50):
    '''check if static cables hit anchor buffer, if so reroute cables around anchors'''
    # make anchor buffers with 50m radius
    if anchor:
        anchor_buffs = []
        for anch in proj.anchorList.values():
            anchor_buffs.append(anch.makeBuffer())
      
    # make cable linestrings including joint locs and static cable (no dynamic cable for simplicity)
    if cable:
        cable_line = {}
        for name, cab in proj.cableList.items():
            cable_line[name] = cab.makeLine(include_dc=False)

    # Function to calculate angle of a point relative to a center
    def angle(pt):
        return np.arctan2(pt[1] - center[1], pt[0] - center[0])    

    # Loop through each cable linestring and anchor buffer
    for name,cab in cable_line.items():
        for anch in anchor_buffs:
            if cab.intersects(anch):
                # Get the start and end of the detour (the two closest points to the buffer)
                new_points = []
                # make additional points on the line on either side of anchor
                dist_to_anch = cab.line_locate_point(anch.centroid)
                if dist_to_anch > 100:
                    new_points.append(cab.interpolate(dist_to_anch - 100))
                if cab.length - dist_to_anch > 100:
                    new_points.append(cab.interpolate(dist_to_anch + 100))

                # pull out the coordinates of the first new point
                start = np.array(new_points[0].coords[-1])

                # Get buffer center and radius
                center = np.array(anch.centroid.coords[0])
                radius = anch.boundary.distance(sh.Point(center))+padding

                # Calculate angle for start point relative to center
                angle_start = angle(start)

                # Generate point along the arc (detour)
                arc_point = [center[0] + radius * np.cos(angle_start+np.pi/2), center[1] + radius * np.sin(angle_start+np.pi/2)]

                # determine relative positions of new routing points among other routing points
                rel_dist = []
                orig_coords = []
                for i,x in enumerate(proj.cableList[name].subcomponents[2].x):
                    y = proj.cableList[name].subcomponents[2].y[i]
                    rel_dist.append(cab.line_locate_point(sh.Point([x,y])))
                    orig_coords.append([x,y])
                all_dists = np.hstack((rel_dist,dist_to_anch-100, dist_to_anch+100, dist_to_anch))
                all_points = np.vstack((orig_coords,[coord.coords[0] for coord in new_points],arc_point))
                sorted_idxs = np.argsort(all_dists)
                final_points = all_points[sorted_idxs]

                # add new routing point in cable object
                proj.cableList[name].subcomponents[2].updateRouting(final_points)





            
def configureAdjuster(mooring, adjuster=None, method='horizontal',
                      i_line=0, span=None, project=None, target=None):
    '''Configures adjuster function for mooring object
    
    mooring : FAModel Mooring object
    adjuster : function, optional
        Function to adjust the mooring object
    method : str
        'horizontal' or 'pretension' ; method of adjusting the mooring
    i_line : int
        Line section number (0 is closest to end A) to adjust the length of 
    span : float
        Horizontal distance from fairlead to anchor (or fairlead to fairlead for shared lines)
    project : FAModel project class 
        Project class this mooring object is associated with. Required if adjuster function provided and 
        target is not provided, or method=pretension
    target : target value(s) for method - either pretension value or horizontal force value in x and y
    
    '''
    #calculate target pretension or horizontal tension if none provided
    if adjuster != None and target == None:
        targetdd = deepcopy(mooring.dd)
        if project == None:
            raise Exception('Project class instance needs to be provided to determine target')
        targetdd['zAnchor'] = -project.depth
        if not mooring.shared:
            mooring.rA[2] = -project.depth
        mooring.createSubsystem(dd=targetdd)
        
        if method == 'horizontal':
            mooring.target = np.linalg.norm(mooring.ss.fB_L[:2])
        elif method =='pretension':
            mooring.target = np.linalg.norm(mooring.ss.fB_L)
        else:
            raise Exception('Invalid adjustment method. Must be pretension or horizontal')
        # return mooring depth to accurate val
        if not mooring.shared:
            depth = project.getDepthAtLocation(mooring.rA[0], mooring.rA[1])
            mooring.rA[2] = -depth
    else:
        mooring.target = target
    
    if adjuster!= None:
        mooring.adjuster = adjuster 
        mooring.i_line = i_line 
        
        # check if method is 'pretension' then save slope
        if method == 'pretension':
            if project == None:
                raise Exception('Project class instance needs to be provided to determine slope')
            
            if mooring.dd:
                
                #calculate mooring slope using base depth
                #**** this assumes that the mooring system is designed for the base depth*****
                mooring.slope = (project.depth+mooring.rB[2]) / mooring.dd['span']
                
            else:   
                if span:
                    mooring.slope = (project.depth+mooring.rB[2]) / span
                else:
                    raise Exception('Span required to perform adjustment')
                    
    return(mooring)

def adjustMooring(mooring, method = 'horizontal', r=[0,0,0], project=None, target=1e6,
                       i_line = 0, slope = 0.58, display=False ):
    '''Custom function to adjust a mooring, called by
    Mooring.adjust. Fairlead point should have already
    been adjusted.
    
    There are two methods: "pretension" geometrically adjusts the anchor point and matches pretension, intended for taut moorings.
    "horizontal" leaves anchor point in the same position and matches the horizontal forces, intended for catenary and semi-taut moorings
    
    Parameters
    ----------
    mooring : FAModel Mooring object
    r : array
        platform center location
    project : FAModel Project object this is a part of. 
        Optional, default is None. This is a required input for the "pretension" option to correctly move the anchor position
    target_pretension : float
        Total pretension OR horizontal force in N to target for the mooring line 
    i_line : int
        Index of line section to adjust
    slope: float
        depth over span for baseline case (to match same geometric angle for 'pretension' option)
    
        '''
    from moorpy.helpers import dsolve2
    ss = mooring.ss  # shorthand for the mooring's subsystem

    if method == 'pretension':
        
        # Find anchor location based on desired relation
        fairlead_rad = mooring.rad_fair
        fairlead_z = mooring.z_fair
        
        fairlead = ss.rB # fairlead point (already updated)
        
        #unit direction vector towards ORIGNAL anchor in x,y plane, and inputted slope as the z component
        xydist = mooring.span # np.linalg.norm([ss.rA[0] - ss.rB[0],ss.rA[1] - ss.rB[1]])
        phi = np.pi/2 - np.radians(mooring.heading)
        direction = np.array([np.cos(phi), np.sin(phi), -slope]) # np.array([(ss.rA[0] - ss.rB[0])/xydist, (ss.rA[1] - ss.rB[1])/xydist, -slope])
        
        #use project class to find new anchor interesection point, maintaining original line heading
        if project:
            r_anch = project.seabedIntersect(fairlead, direction)  # seabed intersection  
        else:
            print('Project must be inputted for the pretension method')
            return 
        #update mooring properties
        if display:
            print('R_anch new ', r_anch)
        mooring.dd['zAnchor'] = r_anch[2]
        mooring.z_anch = mooring.dd['zAnchor']
        mooring.rad_anch = np.linalg.norm(r_anch[:2]-r[:2])
        span = mooring.rad_anch - fairlead_rad
        mooring.setEndPosition(r_anch, 'a', sink=True)  # set the anchor position

        #move anchor attachments
        for i,att in enumerate(mooring.attached_to):
            iend = mooring.rA if i == 0 else mooring.rB
            if type(att).__name__ in 'Anchor':
                # this is an anchor, move anchor location
                if project:
                    project.updateAnchor(att) 
                else: 
                    att.r = iend
                if att.mpAnchor:
                    att.mpAnchor.r = att.r
                    
        # Estimate the correct line length to start with based on % of total length
        L_tot = sum([line.L for line in ss.lineList])
        initial_L_ratio = ss.lineList[i_line].L/L_tot
        ss.lineList[i_line].setL(np.linalg.norm(mooring.rB - mooring.rA)*initial_L_ratio)
            
        # Next we could adjust the line length/tension (if there's a subsystem)
           
        def eval_func(X, args):
            '''Tension evaluation function for different line lengths'''
            ss.lineList[i_line].L = X[0]  # set the first line section's length
            ss.staticSolve(tol=0.0001)  # solve the equilibrium of the subsystem
            return np.array([ss.TB]), dict(status=1), False  # return the end tension

        # run dsolve2 solver to solve for the line length that matches the initial tension
        X0 = [ss.lineList[i_line].L]  # start with the current section length
        if display:
            L_final, T_final, _ = dsolve2(eval_func, X0, Ytarget=[target], 
                                  Xmin=[1], Xmax=[1.1*np.linalg.norm(ss.rB-ss.rA)],
                                  dX_last=[1], tol=[0.01], maxIter=50, stepfac=4, display=5)
        else:
            L_final, T_final, _ = dsolve2(eval_func, X0, Ytarget=[target], 
                                  Xmin=[1], Xmax=[1.1*np.linalg.norm(ss.rB-ss.rA)],
                                  dX_last=[1], tol=[0.01], maxIter=50, stepfac=4)
        ss.lineList[i_line].L = L_final[0]
        sec = mooring.getSubcomponent(i_line)
        sec['L'] = L_final[0]
        mooring.dd['span'] = span
        mooring.span = span
            
    elif method == 'horizontal':
        def func_TH_L(X, args):
            '''Apply specified section L, return the horizontal pretension error.'''
            ss.lineList[i_line].setL(X[0])
            ss.staticSolve()
            
            #Fx is the horizontal pretension
            Fx = np.linalg.norm([ss.fB_L[0], ss.fB_L[1]])
            
            return np.array([Fx - target]), dict(status=1) , False
            
        X0 = [ss.lineList[i_line].L]
        if display:
            x, y, info = dsolve2(func_TH_L, X0,  tol=[0.01], 
                                 args=dict(direction='horizontal'), 
                                 Xmin=[10], Xmax=[2000], dX_last=[10], 
                                 maxIter=50, stepfac=4, display = 5)
        else:
            x, y, info = dsolve2(func_TH_L, X0,  tol=[0.01], 
                                 args=dict(direction='horizontal'), 
                                 Xmin=[10], Xmax=[2000], dX_last=[10], 
                                 maxIter=50, stepfac=4)

    else:
        print('Invalid method. Must be either pretension or horizontal')
        
def yamlList(in_list):
    '''
    Function to convert a list to a ruaml list class with flow style set.
    Helpful for unloading ontologies in a readable manner.

    Parameters
    ----------
    in_list : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    yaml_list = ruamel.yaml.comments.CommentedSeq(in_list)
    yaml_list.fa.set_flow_style()
    return(yaml_list)

def cleanDataTypes(info, convert_lists=True):
    '''
    cleans data types in a dictionary to be yaml-writeable data types, usually for 
    dumping a yaml file i.e. in project.unload()
    
    Will convert lists to a ruamel.yaml list object with flow style set for ease of reading
    when unloading an ontology
    
    Example: will convert a numpy float64 to regular float

    Parameters
    ----------
    info : dictionary
        Dictionary (can be nested dictionary & can contain lists) that needs to convert
        data types to yaml-writeable data types.

    Returns
    -------
    info : dictionary
        Dictionary with yaml-writeable data types

    '''
    
    def fixType(val):
        if 'str' in type(val).__name__:
            valnew = str(val)
        elif 'float' in type(val).__name__:
            valnew = float(val)
        elif 'int' in type(val).__name__ and not 'Joint' in type(val).__name__:
            valnew = int(val)
        else:
            valnew = val

        return(valnew)
    
    def gothroughdict(dat):
        for key,value in dat.items():
            if isinstance(value,dict):
                value = gothroughdict(value)
            elif isinstance(value,list):
                value = gothroughlist(value)
            elif 'array' in type(value).__name__:
                value = value.tolist()
                value = gothroughlist(value)
            dat[key] = fixType(value)
        return(dat)
                
    def gothroughlist(dat):
        bottom_list = True
        for i,value in enumerate(dat):
            if isinstance(value,dict):
                value = gothroughdict(value)
                bottom_list = False
            elif isinstance(value,list):
                value = gothroughlist(value)
                bottom_list = False
            elif 'array' in type(value).__name__:
                value = value.tolist()
                value = gothroughlist(value)
                bottom_list = False
                
            dat[i] = fixType(value)
        if bottom_list and convert_lists:
            dat = yamlList(dat)
            
        return(dat)
    
    # recursively go through and clean data types for everything in this dictionary
    info = gothroughdict(info) 
    # return cleaned dictionary           
    return(info)


def createRAFTDict(project):
    from famodel.turbine.turbine import Turbine
    # Create a RAFT dictionary from a project class to create RAFT model
    rd = {'array':{'keys':['ID', 'turbineID', 'platformID', 'mooringID', 'x_location', 'y_location', 'heading_adjust'],
                   'data':[]}}
    turb = 0
    for pf in project.platformList.values():
        for att in pf.attachments.values():
            if isinstance(att['obj'],Turbine):
                turb = att['obj'].dd['type']
                break
        rd['array']['data'].append([pf.id, turb, pf.dd['type'], 0, pf.r[0], pf.r[1],np.degrees(pf.phi)])
        rd['site'] = {'water_depth':project.depth,'rho_water':project.rho_water,'rho_air':project.rho_air,'mu_air':project.mu_air}
        rd['site']['shearExp'] = .12
        
    rd['turbines'] = project.turbineTypes
    rd['platforms'] = project.platformTypes

    return rd


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

def updateYAML_array(fname,outx,outy,turbID,pfID,moorID,hadjust,newFile):
    '''
    Write turbines and locations to yaml file. Recommend using a different file name for the output than 
    the input yaml file, because the array table section will be written out in a manner that is not as readable
    
    Parameters
    ----------
    fname : str
        filename of yaml to read from
    outx : array
        1D array of x coordinates for platform locations
    outy : array
        1D array of y coordinates for platform locations
    turbID : array
        1D array of ID number of turbine for each platform, referencing turbine listing in ontology yaml
    pfID : array
        1D array of ID number of platform type for each platform, referencing platform listing in ontology yaml
    moorID : str or int
        1D array of ID of mooring system for each platform, referencing mooring system in ontology yaml
    hadjust : array
        1D array of angle rotation for each platform 
    newFile : str
        New file to write yaml to
    '''
    import ruamel.yaml
    yaml = ruamel.yaml.YAML()
    
    # read in yaml file 
    with open(fname) as fp:
        data = yaml.load(fp)
        
    # add rows for all platforms with general info
    data['array']['data'] = [] # remove any existing rows in the array table
    for i in range(0,len(outx)):
        data['array']['data'].append(['fowt'+str(i),turbID[i],pfID[i],moorID[i],float(outx[i]),float(outy[i]),hadjust[i]])
    
    # write to yaml file
    with open(newFile,'w') as f:    
        yaml.dump(data,f)

def updateYAML_MooringConfig(fname,ms,newfile):
    '''
    Update a yaml file with mooring configuration and mooring line type info from a moorpy system

    Parameters
    ----------
    fname : str
        YAML file to read in
    ms : object
        MoorPy system
    newfile : str
        YAML file to write to

    Returns
    -------
    None.

    '''
    import ruamel.yaml 
    yaml = ruamel.yaml.YAML()
    from moorpy.subsystem import Subsystem
    from moorpy.helpers import lines2ss
    
    # read in yaml file
    with open(fname) as fp:
        data = yaml.load(fp)
        
    # fill in mooring line types info
    for mtype in ms.lineTypes:
        data['mooring_line_types'][mtype] = ms.lineTypes[mtype]
    
    for i,line in enumerate(ms.lineList):
        if not isinstance(line,Subsystem):
            # convert to subsystem
            lines2ss(ms)
        types = []
        lengths = []
        connType = []
        for seg in line:
            types.append(seg.type['name'])
            lengths.append(seg.L)
        for j,conn in enumerate(line.pointList):
            connType.append({})
            if conn.m != 0:
                connType.append({'m':conn.m})
                connType[-1]['v'] = conn.v 
                connType[-1]['Cd'] = conn.cd
        
        # 

            
            
    

def updateYAML_mooring(fname,ms,newfile):
    '''
    Update a yaml file with mooring line information and platform locations from a moorpy system.
    Does not support cables currently.

    Parameters
    ----------
    fname : str
        YAML file to read in
    ms : object
        MoorPy system
    newfile : str
        YAML file to write to

    Returns
    -------
    None.

    '''
    import ruamel.yaml 
    yaml = ruamel.yaml.YAML()
    
    # read in yaml file
    with open(fname) as fp:
        data = yaml.load(fp)
        
    # fill in mooring line types info
    for mtype in ms.lineTypes:
        data['mooring_line_types'][mtype] = ms.lineTypes[mtype]
    
    # convert to subsystems if needed
    
    
    # fill in mooring systems info and parse for similar mooring systems (same line types, lengths, spans, # lines, headings)
    # msys = []
    # for pf in ms.bodyList:
    #     psys = []
    #     for pidx in pf.attachedP:
    #         lineDetails = []
    #         for lidx in ms.pointList[pidx-1].attached:
    #             ss = ms.lineList[lidx-1]
    #             lineDetails.append(ss.span) # span
    #             # calc heading
    #             heading = np.pi/2 - np.atan2((ss.rB[1]-ss.rA[1]),(ss.rB[0]-ss.rA[0]))
    #             lineDetails.append(heading) # heading
    #             # get segment details
    #             segType = []
    #             segL = []
    #             for seg in ss.lineList:
    #                 segType.append(seg.type['name'])
    #                 segL.append(seg.L)
    #             lineDetails.extend(segType,segL) # segment details
    #         # append line details to system for that platform
    #         psys.append(lineDetails)
    #     # append platfrom line system to mooring system
    #     msys.append(psys)
    '''spans = []
    headings = []
    segTypes = []
    segLs = []
    for pf in ms.bodyList:
        spansP = []
        headingsP = []
        segTypesP = []
        segLsP = []
        for pidx in pf.attachedP:
            for lidx in ms.pointList[pidx-1].attached: # should only be one line
                ss = ms.lineList[lidx-1]
                spansP.append(ss.span)
                # calc heading
                heading = np.pi/2 - np.atan2((ss.rB[1]-ss.rA[1]),(ss.rB[0]-ss.rA[0]))
                headingsP.append(heading) # heading
                # get segment details
                segType = []
                segL = []
                for seg in ss.lineList:
                    segType.append(seg.type['name'])
                    segL.append(seg.L)
                segTypesP.append(segType)
                segLsP.append(segL)
        # append line details to system for that platform
        spans.append(spansP)
        headings.append(headingsP)
        segTypes.append(segTypesP)
        segLs.append(segLsP)
    
    # first find where lengths of spans are the same
    lenspans = [len(span) for span in spans]
    uspanL = np.unique(lenspans)
    ind = []
    for i in range(0,len(uspanL)):
        ind.append(np.where(spans == uspanL[i])[0])
    # then find where spans are the same    
    
    spinds = []
    for ix in ind:
        for jx in ix:
            spinds_in = []
            for k in range(0,len(spans[jx])):
                spinds_in.append(np.where(spans[ix][k]==spans[jx][k])[0])
            if spinds_in[-1]'''
                

            
        
    
        
    # add rows for all platforms with general info
    data['array']['data'] = []