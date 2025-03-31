
import numpy as np
import time
import yaml
import os
import re
from copy import deepcopy
from famodel.cables.cable_properties import getCableProps, getBuoyProps, loadCableProps,loadBuoyProps



def cart2pol(x, y):
    rho = np.hypot(x, y)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)


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
            #breakpoint()
            m_headings[i] = 2*np.pi + mh
    ang_diff = m_headings - c_heading
    inds_to_fix = np.where([round(abs(angd),8)<round(rad_buff,8) for angd in ang_diff])[0]
    if len(inds_to_fix)>0:
        return([m_headings[ind] for ind in inds_to_fix])
    else:
        return([])
    
        
def head_adjust(att,heading,rad_buff=np.radians(30),endA_dir=1):
    '''
    function to adjust heading of cable based on angle buffer from mooring lines

    Parameters
    ----------
    att : list
        list of objects to attach to. 1 object if only concerned about the attached object associated with that side
    heading : float
        Cable heading at attachment to att in radians
    rad_buff : float
        Buffer angle in radians
    endA_dir : float, optional
        Either 1 or -1, controls sign of new heading for end B. Only altered to -1 if dynamic
        cable from end A will get close to end B moorings. Default is 1.

    Returns
    -------
    headnew : float
        New cable heading

    '''
    if heading<0:
        headnew = np.pi*2 + heading
    else:
        headnew = heading
    #breakpoint()
    attheadings = []
    flipheads = False # whether to flip headings ( for if you are looking at mooring headings of platform on the other end)
    for at in att:
        if flipheads:
            atmh = at.mooring_headings + at.phi + np.pi
            for j,a in enumerate(atmh):
                # keep everything under 2pi angle
                if a>2*np.pi:
                    atmh[j] = a-2*np.pi
        else:
            atmh = at.mooring_headings + at.phi
        attheadings.extend(np.pi/2 - atmh)
        flipheads = True

    interfere_h = check_headings(attheadings,headnew,rad_buff)
    # if the headings interfere, adjust them by angle buffer
    for mhead in interfere_h:
        ang_diff_dir = np.sign(headnew - mhead)
        headnew = mhead + rad_buff*endA_dir*ang_diff_dir #headnew + np.sign(ang_diff)*(rad_buff - abs(ang_diff))*endA_dir
        interfere_hi = check_headings(attheadings,headnew,rad_buff)
        
        for i in interfere_hi:
            # try rotating other way
            headnew = mhead - rad_buff*endA_dir*ang_diff_dir
            # re-check offsets
            interfere_hij = check_headings(attheadings,headnew,rad_buff)
            if not interfere_hij:
                return(headnew)
            else:
                # cut buffer in half and try again
                newbuff = rad_buff/2
                headnew = mhead - newbuff*endA_dir*ang_diff_dir
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
        dd['d_nom'] = mProps['d_nom']
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
    for k in range(0,len(lineConfigs[lcID]['sections'])): # loop through each section in the line
    
        lc = lineConfigs[lcID]['sections'][k] # set location for code clarity later
        # determine if it's a line type or a connector listed
        if 'type' in lc or 'mooringFamily' in lc: 
            # this is a line
            if lineLast: # previous item in list was a line (or this is the first item in a list)
                # no connector was specified for before this line - add an empty connector
                c_config.append({})                        
            # set line information
            lt = MooringProps(lc, proj.lineTypes, proj.rho_water, proj.g)                                             
            # lt = self.lineTypes[lc['type']] # set location for code clarity and brevity later
            # set up sub-dictionaries that will contain info on the line type
            m_config['sections'].append({'type':lt})# {'name':str(ct)+'_'+lc['type'],'d_nom':lt['d_nom'],'material':lt['material'],'d_vol':lt['d_vol'],'m':lt['m'],'EA':float(lt['EA'])}})
            m_config['sections'][ct]['type']['name'] = str(ct)+'_'+str(lt['name'])
            # make EA a float not a string
            m_config['sections'][ct]['type']['EA'] = float(lt['EA'])  
            # set line length
            m_config['sections'][ct]['L'] = lc['length']
            # update counter for line types 
            ct = ct + 1
            # update line last boolean
            lineLast = 1
            
        elif 'connectorType' in lc:
            # this is a connector
            if lineLast == 0:
                # last item in list was a connector
                raise Exception(f"Two connectors were specified in a row for line configuration '{lcID}', please remove one of the connectors")
            else:
                # last item in list was a line
                c_config.append(connectorTypes[lc['connectorType']]) # add connector to list
                c_config[-1]['type'] = lc['connectorType']
                # update lineLast boolean
                lineLast = 0
        else:
            # not a connector or a line
            raise Exception(f"Please make sure that all section entries for line configuration '{lcID}' are either line sections (which must have a 'type' key) or connectors (which must have a 'connectorType' key")

    # check if line is a shared symmetrical configuration
    if 'symmetric' in lineConfigs[lcID] and lineConfigs[lcID]['symmetric']:
        if not lineLast: # check if last item in line config list was a connector
            for ii in range(0,ct):
                # set mooring configuration 
                m_config['sections'].append(m_config['sections'][-1-2*ii])
                # set connector (since it's mirrored, connector B becomes connector A)
                c_config.append(c_config[-2-2*ii])
        else: # double the length of the end line
            m_config['sections'][-1]['L'] = m_config['sections'][-1]['L']*2
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
    m_config['zAnchor'] = -proj.depth 
    m_config['span'] = lineConfigs[lcID]['span']
    m_config['name'] = lcID
    # add fairlead radius and depth to dictionary
    m_config['rad_fair'] = proj.platformList[pfID].rFair
    m_config['z_fair'] = proj.platformList[pfID].zFair
    
    m_config['connectors'] = c_config  # add connectors section to the mooring dict
    
    return(m_config) #, c_config)

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
        else:
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
    ad['type'] = proj.anchorTypes[lineAnch]['type']
    ad['name'] = lineAnch
    
    return(ad)

def cleanDataTypes(info):
    '''
    cleans data types in a dictionary to be yaml-writeable data types, usually for 
    dumping a yaml file i.e. in project.unload()
    
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
        #print(valnew)
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
        for i,value in enumerate(dat):
            if isinstance(value,dict):
                value = gothroughdict(value)
            elif isinstance(value,list):
                value = gothroughlist(value)
            elif 'array' in type(value).__name__:
                value = value.tolist()
                value = gothroughlist(value)
            dat[i] = fixType(value)
        return(dat)
    
    # recursively go through and clean data types for everything in this dictionary
    info = gothroughdict(info) 
    # return cleaned dictionary           
    return(info)

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
