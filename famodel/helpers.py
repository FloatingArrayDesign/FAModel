
import numpy as np
import time
import yaml
import os
import re



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


def check_headings(m_headings,c_heading,rad_buff):
    breakpoint()
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
        Cable heading at attachment to att.
    rad_buff : float
        Buffer angle in radians

    Returns
    -------
    heading : float
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
    # if hasattr(att,'mooring_headings'):                       
    #     # collect list of interfering mooring headings
    #     attheadings = np.pi/2 - (att.mooring_headings + att.phi)
    interfere_h = check_headings(attheadings,headnew,rad_buff)
    # if the headings interfere, adjust them by angle buffer
    # if interfere_h:
        # breakpoint()
        #ang_diffs = [headnew-mhead for mhead in interfere_h] # difference between cable heading and mooring headings for mooring headings that interfere with buffer
        # for ang_diff in ang_diffs:
        #     headnew = headnew + np.sign(ang_diff)*(rad_buff - abs(ang_diff))*endA_dir
        #     interfere_hi = check_headings(attheadings,headnew,rad_buff)
            
        #     for i in interfere_hi:
        #         # try rotating 30 degrees other way
        #         headnew = heading - np.sign(ang_diff)*rad_buff*endA_dir
        #         # re-check offsets
        #         interfere_hij = check_headings(attheadings,headnew,rad_buff)
        #         if not interfere_hij:
        #             return(headnew)
        #         else:
        #             # cut buffer in half and try again
        #             newbuff = rad_buff/2
        #             headnew = heading - np.sign(ang_diff)*newbuff
        #             return(headnew)
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
