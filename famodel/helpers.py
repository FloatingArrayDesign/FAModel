
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

