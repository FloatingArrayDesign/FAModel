
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

