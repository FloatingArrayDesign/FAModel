"""A set of functions for processing seabed information for a Project."""

import matplotlib.pyplot as plt
import numpy as np




def processASC(filename, lat, lon):
    '''Calculate anchor holding capacity based on specified anchor and soil
    information.
    
    Parameters
    ----------
    filename : string
        GEBCO or similar filename ASC file..
    lat : float
        lattitude of reference point to use for array y grid
    long : float
        lattitude of reference point to use for array x grid

    Returns
    -------
    Xs : array
        x values of grid points [m]
    Ys : array
        y values of grid points [m]
    Zs?  : 2D array
        depth grid (positive down?) [m]
    '''

    # stein's code for the conversion here 


    return Xs, Ys, Zs?

