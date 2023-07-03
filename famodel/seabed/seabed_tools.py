"""A set of functions for processing seabed information for a Project."""


import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import geopy.distance



def convertLatLong2Meters(zerozero, lats, longs):
    '''Convert a list of latitude and longitude coordinates into
    x-y positions relative to a project reference point.
    
    Parameters
    ----------
    zerozero : tuple
        A tuple or list of two values, x and y, of the project 
        reference point
    lats : array
        array of latitude coordinates (y positions)
    longs : array
        array of longitude coordinates (x positions)

    Returns
    -------
    Xs : array
        x values of grid points (from longitudes) [m]
    Ys : array
        y values of grid points (from latitudes) [m]
    '''

    Xs = np.zeros(len(longs))
    Ys = np.zeros(len(lats))
    for i in range(len(longs)):
        if longs[i] < zerozero[1]:
            sign = -1
        else:
            sign = 1
        Xs[i] = geopy.distance.distance(zerozero, (zerozero[0], longs[i])).km*1000*sign
    for i in range(len(lats)):
        if lats[i] < zerozero[0]:
            sign = -1
        else:
            sign = 1
        Ys[i] = geopy.distance.distance(zerozero, (lats[i], zerozero[1])).km*1000*sign

    return Xs, Ys


def processASC(gebcofilename, lat, lon, outfilename=""):
    '''Process an ASC file of bathymetry information and convert into
    a rectangular bathymetry grid in units of m relative to the 
    project reference point.
    
    Parameters
    ----------
    filename : string
        GEBCO or similar filename ASC file..
    lat : float
        lattitude of reference point to use for array y grid
    long : float
        lattitude of reference point to use for array x grid
    outfilename : string, optional
        If provided, writes a MoorDyn/MoorPy style bathymetry file

    Returns
    -------
    Xs : array
        x values of grid points [m]
    Ys : array
        y values of grid points [m]
    depths  : 2D array
        water depth grid (positive down) [m]
    lats : array
        array of latitude coordinates (y positions)
    longs : array
        array of longitude coordinates (x positions)
    '''

    depths = -np.loadtxt(gebcofilename, skiprows=6)
    nrows = len(depths)
    ncols = len(depths[0])

    newlinedata = []
    with open(gebcofilename) as f:
        lines = f.readlines()
    for i,line in enumerate(lines):
        if i < 6: newlinedata.append(line.split())
        if i==2: xllcorner = float(line.split()[1])
        if i==3: yllcorner = float(line.split()[1])
        if i==4: cellsize = float(line.split()[1])

    longs = np.linspace(xllcorner, xllcorner+ncols*cellsize, ncols)
    lats  = np.linspace(yllcorner, yllcorner+nrows*cellsize, nrows)

    zerozero = (lat, lon)  # lattitude and longitude of reference point (grid origin)

    Xs, Ys = convertLatLong2Meters(zerozero, lats, longs)

    # assuming that we don't need to change the depth values based on the curvature of the earth
    # assuming that we don't need to adjust the x/y values based on arcseconds due to curvature
    
    # ----- save a MoorDyn/MoorPy-style bathymetry file if requested -----
    
    if len(outfilename) > 0:

        f = open(os.path.join(os.getcwd(), outfilename), 'w')
        f.write('--- MoorPy Bathymetry Input File ---\n')
        f.write(f'nGridX {ncols}\n')
        f.write(f'nGridY {nrows}\n')
        f.write(f'      ')
        for ix in range(len(Xs)):
            f.write(f'{Xs[ix]:.2f} ')
        f.write('\n')
        for iy in range(len(Ys)):
            f.write(f'{Ys[iy]:.2f} ')
            for id in range(len(depths[iy])):
                iy2 = len(depths) - iy-1
                f.write(f'{depths[iy2,id]} ')
            f.write('\n')

        f.close()

    return Xs, Ys, depths, lats, longs


def getCoast(Xbath, Ybath, depths):
    '''Gets the x and y coordinates of the coastline from the bathymetry
    data to be used for plotting.
        
    Parameters
    ----------
    Xbath : array
        x values of bathymetry grid points [m]
    Ybath : array
        y values of bathymetry grid points [m]
    depths  : 2D array
        water depth grid (positive down) [m]

    Returns
    -------
    xcoast : array
        x values of coastal grid points [m]
    ycoast : array
        y values of coastal grid points [m]
    '''
    
    xcoast = np.zeros(len(Ybath))
    ycoast = np.zeros(len(Ybath))
    for i in range(len(depths)):
        ixc = np.argmin(np.abs(depths[i]))
        iyc = len(depths) - i-1
        xcoast[i] = Xbath[ixc]
        ycoast[i] = Ybath[iyc]
    
    return xcoast, ycoast


def processBoundary(filename, lat, lon):
    '''Reads boundary information from a CSV file and stores the boundary 
    coordinate list in a set of arrays. This function can be extended to
    deal with multiple boundary sets.
        
    Parameters
    ----------
    filename : string
        Filename containing columns of x and y coordinates of boundary.
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
    '''
    
    zerozero = (lat, lon)  # lattitude and longitude of reference point (grid origin)
    
    delin = pd.read_csv(filename)
    longs = np.array(delin['X_UTM10'])
    lats = np.array(delin['Y_UTM10'])
    
    Xs, Ys = convertLatLong2Meters(zerozero, lats, longs)
    
    return Xs, Ys


def getPlotBounds(latsorlongs_boundary, zerozero, long=True):
    '''Gets the x and y bounds to be used in MoorPy.System.plot() so 
    that the center of the matplotlib plot will center around a 
    reference point for ease of viewing
        
    Parameters
    ----------
    latsorlongs_boundary : array
        An array of latitude or longitude coordinates
    zerozero : tuple
        latitude and longitude of reference point
    long : bool, optional
        flag for whether latitudes or longitudes are being used

    Returns
    -------
    xbmin : float
        x (or y) value to set minimum plotting bounds relative to 
        the project reference point [m]
    xbmax : float
        x (or y) value to set maximum plotting bounds relative to 
        the project reference point [m]
    '''
    
    if long:
        il = 1
    else:
        il = 0
    
    xmed = latsorlongs_boundary[int((len(latsorlongs_boundary) + 1) / 2)]
    xmin = latsorlongs_boundary[0]
    xmax = latsorlongs_boundary[-1]
    
    newxmin = xmin + (zerozero[il]-xmed)
    newxmax = xmax + (zerozero[il]-xmed)
    
    # convert lats/longs into m relative to project reference point
    if newxmin < zerozero[il]:
        sign = -1
    else:
        sign = 1
    if long:
        xbmin = geopy.distance.distance(zerozero, (zerozero[il-1], newxmin)).km*1000*sign
    else:
        xbmin = geopy.distance.distance(zerozero, (newxmin, zerozero[il-1])).km*1000*sign

    # convert lats/longs into m relative to project reference point
    if newxmax < zerozero[il]:
        sign = -1
    else:
        sign = 1
    if long:
        xbmax = geopy.distance.distance(zerozero, (zerozero[il-1], newxmax)).km*1000*sign
    else:
        xbmax = geopy.distance.distance(zerozero, (newxmax, zerozero[il-1])).km*1000*sign
    
    return xbmin, xbmax

