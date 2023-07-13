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
        j = 0
        while np.sign(depths[i,j])==1:
            ixc = j
            j += 1
        #ixc = np.argmin(np.abs(depths[i]))
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


def getInterpNums(xlist, xin, istart=0):  # should turn into function in helpers
    '''
    Paramaters
    ----------
    xlist : array
        list of x values
    xin : float
        x value to be interpolated
    istart : int
        first lower index to try
    
    Returns
    -------
    i : int
        lower index to interpolate from
    fout : float
        fraction to return   such that y* = y[i] + fout*(y[i+1]-y[i])
    '''
    
    nx = len(xlist)
  
    if xin <= xlist[0]:  #  below lowest data point
        i = 0
        fout = 0.0
  
    elif xlist[-1] <= xin:  # above highest data point
        i = nx-1
        fout = 0.0
  
    else:  # within the data range
 
        # if istart is below the actual value, start with it instead of 
        # starting at 0 to save time, but make sure it doesn't overstep the array
        if xlist[min(istart,nx)] < xin:
            i1 = istart
        else:
            i1 = 0

        for i in range(i1, nx-1):
            if xlist[i+1] > xin:
                fout = (xin - xlist[i] )/( xlist[i+1] - xlist[i] )
                break
    
    return i, fout


def getDepthFromBathymetry(x, y, bathGrid_Xs, bathGrid_Ys, bathGrid):   #BathymetryGrid, BathGrid_Xs, BathGrid_Ys, LineX, LineY, depth, nvec)
    ''' interpolates local seabed depth and normal vector
    
    Parameters
    ----------
    x, y : float
        x and y coordinates to find depth and slope at [m]
    
    Returns
    -------        
    depth : float
        local seabed depth (positive down) [m]
    nvec : array of size 3
        local seabed surface normal vector (positive out) 
    '''

    # get interpolation indices and fractions for the relevant grid panel
    ix0, fx = getInterpNums(bathGrid_Xs, x)
    iy0, fy = getInterpNums(bathGrid_Ys, y)


    # handle end case conditions
    if fx == 0:
        ix1 = ix0
    else:
        ix1 = min(ix0+1, bathGrid.shape[1])  # don't overstep bounds
    
    if fy == 0:
        iy1 = iy0
    else:
        iy1 = min(iy0+1, bathGrid.shape[0])  # don't overstep bounds
    

    # get corner points of the panel
    c00 = bathGrid[iy0, ix0]
    c01 = bathGrid[iy1, ix0]
    c10 = bathGrid[iy0, ix1]
    c11 = bathGrid[iy1, ix1]

    # get interpolated points and local value
    cx0    = c00 *(1.0-fx) + c10 *fx
    cx1    = c01 *(1.0-fx) + c11 *fx
    c0y    = c00 *(1.0-fy) + c01 *fy
    c1y    = c10 *(1.0-fy) + c11 *fy
    depth  = cx0 *(1.0-fy) + cx1 *fy

    # get local slope
    dx = bathGrid_Xs[ix1] - bathGrid_Xs[ix0]
    dy = bathGrid_Ys[iy1] - bathGrid_Ys[iy0]
    
    if dx > 0.0:
        dc_dx = (c1y-c0y)/dx
    else:
        dc_dx = 0.0  # maybe this should raise an error
    
    if dx > 0.0:
        dc_dy = (cx1-cx0)/dy
    else:
        dc_dy = 0.0  # maybe this should raise an error
    
    nvec = np.array([dc_dx, dc_dy, 1.0])/np.linalg.norm([dc_dx, dc_dy, 1.0])  # compute unit vector      

    return depth, nvec


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

