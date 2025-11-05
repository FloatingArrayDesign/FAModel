"""A set of functions for processing seabed information for a Project."""


import os
import matplotlib.pyplot as plt
import numpy as np






def readBathymetryFile(filename, dtype=float):

    with open(filename, 'r') as f:
        # skip the header
        line = next(f)
        # collect the number of grid values in the x and y directions from the second and third lines
        line = next(f)
        nGridX = int(line.split()[1])
        line = next(f)
        nGridY = int(line.split()[1])
        # allocate the Xs, Ys, and main bathymetry grid arrays
        bathGrid_Xs = np.zeros(nGridX)
        bathGrid_Ys = np.zeros(nGridY)
        bathGrid = np.zeros([nGridY, nGridX], dtype=dtype)  # MH swapped order June 30
        # read in the fourth line to the Xs array
        line = next(f)
        bathGrid_Xs = [float(line.split()[i]) for i in range(nGridX)]
        strlist = []
        # read in the remaining lines in the file into the Ys array (first entry) and the main bathymetry grid
        for i in range(nGridY):
            line = next(f)
            entries = line.split()
            bathGrid_Ys[i] = entries[0]
            if dtype==float:
                bathGrid[i,:] = entries[1:]
            if dtype==str:
                strlist.append(entries[1:])
        if dtype==str:
            bathGrid = np.array(strlist)
    
    return bathGrid_Xs, bathGrid_Ys, bathGrid

            
def getSoilTypes(filename):
    '''function to read in a preliminary input text file format of soil type information'''

    soilProps = {}

    f = open(filename, 'r')
    
    for line in f:
        if line.count('---') > 0 and (line.upper().count('SOIL TYPES') > 0):
            line = next(f) # skip this header line, plus channel names and units lines
            var_names = line.split()
            line = next(f)
            line = next(f)
            while line.count('---') == 0:
                entries = line.split()
                soilProps[entries[0]] = {}
                for iv,var in enumerate(var_names[1:]):
                    # convert entries to strings unless there is 
                    if entries[iv+1] == '-':
                        soilProps[entries[0]][var] = [0]
                    else:
                        soilProps[entries[0]][var] = [float(entries[iv+1])]
                line = next(f)
    
    f.close()

    return soilProps


def processBoundary(filename, lat, lon,meters=True):
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
    
    import pandas as pd
    
    zerozero = (lat, lon)  # lattitude and longitude of reference point (grid origin)
    
    delin = pd.read_csv(filename)
    longs = np.array(delin['X_UTM10'])
    lats = np.array(delin['Y_UTM10'])
    
    if meters:
        Xs = longs
        Ys = lats
    #else:
        #Xs, Ys = convertLatLong2Meters(zerozero, lats, longs)
    
    return Xs, Ys









def resampleGrid(x_new, y_new, x_old, y_old, grid_values):
    '''Interpolate an existing array of values on a rectangular grid to a new
    rectangular grid.
    
    Parameters
    ----------
    x_new : list
        x values of the new grid to interpolate to
    y_new : list
        y values of the new grid to interpolate to
    x_old : list
        x values of the original grid
    y_old : list
        x values of the original grid
    grid_values : 2D array
        The values on the old grid to be interpolated from (dimensions must
        match the length of y_old and x_old, in that order).
    
    Returns
    -------
    grid_values_new : 2D array
        Interpolated grid values on y_new and x_new grid lines.
    '''
    
    grid_values_new = np.zeros([len(y_new), len(x_new)])
    
    for i in range(len(y_new)):
        for j in range(len(x_new)):
            grid_values_new[i,j], _ = getDepthFromBathymetry(x_new[j], y_new[i],
                                                   x_old, y_old, grid_values)
    
    return grid_values_new


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
    
    if np.isnan(xin):
        raise Exception('xin value is NaN.')
    
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


def interpFromGrid(x, y, grid_x, grid_y, values):
    '''Interpolate from a rectangular grid of values.'''

    # get interpolation indices and fractions for the relevant grid panel
    ix0, fx = getInterpNums(grid_x, x)
    iy0, fy = getInterpNums(grid_y, y)

    # handle end case conditions
    if fx == 0:
        ix1 = ix0
    else:
        ix1 = min(ix0+1, values.shape[1])  # don't overstep bounds
    
    if fy == 0:
        iy1 = iy0
    else:
        iy1 = min(iy0+1, values.shape[0])  # don't overstep bounds
    
    # get corner points of the panel
    c00 = values[iy0, ix0]
    c01 = values[iy1, ix0]
    c10 = values[iy0, ix1]
    c11 = values[iy1, ix1]

    # get interpolated points and local value
    cx0    = c00 *(1.0-fx) + c10 *fx
    cx1    = c01 *(1.0-fx) + c11 *fx
    c0y    = c00 *(1.0-fy) + c01 *fy
    c1y    = c10 *(1.0-fy) + c11 *fy
    value  = cx0 *(1.0-fy) + cx1 *fy

    # get local slope
    dx = grid_x[ix1] - grid_x[ix0]
    dy = grid_y[iy1] - grid_y[iy0]
    
    # deal with being on an edge or a zero-width grid increment
    if dx > 0.0:
        dc_dx = (c1y-c0y)/dx
    else:
        dc_dx = c0y*0  # maybe this should raise an error
    
    if dy > 0.0:
        dc_dy = (cx1-cx0)/dy
    else:
        dc_dy = cx0*0  # maybe this should raise an error
    
    # return the interpolated value, the derivatives, and the grid indices
    return value, dc_dx, dc_dy, ix0, iy0



def getDepthFromBathymetry(x, y, grid_x, grid_y, grid_depth, index=False):
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
    index : bool, optional
        If True, will also retun ix and iy - the indices of the intersected
        grid panel.
    '''
    
    # Call general function for 2d interpolation
    depth, dc_dx, dc_dy, ix0, iy0 = interpFromGrid(x, y, grid_x, grid_y, grid_depth)
    
    # Compute unit vector of the seabed panel
    nvec = np.array([dc_dx, dc_dy, 1.0])/np.linalg.norm([dc_dx, dc_dy, 1.0])
    
    if index:
        return depth, nvec, ix0, iy0
    else:
        return depth, nvec








if __name__ == '__main__':
    
    centroid = (40.928, -124.708)  #humboldt    
    xs = np.arange(-30000,30001,400)
    ys = np.arange(-40000,40001,400)
    
    xs, ys, depths = processGeotiff('humboldt.tif', centroid[0], centroid[1], xs=xs, ys=ys, outfilename='test output.txt')
    
    import moorpy as mp
    ms = mp.System(depth=np.max(depths), bathymetry='test output.txt')
    ms.initialize()
    ms.plot(hidebox=True, args_bath={'cmap':'viridis'})
    '''
    # try converting to a different grid
    x_new = np.arange(-20000, 20001, 800)
    y_new = np.arange(-20000, 20001, 800)
    depths_new = resampleGrid(x_new, y_new, xs, ys, depths)
    '''
    plt.show()
