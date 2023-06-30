"""A set of functions for processing seabed information for a Project."""


import matplotlib.pyplot as plt
import numpy as np
import geopy.distance


def processASC(gebcofilename, lat, lon, outfilename=""):
    '''Process an ASC file of bathymetry information and convert into
    a rectangular bathymetry grid in units of me relative to the 
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
        If provided, writes a MoorDyn/MoorPy style bathymetry file.

    Returns
    -------
    Xs : array
        x values of grid points [m]
    Ys : array
        y values of grid points [m]
    Zs?  : 2D array
        depth grid (positive down?) [m]
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


    # ----- save a MoorDyn/MoorPy-style bathymetry file if requested -----
    
    # assuming that we don't need to change the depth values based on the curvature of the earth
    # assuming (since we only need plots right now) that we don't need to adjust the x/y values based on arcseconds due to curvature
    if len(outfilename) > 0:

        f = open(os.path.join(os.getcwd(), outfilename'), 'w')
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
                f.write(f'{depths[iy,id]} ')
            f.write('\n')

        f.close()

    return Xs, Ys, depths


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

    # >>> placeholder code to be updated <<<
    
    delin = pd.read_csv("humboldt.csv")
    longs = np.array(delin['X_UTM10'][0:75])
    lats = np.array(delin['Y_UTM10'][0:75])
    

    Xs = np.zeros(len(longs))
    Ys = np.zeros(len(lats))
    for i in range(len(longs)):
        if longs[i] < zerozero[1]:
            sign = -1
        else:
            sign = 1
        Xs[i] = geopy.distance.distance(zerozero, (zerozero[0], longs[i])).km*1000*sign
    for i in range(len(lats_ne)):
        if lats_ne[i] < zerozero[0]:
            sign = -1
        else:
            sign = 1
        Ys[i] = geopy.distance.distance(zerozero, (lats[i], zerozero[1])).km*1000*sign

    return Xs, Ys