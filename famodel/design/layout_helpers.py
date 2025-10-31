import moorpy as mp
from moorpy.helpers import dsolve2
import numpy as np
import matplotlib.pyplot as plt
from shapely import Point, Polygon
from numpy import random
from copy import deepcopy
import time


from famodel.mooring.mooring import Mooring
from famodel.seabed.seabed_tools import getDepthFromBathymetry
from famodel.project import Project
from famodel.design.fadsolvers import dsolve2


def create_initial_layout(lease_xs, lease_ys, ms, grid_x, grid_y, grid_depth, update_ms=True, display=0):
    '''
    The first iteration of a layout generator function based off of Katherine's previous work in summer 2023.
    The idea is to come up with turbine locations within a lease area boundary that can be oriented various directions, 
    not overlap other mooring systems, and not extend outside of a lease area boundary.
    In reality, I'd imagine this function would become obsolete, as we could populate a lease area with random points and
    then optimize their positions, but the capabilities within this function can be used as a starting point to 
    incorporate into the optimization process.

    Right now, it loops through a "grid" of x and y positions, spaced relatively close together and calculates the anchor positions 
    of that turbine (using an adjustable mooring system orientation) by extending or shortening the mooring line until it
    contacts the seabed bathymetry along that path. Using these new anchor positions, the function then checks each turbine 
    position on whether 1) it, and the anchor x/y positions are within the lease area bounds and 2) the given triangle that 
    connects the anchor points overlaps any other existing triangles (i.e., footprints). If it satisfies those two criteria, 
    then it appends that turbine x/y position to the list of valid points.

    Parameters
    ----------
    lease_xs : float, array
        The x coordinates of coordinate pairs defining a layout boundary, relative to a certain point (usually a centroid)
    lease_ys : float, array
        The y coordinates of coordinate pairs defining a layout boundary, relative to a certain point (usually a centroid)
    ms : MoorPy System object
        A MoorPy System object defining a mooring system
    grid_x : float, array
        The x coordinates of coordinate pairs defining a bathymetry grid
    grid_y : float, array
        The y coordinates of coordinate pairs defining a bathymetry grid
    grid_depth: float, 2D matrix
        The depth (z coordinates) of the grid defined by grid_x and grid_y

    Returns
    -------
    xs : float, array
        A list of the x coordinates of turbine locations within the array
    ys : float, array
        A list of the y coordinates of turbine locations within the array
    footprintList : list, Polygon objects
        A list of shapely Polygon objects of each mooring system footprint based on anchor locations
    '''

    coords = []

    area = Polygon([(lease_xs[i], lease_ys[i]) for i in range(len(lease_xs))])

    # Brainstorm different initialization approaches
    # - Placing one at a starting point and filling in from there
    # - Choosing a predetermined number of turbine and making them fit
    # - Placing anchors and working backwards
    # not sure which one is the best right now; will stick with the first one of choosing a starting point and filling in around

    # Placing approaches:
    # make a very fine xlocs and ylocs grid
    # loop through and find the first point in the lease area that has all 3 anchors in the lease area
    # - the next point that's tested will obviously be way too close to the first point, but this will allow for better placement
    # make an "orientation" variable in case we want to switch the orientations (ms.transform) (typically either 180 or 0) (later todo item)

    xlocs = np.arange(np.min(lease_xs), np.max(lease_xs)+1, 1000)     # set a really small spacing between
    ylocs = np.arange(np.min(lease_ys), np.max(lease_ys)+1, 1000)

    # if you want to change the "starting" point, you will need to rearrange the xlocs and ylocs variables
    # something like "xlocs = np.hstack([(xlocs[int(len(xlocs))/2):], xlocs[:int(len(xlocs))/2)]])"

    # count how many anchor points there are
    anchors = [point.number for point in ms.pointList if point.r[2]==-ms.depth]        # <<<<< might need to change this assumption later on checking if it's on the seabed
    fairleads = ms.bodyList[0].attachedP

    # initialize a couple storage/bool variables
    invalid = False
    footprintList = []
    msList = []
    counter = 0

    # placeholder to deal with the mooring system orientation
    orientation = -180

    # loop through the xlocs and ylocs variables to test x/y positions to place turbines
    for ix in range(len(xlocs)):
        for iy in range(len(ylocs)):
            anchorGlobalTempList = []
            anchorLocalTempList = []
            
            # set the x/y position of a point to test
            point = [xlocs[ix], ylocs[iy]]      
            
            # orient the mooring system around that point by a certain amount
            orientation = 0
            #orientation += 180
            #orientation = random.choice([0,90,180,270])
            ms.transform(rot=orientation)

            # reinitialize the mooring system after reorientation
            ms.initialize()
            ms.solveEquilibrium()

            # loop through the anchors in the mooring system and evaluate whether they meet the criteria
            for i,anchornum in enumerate(anchors):

                old_anchor_point = ms.pointList[anchornum-1].r + np.array([point[0], point[1], 0])
                fairlead_point = ms.pointList[fairleads[i]-1].r + np.array([point[0], point[1], 0])
                # update the anchor point based on how close it is to the bathymetry
                new_anchor_point = getUpdatedAnchorPosition(old_anchor_point, fairlead_point, grid_x, grid_y, grid_depth)

                # check to make sure the updated anchor point is within the lease area boundary
                if not area.contains(Point(new_anchor_point[0], new_anchor_point[1])):
                    invalid = True
                
                # save the anchor point for later, regardless of whether it fits or not
                anchorGlobalTempList.append(new_anchor_point)
                anchorLocalTempList.append(new_anchor_point - np.array([point[0], point[1], 0]))
            
            
            # create new lists/polygons of using the anchor positions of this one turbine point
            #anchorList.append([anchor_point for anchor_point in anchorGlobalTempList])

            # create a shapely polygon made up of the anchor points
            new_boundary = Polygon( [(anchor_point[0], anchor_point[1]) for anchor_point in anchorGlobalTempList] )
            
            # check to make sure that the newly created polygon does not intersect any other polygons
            for moor_sys in footprintList:
                if moor_sys.intersects(new_boundary):
                    invalid = True
            

            


            # if all checks pass, then include this point to the list of coordinates of the farm and include the boundary polygon to reference later
            if invalid==False:

                # save the point to a list of "coords"
                coords.append(point)
                if display > 0: print(f"Appending Point ({point[0]:6.1f}, {point[1]:6.1f}) to the coords list")

                # add to the counter for the number of turbines that meet criteria
                counter += 1
                if display > 0: print(f'nTurbines = {counter}')

                # save the polygon footprint of the anchor points
                footprintList.append(Polygon( [(anchor_point[0], anchor_point[1]) for anchor_point in anchorGlobalTempList] ) )


                # if you want to update the mooring system line lengths to match pretension
                if update_ms:
                    if display > 0: print(f"Updating the mooring system at Point ({point[0]:6.1f}, {point[1]:6.1f}) ")
                    
                    # create a copy of the MoorPy System and of the anchor point list
                    mscopy = deepcopy(ms)

                    # adjust the MoorPy System line lengths to keep the same pretension (calls another internal function)
                    #ms_new = adjustMS4Pretension(mscopy, anchorLocalTempList)
                    ms_new = adjustMS4Bath(mscopy, point, grid_x, grid_y, grid_depth, display=display)

                else:
                    ms_new = deepcopy(ms)
                
                # save the adjusted (or not adjusted mooring system)
                msList.append(ms_new)
            


            # reset the invalid flag variable in case it was changed to true
            invalid = False
    

    # extract the x and y variables from the list of points
    xs = [xy[0] for xy in coords]
    ys = [xy[1] for xy in coords]

    return xs, ys, footprintList, msList



def create_layout(bound_xs, bound_ys, subsystem, grid_x, grid_y, grid_depth,
                  spacing_x, spacing_y, headings=[60, 180, 300]):
    '''
    Create a rectangular grid layout.
    
    Parameters
    ----------
    lease_xs : float, array
        The x coordinates of coordinate pairs defining a layout boundary, relative to a certain point (usually a centroid)
    lease_ys : float, array
        The y coordinates of coordinate pairs defining a layout boundary, relative to a certain point (usually a centroid)
    subsystem : MoorPy Subsystem object
        A MoorPy Subsystem object defining the mooring configuration to be used.
    grid_x : float, array
        The x coordinates of coordinate pairs defining a bathymetry grid
    grid_y : float, array
        The y coordinates of coordinate pairs defining a bathymetry grid
    grid_depth: float, 2D matrix
        The depth (z coordinates) of the grid defined by grid_x and grid_y
    spacing_x : float
        The x spacing between turbines [m].
    spacing_y : float
        The y spacing between turbines [m].

    Returns
    -------
    xs : float, array
        A list of the x coordinates of turbine locations within the array
    ys : float, array
        A list of the y coordinates of turbine locations within the array
    footprintList : list, Polygon objects
        A list of shapely Polygon objects of each mooring system footprint based on anchor locations
    '''
    

    # make sure the subsystem is initialized
    subsystem.initialize()
    
    # save dimensions from the subsystem
    rad_anch = np.linalg.norm(subsystem.rA[:2])
    rad_fair = np.linalg.norm(subsystem.rB[:2])
    z_anch   = subsystem.rA[2]
    z_fair   = subsystem.rB[2]
    
    # initialize some lists
    coords = []
    mooringList = []
    footprintList = []
    
    # make the bounds into a shapely Polygon
    area = Polygon([(bound_xs[i], bound_ys[i]) for i in range(len(bound_xs))])

    # Grid of turbine locations (only those in the boundaries will be kept)
    xlocs = np.arange(np.min(lease_xs), np.max(lease_xs)+1, spacing_x)
    ylocs = np.arange(np.min(lease_ys), np.max(lease_ys)+1, spacing_y)

    mooring_count = 0

    # loop through the xlocs and ylocs variables to test x/y positions to place turbines
    for ix in range(len(xlocs)):
        for iy in range(len(ylocs)):
            
            valid = True  # flag for whether the turbine position satisfies requirements
            
            # set the x/y position of a point to test
            point = [xlocs[ix], ylocs[iy]]
            
            # make sure the turbine location is in the boundary
            if not area.contains(Point(point)):
                valid = False

            # assume "orientation" is always 0
            # initialize a list
            anchorlist = []
            ssList = []

            # at the current grid point, set the anchor and fairlead points of the subsystem using a list of line heading angles and adjust for bathymetry
            for ang in headings:
                
                if not valid:
                    break
                
                th = np.radians(ang)
                
                # set the local anchor and fairlead points
                r_anch = np.hstack([rad_anch*np.array([np.cos(th), np.sin(th)])+point, z_anch])
                r_fair = np.hstack([rad_fair*np.array([np.cos(th), np.sin(th)])+point, z_fair])
                
                mooring_count += 1
                print(f"Mooring count is {mooring_count}.")
                
                ss = deepcopy(subsystem)        # make a copy from the original since we'll be iterating on this object

                # set the anchor and fairlead points of the subsystem
                #subsystem_copy.pointList[0].setPosition(r_anch)
                ss.setEndPosition(r_anch, endB=0)
                #subsystem_copy.pointList[-1].setPosition(r_fair)
                ss.setEndPosition(r_fair, endB=1)
                ss.staticSolve()


                # adjust subsystem for bathymetry (adjusting anchor points and line lengths)
                adjustSS4Bath(ss, grid_x, grid_y, grid_depth, display=0)

                new_anchor_point = ss.rA
                anchorlist.append(new_anchor_point)
                
                ssList.append(ss)  # add it to a temporary list for just this turbine

                # if any new anchor point is outside the bounds of the Polygon area, then this point is invalid
                if not area.contains(Point(new_anchor_point)):
                    valid = False

            # if not valid, skip the rest of this point in the for loop
            if not valid:
                continue

            # after checking all new anchor points for each line heading, check to make sure the new footprint doesn't overlap with any others
            new_footprint = Polygon( [(anchor_point[0], anchor_point[1]) for anchor_point in anchorlist] )
        
            # check to make sure that the newly created polygon does not intersect any other polygons
            for footprint in footprintList:
                if footprint.intersects(new_footprint):
                    valid = False

            # make the moorings and add to the master lists if valid
            if valid:

                for i, ss in enumerate(ssList):
                    mooringList.append(Mooring(subsystem=ss, rA=ss.rA, rB=ss.rB))
                
                coords.append(point)

                footprintList.append( Polygon( [(anchor_point[0], anchor_point[1]) for anchor_point in anchorlist] ) )
    

    return np.array(coords), mooringList, footprintList


def create_rotated_layout(bound_xs, bound_ys, spacing_x, spacing_y, grid_rotang, grid_skew_x, grid_skew_y, grid_trans_x, grid_trans_y, fullMPsystem = True, ms = None, rad_anch = None, rotations=None, center=None):
        '''
        Create a rectangular grid layout.
        
        Parameters
        ----------
        bound_xs : list
            The x coordinates of coordinate pairs defining a layout boundary, relative to a certain point (usually a centroid)
        bound_ys : list
            The y coordinates of coordinate pairs defining a layout boundary, relative to a certain point (usually a centroid)
        spacing_x : float
            The x spacing between turbines [m].
        spacing_y : float
            The y spacing between turbines [m].
        grid_rotang : float
            Rotation of y axis in deg (0 deg is due North, 90 deg is due West)
        grid_skew_x : float
            Angle of parallelogram between adjacent rows in deg
        grid_skew_y : float
            Angle of parallelogram between adjacent columns in deg
        grid_trans_x : float
            x offset to add to all turbine positions
        grid_trans_y : float
            y offset to add to all turbine positions   
        fullMPsystem : bool
            if True, create/rotation full moorpy systems (slower). if False, use rad_anch as circular buffer zone
        ms : MoorPy system
            mooring system to rotate, need to input if fullMPsystem = True
        rad_anch : float
            mooring anchoring radius. need to input if fullMPsystem = False
        rotations: list
            list of two mooring orientations in deg relative to the rotated y axis (used for every other row). not used if fullMPsystem = False
        center: list
            the coordinate of the center of the layout. Default: the midpoint of the x and y bounds
        Returns
        -------
        x_coords : array
            array of turbine x coordinates
        y_coords : array
            array of turbine y coordinates
        moorings : list
            list of mooring systems 
        area : shapely polygon
            polygon of boundary 
        '''
        
        #boundary of area
        area = Polygon([(bound_xs[i], bound_ys[i]) for i in range(len(bound_xs))])
    
        # Shear transformation in X
        # Calculate trigonometric values
        cos_theta = np.cos(np.radians(-grid_rotang))
        sin_theta = np.sin(np.radians(-grid_rotang))
        tan_phi_x = np.tan(np.radians(grid_skew_x))
        tan_phi_y = np.tan(np.radians(grid_skew_y))

        # Compute combined rotation and skew transformation matrix
        transformation_matrix = np.array([[cos_theta - sin_theta*tan_phi_y, cos_theta*tan_phi_x - sin_theta],
                                          [sin_theta + cos_theta*tan_phi_y, sin_theta*tan_phi_x + cos_theta]])
        
        # Generate points in the local coordinate system
        points = []
        moorings = []
        labels_list = [] # list with grid labels, so that each point now to which horizontal or vertical line it belongs
        break_flag = False
        # LOCAL COORDINATE SYSTEM WITH (0,0) LEASE AREA CENTROID
        # Therefore, +/- self.boundary_centroid_y/x cover the entire area
        # Loop through y values within the boundary_centroid_y range with grid_spacing_y increments
        iy = 0
        
        ywidth = np.max(bound_ys) - np.min(bound_ys)
        xwidth = np.max(bound_xs) - np.min(bound_xs)
        if center==None:
            ycenter = (np.max(bound_ys) + np.min(bound_ys))/2
            xcenter = (np.max(bound_xs) + np.min(bound_xs))/2
        else:
            xcenter = center[0]
            ycenter = center[1]
        
        for y in np.arange(np.min(bound_ys) - ywidth*1.0, np.max(bound_ys) + ywidth*1.0, spacing_y):    # extending by 1.0*width in x and y to make sure rotations include everything
            # Loop through x values within the boundary_centroid_x range with grid_spacing_x increments
            ix = 0
            for x in np.arange(np.min(bound_xs) - xwidth*1.0, np.max(bound_xs) + xwidth*1.0, spacing_x):
                # Apply transformation matrix to x, y coordinates
                local_x, local_y = np.dot(transformation_matrix, [x - xcenter, y - ycenter])
                # Add grid translation offsets to local coordinates
                local_x += grid_trans_x
                local_y += grid_trans_y
                
                
                # Create a Point object representing the transformed coordinates
                # Transform back into global coordinate system with by adding centroid to local coordinates
                #point = Point(local_x + np.min(bound_xs), local_y + np.min(bound_ys))
                point = Point(local_x + xcenter, local_y + ycenter)
                
                if fullMPsystem:
                    
                    if ms == None:
                        raise ValueError('NEED TO INPUT MOORPY SYSTEM')
                         
                    # Check if the point lies within the specified shape (boundary_sh_int)
                    
                    # deep copy of mooring system to apply translations and rotation
                    mss = deepcopy(ms)
                    
                    # select every other column for rotation and add to farm rotation (mooring rotation is relative to y')
                    rot = rotations[iy % 2] + grid_rotang
                    
                    mss.transform(trans = [point.x, point.y], rot = -rot) #moorpy rotation convention is opposite
                    mss.initialize()
                    mss.solveEquilibrium()
                    
                    contained = True
                    for l in mss.lineList:
                      anchor = Point(l.rA[0], l.rA[1])
                      
                      if not area.contains_properly(anchor):
                          contained = False
                
                else:
                    if rad_anch == None:
                        raise ValueError('NEED TO INPUT RAD_ANCH')
                    buff = point.buffer(rad_anch)                    
                    contained = True
                        
                    if not area.contains_properly(buff):
                        contained = False
                    
                if contained:
                    # If the point is within the shape, append it to the list of points
                    points.append(point)
                    
                    if fullMPsystem:
                        moorings.append(mss)
                    # Save grid label
                    labels_list.append([ix,iy-1]) # y -1 so that labels are again starting at 0
                    # If the number of points collected reaches the desired threshold (nt), set break_flag to True and exit the loop
                    # if len(points) >= self.nt:
                    #     break_flag = True
                    #     break
                ix += 1    
            iy += 1
            # If break_flag is True, exit the outer loop as well
            if break_flag:
                break
             
        x_coords = np.array([point.x for point in points])#/1000 
        y_coords = np.array([point.y for point in points])#/1000   

        return(x_coords, y_coords, moorings, area)



def getUpdatedAnchorPosition(old_anchor_point, fairlead_point, grid_x, grid_y, grid_depth, ratio=1000):
    '''
    Compute a new anchor position for a taut mooring line by looking along the
    a line from old anchor to fairlead and seeing where it intersects the seabed.
    
    Paramaters
    ----------
    old_anchor_point : float, array
        list of a xyz coordinate of an anchor point
    fairlead_point : float, array
        list of a xyz coordinate of a fairlead point
    grid_x : float, array
        The x coordinates of coordinate pairs defining a bathymetry grid
    grid_y : float, array
        The y coordinates of coordinate pairs defining a bathymetry grid
    grid_depth: float, 2D matrix
        The depth (z coordinates) of the grid defined by grid_x and grid_y
    ratio: int or float (optional)
        the value of how far to extend a mooring line until it intersects the bathymetry grid plane
    
    Returns
    -------
    new_anchor_point : float, array
        list of a xyz coordinate of the updated anchor point so that it intersects the local bathymetry grid plane
    '''
    
    # calculate the actual depth based on bathymetry of the x/y coordinates of the anchor
    x = old_anchor_point[0]
    y = old_anchor_point[1]
    #depth, nvec, ix0, iy0 = getDepthFromBathymetry(x, y, grid_x, grid_y, grid_depth)  # needed to adjust gDFB function <<<<<< can change later
    depth, nvec = getDepthFromBathymetry(x, y, grid_x, grid_y, grid_depth)  # needed to adjust gDFB function <<<<<< can change later
    
    # create points of a line that connect the fairlead to the anchor
    p0 = fairlead_point
    p1 = old_anchor_point
    
    '''
    # but adjust the "anchor" point to way below the bathymetry, if it is found that the initial anchor position is above the bathymetry
    if p1[2] > -depth:
        p1 = np.array([p0[0]+ratio*(p1[0]-p0[0]), p0[1]+ratio*(p1[1]-p0[1]), p0[2]+ratio*(p1[2]-p0[2])])
    '''
    
    # Find the intersection point between the mooring Line (assumed straight)
    # and the bathymetry grid panel 
    u = p1 - p0     # vector from fairlead to original anchor
    w = np.array([x, y, -depth]) - p0  # vector from fairlead to a point on the grid panel of the original anchor
    
    fac = np.dot(nvec, w) / np.dot(nvec, u)  # fraction along u where it crosses the seabed (can be greater than 1)
    
    new_anchor_point = p0 + u*fac

    return new_anchor_point



"""
def getInterpNums(xlist, xin, istart=0):  # should turn into function in helpers
    '''
    Paramaters
    ----------
    xlist : array
        list of x values
    xin : float
        x value to be interpolated
    istart : int (optional)
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

def getDepthFromBathymetry(x, y, bathGrid_Xs, bathGrid_Ys, bathGrid, point_on_plane=False):   #BathymetryGrid, BathGrid_Xs, BathGrid_Ys, LineX, LineY, depth, nvec)
    ''' interpolates local seabed depth and normal vector
    
    Parameters
    ----------
    x, y : float
        x and y coordinates to find depth and slope at [m]
    bathGrid_Xs, bathGrid_Ys: float, array
        The x and y coordinates defining a bathymetry grid
    bathGrid: float, 2D matrix
        The depth (z coordinates) of the grid defined by bathGrid_Xs and bathGrid_Ys
    point_on_plane: bool (optional):
        determines whether to return the indices that go with the bathGrid arrays to return a point on the bathymetry grid plane
    
    Returns
    -------        
    depth : float
        local seabed depth (positive down) [m]
    nvec : array of size 3
        local seabed surface normal vector (positive out) 
    ix0 : int
        index of the point on the bathymetry grid plane that goes with bathGrid_Xs
    iy0 : int
        index of the point on the bathymetry grid plane that goes with bathGrid_Xs
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

    if not point_on_plane:
        return depth, nvec
    else:
        return depth, nvec, ix0, iy0
"""



def adjustMS4Bath(ms, ms_xy, grid_x, grid_y, grid_depth, iLine=-1, nLines_in_ms=3, nLines_in_line=3, display=0, extend=True):
    '''Function that updates a MoorPy System object's anchor positions in response to bathymetry and then updates
    the line lengths to keep the same pretension that was there before the bathymetry adjustments

    Parameters
    ----------
    ms : MoorPy System object
        A MoorPy System object defining a mooring system
    ms_xy: float, array
        The 2D x/y position of the system's coordinate system relative to a reference point (e.g., a centroid)
        to reference the proper bathymetry location
    grid_x : float, array
        The x coordinates of coordinate pairs defining a bathymetry grid
    grid_y : float, array
        The y coordinates of coordinate pairs defining a bathymetry grid
    grid_depth: float, 2D matrix
        The depth (z coordinates) of the grid defined by grid_x and grid_y
    iLine: int, optional
        the index of the line object that is to be adjusted (among the indices of line objects from one anchor to one fairlead, like a subsystem)
    nLines: int, optional
        the number of mooring lines that surround the MoorPy Body
    display: int, optional
        an option for print statement outputting
    extend: boolean, optional
        True for updating anchor positions to bathymetry along the vector of the mooring line, or False for dropping/lifting the anchor at the same x/y position

    Returns
    -------
    ms: MoorPy System object
        The updated MoorPy System object with new anchor positions and line lengths that match the initial pretensions
    '''

    # NOTE: This function can probably be put in system.py as a method since it adjusts a System object
    if np.any([isinstance(line, mp.Subsystem) for line in ms.lineList]):
        subsystem_flag = True
    else:
        subsystem_flag = False

    ### COLLECT INFORMATION ABOUT THE INPUT MOORING SYSTEM (OR SUBSYSTEMS(S)) ###

    T_init_list = [np.linalg.norm(ss.fB_L[0:2]) for ss in ms.lineList]

    # collect point numbers for all anchor points
    anchors = [point.number for point in ms.pointList if point.type==1 and point.number not in ms.bodyList[0].attachedP]
    # collect point numbers for all anchor points if there are any subsystems in the lineList
    #anchors_subsystem = [anchornum for anchornum in anchors for line in ms.lineList for linenum in ms.pointList[anchornum-1].attached if isinstance(line, mp.Subsystem) and linenum==line.number]
    # split anchor point numbers up based on whether they are attached to a subsystem or just a line
    #anchors_lines = list(set(anchors).difference(anchors_subsystem))
    
    # collect point numbers for all "fairleads"
    # (fairleads are defined as the points attached to the body where "upper_points" are the points that the lines that are to be adjusted are attached to at the top)
    if not subsystem_flag:
        iLines = np.arange(iLine, 1e3, nLines_in_line, dtype=int)[:nLines_in_ms]      # create a list of the indices of all lines in a mooring system to vary (doesn't always need to be line connected to the fairlead)
        upper_points = np.sort([point.number for point in ms.pointList for iL in iLines if all(point.r==ms.lineList[iL].rB)])   # collect the numbers of the points where the lines of interest are attached to at the top
    
    fairleads = [point.number for point in ms.pointList if point.type==1 and point.number in ms.bodyList[0].attachedP]      # collect the numbers of the points that are fairleads
    # collect point numbers for all "upper_points" if there are any subsystems in the list
    #upper_points_subsystem = [fairleadnum for fairleadnum in upper_points for line in ms.lineList for linenum in ms.pointList[fairleadnum-1].attached if isinstance(line, mp.Subsystem) and linenum==line.number]
    # split the upper_points list up based on whether they are attached to a subsystem or not
    #upper_points_lines = np.sort(list(set(upper_points).difference(upper_points_subsystem)))

    if not subsystem_flag:
        # collect line numbers that are attached to the points of interest
        #lower_lines = [line.number for line in ms.lineList if not isinstance(line, mp.Subsystem) for point in ms.pointList if point.number in anchors if all(line.rA==point.r)]
        upper_lines = [line.number for line in ms.lineList if not isinstance(line, mp.Subsystem) for point in ms.pointList if point.number in upper_points if all(line.rB==point.r)]
        fairleads_lines = [line.number for line in ms.lineList if not isinstance(line, mp.Subsystem) for point in ms.pointList if point.number in fairleads if all(line.rB==point.r)]
        # collect the upper tensions of each line attached to the points of interest
        upper_lines_TB = [ms.lineList[linenum-1].TB for linenum in upper_lines]
        fairleads_lines_TB = [ms.lineList[linenum-1].TB for linenum in fairleads_lines]

    # separate the subsystem objects from the rest to use later, separately from the Line objects
    subsystems = [line for line in ms.lineList if isinstance(line, mp.Subsystem)]

    ### CALCULATE AND SET NEW ANCHOR POSITIONS FOR ONLY LINE OBJECTS ###
    for i,anchornum in enumerate(anchors):
        anchor_point_local = ms.pointList[anchornum-1].r
        anchor_point_global = anchor_point_local + np.array([ms_xy[0], ms_xy[1], 0])
        fairlead_point_local = ms.pointList[fairleads[i]-1].r
        fairlead_point_global = fairlead_point_local + np.array([ms_xy[0], ms_xy[1], 0])

        if extend:  # if you wish to "extend" or "retract" the anchor point along the vector of the mooring line
            new_anchor_point_global = getUpdatedAnchorPosition(anchor_point_global, fairlead_point_global, grid_x, grid_y, grid_depth)
            if new_anchor_point_global[2] < anchor_point_global[2]:
                if display > 0: print("'Extending' the anchor point to the bathymetry, along the vector of the mooring line")
            elif new_anchor_point_global[2] > anchor_point_global[2]:
                if display > 0: print("'Retracting' the anchor point to the bathymetry, along the vector of the mooring line")
            else:
                if display > 0: print("No change in the anchor depth")
        else:       # if you wish to "drop" or "lift" the anchor point at the same x/y position
            new_depth, _ = ms.getDepthFromBathymetry(anchor_point_global[0], anchor_point_global[1])
            new_anchor_point_global = np.array([anchor_point_global[0], anchor_point_global[1], -new_depth])
            if new_anchor_point_global[2] < anchor_point_global[2]:
                if display > 0: print("'Dropping' the anchor point to the bathymetry, at the same x/y position")
            elif new_anchor_point_global[2] > anchor_point_global[2]:
                if display > 0: print("'Lifting' the anchor point to the bathymetry, at the same x/y position")
            else:
                if display > 0: print("No change in the anchor depth")
        
        new_anchor_point_local = new_anchor_point_global - np.array([ms_xy[0], ms_xy[1], 0])
        ms.pointList[anchornum-1].setPosition(new_anchor_point_local)
        # setPosition sets the point.r value to the input, and also updates the end position of the line object
        # setPosition also doesn't allow the input position to be less than ms.depth (which shouldn't matter if the input ms to this function is already at seabedMod=2)
        if subsystem_flag:
            ms.lineList[i].setEndPosition(new_anchor_point_local, endB=0)
            ms.lineList[i].depth = -new_anchor_point_local[2]
            ms.lineList[i].pointList[0].setPosition(np.array([ms.lineList[i].pointList[0].r[0], 0, new_anchor_point_local[2]]))

    # resolve for equilibrium
    ms.solveEquilibrium()

    if subsystem_flag: 
        for i,ss in enumerate(subsystems):
            L = adjustSS4Pretension(ss, i_line=1, T_init=T_init_list[i], horizontal=True, display=3, tol=0.001) 
            ss.lineList[1].setL(L)
            ms.solveEquilibrium()

    else:
        ## update line lengths to match pretension ##
        def eval_func(X, args):
            '''Tension evaluation function for different line lengths'''
            L = X[0]                                # extract the solver variable
            # set args variables
            ms_copy = args['ms']        #ms_copy = deepcopy(args['ms'])
            iLineX = args['iLineX']
            iLineFair = args['iLineFair']
            # set System variables and solve for new tension
            ms_copy.lineList[iLineX].L = L
            ms_copy.solveEquilibrium()
            T = np.linalg.norm(ms_copy.lineList[iLineFair].fB)
            return np.array([T]), dict(status=1), False

        #upper_lines_byL = [ upper_lines[i] for i in np.flip(np.argsort([ms.lineList[ul-1].L for ul in upper_lines])) ]
        #fairleads_lines_byL = [ fairleads_lines[i] for i in np.flip(np.argsort([ms.lineList[ul-1].L for ul in upper_lines])) ]

        # loop through the upper lines and run dsolve2 solver to solve for the line length that matches that initial tension
        for i,upper_linenum in enumerate(upper_lines):
            # set initial variables
            T_init = fairleads_lines_TB[i]      #T_init = upper_lines_TB[i]
            EA = ms.lineList[upper_linenum-1].type['EA']
            L_init = ms.lineList[upper_linenum-1].L
            X0 = [ 10 ]                         #X0 = [ L_init/(T_init/EA+1) ]      # setting to start at 10 to start from really taut and extend to longer, always
            if display > 0: print(f"    Updating Line {upper_linenum} length to match pretension")
            
            # run dsolve2 to solve for the line length that produces the same initial tension
            L_final, T_final, _ = dsolve2(eval_func, X0, Ytarget=[T_init], args=dict(ms=ms, iLineX=upper_linenum-1, iLineFair=fairleads_lines[i]-1), maxIter=200, stepfac=4, display=display, tol=1e-4)
            # has the option to solve for an intermediate line length that results in the same tension on the fairlead line (different line)

            # set the new line length into the ms System
            ms.lineList[upper_linenum-1].L = L_final[0]
            if display > 0: print(f'    L0 = {X0[0]:6.1f}, LF = {L_final[0]:6.1f}')
            if display > 0: print(f'    T0 = {T_init:8.2e}, TF = {T_final[0]:8.2e}')
    
    
    ms.solveEquilibrium()
    
    return ms


def adjustSS4Pretension(ssin, i_line=0, T_init=None, horizontal=False, display=0, stepfac=10, tol=0.01):

    ss = deepcopy(ssin)

    if T_init==None:
        if horizontal:
            T_init = np.linalg.norm(ss.fB_L[0:2])
        else:
            T_init = ss.TB  # save the initial pretension
    
    # can update the subsystem initially if need be (Subsystem.staticSolve is the equivalent to System.solveEquilibrium)
    #ss.staticSolve()
    #T0 = ss.TB
    
    # update line lengths to match pretension
    def eval_func(X, args):
        '''Tension evaluation function for different line lengths'''
        L = X[0]
        ss.lineList[i_line].L = L
        ss.staticSolve()
        if horizontal:
            T = np.linalg.norm(ss.fB_L[0:2])
        else:
            T = ss.TB 
        return np.array([T]), dict(status=1), False

    # run dsolve2 solver to solve for the upper line length that matches the initial tension
    L_init = ss.lineList[i_line].L
    if display > 0: print(f"    Updating Subsystem {ss.number}'s Line {ss.lineList[i_line].number} length to match pretension")
    
    #X0 = [L_init]
    X0 = [10]
    
    # run dsolve2 to solve for the line length that sets the same pretension
    L_final, T_final, _ = dsolve2(eval_func, X0, Ytarget=[T_init], 
                                  Xmin=[1], Xmax=[1.1*np.linalg.norm(ss.rB-ss.rA)],
                                  dX_last=[1], tol=[tol], 
                                  maxIter=200, stepfac=stepfac, display=display)
    
    #ss.lineList[i_line].setL(L_final[0])  # assign the solved_for length
    if display > 0: print(f'    L_init = {L_init:6.1f}, LF = {L_final[0]:6.1f}')
    if display > 0: print(f'    T_init = {T_init:8.2e}, TF = {T_final[0]:8.2e}')
    
    #ss.staticSolve()    # reset the subsystem

    return L_final[0]
    #return ss





def adjustMooring(mooring, layout, r_fair, r_anch, adjust={}):
    '''Adjust a Mooring object for change in layout considering the seabed,
    which is contained in Project object that is also passed in.
    The Mooring adjustment should work regardless of whether the mooring
    is only 2D or also includes a 3D representation via MoorPy Subsystem.
    
    When a subsystem is involved, a dictionary can be past via 'adjust'
    to ask for the pretension to be adjusted to a desired value.
    
    Parameters
    ----------
    mooring : Mooring object
        The Mooring to be adjusted.
    layout : Layout object
        An object of the Layout class that contains seabed information.
    r_fair : float, array
        Absolute xyz coordinate of a fairlead point [m].
    r_anch : float, array
        Absolute xyz coordinate of the anchor point (guess to be adjusted) [m].
    adjust : dict
        Dictionary specifying a method of adjusting the mooring to maintain a
        desired characteristic. Currently only pretension is supported:
        {'pretension' : {'target' : XXX N, 'i': line index to adjust}}.
    '''
    
    # Update the anchor position if it isn't already on the seabed
    if not np.isclose(r_anch[2], layout.getDepthAtLocation(*r_anch[:2])):
        r_anch = layout.getUpdatedAnchorPosition(r_fair, r_anch)
    
    # Set the mooring end positions (this will update any Subsystem too)
    mooring.setEndPosition(r_anch, 'a')
    mooring.setEndPosition(r_fair, 'b')

    # If requested, update the line lengths to maintain pretension
    if mooring.subsystem and 'pretension' in adjust:
        target = adjust['pretension']['target']
        i_line = int(adjust['pretension']['i'])

        ss = mooring.subsystem  # shorthand
        
        # update line lengths to match pretension
        def eval_func(X, args):
            '''Tension evaluation function for different line lengths'''
            ss.lineList[i_line].L = X[0]  # set specified line section's length
            ss.staticSolve(tol=0.0001)    # solve the equilibrium of the subsystem
            return np.array([ss.TB]), dict(status=1), False   # return the end tension

        # run dsolve2 solver to solve for the upper line length that matches the initial tension
        X0 = [ss.lineList[i_line].L]  # initial value is current section length
        
        L_final, T_final, _ = dsolve2(eval_func, X0, Ytarget=[target], 
                                      Xmin=[1], Xmax=[1.1*np.linalg.norm(ss.rB-ss.rA)],
                                      dX_last=[1], tol=[0.1], 
                                      maxIter=200, stepfac=4, display=2)
        
        ss.lineList[i_line].L = L_final[0]
        
        print(f"Adjusted mooring to pretension of {T_final[0]:.0f} N and {L_final[0]:.2f}-m section length.")


def makeMooringListN(subsystem0, N):
    '''Simple function for making a mooringList of N mooring objects, by
    by duplication one provided subsystem. They can be positioned later.
    '''
    
    #mooringList = []
    # Initialize empty list
    #mooringList = [None] * N
    mooringList = {}
    
    for i in range(N):
    
        # Make a copy from the original
        ss = deepcopy(subsystem0)  
                      
        #mooringList.append(Mooring(subsystem=ss,id=i))
        mooringList[i] = Mooring(subsystem=ss,id=i) 
        #mooringList[i].rA = ss.rA
        #mooringList[i].rB = ss.rB
        # Make a new mooring object to hold the copied subsystem
        #mooringList.append(Mooring(subsystem=ss, rA=ss.rA, rB=ss.rB))

    return mooringList


def getLower(A):
    '''Return a vector of the serialized lower-triangular elements of matrix A'''
    return A[np.tril_indices_from(A , k=-1)]

if __name__ == '__main__':
    

    # initialize the area bounds
    lease_xs = np.array([ 2220.61790941,  2220.61787966,  3420.61793096,  3420.61791961,
                         3420.61801382,  3420.61803977,  3420.61807596,  3420.61811471,
                         3420.61822821,  4620.61806679,  4620.61816982,  4620.6182304 ,
                         4620.61832506,  4620.61827356,  4620.61840937,  4620.61846802,
                         5820.61842928,  5820.61837593,  5820.61846352,  5820.61851395,
                         5820.61860812,  7020.61852522,  7020.61862679,  7020.61856713,
                         5820.61872489,  5820.61881681,  5820.61883604,  4620.61889509,
                         4620.61904038,  4620.61903256,  3420.619089  ,  3420.61917709,
                         3420.61920431,  2220.61926344,  2220.61944036,  2220.61939526,
                         1020.61946852,  -179.38041543,  -179.38053362,  -179.38049454,
                         -179.38063026,  -179.38059316,  -179.38071736,  -179.38082411,
                         -179.38081907,  -179.38086534,  -179.38087289,  -179.38099606,
                         -179.38098055, -1379.38097106, -1379.38100559, -1379.38111504,
                         -2579.38099329, -2579.38099701, -2579.38112028, -2579.38119101,
                         -3779.38113702, -3779.38108251, -3779.38118818, -4979.38109205,
                         -4979.38117059, -4979.38120063, -6179.38111446, -6179.38120826,
                         -7379.38107666, -7379.38119559, -8579.38111169, -8579.38121027,
                         -7379.38119975, -7379.38133324, -7379.38133537, -7379.38136608,
                         -7379.3815011 , -7379.38148217, -6179.38158388, -6179.38169844,
                         -4979.38166381, -4979.38179166, -3779.38188846, -2579.38201223,
                         -1379.38197896, -1379.38206961,  -179.38223745,  -179.38215968,
                         1020.61770826,  2220.61761805,  2220.61770748,  2220.61772366,
                         2220.6178225 ,  2220.61783243,  2220.61790941])

    lease_ys = np.array([ 10997.93747912,   9797.9374786 ,   9797.93762135,   8597.93777362,
                     7397.93779263,   6197.93779256,   4997.9378827 ,   3797.93806299,
                     2597.93810968,   2597.93821722,   1397.93821124,    197.93829277,
                     -1002.06153985,  -2202.06150557,  -3402.0613888 ,  -4602.0614069 ,
                     -4602.06132321,  -5802.06115355,  -7002.06101498,  -8202.06101547,
                     -9402.06093466,  -9402.06087921, -10602.06081225, -11802.06066334,
                     -11802.06074832, -13002.06064456, -14202.06056929, -14202.06071115,
                     -15402.06052442, -16602.06047346, -16602.0605191 , -17802.06042198,
                     -19002.060349  , -19002.06047655, -20202.06042422, -21402.06027959,
                     -21402.06033391, -21402.06033968, -20202.06050384, -19002.06057741,
                     -17802.06066473, -16602.06066129, -15402.06078244, -14202.06092037,
                     -13002.06096561, -11802.06102581, -10602.06121353,  -9402.06119237,
                     -8202.0613002 ,  -8202.06134893,  -7002.06145163,  -5802.06156475,
                     -5802.06157941,  -4602.06171868,  -3402.06175304,  -2202.06190779,
                     -2202.0618502 ,  -1002.06194537,    197.93795779,    197.93795451,
                     1397.93778023,   2597.93771457,   2597.93761946,   3797.93767253,
                     3797.93747917,   4997.93739719,   4997.93747426,   6197.9373365 ,
                     6197.93731811,   7397.93724719,   8597.93717933,   9797.93711807,
                     10997.93706522,  12197.93690464,  12197.93699804,  13397.93696353,
                     13397.93699185,  14597.93692719,  14597.93693641,  14597.93707605,
                     14597.93701544,  15797.93706205,  15797.93707378,  16997.93692194,
                     16997.93704183,  16997.93701933,  15797.93709866,  14597.93727451,
                     13397.93732149,  12197.93746482,  10997.93747912])
    
    # initialize dummy mooring system to use to organize turbines within a layout
    ms_type = 1
    
    if ms_type==1:
        ms = mp.System(file='sample_deep.txt')

    elif ms_type==2:
        depth     = 200                             # water depth [m]
        angles    = np.radians([60, 300])      # line headings list [rad]
        rAnchor   = 600                            # anchor radius/spacing [m]
        zFair     = -21                             # fairlead z elevation [m]
        rFair     = 20                              # fairlead radius [m]
        lineLength= 650                            # line unstretched length [m]
        typeName  = "chain1"                        # identifier string for the line type

        ms = mp.System(depth=depth)

        # add a line type
        ms.setLineType(dnommm=120, material='chain', name=typeName)  # this would be 120 mm chain

        # Add a free, body at [0,0,0] to the system (including some properties to make it hydrostatically stiff)
        ms.addBody(1, np.zeros(6), m=1e6, v=1e3, rM=100, AWP=1e3)

        # For each line heading, set the anchor point, the fairlead point, and the line itself
        for i, angle in enumerate(angles):

            # create end Points for the line
            ms.addPoint(1, [rAnchor*np.cos(angle), rAnchor*np.sin(angle), -depth])   # create anchor point (type 0, fixed)
            ms.addPoint(1, [  rFair*np.cos(angle),   rFair*np.sin(angle),  zFair])   # create fairlead point (type 0, fixed)
            
            # attach the fairlead Point to the Body (so it's fixed to the Body rather than the ground)
            ms.bodyList[0].attachPoint(2*i+2, [rFair*np.cos(angle), rFair*np.sin(angle), zFair]) 

            # add a Line going between the anchor and fairlead Points
            ms.addLine(lineLength, typeName, pointA=2*i+1, pointB=2*i+2)
        
        # ----- Now add a SubSystem line! -----
        ss = mp.Subsystem(mooringSys=ms, depth=depth, spacing=rAnchor, rBfair=[10,0,-20])

        # set up the line types
        ms.setLineType(180, 'chain', name='one')
        ms.setLineType( 50, 'chain', name='two')

        # set up the lines and points and stuff
        ls = [350, 300]
        ts = ['one', 'two']
        ss.makeGeneric(lengths=ls, types=ts)

        # add points that the subSystem will attach to...
        ms.addPoint(1, [-rAnchor, 100, -depth])   # Point 5 - create anchor point (type 0, fixed)
        ms.addPoint(1, [ -rFair ,   0,  zFair])   # Point 6 - create fairlead point (type 0, fixed)
        ms.bodyList[0].attachPoint(6, [-rFair, 0, zFair])  # attach the fairlead Point to the Body 

        #ms.addLine(length, type, pointA, pointB)

        # string the Subsystem between the two points!
        ms.lineList.append(ss)  # add the SubSystem to the System's lineList
        ss.number = 3
        ms.pointList[4].attachLine(3, 0)  # attach it to the respective points
        ms.pointList[5].attachLine(3, 1)  # attach it to the respective points

        #ms.bodyList[0].type = 1

    

    ms.initialize()                                             # make sure everything's connected
    ms.solveEquilibrium()                                       # equilibrate
    





    # create a subsystem
    ss = mp.Subsystem(depth=2000, spacing=2000, rBfair=[10,0,-20])

    # set up the line types
    ss.setLineType(180, 'polyester', name='one')

    # set up the lines and points and stuff
    lengths = [2000]
    types = ['one']
    ss.makeGeneric(lengths, types)

    # plotting examples
    ss.setEndPosition([-2000  , 0,-2000], endB=0)
    ss.setEndPosition([-10, 0,  -20], endB=1)
    ss.staticSolve()

    #ss.pointList[0].setPosition(np.array([-2000, 0, -2000]))
    #ss.pointList[-1].setPosition(np.array([-10, 0, -20]))
    #ss.solveEquilibrium()






    # initialize dummy bathymetry variables to check anchor depths
    grid_x = np.array([-65000, 65000])
    grid_y = np.array([-65000, 65000])
    grid_depth = np.array([[2367, -211],
                            [3111, -338]])
    
    start_time = time.time()

    coords, mooringList, footprintList = create_layout(lease_xs, lease_ys, ss, grid_x, grid_y, grid_depth,
                                                       spacing_x=8400, spacing_y=8600)
    
    end_time = time.time() - start_time
    print(end_time, end_time/60)
    # create a layout of turbine positions
    #xs, ys, footprintList, msList = create_initial_layout(lease_xs, lease_ys, ms, grid_x, grid_y, grid_depth, display=1)

    # plot the result
    fig, ax = plt.subplots(1,1)
    ax.plot(coords[:,0], coords[:,1], color='k', marker='o', linestyle='')
    ax.plot(lease_xs, lease_ys, color='r')
    for polygon in footprintList:
        x, y = polygon.exterior.coords.xy
        ax.plot(x, y, color='b', alpha=0.5)
    ax.set_aspect('equal')


    # try making a Project with the above
    
    project = Project()
    
    project.loadBathymetry('bathymetry_sample.txt')

    project.mooringList = mooringList

    project.plot3d()

    plt.show()


    a = 2




    # Next Steps:
    # - be able to adjust the starting point to see if there are any other arrangements that can fit more turbines
    # - be able to adjust the "bounds" on each turbine (i.e., change from a circle around each turbine to maybe a triangle)
    # - be able to account for bathymetry for each anchor point (will likely need FAModel integration, as this is already set up to do this)



def convertm2km(coords):
    above_1000 = any(coord > 1000 for coord in coords)
    if above_1000:
        coords = [coord / 1000 for coord in coords]
    return coords











