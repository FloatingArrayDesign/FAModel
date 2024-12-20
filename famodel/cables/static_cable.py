# class for a static subsea power cable

import numpy as np
from copy import deepcopy
from moorpy.subsystem import Subsystem
from moorpy import helpers
from famodel.mooring.connector import Connector, Section
from famodel.famodel_base import Edge

class StaticCable(Edge):
    '''
    Class for a static power cable. It inherits from Edge.
    A StaticCable object will likely be within a SubseaCable object, which could
    also include a DynamicCable.
    Both ends will be attached to subsea joints.
    The main thing done by this class is describing the cable routing on the
    seabed.
    '''
    
    def __init__(self, id, dd=None, subsystem=None, anchor=None, rA=[0,0,0], rB=[0,0,0],
                 rad_anch=500, rad_fair=58, z_anch=-100, z_fair=-14, 
                 rho=1025, g=9.81, L=2000,A=0,conductorSize=None, 
                 type='static',voltage=66,powerRating=None,cable_type=None,routing=[],
                 burial=None,zAnchor=None,**kwargs):
        '''
        Parameters
        ----------
 
        '''
        Edge.__init__(self, id)  # initialize Edge base class
        # Design description dictionary for this dynamic cable
        self.dd = dd
        
        self.n_sec = 1
        
        # Store the cable type properties dict here for easy access (temporary - may be an inconsistent coding choice)
        self.cableType = self.makeCableType(self.dd['cable_type'])  # Process/check it into a new dict

        # cable length
        self.L = L

        # end point absolute coordinates, to be set later
        self.rA = rA
        self.rB = rB
        
        # relevant site info
        self.rho = rho
        self.g = g
        
        # Vectors of optional vertex points along the cable route. 
        # Nonzero radius wraps around the vertex at that radius.
        self.coordinates = routing
        if routing:
            # call function to update routing variables
            self.updateRouting()
        
        # Dictionaries for addition information
        self.loads = {}
        self.reliability = {}
        self.cost = {}
        self.failure_probability = {}
    
    
    def getLength(self):
        '''Compute the length of the cable based on the end point locations 
        (rA, rB) and any routing information (self.x, self.y. self.rad).'''
        # get some calculations from Stein's arc-shaped work on IProTech geometry
        # >>> calcs here <<<
        def get_length_between_points(coord1,coord2):
            valx = coord1[0] - coord2[0]
            valy = coord1[1] - coord2[1]
            return(np.hypot(valx,valy))
        
        L = 0
        
        if not hasattr(self,'x') or not self.x: 
            # there is not routing here, assume straight line between ends
            L += get_length_between_points(self.rB,self.rA)
            
        else:
            # get L from rA to first point
            L += get_length_between_points([self.x[0],self.y[0]],self.rA[0:2])
            # get L from rB to last point
            L += get_length_between_points([self.x[-1],self.y[-1]],self.rB[0:2])
            
            # get L between coordinate points
            if len(self.x)>1:
                for i in range(len(self.x)-1):                
                    L += get_length_between_points([self.x[i+1],self.y[i+1]],
                                                   [self.x[i],self.y[i]])
                
        self.L = L                            
        
        return self.L
    
    def makeCableType(self,cabDict):
        '''
        Processes dictionary info to make cableType dictionary

        Parameters
        ----------
        cabDict : dict
            Dictionary of information on cable type to be processed

        Returns
        -------
        cableType : dict
            Uniform cable dictionary

        '''
        # fix later, for now just use cabDict
        cableType=cabDict
        return(cableType)
    
    
    def initializeRoute(self):
        # Initially assume that the static portion of the cable goes straight
        # between the two ends of the dynamic cables.
        
        '''  >>> to be updated >>>
        # x y coordinates of the FOWTs that the cable runs between
        rtA = self.system.coords[self.AttachA-1]
        rtB = self.system.coords[self.AttachB-1]
        
        # set cable end points
        self.rA[0] = rtA[0] + self.DynCableA.span*np.cos(self.headingA)  # x
        self.rA[1] = rtA[1] + self.DynCableA.span*np.sin(self.headingA)  # y
        self.rB[0] = rtB[0] + self.DynCableB.span*np.cos(self.headingB)  # x
        self.rB[1] = rtB[1] + self.DynCableB.span*np.sin(self.headingB)  # y
        '''
        
        # fill in x, y, r lists for just straight connection between end points
        self.x = np.array([self.rA[0], self.rB[0]])
        self.y = np.array([self.rA[1], self.rB[1]])
        self.r = np.array([0, 0])
        
        # (fancier routing could be applied later by layout design methods)


    def resolveRoute(self):
        '''Takes established cable route points, considers seabed
        bathymetry and embedment depth, and computes points along
        the path along with its length.'''
        
        # coordinates for plotting
        self.xs = []
        self.ys = []
        
        # figure out x and y coordinates
        # add points along each stretch
        for i in range(len(self.x))-1:
            n = 3 # number of segments along a stretch
            np.append(self.xs, np.linspace(self.x[i], self.x[i+1], n))
            np.append(self.ys, np.linspace(self.y[i], self.y[i+1], n))
        # add last point
        np.append(self.xs, self.x[-1])
        np.append(self.ys, self.y[-1])
        # >>> Stein to add radii around points <<<
        
        # get z coordinates along seabed
        self.zs = self.system.projectAlongSeabed(self.xs, self.ys)
        
        # calculate cable length (discretized approach for now)
    
    def updateRouting(self,coords=None):
        '''
        Updates routing parameters based on passed in coordinates or existing coordinates

        Parameters
        ----------
        coords : list, optional
            List of xy coordinates (and potentially burial and radius) for routing

        Returns
        -------
        None.

        '''
        if coords:
            self.coordinates = coords
           
        if self.coordinates:
            self.x = [coord[0] for coord in self.coordinates]  # cable route vertex global x coordinate [m]
            self.y = [coord[1] for coord in self.coordinates]  # cable route vertex global y coordinate [m]
            # Check if radius available
            if len(self.coordinates[0]) <= 2:
                self.r = []  # cable route vertex corner radius [m]
                self.burial = []
            elif len(self.coordinates[0]) == 3:
                self.burial = [coord[2] for coord in self.coordinates] # cable route depth at coordinates
                self.r = []
            elif len(self.coordinates[0]) == 4:
                self.burial = [coord[2] for coord in self.coordinates] # cable route depth at coordinates
                self.r = [coord[3] for coord in self.coordinates]  # cable route vertex corner radius [m]
                
        # update static cable length
        self.getLength()
        # update overall cable length
        overall_cable = self.part_of
        if overall_cable:
            overall_cable.getL()
        
            
        
    
    def checkCableExclusions(self, cable):
        '''Checks whether a cable crosses over any exclusions
        or other out of bounds areas.
        '''

        # select cable
        
        # check its path against any exclusion areas or boundaries
        
        # make a list of any exclusion/nodes that it is too close to
        
        return 0 #score, list_of_violations
    
    
    
    def getCost(self):
        
        # >>> UPDATE <<<
        self.getLength()
        cost = 0
        if 'cost' in self.dd['cable_type']:
            cost += self.dd['cable_type']['cost']*self.L
            
        self.cost['cable'] = cost
        # sum up the costs in the dictionary and return
        return sum(self.cost.values()) 
    
    
    # maybe cable routing methods would go here (currently drafted in cable.py)
