# class for a mooring line in a floating array

import numpy as np

class Mooring():
    '''
    Class for a floating array mooring line (anchored or shared).
    The idea is to have a common object that can be used for 2D
    layout optimization as well as (optionally) 3D design.
    Work in progress. Eventually will inherit from Edge.
    '''
    
    def __init__(self, subsystem=None, rA=[0,0,0], rB=[0,0,0],
                 rad_anch=500, rad_fair=58, z_anch=-100, z_fair=-14):
        '''
        Initialize an empty object for a mooring line.
        Eventually this will fully set one up from ontology inputs.
        '''
        
        self.subsystem = subsystem  # The MoorPy subsystem that corresponds to the mooring line
        
        # end point absolute coordinates, to be set later
        self.rA = rA
        self.rB = rB
        self.heading = 0
        
        # relative positions
        self.rad_anch = rad_anch
        self.rad_fair = rad_fair
        self.z_anch   = z_anch  
        self.z_fair   = z_fair
        
        self.cost = {}
    
    
    def adjust(self, r_center=None, heading=None, project=None, degrees=False, **kwargs):
        '''Adjusts mooring position based on changed platform location or
        heading [deg]. It will call a custom "designer" function if one is
        provided. Otherwise it will just update the end positions.
        
        Parameters
        ----------
        r_center
            The x, y coordinates of the platform (undisplaced) [m].
        heading : float
            The absolute heading of the mooring line [deg].
        project : FAModel Project, optional
            A Project-type object for site-specific information used in custom
            mooring line adjustment functions (mooring.designer).
        **kwargs : dict
            Additional arguments passed through to the designer function.
        '''
        
        # Adjust heading if provided
        if not heading == None:
            if degrees:
                self.heading = np.radians(heading)
            else:
                self.heading = heading
            
        # heading 2D unit vector
        u = np.array([np.cos(self.heading), np.sin(self.heading)])
        
        r_center = np.array(r_center)[:2]
        
        # Set the updated fairlead location
        self.setEndPosition(np.hstack([r_center + self.rad_fair*u, self.z_fair]), 'b')
        
        # run custom function to update the mooring design (and anchor position)
        # this would also szie the anchor maybe?
        if self.designer:
            self.designer(self, project, r_center, u, **kwargs)
        
        else: # otherwise just set the anchor position based on a set spacing
            self.setEndPosition(np.hstack([r_center + self.rad_anch*u, self.z_anch]), 'a', sink=True)
        
    
    
    def setEndPosition(self, r, end, sink=False):
        '''Set the position of an end of the mooring.
        
        Parameters
        ----------
        r : list
            Cordinates to set the end at [m].
        end
            Which end of the edge is being positioned, 'a' or 'b'.
        sink : bool
            If true, and if there is a subsystem, the z position will be on the seabed.
        '''
        
        if end in ['a', 'A', 0]:
            self.rA = np.array(r)
            
            if self.subsystem:
                self.subsystem.setEndPosition(self.rA, False, sink=sink)
            
        elif end in ['b', 'B', 1]:
            self.rB = np.array(r)
            
            if self.subsystem:
                self.subsystem.setEndPosition(self.rB, True, sink=sink)
                
        else:
            raise Exception('End A or B must be specified with either the letter, 0/1, or False/True.')
    
    
    def getCost(self):
        
        # mooring line cost
        line_cost = 0
        if self.subsystem:
            for line in self.subsystem.lineList:
                line_cost += line.getCost()
        else:
            pass # could have some proxy cost calculation
        
        self.cost['line'] = line_cost
        
        # anchor cost  (it should already be computed)
        
        # sum up the costs in the dictionary and return
        return sum(self.cost.values()) 
        

class Platform():
    '''
    Class for a mooring floating platform.
    Eventually will inherit from Node.
    '''
    
    def __init__(self, r=[0,0], heading=0, mooring_headings=[60,180,300]):
        '''
        
        Parameters
        ----------
        r 
            x and y coordinates [m].
        theta, float (optional)
            The heading of the object [deg].
        mooring_headings (optional)
            relative headings of mooring lines [deg].
        '''
        
        self.r = np.array(r)  # coordinate of platform
        
        self.theta = np.radians(heading)  # heading offset of platform
        
        self.mooring_headings = np.radians(mooring_headings) # headings of mooring lines [rad]
        
        self.n_mooring = len(mooring_headings) # number of mooring lines
        
        # self.anchor_rads   = np.zeros(self.n_mooring)      # anchoring radius of each mooring [m]
        # self.anchor_coords = np.zeros([self.n_mooring, 2]) # coordinates of each anchor [m]
        
        self.mooringList = []  # to be filled by references to Mooring objects
    
    
    def setPosition(self, r, heading=None):
        '''
        Set the position/orientation of the platform as well as the associated
        anchor points.
        
        Parameters
        ----------
        r : list
            x and y coordinates to position the node at [m].
        heading, float (optional)
            The heading of the object [deg].
        '''
        
        '''
        # Future approach could be
        # first call the Node method to take care of the platform and what's directly attached
        Node.setPosition(self, r, heading=heading)
        # then also adjust the anchor points
        '''
        
        # Store updated position and orientation
        self.r = np.array(r)
        
        if not heading == None:
            self.theta = np.radians(heading)
        
        # Get rotation matrix...
        self.R = np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]])
        
        # Update the position of any Moorings
        for i, mooring in enumerate(self.mooringList):
        
            # heading of the mooring line
            heading_i = self.mooring_headings[i] + self.theta
            
            # adjust the whole Mooring
            mooring.adjust(self.r, heading=heading_i)
        
