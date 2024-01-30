# class for a floating platform

import numpy as np


class Platform():
    '''
    Class for a moored floating platform with methods for requesting
    adjustments of mooring lines or power cables when the platform is 
    repositioned.
    Eventually will inherit from Node.
    '''
    
    def __init__(self, r=[0,0], heading=0, mooring_headings=[60,180,300]):
        '''
        
        Parameters
        ----------
        r 
            x and y coordinates [m].
        phi, float (optional)
            The heading of the object [deg].
        mooring_headings (optional)
            relative headings of mooring lines [deg].
        '''
        
        # Design description dictionary for this Platform
        self.dd = {}
        
        # Platform position and orientation
        self.r = np.array(r)  # x, y coordinates of platform [m]
        self.phi = np.radians(heading)  # heading offset of platform [rad]
        
        self.mooring_headings = np.radians(mooring_headings) # headings of mooring lines [rad]
        
        self.n_mooring = len(mooring_headings) # number of mooring lines
        
        # self.anchor_rads   = np.zeros(self.n_mooring)      # anchoring radius of each mooring [m]
        # self.anchor_coords = np.zeros([self.n_mooring, 2]) # coordinates of each anchor [m]
        
        self.mooringList = []  # to be filled by references to Mooring objects
        
        # Dictionaries for addition information
        self.loads = {}
        self.reliability = {}
        self.cost = {}
    
    
    def setPosition(self, r, heading=None, degrees=False):
        '''
        Set the position/orientation of the platform as well as the associated
        anchor points.
        
        Parameters
        ----------
        r : list
            x and y coordinates to position the node at [m].
        heading, float (optional)
            The heading of the platform [deg or rad] depending on
            degrees parameter (True or False).
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
            self.phi = np.radians(heading)
        
        # Get 2D rotation matrix
        self.R = np.array([[np.cos(phi), -np.sin(phi)],[np.sin(phi), np.cos(phi)]])
        
        # Update the position of any Moorings
        for i, mooring in enumerate(self.mooringList):
        
            # Heading of the mooring line
            heading_i = self.mooring_headings[i] + self.phi
            
            # Reposition the whole Mooring
            mooring.reposition(self.r, heading=heading_i)
        
