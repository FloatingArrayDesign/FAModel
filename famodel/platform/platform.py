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
        
    def mooringSystem(self,rotateBool=0,mList=None):
        '''
        Create a MoorPy system for the platform based on the mooring subsystems provided in 
        the mooringList of Mooring classes and the mooring headings

        Parameters
        ----------
        rotateBool : boolean
            controls whether the list of subsystems needs rotation at the mooring headings (1) or if it is already rotated (0)
            Default is 0 (false)
        mList : list
            List of mooring subsystems on the platform
            Default is None

        Returns
        -------
        None

        '''
        import moorpy as mp
        
        #check if the subsystems were passed in from the function call
        if mList:
            if self.mooringList:
                print('There is a mooring list already in the platform class instance. This is being overwritten because a mooring list was specified in the function call.')
            #if so, set the platfrom mooringList to the passed-in list
            self.mooringList = mList
        
        #create new MoorPy system and set its depth
        self.ms = mp.System(depth=self.mooringList[0].subsystem.depth)
        
        if rotateBool:
            #rotation
            self.setPosition(self.r)
            # for i in range(0,len(self.mooringList)):
                # ssloc = self.mooringList[i].subsystem
                # #get anchor and fairlead radius
                # anch_rad = np.sqrt(ssloc.rA[0]**2+ssloc.rA[1]**2)
                # fair_rad = np.sqrt(ssloc.rB[0]**2+ssloc.rB[1]**2)
                # #set anchor and fairlead new location based on headings
                # ssloc.setEndPosition([anch_rad*np.cos(self.phi[i]), anch_rad*np.sin(self.phi[i]), -ssloc.depth], endB=0)
                # ssloc.setEndPosition([fair_rad*np.cos(self.phi[i]), fair_rad*np.sin(self.phi[i]), ssloc.rB[2]], endB=1)
        
        #make mooring system from subsystems
        for i in range(0,len(self.mooringList)):
            ssloc = self.mooringList[i].subsystem
            #add subsystem as a line to the linelist
            self.ms.lineList.append(ssloc)
            ssloc.number = i
            #add anchor point as a fixed point
            self.ms.addPoint(1,ssloc.rA)
            #attach subsystem line to the anchor point
            self.ms.pointList[-1].attachLine(i,0)
            #add fairlead point as a coupled point
            self.ms.addPoint(-1,ssloc.rB)
            #attach subsystem line to the fairlead point
            self.ms.pointList[-1].attachLine(i,1)
        
        #initialize and plot
        self.ms.initialize()
        self.ms.solveEquilibrium()
        fig,ax = self.ms.plot()
        
        