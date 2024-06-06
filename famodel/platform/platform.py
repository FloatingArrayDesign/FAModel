# class for a floating platform

import numpy as np
from famodel.famodel_base import Node
from famodel.mooring.mooring import Mooring


class Platform(Node):
    '''
    Class for a moored floating platform with methods for requesting
    adjustments of mooring lines or power cables when the platform is 
    repositioned.
    Eventually will inherit from Node.
    '''
    
    def __init__(self, id, r=[0,0], heading=0, mooring_headings=[60,180,300],rFair=None,zFair=None):
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
        # Initialize as a node
        Node.__init__(self,id)
        
        # Design description dictionary for this Platform
        self.dd = {}
        
        # Platform position and orientation
        self.r = np.array(r)  # x, y coordinates of platform [m]
        self.phi = np.radians(heading)  # heading offset of platform [rad]
        self.rFair = rFair
        self.zFair = zFair
        
        self.body = None # body object in MoorPy associated with the platform
        
        self.mooring_headings = [np.radians(mooring_headings)] # headings of mooring lines [rad]
        
        self.n_mooring = len(mooring_headings) # number of mooring lines
        
        self.endB = {} # dictionary with key as mooring object names and values of booleans (one for each mooring line) describing whether platform is connected to end B of the line (important for shared lines)
        # self.anchor_rads   = np.zeros(self.n_mooring)      # anchoring radius of each mooring [m]
        # self.anchor_coords = np.zeros([self.n_mooring, 2]) # coordinates of each anchor [m]
        
        # >>> replace these with the Node attachments dict? <<<
        self.mooringList = {}  # dictionary to be filled by references to Mooring objects
        # self.anchorList = {} # dictionary of references to anchor objects connected to this platform
        
        # Dictionaries for addition information
        self.loads = {}
        self.reliability = {}
        self.cost = {}
        self.failure_probability = {}
    
    
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
            if degrees:
                self.phi = np.radians(heading)
            else:
                self.phi = heading
        
        # Get 2D rotation matrix
        self.R = np.array([[np.cos(self.phi), -np.sin(self.phi)],[np.sin(self.phi), np.cos(self.phi)]])
        
        # Update the position of any Moorings
        for i, mooring in enumerate(self.attachments):
            if isinstance(self.attachments[mooring]['obj'], Mooring): 
        
                # Heading of the mooring line
                heading_i = self.mooring_headings[i] + self.phi
                
                # Reposition the whole Mooring
                self.attachments[mooring]['obj'].reposition(self.r, heading=heading_i)
        
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
        
        # check if the subsystems were passed in from the function call
        if not mList:
            mList = []
            for i in self.attachments:
                if isinstance(self.attachments[i]['obj'],Mooring):
                    mList.append(self.attachments[i]['obj'])
            
        
        # create new MoorPy system and set its depth
        self.ms = mp.System(depth=mList[0].subsystem.depth)
        
        if rotateBool:
            # rotation
            self.setPosition(self.r)
        
        # make mooring system from subsystems
        for i in self.attachments:
        
            # only process moorings that have subsystems for now
            
            if type(self.attachments[i]['ref']) == Mooring:
                mooring = self.attachments[i]['ref']
                ssloc = mooring.subsystem
                
                if ssloc:  # only proceed it's not None
                    # add subsystem as a line to the linelist
                    self.ms.lineList.append(ssloc)
                    ssloc.number = i
                    # check whether a moorpy anchor object exists for this mooring line
                    if not mooring.anchor.mpAnchor:
                        # create anchor moorpy object
                        mooring.anchor.makeMoorPyAnchor()
                    # add anchor point from anchor class and fairlead point adjusted to include location offsets, attach subsystem
                    self.ms.pointList.append(mooring.anchor.mpAnchor) # anchor
                    # attach subsystem line to the anchor point
                    self.ms.pointList[-1].attachLine(i,0)
                    # add fairlead point as a coupled point
                    self.ms.addPoint(-1,ssloc.rB)
                    # attach subsystem line to the fairlead point
                    self.ms.pointList[-1].attachLine(i,1)
        
        # initialize and plot
        self.ms.initialize()
        self.ms.solveEquilibrium()
        fig,ax = self.ms.plot()
        
