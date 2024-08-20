# class for a floating platform

import numpy as np
from famodel.famodel_base import Node
from famodel.mooring.mooring import Mooring
from famodel.turbine.turbine import Turbine
import matplotlib.pyplot as plt
from copy import deepcopy
from famodel.cables.cable import Cable
from famodel.anchors.anchor import Anchor

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
        
        self.mooring_headings = list(np.radians(mooring_headings)) # headings of mooring lines [rad]
        
        self.n_mooring = len(mooring_headings) # number of mooring lines
        
        self.endB = {} # dictionary with key as mooring object names and values of booleans (one for each mooring line) describing whether platform is connected to end B of the line (important for shared lines)
        # self.anchor_rads   = np.zeros(self.n_mooring)      # anchoring radius of each mooring [m]
        # self.anchor_coords = np.zeros([self.n_mooring, 2]) # coordinates of each anchor [m]
        
        self.rc = None # optional [row,col] information in array (useful in updateUniformArray etc.)
        # >>> replace these with the Node attachments dict? <<<
        # self.mooringList = {}  # dictionary to be filled by references to Mooring objects
        # self.anchorList = {} # dictionary of references to anchor objects connected to this platform
        
        # Dictionaries for addition information
        self.envelopes = {}  # 2D motion envelope, buffers, etc. Each entry is a dict with x,y or shape
        self.loads = {}
        self.reliability = {}
        self.cost = {}
        self.failure_probability = {}
    
    
    def setPosition(self, r, heading=None, degrees=False,project=None):
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
        
        # correction for xy coords
        corr = np.radians(90)
        
        # Get 2D rotation matrix
        self.R = np.array([[np.cos(corr-self.phi), -np.sin(corr-self.phi)],[np.sin(corr-self.phi), np.cos(corr-self.phi)]])
        
        # Update the position of any Moorings
        count = 0 # mooring counter (there are some attachments that aren't moorings)
        for i, att in enumerate(self.attachments):
            if isinstance(self.attachments[att]['obj'], Mooring): 
        
                # Heading of the mooring line
                heading_i = self.mooring_headings[count] + self.phi
                # Reposition the whole Mooring
                self.attachments[att]['obj'].reposition(r_center=self.r, heading=heading_i,project=project)
                                
                count += 1
                
            if isinstance(self.attachments[att]['obj'], Cable):
                
                cab = self.attachments[att]['obj']
                
                # update headings stored in subcomponents
                headings = [cab.subcomponents[0].headingA + self.phi, cab.subcomponents[-1].headingB + self.phi]
                
                # reposition the cable
                cab.reposition(headings=headings)
        
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
        
        
    def getWatchCircle(self, plot=0, ang_spacing=45, RNAheight=150,
                       shapes=True,Fth=None,SFs=True):
        '''
        Compute watch circle of platform based on rated thrust.
        
        Parameters
        ----------
        ang_spacing : float
            Angle increment to evaluate watch circle at [deg].
        plot : bool
            Plots the shape of the watch circle
        RNAheight : float
            Height of the rotor-nacelle assembly
        shapes : bool
            Whether or not to create shapely objects
        Fth : float
            Thrust force
        SFs : bool
            WHether or not to calculate safety factors etc for the line
        Returns
        -------
        x: list of x-coordinates for watch circle
        y: list of y-coordinates for watch circle
        maxVals: dictionary of minimum safety factors for line tension, cable tension and cable curvature, and the minimum sag of cables

        '''
        self.body.type = -1
        ms = self.body.sys  # work with the mooring system the body is part of
        ms.initialize()
        
        x = []
        y = []
        
        moorings = [] # list of mooring lines attached
        cables = [] # list of cables attached
        
        # find turbines, cables, and mooorings attached to platform
        for i in self.attachments:
            if isinstance(self.attachments[i]['obj'],Turbine):
                turbine = self.attachments[i]['obj']
            elif isinstance(self.attachments[i]['obj'],Mooring):
                moorings.append(self.attachments[i]['obj'])
            elif isinstance(self.attachments[i]['obj'],Cable):
                # find cable subcomponent attached to this cable
                if self.attachments[i]['end'] =='a' or 'A':
                    cables.append(self.attachments[i]['obj'].subcomponents[0])
                elif self.attachments[i]['end'] == 'b' or 'B':
                    cables.append(self.attachments[i]['obj'].subcomponents[-1])
        
        if Fth:
            thrust = Fth
        else:
            try:
                # create rotor
                turbine.makeRotor()
                # get thrust curve
                turbine.calcThrustForces()
            except Exception as e:
                print(e)
                print('Could not get thrust forces from RAFT, using IEA 15 MW turbine thrust as default')
                thrust = 1.95e6
        
        # btenMax = np.zeros((len(moorings),1))
        # atenMax = np.zeros((len(moorings),1))
        # CbtenMax = np.zeros((len(cables),1))
        # CatenMax = np.zeros((len(cables),1))
        minSag = [None]*len(cables)
        minCurvSF = [None]*len(cables)
        CminTenSF = [None]*len(cables)
        minTenSF = [None]*len(moorings)
        
        for ang in range(0, 360+ang_spacing, ang_spacing):
            print(ang)
            fx = thrust*np.cos(np.radians(ang))
            fy = thrust*np.sin(np.radians(ang))
            
            self.body.f6Ext = np.array([fx, fy, 0, fy*RNAheight, fx*RNAheight, 0])       # apply an external force on the body [N]                       
            
            ms.solveEquilibrium3(DOFtype='both')                        # equilibrate (specify 'both' to enable both free and coupled DOFs)
            
            if SFs:
                # get tensions on mooring line
                for j,moor in enumerate(moorings): 
                    MBLA = float(moor.ss.lineList[0].type['MBL'])
                    MBLB = float(moor.ss.lineList[-1].type['MBL'])
                    # print(MBLA,MBLB,moor.ss.TA,moor.ss.TB,MBLA/moor.ss.TA,MBLB/moor.ss.TB,abs(MBLA/moor.ss.TA),abs(MBLB/moor.ss.TB))
                    MTSF = min([abs(MBLA/moor.ss.TA),abs(MBLB/moor.ss.TB)])
                    # atenMax[j], btenMax[j] = moor.updateTensions()
                    if not minTenSF[j] or minTenSF[j]>MTSF:
                        minTenSF[j] = deepcopy(MTSF)
                
                # get tensions, sag, and curvature on cable
                for j,cab in enumerate(cables):
                    MBLA = cab.ss.lineList[0].type['MBL']
                    MBLB = cab.ss.lineList[-1].type['MBL']
                    CMTSF = min([abs(MBLA/cab.ss.TA),abs(MBLB/cab.ss.TB)])
                    if not CminTenSF[j] or CminTenSF[j]>CMTSF:
                        CminTenSF[j] = deepcopy(CMTSF)
                    # CatenMax[j], CbtenMax[j] = cab.updateTensions()
                    cab.ss.calcCurvature()
                    mCSF = cab.ss.getMinCurvSF()
                    if not minCurvSF[j] or minCurvSF[j]>mCSF:
                        minCurvSF[j] = mCSF
                    # determine number of buoyancy sections
                    nb = len(cab.dd['buoyancy_sections'])
                    m_s = []
                    for k in range(0,nb):
                        m_s.append(cab.ss.getSag(2*k))
                    mS = min(m_s)
                    if not minSag[j] or minSag[j]<mS:
                        minSag[j] = deepcopy(mS)
        
            x.append(self.body.r6[0])       
            y.append(self.body.r6[1])
        
        # Convert to np array and save in object envelope
        x = np.array(x)
        y = np.array(y)
        self.envelopes['mean'] = dict(x=np.array(x), y=np.array(y))
        
        if shapes:  # want to *optionally* make a shapely polygon
            from shapely import Polygon
            self.envelopes['mean']['shape'] = Polygon(list(zip(x,y)))
         

        if plot:
            plt.figure()
            plt.plot(x,y)
            
        # restore platform to equilibrium position 
        self.body.f6Ext = np.array([0, 0, 0, 0, 0, 0])
        ms.solveEquilibrium3(DOFtype='both')
                
        if SFs:
            maxVals = {'minTenSF':minTenSF,'minTenSF_cable':CminTenSF,'minCurvSF':minCurvSF,'minSag':minSag}# np.vstack((minTenSF,CminTenSF,minCurvSF,minSag))    
            return(x,y,maxVals)
        else:
            return(x,y)
    
    def getMoorings(self):
        '''
        Function to get mooring lines connected to this platform

        Returns
        -------
        mooringList: dict
            dictionary of mooring objects connected to the platform

        '''
        mooringList = {}
        for i in self.attachments:
            if isinstance(self.attachments[i]['obj'],Mooring):
                mooringList[self.attachments[i]['id']] = self.attachments[i]['obj']
        
        return(mooringList)
                
    def getCables(self):
        '''
        Function to get cables connected to this platform

        Returns
        -------
        cableList: dict
            dictionary of cable objects connected to the platform

        '''
        cableList = {}
        for att in self.attachments:
            if isinstance(self.attachments[att]['obj'],Cable):
                cableList[self.attachments[att]['id']] = self.attachments[att]['obj']
                
        return(cableList)
    
    def getAnchors(self):
        '''
        Function to find anchors associated with this platform
        
        Returns
        -------
        anchorList: dict
            dictionary of cable objects connected to the platform

        '''
        anchorList = {}
        mList = self.getMoorings()
        for moor in mList.values():
            for att in moor.attached_to:
                if isinstance(att,Anchor):
                    anchorList[att.id] = att
                
        return(anchorList)
        

    def getbufferzones(self, buffer_rad=50):     
        ''' Function to calculate buffer zones around mooring lines and anchors.
        
        Parameters
        ----------
        buffer_rad: Radius of buffer zones in m
                
        '''
        from shapely.geometry import LineString, MultiLineString, Polygon
        from shapely.ops import unary_union
        
        # get anchor objects connected to platform
        anchorList = self.getAnchors()
        
        # Create LineString geometries and buffer them
        buffer_group = []
        for anch in anchorList:
            # im = 3*i + j  # global index of mooring/anchor
            line = LineString([self.r, anchorList[anch].r[:2]])
            #line = LineString([self.turb_coords[i,:], self.anchor_coords[im,:]])
            buffer = line.buffer(buffer_rad)
            buffer_group.append(buffer)
        
        # Combine the buffered lines connected to the same turbine into one polygon
        polygon = unary_union(buffer_group)  # Combine buffers for each turbine
        self.envelopes['buffer_zones'] = {}
        if isinstance(polygon, MultiLineString):
            # Convert MultiLineString to Polygon
            self.envelopes['buffer_zones']['shape'] = Polygon(polygon)
        else:
            self.envelopes['buffer_zones']['shape'] = polygon
        # get coords of object
        xx,yy = self.envelopes['buffer_zones']['shape'].exterior.coords.xy
        x = xx.tolist()
        y = yy.tolist()
        self.envelopes['buffer_zones']['x'] = x
        self.envelopes['buffer_zones']['y'] = y

        return  self.envelopes['buffer_zones']

