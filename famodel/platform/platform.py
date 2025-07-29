# class for a floating platform

import numpy as np
from famodel.famodel_base import Node
from famodel.mooring.mooring import Mooring
from famodel.turbine.turbine import Turbine
import matplotlib.pyplot as plt
from copy import deepcopy
from famodel.cables.cable import Cable
from famodel.anchors.anchor import Anchor
from famodel.cables.cable import DynamicCable

class Platform(Node):
    '''
    Class for a moored floating platform with methods for requesting
    adjustments of mooring lines or power cables when the platform is 
    repositioned.
    Eventually will inherit from Node.
    '''
    
    def __init__(self, id, r=[0,0,0], heading=0, mooring_headings=[60,180,300],rFair=None,zFair=None):
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
        
        self.entity = None # describes what type of platform this is/what its topside carries (for floating wind turbine, entity = 'FOWT', for substation, entity = 'OSS')
        
        self.rc = None # optional [row,col] information in array (useful in updateUniformArray etc.)

        
        # Dictionaries for addition information
        self.envelopes = {}  # 2D motion envelope, buffers, etc. Each entry is a dict with x,y or shape
        self.loads = {}
        self.reliability = {}
        self.cost = {}
        self.failure_probability = {}
        self.raftResults = {}
    
    
    def setPosition(self, r, heading=None, degrees=False,project=None):
        '''
        Set the position/orientation of the platform as well as the associated
        anchor points. 
        
        "Note: must only be used for a platform that's only attached with anchored lines"

        Parameters
        ----------
        r : list
            x and y coordinates to position the node at [m].
        heading, float (optional)
            The heading of the platform [deg or rad] depending on
            degrees parameter (True or False) in compass direction
        '''
        

        # first call the Node method to take care of the platform and what's directly attached
        if heading: # save compass heading in radians
            if degrees == True:
                self.phi = np.radians(heading)
            else:
                self.phi = heading
        # send in cartesian heading to node.setPosition (+ rotations CCW here)    
        Node.setPosition(self, r, theta=-self.phi)
        # then also adjust the anchor points

        
        # Update the position of any Moorings
        count = 0 # mooring counter (there are some attachments that aren't moorings)
        for i, att in enumerate(self.attachments):
            if isinstance(self.attachments[att]['obj'], Mooring): 
        
                # Heading of the mooring line
                heading_i = self.mooring_headings[count] + self.phi
                # Reposition the whole Mooring if it is an anchored line
                if not self.attachments[att]['obj'].shared:
                    self.attachments[att]['obj'].reposition(r_center=self.r, heading=heading_i,project=project)
                
                count += 1
                
            if isinstance(self.attachments[att]['obj'], Cable):
                
                cab = self.attachments[att]['obj']
                
                # update heading stored in subcomponent for attached end
                pf_phis = [cab.attached_to[0].phi, cab.attached_to[1].phi]
                headings = [cab.subcomponents[0].headingA + pf_phis[0], cab.subcomponents[-1].headingB + pf_phis[1]]
                
                # reposition the cable
                cab.reposition(headings=headings,project=project)
    
    
    def mooringSystem(self,rotateBool=0,mList=None,bodyInfo=None, project=None):
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
            
        if project and len(project.grid_depth) > 1:
            # calculate the maximum anchor spacing
            anchor_spacings = [np.linalg.norm(mooring.rA[0:2] - self.r[:2]) for mooring in mList]
            # get the bathymetry range that is related to this platform
            margin = 1.2
            small_grid_x = np.linspace((self.r[0] - np.max(anchor_spacings)*margin), (self.r[0] + np.max(anchor_spacings)*margin), 10)
            small_grid_y = np.linspace((self.r[1] - np.max(anchor_spacings)*margin), (self.r[1] + np.max(anchor_spacings)*margin), 10)
            # interpolate the global bathymetry
            small_grid_depths = np.zeros([len(small_grid_y), len(small_grid_x)])
            for i,x in enumerate(small_grid_x):
                for j,y in enumerate(small_grid_y):
                    small_grid_depths[j,i] = project.getDepthAtLocation(x, y)
            #self.ms = mp.System(bathymetry=dict(x=small_grid_x, y=small_grid_y, depth=small_grid_depths))
            self.ms = mp.System(bathymetry=dict(x=project.grid_x, y=project.grid_y, depth=project.grid_depth))
        
        elif mList[0].ss:
            # create new MoorPy system and set its depth
            self.ms = mp.System(depth=mList[0].ss.depth)
        else:
            self.ms = mp.System(depth=mList[0].z_anch)
        
        r6 = [self.r[0],self.r[1],self.r[2],0,0,0]
        # create body
        if bodyInfo:
            self.ms.addBody(0,r6,m=bodyInfo['m'],v=bodyInfo['v'],rCG=np.array(bodyInfo['rCG']),rM=np.array(bodyInfo['rM']),AWP=bodyInfo['AWP'])
        else:
            self.ms.addBody(0,r6,m=19911423.956678286,rCG=np.array([ 1.49820657e-15,  1.49820657e-15, -2.54122031e+00]),v=19480.104108645974,rM=np.array([2.24104273e-15, 1.49402849e-15, 1.19971829e+01]),AWP=446.69520543229874)
        
        if rotateBool:
            # rotation
            self.setPosition(self.r)
        
        # make mooring system from subsystems
        for i,attID in enumerate(self.attachments):
        
            # only process moorings that have subsystems for now
            
            if type(self.attachments[attID]['obj']) == Mooring:
                mooring = self.attachments[attID]['obj']
                if mooring.ss:
                    ssloc = mooring.ss
                else:
                    ssloc = mooring.createSubsystem()
                
                if ssloc:  # only proceed it's not None
                    '''
                    # add subsystem as a line to the linelist
                    self.ms.lineList.append(ssloc)
                    ssloc.number = i+1
                    '''
                    for att in mooring.attached_to:
                        if isinstance(att,Anchor):
                            # check whether a moorpy anchor object exists for this mooring line
                            # if not att.mpAnchor:
                            # create anchor moorpy object
                            att.makeMoorPyAnchor(self.ms)
                            # else:
                            #     # add anchor point from anchor class and fairlead point adjusted to include location offsets, attach subsystem
                            #     self.ms.pointList.append(att.mpAnchor) # anchor
                            # attach subsystem line to the anchor point
                            self.ms.pointList[-1].attachLine(i,0)
                            # add fairlead point as a coupled point
                            self.ms.addPoint(1,ssloc.rB)
                            # attach subsystem line to the fairlead point
                            self.ms.pointList[-1].attachLine(i,1)
                            # attach fairlead point to body
                            self.ms.bodyList[0].attachPoint(len(self.ms.pointList),self.ms.pointList[-1].r-np.append(self.r[:2], [0]))
        # initialize and plot
        self.ms.initialize()
        self.ms.solveEquilibrium()
        #fig,ax = self.ms.plot()
        
        return(self.ms)
        
        
    def getWatchCircle(self, plot=0, ang_spacing=45, RNAheight=150,
                       shapes=True,Fth=None,SFs=True,ms=None):
        '''
        Compute watch circle of platform, as well as mooring and cable tension safety factors and 
        cable sag safety factors based on rated thrust.
        
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
        if not ms:
            self.body.type = -1
            ms = self.body.sys  # work with the mooring system the body is part of
            body = self.body
        else:
            body = ms.bodyList[0]
            body.type = -1
        ms.initialize()
        ms.solveEquilibrium()

        x = []
        y = []
        
        moorings = [] # list of mooring lines attached
        cables = [] # list of cables attached
        dcs = []
        
        # find turbines, cables, and mooorings attached to platform
        moorings = self.getMoorings().values()
        cables = self.getCables().values()
        for i,cab in enumerate(cables):
            dcs.extend([sub for sub in cab.subcomponents if hasattr(sub,'ss')])
        anchors = self.getAnchors().values()
        for i in self.attachments:
            if isinstance(self.attachments[i]['obj'],Turbine):
                turbine = self.attachments[i]['obj']
        #     elif isinstance(self.attachments[i]['obj'],Mooring):
        #         moorings.append(self.attachments[i]['obj'])
        #     elif isinstance(self.attachments[i]['obj'],Cable):
        #         # find cable subcomponent attached to this cable
        #         if self.attachments[i]['end'] =='a' or 'A':
        #             cables.append(self.attachments[i]['obj'].subcomponents[0])
        #         elif self.attachments[i]['end'] == 'b' or 'B':
        #             cables.append(self.attachments[i]['obj'].subcomponents[-1])
        
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
        minSag = [None]*len(dcs)
        minCurvSF = [None]*len(dcs)
        CminTenSF = [None]*len(dcs)
        minTenSF = [None]*len(moorings)
        F = [None]*len(moorings)
        for ang in range(0, 360+ang_spacing, ang_spacing):
            print('Analyzing offset at angle ',ang)
            fx = thrust*np.cos(np.radians(ang))
            fy = thrust*np.sin(np.radians(ang))
            
            body.f6Ext = np.array([fx, fy, 0, fy*RNAheight, fx*RNAheight, 0])       # apply an external force on the body [N]                       
            
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
                        if not moor.shared:
                            if self.attachments[moor.id]['end'] == 'a':
                                # anchor attached to end B
                                F[j] = moor.ss.fB
                            else:
                                F[j] = moor.ss.fA
                
                # get tensions, sag, and curvature on cable
                for j,cab in enumerate(dcs):
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
                    if 'buoyancy_sections' in cab.dd and cab.dd['buoyancy_sections']:
                        nb = len(cab.dd['buoyancy_sections'])
                        m_s = []
                        for k in range(0,nb):
                            m_s.append(cab.ss.getSag(2*k))
                        mS = min(m_s)
                        if not minSag[j] or minSag[j]<mS:
                            minSag[j] = deepcopy(mS)
        
            x.append(body.r6[0])       
            y.append(body.r6[1])
        
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
        body.f6Ext = np.array([0, 0, 0, 0, 0, 0])
        ms.solveEquilibrium3(DOFtype='both')
        
                
        if SFs:
            # save anchor loads
            for j,moor in enumerate(moorings):
                for att3 in moor.attached_to:
                    if isinstance(att3,Anchor):
                        att3.loads['Hm'] = np.sqrt(F[j][0]**2+F[j][1]**2)
                        att3.loads['Vm'] = F[j][2]
                        att3.loads['thetam'] = np.degrees(np.arctan(att3.loads['Vm']/att3.loads['Hm'])) #[deg]
                        att3.loads['mudline_load_type'] = 'max'
                    
            maxVals = {'minTenSF':minTenSF,'minTenSF_cable':CminTenSF,'minCurvSF':minCurvSF,'minSag':minSag,'maxF':F}# np.vstack((minTenSF,CminTenSF,minCurvSF,minSag))    
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
        

    def getBufferZones(self, buffer_rad=50,buffType=['Anchor','Mooring','Platform']):     
        ''' Function to calculate buffer zones around mooring lines and anchors.
        
        Parameters
        ----------
        buffer_rad: Radius of buffer zones in m
                
        '''
        from shapely.geometry import LineString, MultiLineString, Polygon, Point
        from shapely.ops import unary_union
        
        # get anchor objects connected to platform
        anchorList = self.getAnchors()
        moorList = self.getMoorings()
        
        # Create LineString geometries and buffer them
        buffer_group = []
        if 'Anchor' in buffType:
            for anch in anchorList:
                # im = 3*i + j  # global index of mooring/anchor
                #line = LineString([self.r, anchorList[anch].r[:2]])
                ##line = LineString([self.turb_coords[i,:], self.anchor_coords[im,:]])
                point = Point(anchorList[anch].r[0],anchorList[anch].r[1])
                buffer = point.buffer(buffer_rad)
                #buffer = line.buffer(buffer_rad)
                buffer_group.append(buffer)
        if 'Mooring' in buffType:
            for moor in moorList:
                if 'mean' in moor.envelopes:
                    poly = moor.envelopes['mean']['shape']
                    buffer = poly.buffer(buffer_rad)
                else:
                    line = LineString(moor.rA[:2],moor.rB[:2])
                    buffer = line.buffer(buffer_rad)
                    
                buffer_group.append(buffer)
        if 'Platform' in buffType:
            if 'mean' in self.envelopes:
                poly = self.envelopes['mean']['shape']
                buffer = poly.buffer(buffer_rad)
            else:
                point = Point(self.r[0],self.r[1])
                buffer = point.buffer(buffer_rad*3)
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

