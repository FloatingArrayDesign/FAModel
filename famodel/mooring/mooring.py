# class for a mooring line

import numpy as np
from copy import deepcopy
from moorpy.subsystem import Subsystem
from moorpy import helpers
from famodel.mooring.connector import Connector, Section
from famodel.famodel_base import Edge

class Mooring(Edge):
    '''
    Class for a floating array mooring line (anchored or shared).
    The idea is to have a common object that can be used for 2D
    layout optimization as well as (optionally) 3D design.
    Work in progress. Eventually will inherit from Edge.
    '''
    
    def __init__(self, dd=None, subsystem=None, anchor=None,
                 rho=1025, g=9.81,id=None,shared=0):
        '''
        Parameters
        ----------
        dd: dictionary
            Design dictionary that contains all information on a mooring line needed to create a MoorPy subsystem
            Layout: {
                     sections:
                         {
                             0 
                                 {
                                  type:
                                      {
                                         name, d_nom, material, d_vol, m, EA, EAd, EAd_Lm, MBL, cost, weight
                                      }
                                  L
                                 }
                         }
                     connectors:
                     [
                        ...
                     ]
                     span
                     zAnchor
                     z_fair
                     rad_fair
                     EndPositions:
                                  {
                                    endA, endB
                                  }
                    }
        Initialize an empty object for a mooring line.
        Eventually this will fully set one up from ontology inputs.
        
        >>> This init method is a placeholder that currently may need
        some additional manual setup of the mooring object after it is
        called. <<<
        
        '''    
        Edge.__init__(self, id)  # initialize Edge base class
        # Design description dictionary for this Mooring
        self.dd = dd

        # let's turn the dd into something that holds subdict objects of connectors and sections
        if self.dd:
            # >>> list of sections ?  And optional of section types (e,g, chian, poly) 
            # and dict of scaling props (would be used by linedesign) ?
            self.n_sec = len(self.dd['sections'])
            
            # Turn what's in dd and turn it into Sections and Connectors
            for i, con in enumerate(self.dd['connectors']):
                if con and 'type' in con:
                    Cid = con['type']+str(i)
                else:
                    Cid = 'Conn'+str(i)
                self.dd['connectors'][i] = Connector(Cid,**self.dd['connectors'][i])
            
            for i, sec in enumerate(self.dd['sections']):
                self.dd['sections'][i] = Section('Section'+str(i),**self.dd['sections'][i])
                #self.dd['connectors'][i  ].attach(self.dd['sections'][i], end=0)
                #self.dd['connectors'][i+1].attach(self.dd['sections'][i], end=1)
            
            # >>> add some error checks for the correct lengths <<<
            
            # Connect them and store them in self(Edge).subcomponents!
            subcons = []  # temporary list of node-edge-node... to pass to the function
            for i in range(self.n_sec):
                subcons.append(self.dd['connectors'][i])
                subcons.append(self.dd['sections'][i])
            subcons.append(self.dd['connectors'][-1])
            self.addSubcomponents(subcons)  # Edge method to connect and store em
            
            # Indices of connectors and sections in self.subcomponents list
            self.i_con = list(range(0, 2*self.n_sec+1, 2))
            self.i_sec = list(range(1, 2*self.n_sec+1, 2))
        
        # MoorPy subsystem that corresponds to the mooring line
        self.ss = subsystem
        self.ss_mod = None
        # workaround for users who are making mooring objects based on pre-existing subsystems
        if self.ss:
            self.dd = {}
            if not 'zAnchor' in self.dd:
                self.dd['zAnchor'] = -self.ss.depth
            if not 'span' in self.dd:
                self.dd['span'] = self.ss.span
            if not 'rad_fair' in self.dd:
                self.dd['rad_fair'] = self.ss.rad_fair
            if not 'z_fair' in self.dd:
                self.dd['z_fair'] = self.ss.z_fair
        
        
        # relative positions (variables could be renamed)
        self.rad_anch = self.dd['span'] + self.dd['rad_fair']
        self.rad_fair = self.dd['rad_fair']
        self.z_anch   = self.dd['zAnchor']
        self.z_fair   = self.dd['z_fair']
        
        # end point absolute coordinates, to be set later
        self.rA = np.array([-self.rad_anch, 0, self.z_anch])
        self.rB = np.array([-self.rad_fair, 0, self.z_fair])
        self.heading = 270  # compass heading from B to A [deg]
        
        self.adjuster = None  # custom function that can adjust the mooring
        
        self.shared = shared # int for if the mooring line is a shared line (1) or part of a shared line (2)
        self.symmetric = False # boolean for if the mooring line is a symmetric shared line
        
        # relevant site info
        self.rho = rho
        self.g = g
        
        # Dictionaries for addition information
        self.envelopes = {}  # 2D motion envelope, buffers, etc.
        self.loads = {}
        self.reliability = {}
        self.cost = {}
        self.failure_probability = {}
    
    
    def update(self, dd=None):
        '''Update the Mooring object based on the current state of the design
        dictionary (self.dd) or, if passed in, a new design dictionary (dd).
        '''
        
        if not dd == None:  # if dd passed in
            self.dd.update(dd)  # move contents of dd into Mooring.dd
        
        # Update section lengths and types
        for i, sec in enumerate(dd['sections']):
            if self.ss:
                self.ss.lineList[i].setL(sec[i]['L'])
                self.ss.lineTypes[i] = sec[i]['type']
        
        #TODO: update any other things (connectors, positions...)
        
        
    def setSectionLength(self, L, i):
        '''Sets length of section, including in the subdsystem if there is
        one.'''
        
        self.dd['sections'][i]['L'] = L  # set length in dd (which is also Section/subcomponent)
        
        if self.ss:  # is Subsystem exists, adjust length there too
            self.ss.lineList[i].setL(L)
    
    """
    def setSectionDiameter(self, d, i):
        '''Adjusts diameter of section (including in the subdsystem if there is
        one, since it should point to the lineType dict in the Mooring.'''
        
        # just adjust the dict? ss.lineTypes should already point to it
        self.sectionType[i].update( setLineType( ... d))
    """
    
    def setSectionType(self, lineType, i):
        '''Sets lineType of section, including in the subdsystem 
        if there is one.'''
        
        # set type dict in dd (which is also Section/subcomponent)
        self.dd['sections'][i]['type'] = lineType  
        
        if self.ss:  # is Subsystem exists, adjust length there too
            self.ss.lineTypes[i] = lineType
    
    
    
    
    def reposition(self, r_center=None, heading=None, project=None, degrees=False, **kwargs):
        '''Adjusts mooring position based on changed platform location or
        heading. It can call a custom "adjuster" function if one is
        provided. Otherwise it will just update the end positions.
        
        Parameters
        ----------
        r_center : list or nested list
            The x, y coordinates of the platform(s) (undisplaced) [m]. If shared mooring, must be a list of lists, with the
            coordinates for each platform connection. In this case, end A platform connection is the first entry and end B 
            platform connection is the second entry.
        heading : float
            The absolute heading compass direction of the mooring line from end B
            [deg or rad] depending on degrees parameter (True or False). Must account for the platform heading as well.
        project : FAModel Project, optional
            A Project-type object for site-specific information used in custom
            mooring line adjustment functions (mooring.adjuster).
        **kwargs : dict
            Additional arguments passed through to the designer function.
        '''
        
        # Adjust heading if provided
        if not heading == None:
            if degrees:
                self.heading = heading
            else:
                self.heading = np.degrees(heading)

        phi = np.radians(90-self.heading) # heading in x-y radian convention [rad] (CONVERTS FROM COMPASS HEADINGS!)
        # heading 2D unit vector
        u = np.array([np.cos(phi), np.sin(phi)])
        
        if self.shared == 1:
            r_centerB = np.array(r_center)[1][:2]
            r_centerA = np.array(r_center)[0][:2]
        else:
            r_centerB = np.array(r_center)[:2]
            
        # Set the updated end B location
        
        self.setEndPosition(np.hstack([r_centerB + self.rad_fair*u, self.z_fair]), 'b')
        
        # Run custom function to update the mooring design (and anchor position)
        # this would also szie the anchor maybe?
        if self.adjuster:
            self.adjuster(self, project, r_centerB, u, **kwargs)
            
        elif self.shared == 1: # set position of end A at platform end A
            self.setEndPosition(np.hstack([r_centerA - self.rad_fair*u, self.z_fair]),'a')
        
        else: # otherwise just set the anchor position based on a set spacing
            self.setEndPosition(np.hstack([r_centerB + self.rad_anch*u, self.z_anch]), 'a', sink=True)
        
    
    
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
            
            if self.ss:
                self.ss.setEndPosition(self.rA, False, sink=sink)
            
        elif end in ['b', 'B', 1]:
            self.rB = np.array(r)
            
            if self.ss:
                self.ss.setEndPosition(self.rB, True, sink=sink)
                
        else:
            raise Exception('End A or B must be specified with either the letter, 0/1, or False/True.')
    
    
    def getCost(self):
        
        # mooring line cost
        line_cost = 0
        if self.ss:
            for line in self.ss.lineList:
                line_cost += line.getCost()
        else:
            pass # could have some proxy cost calculation
        
        self.cost['line'] = line_cost
        
        # anchor cost  (it should already be computed)
        
        # sum up the costs in the dictionary and return
        return sum(self.cost.values())
    
    def updateTensions(self):
        ''' Gets tensions from subsystem and updates the max tensions dictionary if it is larger than a previous tension
        '''
        if not 'TAmax' in self.loads:
            self.loads['TAmax'] = 0
        if not 'TBmax' in self.loads:
            self.loads['TBmax'] = 0
        # get anchor tensions
        if abs(self.ss.TA) > self.loads['TAmax']:
            self.loads['TAmax'] = deepcopy(self.ss.TA)
        # get TB tensions
        if abs(self.ss.TB) > self.loads['TBmax']:
            self.loads['TBmax'] = deepcopy(self.ss.TB)
            
        return(self.loads['TAmax'],self.loads['TBmax'])
    
    
    def createSubsystem(self, case=0,pristine=True,dd=None):
        ''' Create a subsystem for a line configuration from the design dictionary
        
        Parameters
        ----------
        pristine : bool, optional
            0/False: modified line (marine growth, corrosion, etc) 1/True: pristine line (default)
        case : int, optional
            Selector shared/suspended cases:
                - 0 (default): end A is on the seabed
                - 1: assembly is suspended and end A is at another floating system
                - 2: the assembly is suspended and assumed symmetric, end A is the midpoint
        dd : dict, optional
            Dictionary describing the design
        flag : whether or not the warning for a subsystem already existing has been printed (to prevent mass printing to screen
                                                                                             for big projects)
        '''
        # set design dictionary as self.dd if none given, same with connectorList
        if not dd:
            dd = self.dd
        ss=Subsystem(depth=-dd['zAnchor'], rho=self.rho, g=self.g, 
                          span=dd['span'], rad_fair=self.rad_fair,
                          z_fair=self.z_fair)
        lengths = []
        types = []
        # run through each line section and collect the length and type
        for i, sec in enumerate(dd['sections']):
            lengths.append(sec['L'])
            # points to existing type dict in self.dd for now
            types.append(sec['type']) # list of type names
            #types.append(sec['type']['name']) # list of type names
            #self.ss.lineTypes[i] = sec['type']  

        
        # make the lines and set the points 
        ss.makeGeneric(lengths, types, 
            connectors=[dd['connectors'][ic+1] for ic in range(len(dd['connectors'])-2)], 
            suspended=case)
        ss.setEndPosition(self.rA,endB=0)
        ss.setEndPosition(self.rB,endB=1)
        
        # note: next bit has similar code/function as Connector.makeMoorPyConnector <<<
        
        # add in connector info to subsystem points
        if case == 0: # has an anchor - need to ignore connection for first point because anchor is a point itself so can't have a point attached to a point
            startNum = 1
        else: # no anchor - need to include all connections
            startNum = 0 

        for i in range(startNum,len(ss.pointList)):                               
            point = ss.pointList[i]
            point.CdA = dd['connectors'][i]['CdA']
        # solve the system
        ss.initialize()
        ss.staticSolve()
        
        # save ss to the correct Mooring variable
        if pristine:
            # save to ss
            self.ss = ss
            return(self.ss)
        else:
            # save to modified ss (may have marine growth, corrosion, etc)
            self.ss_mod = ss
            return(self.ss_mod)
        
    def mirror(self,create_subsystem=True):
        ''' Mirrors a half design dictionary. Useful for symmetrical shared mooring lines where 
        only half of the line is provided 

        Parameters
        ----------
        create_subsystem : bool, optional
            Controls whether to create a new subsystem after mirroring dd.
            Default is True
        '''
        # find out if the connector at end A (center of line) is empty 
        if self.dd['connectors'][0]:
            # do not double the middle section length
            doubleL = False
        else:
            # double the middle section length (remove the empty connector essentially)
            doubleL = True
        from copy import deepcopy
        # sections list is currently reversed (second half of full list)
        addSections = deepcopy(self.dd['sections'])
        addConns = deepcopy(self.dd['connectors'])
        # reverse to get first half of full list
        self.dd['sections'].reverse()
        self.dd['connectors'].reverse()
        # combine middle mirrored sections into one by removing from add list and doubling length
        if doubleL:
            # remove mirror section and connector in middle, double length of middle section
            addSections.pop(0)
            addConns.pop(0)
            self.dd['sections'][-1]['L'] = self.dd['sections'][-1]['L']*2

        # combine addSections and current reversed sectiosn list
        self.dd['sections'].extend(addSections)
        self.dd['connectors'].extend(addConns)
        
        # Connect them and store them in self(Edge).subcomponents!
        subcons = []  # temporary list of node-edge-node... to pass to the function
        for i in range(len(self.dd['sections'])):
            subcons.append(self.dd['connectors'][i])
            subcons.append(self.dd['sections'][i])
        subcons.append(self.dd['connectors'][-1])
        self.addSubcomponents(subcons)  # Edge method to connect and store em
        
        self.n_sec = len(self.dd['sections'])
        # Indices of connectors and sections in self.subcomponents list
        self.i_con = list(range(0, 2*self.n_sec+1, 2))
        self.i_sec = list(range(1, 2*self.n_sec+1, 2))
    

        # creat subsystem if asked
        if create_subsystem:
            self.createSubsystem(case=1)



        
    
    """
    # rough method ideas...maybe not necessary or can use existing dict methods
    def ddApply():
        '''Applies the current contents of self.dd to update all the parts
        of the Mooring object, including the subsystem if enabled.'''
    
    def ddLoad(self, new_dd):
        '''Takes in a new self.dd then applies it?'''
        
        self.dd = new_dd
        self.ddApply()
    
    def ddSave(self):
        '''Takes the current design info stored in the self attributes
        and writes the values into self.dd. Also returns it, for use in
        writing output yamls...'''
        
        ... save things into self.dd...
        
        return self.dd
    """
    
    def addMarineGrowth(self, mgDict, project=None, idx=None):
        '''Re-creates sections part of design dictionary to account for marine 
        growth on the subystem, then calls createSubsystem() to recreate the line

        Parameters
        ----------
        mgDict : dictionary
            Provides marine growth thicknesses and the associated depth ranges
            {
                th : list with 3 values in each entry - thickness, range lower z-cutoff, range higher z-cutoff [m]
                    *In order from sea floor to surface*
                    example, if depth is 200 m: - [0.00,-200,-100]
                                                - [0.05,-100,-80]
                                                - [0.10,-80,-40]
                                                - [0.20,-40,0]
                rho : list of densities for each thickness, or one density for all thicknesses, [kg/m^3] (optional - default is 1325 kg/m^3)
                }
        project : object, optional
            A FAModel project object, with the mooringListPristine used as the basis
            to build the marine growth model, necessary if the addMarineGrowth method
            will be called in a loop (or even just multiple times) to improve accuracy 
            of change depths, which may decrease significantly after solveEquilibrium() 
            for the moorpy model. The default is None.
        idx : tuple, optional
            A key for the pristineMooringList in the project object that is associated
            with this mooring object. Since the pristineMooringList is a deepcopy of the 
            project mooringList, the mooring objects are not the same and therefore if the 
            project object is provided in the method call, the index must also be provided.

        Returns
        -------
        changePoints : list
            List of point indices in the moorpy subsystem that are at the changeDepth
        changeDepths : list
            List of cutoff depths the changePoints should be located at

        '''
        # set location of reference mooring object
        # if project: # use pristine line
        #     oldLine = project.mooringListPristine[idx]
        # else: # use current mooring object
        #     oldLine = self
        
        # create a reference subsystem if it doesn't already exist
        if not self.ss:
            self.createSubsystem()     
        oldLine = self.ss
        # set up variables
        LTypes = [] # list of line types for new lines (types listed are from reference object)
        Lengths = [] # lengths of each section for new line
        Mats = [] # materials list for new line        
        connList = [] # new list of connectors (need to add empty connector objects in between changeDepths)
        LThick = [] # list of mg thicknesses for new lines
        ln_raw = [] # list of line lengths from rA to current split in line (needed to determine length of new sections when there are multiple splits in one line section)
        # set up variables needed to check before/after of changeDepths
        changePoints = []
        changeDepths = [] # index of list that has the corresponding changeDepth
        
        # set first connector
        connList.append(self.dd['connectors'][0])
        # go through each line section
        for i in range(0,len(oldLine.lineList)):
            slthick = [] # mg thicknesses for the section (if rA is above rB, needs to be flipped before being added to full subsystem list LThick)
            slength = [] # line lengths for the section (if rA is above rB, needs to be flipped before being added to full subsystem list)
            schangeDepth = [] # changeDepths for the section (if rA is above rB, needs to be flipped before being added to full subsystem list)
            # set reference subsystem line section location
            ssLine = oldLine.lineList[i]
            # add line material, type to list
            Mats.append(ssLine.type['material'])
            LTypes.append(ssLine.type['name'])
                       
            # check whether rA is above rB (can happen for sections of shared lines)
            if ssLine.rA[2]>ssLine.rB[2]: # set low and high point locations accordingly
                low = ssLine.rB
                high = ssLine.rA
                flip = 1
            else:
                low = ssLine.rA
                high = ssLine.rB
                flip = 0

            th = mgDict['th'] # set location for ease of coding
            # look up what thickness this line section starts at (if lowest point is not on the sea floor, first segment will have a thickness other than the sea floor thickness)
            rAth = 0 # exit while loop when we find thickness at low
            count1 = 0 # counter
            while rAth==0 and count1 <= len(th):
                if flip:
                    if high[2] <= th[count1][2]:
                        LThick.append(th[count1][0])
                        rAth = 1 # exit while loop
                else:
                    if low[2] <= th[count1][2]:
                        LThick.append(th[count1][0])
                        rAth = 1 # exit while loop
                count1 = count1 + 1 # iterate counter
                
            # determine if this line section will be split
            for j in range(0,len(th)): # go through all changeDepths
                if flip:
                    rs = 2
                else:
                    rs = 1
                if th[j][rs]>low[2] and th[j][rs]<high[2]:
                    # line section will be split - add line type, mg thickness, and material to list
                    LTypes.append(ssLine.type['name'])
                    slthick.append(th[j][0])
                    Mats.append(ssLine.type['material'])
                    # add an empty connector object to list for split location
                    connList.append(Connector(str(len(connList))+'_empty'))
                    changePoints.append(len(connList)-1)
                    schangeDepth.append([j,rs])
                    
                    # get length of line between each node
                    lenseg = ssLine.L/ssLine.nNodes
                    
                    old_line = ssLine.getLineCoords(Time=0) # get the coordinates of the line
                    #find length of each new section by finding node at changeDepth
                    for k in range(0, ssLine.nNodes-1): # go through each node in the line
                        if flip: # need to check the node ahead is <= the changeDepth to see which node is split
                            if old_line[2][k+1]<=th[j][rs] and old_line[2][k]>th[j][rs]:
                                nodeD = k+1 # nodeD is closest node below changeDepth
                                xp = old_line[2][::-1] # flip because np.interp doesn't take
                                yp = old_line[1][::-1]
                                fp = old_line[0][::-1]
                        else:
                            if old_line[2][k]<=th[j][rs] and old_line[2][k+1]>th[j][rs]:
                                nodeD = k # the node right below changeDepth depth
                                xp = old_line[2][:]
                                yp = old_line[1][:]
                                fp = old_line[0][:]
                    
                    # interpolate to find x & y coordinates at chosen depth (since node might not be exactly at the right depth)
                    xChange = float(np.interp(th[j][rs], xp, fp))
                    yChange = float(np.interp(th[j][rs], xp, yp))
                    
                    # get the "raw length" of the new lower line (from lower end of section to split point) - if there is multiple splits in one line section this raw length may not be the actual length of the new line
                    if flip: # node numbers start at end A (top) so need to subtract from total line length
                        # raw length = total line length - nodeD*(segment length) + 3d pythagorean theorem (to find length from nodeD to actual cutoff location)
                        ln_raw.append(ssLine.L - lenseg*nodeD + np.sqrt((xChange-old_line[0][nodeD])**2 + (yChange-old_line[1][nodeD])**2 + (th[j][rs]-old_line[2][nodeD])**2))
                    else:
                        # raw length = nodeD*(segment length) + 3d pythagorean theorem 
                        ln_raw.append(lenseg*nodeD + np.sqrt((xChange-old_line[0][nodeD])**2 + (yChange-old_line[1][nodeD])**2 + (th[j][rs]-old_line[2][nodeD])**2))
                    
                    
                    if len(slthick)>1: # line has multiple cuts (upper cut sections have to count the length only from previous nodeD)
                        slength.append(float(ln_raw[-1]-ln_raw[-2]))
                        
                    else: # first split (raw length is actual length)
                        slength.append(float(ln_raw[-1]))
                
            if slthick: # add the length of the top line (from last split to upper end of section) if there was a split in the line
                slength.append(float(ssLine.L-ln_raw[-1]))
                # if rA above rB, reverse the order of the section-level lists (we gathered info from lowest depths up, but this line segment travels from top to bottom)
                if flip:
                    slength.reverse()
                    slthick.reverse()
                    schangeDepth.reverse()
                # Append section-level lists to the subsystem-level lists
                Lengths.extend(slength)
                LThick.extend(slthick)
                changeDepths.extend(schangeDepth)
            else: # line section was not split, add full line length
                Lengths.append(ssLine.L)
                
            # add connector at end of section to list
            connList.append(self.dd['connectors'][i+1])
                
        # Set up list variables for pristine line info
        EA = []
        m = []
        d_ve_old = []
        cd = []
        cdAx = []
                                            
        # create arrays
        d_nom_old = np.zeros((len(LTypes),1))        
        ve_nom_adjust = np.zeros((len(LTypes),1))
        mu_mg = np.zeros((len(LTypes),1))
        rho_mg = np.ones((len(LTypes),1))*1325
        # adjust rho value if alternative provided
        if 'rho' in mgDict:
            if not type(mgDict['rho']) is list:
                # just one density given for all marine growth on the line
                rho_mg = rho_mg*mgDict['rho']/1325
            else: # density given for each thickness of marine growth
                for i in range(0,len(rho_mg)):
                    # look up what thickness number this rho is related to
                    for j in range(0,len(th)):
                        # compare thickness to th list
                        if LThick == th[j][0]:
                            # assign rho_mg based on the rho_mg of the thickness
                            rho_mg[i] = mgDict['rho'][j]                   
                
    
        nd = [] # list of dictionaries for new design dictionary sections part
        
        for j,ltyp in enumerate(LTypes):
            st =  deepcopy(oldLine.lineTypes)
            for k in st:
                if st[k]['name'] == ltyp:
                    ltyp = k
            # add in information for each line type without marine growth
            EA.append(st[ltyp]['EA'])
            m.append(st[ltyp]['m'])
            d_ve_old.append(st[ltyp]['d_vol'])
            # new dictionary for this line type
            nd.append({'type':{}, 'L':{}}) # new design dictionary
            ndt = nd[j]['type']
            
            # load in line props from MoorProps
            opt = helpers.loadLineProps(None)
            
            if 'd_nom' in st[ltyp]:
                d_nom_old[j] = st[ltyp]['d_nom']
                # get ratio between ve and nom diameter normally
                ve_nom_adjust[j] = d_ve_old[j]/d_nom_old[j]
            elif Mats[j] in opt:
                # get ratio between ve and nom diameter from MoorProps yaml                
                ve_nom_adjust[j] = opt[Mats[j]]['dvol_dnom']
            # get cd and cdAx if given, or assign to default value
            if Mats[j] in opt and not 'Cd' in st[ltyp]:
                cd.append(opt[Mats[j]]['Cd'])
            elif 'Cd' in st[ltyp]:
                cd.append(st[ltyp]['Cd'])
            else:
                #print('No Cd given in line type and material not found in MoorProps yaml. Default Cd of 1 will be used.')
                cd.append(2)
            if Mats[j] in opt and not 'CdAx' in st[ltyp]:
                cdAx.append(opt[Mats[j]]['CdAx'])
            elif 'CdAx' in st[ltyp]:
                cdAx.append(st[ltyp]['CdAx'])
            else:
                #print('No CdAx given in line type and material not found in MoorProps yaml. Default CdAx of 0.5 will be used.')
                cdAx.append(0.5)
            
            if LThick[j] == 0:
                nd[j]['type'] = deepcopy(st[ltyp])
                nd[j]['type']['name'] = j
            else:
                # get mu for material
                if Mats[j] == 'chain' or Mats[j] == 'chain_studlink':
                    mu_mg[j] = 2
                else:
                    mu_mg[j] = 1
                
                # re-form dictionaries with marine growth values            
                # calculate nominal diameter
                d_nom_old[j] = d_ve_old[j]/ve_nom_adjust[j] # m
                
                # calculate new line diameter that includes marine growth
                ndt['d_nom'] = float(d_nom_old[j]+2*LThick[j]) #m
                
                # calculate the new mass per meter including marine growth
                growthMass = np.pi/4*(ndt['d_nom']**2-d_nom_old[j]**2)*rho_mg[j]*mu_mg[j] # marine growth mass
                ndt['m'] =  float(growthMass + m[j]) # kg/m (total mass)
                
                # calculate the submerged weight per meter including marine growth
                ndt['w'] = float(growthMass*(1-self.rho/rho_mg[j])*self.g + (m[j]-np.pi/4*d_ve_old[j]**2*self.rho)*self.g) # N/m
                
                # calculate new volume-equivalent diameter (cannot use regular chain/polyester conversion because of added marine growth)
                ndt['d_vol'] = np.sqrt(4*((ndt['m']*self.g-ndt['w'])/self.rho/self.g)/np.pi)
                
                # calculate new increased drag coefficient from marine growth
                # convert cd to cd for nominal diameter, then multiply by inverse of new ve_nom_adjust (ratio of d_nom with mg to d_ve with mg) to return to cd for volume equivalent diameter
                ndt['Cd'] = float(cd[j]*ve_nom_adjust[j]*(ndt['d_nom']/ndt['d_vol']))
                ndt['CdAx'] = float(cdAx[j]*ve_nom_adjust[j]*(ndt['d_nom']/ndt['d_vol']))
                
                # add line details to dictionary
                ndt['material'] = Mats[j]
                ndt['name'] = str(j)
                if 'MBL' in oldLine.lineTypes[ltyp]:
                    ndt['MBL'] = oldLine.lineTypes[ltyp]['MBL']
                if 'cost' in oldLine.lineTypes[ltyp]:
                    ndt['cost'] = oldLine.lineTypes[ltyp]['cost']
                ndt['EA'] = EA[j]
                if 'EAd' in oldLine.lineTypes[ltyp]:
                    ndt['EAd'] = oldLine.lineTypes[ltyp]['EAd']
            # add lengths                 
            nd[j]['L'] = Lengths[j]
        
        # # overwrite design dictionary with new dictionary
        # self.dd['sections'] = nd
        # # overwrite connectorList with new connectorList
        # self.connectorList = connList
        
        # fill out rest of new design dictionary
        nd1 = deepcopy(self.dd)
        nd1['sections'] = nd
        nd1['connectors'] = connList
        
        # call createSubsystem() to make moorpy subsystem with marine growth
        if self.shared:
            self.createSubsystem(case=1,dd=nd1,pristine=0)
        else:
            self.createSubsystem(dd=nd1,pristine=0)
            
        return(changeDepths,changePoints)


    def addCorrosion(self,corrosion_mm=10):
        '''
        Calculates MBL of chain line with corrosion included

        Parameters
        ----------
        corrosion_mm : TYPE, optional
            DESCRIPTION. The default is 10.

        Returns
        -------
        None.

        '''
        for i,sec in enumerate(self.dd['sections']):
            if sec['material']=='chain':
                MBL_cor = sec['MBL']*( (sec['d_nom']-(corrosion_mm/1000))/sec['d_nom'] )**2  # corroded MBL
            else:
                MBL_cor = sec['MBL']


    def getEnvelope(self,ang_spacing=45,SFs=True):
        '''Computes the motion envelope of the Mooring based on the watch 
        circle(s) of what it's attached to. If those aren't already 
        calculated, this method will call the relevant getWatchCircle method.
        
        Parameters
        ----------
        ang_spacing : int
            Spacing between angles to calculate for watch circle
        SFs : bool
            Whether or not to calculate safety factors etc for the lines in the watch circle function
        
        Returns
        x: list of x-coordinates of motion envelope
        y: list of y-coordinates of motion envelope
        -------
        '''
        
        from famodel.platform.platform import Platform
        import shapely as sh
        
        pts = []  # list of watch circle points (on either end)
        
        # Loop through ends A and B
        for i in range(2):
            att = self.attached_to[i] # the object the end is attached to
            
            # Only consider a watch circle if it's attached to a platform
            if isinstance(att, Platform):
                print('Found a platform')
                
                # If no watch circle saved, compute it 
                if not 'mean' in att.envelopes:  
                    att.getWatchCircle(ang_spacing=ang_spacing,SFs=SFs)
                
                # Add each vertex of the watch circle to the list
                for x, y in zip(att.envelopes['mean']['x'], 
                                att.envelopes['mean']['y']):
                    
                    pts.append((x, y))
            
            # Otherwise add a single point (e.g. for an anchor)
            else:
                pts.append((att.r[0], att.r[1]))
                
        # reset platform to equilibrium 
        
        # Convert to shapely polygon, then get convex hull
        poly = sh.Polygon(pts)  # make polygon object
        hull = sh.convex_hull(poly)  # convex hull (wraps around points)
        
        # Convert to np array and save in object envelope
        x = np.array(hull.exterior.xy[0])
        y = np.array(hull.exterior.xy[1])
        self.envelopes['mean'] = dict(x=x, y=y, shape=hull)
        
        return(x,y)

