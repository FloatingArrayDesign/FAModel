# class for a dynamic subsea power cable

import numpy as np
from copy import deepcopy
from moorpy.subsystem import Subsystem
from moorpy import helpers
from famodel.mooring.connector import Connector, Section
from famodel.famodel_base import Edge
from famodel.cables import cable_properties as cp
from famodel.cables.components import Joint

class DynamicCable(Edge):
    '''
    Class for a dynamic power cable. It inherits from Cable(Edge, dict)
    which describes the bare uniform cable before accessories are added.
    A DynamicCable object will likely be within a SubseaCable object, which could
    also include a StaticCable.
    End A of the dynamic cable could attach to a subsea Joint, or in the case
    of a suspended cable it would attach to another Platform.
    '''
    
    def __init__(self, id, dd=None, subsystem=None, rA=[0,0,0], rB=[0,0,0],
                 rad_anch=None, rad_fair=58, z_anch=-100, z_fair=-14, 
                 rho=1025, g=9.81,span=2000,length=2200,A=None, 
                 zJTube=-30,voltage=66,powerRating=None,cable_type=None,
                 headingA=None,headingB=None,buoyancy_sections=None,shared=0,**kwargs):
        '''
        Parameters
        ----------
        
        Unused to delete (?): anchor, A=None, conductorSize=None, type='dynamic',
        zJTube=-30,voltage=66,powerRating=None,cable_type=None, zAnchor
        
        '''
        Edge.__init__(self, id)  # initialize Edge base class
        # Design description dictionary for this dynamic cable
        self.dd = dd
        
        self.n_sec = 1
        
        if 'span' in self.dd:
            self.span = self.dd['span']
        else:
            self.span = None# <<< what about self.dd['span']? TODO: ensure they stay consistent
        self.depth = -z_anch  # <<< may want to make 'depth' an input
        
        # Store the cable type properties dict here for easy access (temporary - may be an inconsistent coding choice)
        self.cableType = self.makeCableType(self.dd['cable_type'])  # Process/check it into a new dict
        # ^^^ curiuos if we can simplify this
        
        # Save some constants for use when computing buoyancy module stuff
        self.d0 = self.cableType['d_vol']  # diameter of bare dynamic cable
        self.m0 = self.cableType['m']      # mass/m of bare dynamic cable
        self.w0 = self.cableType['w']      # weight/m of bare dynamic cable
        
        # Turn what's in dd into a list of buoyancy section info dicts
        if 'buoyancy_sections' in self.dd:
            for i, bs in enumerate(self.dd['buoyancy_sections']):
                for key in ['L_mid', 'module_props', 'N_modules', 'spacing']:
                    if not key in bs:  # make sure no entry is missing
                        raise Exception(f'Required entry {key} not found in buoyancy_sections entry')
                    #could also check for buoyancy section position conflicts/overlaps
                
        
        # MoorPy subsystem that corresponds to the dynamic cable
        self.ss = subsystem
        self.ss_mod = None # modified subsystem that could include marine growth etc
        
        # end point absolute coordinates, to be set later
        self.rA = rA
        self.rB = rB
        
        if 'headingA' in self.dd:
            self.headingA = self.dd['headingA']
            if 'headingB' in self.dd:
                self.headingB = self.dd['headingB']
            else:
                self.headingB = 0
        elif 'headingB' in self.dd:
            self.headingB = self.dd['headingB'] # <<< ??
            self.headingA = 0

        else:
            self.headingA = 0
            self.headingB = 0
        
        
        if 'L' in self.dd:
            self.L = self.dd['L']
        else:
            self.L = length
        # self.dd['lenght'] <<< also/or use this?
        
        # relative positions (variables could be renamed)
        self.rad_fair = rad_fair
        if not rad_anch:
            self.rad_anch = self.dd['span'] + self.rad_fair
        else:
            self.rad_anch = rad_anch
        self.z_anch   = z_anch  
        self.z_fair   = z_fair
        
        self.adjuster = None  # custom function that can adjust the mooring
        
        self.shared = shared # boolean for if a fully suspended cable  (2 = symmetric)
        
        # relevant site info
        self.rho = rho
        self.g = g
        
        # Dictionaries for addition information
        self.loads = {}
        self.reliability = {}
        self.cost = {}
        self.failure_probability = {}
        
    
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
        if 'd' in cableType:
            cableType['d_vol'] = cableType.pop('d')
        return(cableType)
    
    
    def updateSubsystem(self, pristine=True):
        '''Adjusts the subsystem properties when the buoyancy section info changes.
        The contents of self.dd['buoyancy_sections'] should already be changed before
        calling this method. The length (L) could also change.
        
        Parameters
        ----------
        pristine : bool
            True if using unmodified form (default), false if modified 
            (e.g. marine growth), in which case the subsystem is separate.
        '''
        
        # make sure subsystem already exists
        if pristine:           
            if not self.ss:
                raise Exception('Subsystem does not yet exist in this DynamicCable')
            ss = self.ss
        else:
            if not self.ss_mod:
                raise Exception('Subsystem does not yet exist in this DynamicCable')
            ss = self.ss_mod
        
        if not self.L == self.dd['length']:  # ensure consistency (currently redundant)
            breakpoint()
            print("Lengths not consistent in DynamicCable dd and .L")
        
        # make sure the number of buoyancy sections matches the subsystem
        if len(self.dd['buoyancy_sections'])*2 + 1 == self.ss.nLines:  
            cstart = 1  # typical case where buoyancy sections are in middle of cable
        elif len(self.dd['buoyancy_sections'])*2 == self.ss.nLines: 
            cstart = 0  # case where buoyancy section starts right at end A
        else:
            raise Exception("Number of buoyancy sections doesn't match subsystem")
            
            
        
        
        currentL = 0  # the length along the cable as we process the sections
        iLine = 0  # index of the current Line along the subsystem
        
        # Turn what's in dd into a list of buoyancy section info dicts
        for i, bs in enumerate(self.dd['buoyancy_sections']):
        
            # bs contains 'L_mid', 'module_props', 'N_modules', 'spacing'
            
            # call function to calculate equivalent buoyancy 
            Ls, m, w, d_outer = self.calcEquivBuoyancy(bs)
            
            # Save end A and B arc lengths of the buoyancy section
            bs['L_A'] = bs['L_mid'] - Ls/2
            bs['L_B'] = bs['L_mid'] + Ls/2
            
 
            
            # If this buoyancy section isn't at the very start of the cable
            if cstart:  
                
                halfLs = Ls/2  
                # if cstart (end A starts with a buoy section) check if it's the first section of a shared half line (section length is half of actual section length)
                
                
                iLine +=1 
                
                # >>> note: this approach clashes/overlaps with the 'case' approach - should pick one <<<
                
                # Set length of bare cable section before this buoyancy section
                #self.dd['sections'][iLine-1]['length'] = bs['L_mid'] - Ls/2 - currentL
                
                self.ss.lineList[iLine-1].setL(bs['L_mid'] - halfLs - currentL)
                currentL = bs['L_mid'] - halfLs # save the end location of the section
            else:
                if self.shared == 2 and i == 0:
                    # adjust halfLs to be full Ls (half of the total buoyant segment length)
                    halfLs = Ls 
            
            
            # update properties of the corresponding Subsystem Line
            self.ss.lineList[2*i+cstart].setL(Ls)
            #self.dd['sections'][iLine]['length'] = Ls
            
            self.ss.lineTypes[2*i+cstart]['m'] = m
            self.ss.lineTypes[2*i+cstart]['w'] = w
            self.ss.lineTypes[2*i+cstart]['d_vol'] = d_outer
            '''
            self.dd['sections'][iLine]['type']['m'] = m # likely redundant if linked to ss
            self.dd['sections'][iLine]['type']['w'] = w
            self.dd['sections'][iLine]['type']['d_vol'] = d_outer
            '''
            # update end location of the section 
            currentL += Ls 
           
            if i == len(self.dd['buoyancy_sections'])-1:
                # this is the last section - adjust cable length at the end
                L_end = self.L - bs['L_mid'] - halfLs

                self.ss.lineList[-1].setL(L_end)
                #self.dd['sections'][-1]['length'] = L_end
                currentL += L_end
                
            iLine +=1 
        
        print(f"Total DynamicCable.ss length is {sum([l.L0 for l in self.ss.lineList])}")
    
    
    def getBuoyancyOverlaps(self):
        '''Calculate the margin for any buoyancy sections overlapping with
        each other or with the cable ends. The returned values could be used
        as constraints to avoid overlap in an optimization.
        
        Returns
        -------
        array
            List of overlap margins (negative means overlap) of end A, of each
            successive buoyancy section (1-N), and then of end B.
        '''
        
        bss = self.dd['buoyancy_sections']
        
        overlap_margin = np.zeros(len(self.i_buoy)+1)
        
        overlap_margin[0] = bss[0]['L_A']
        
        for i in range(len(self.i_buoy)-1):
            overlap_margin[i] = bss[i+1]['L_A'] - bss[i]['L_B']
        
        overlap_margin[-1] = self.dd['length'] - bss[-1]['L_B']
        
        return overlap_margin
    
    
    def calcEquivBuoyancy(self,bs):
        '''Calculates new cable information that includes equivalent buoyancy of buoys

        Parameters
        ----------
        bs : dictionary
            Dictionary describing buoy information

        Returns
        -------
        L : float
            Length of buoyant cable section
        m : float
            Mass/length of buoyant cable section
        w : float
            Weight per unit length of buoyant cable section
        d_outer : float
            Volumetric diameter of buoyant cable section

        '''
        L = (bs['N_modules'] - 1) * bs['spacing']  # length of section approximating buoyancy section
        
        # compute what diameter of buoyancy module is needed to achieve this buoyancy per unit length
        d_inner = self.d0  # inner diameter for buoyancy module [m]
        rho_buoy = bs['module_props']['density']  # constant density of buoyancy modules [kg/m^3]
        v_buoy = bs['module_props']['volume']*bs['N_modules']  # displaced volume [m^3]
        
        # equivalent outer diameter of buoyancy section
        d_outer = np.sqrt(4*v_buoy/L/np.pi + d_inner**2) 
        
        # mass per meter of spread buoyancy module [kg/m]
        m_buoy = rho_buoy*v_buoy/L
        
        m = self.m0 + m_buoy  # mass per unit length of combined cable + spread buoyancy modules [kg/m]
        w = m*self.g - self.rho*(np.pi/4*(d_outer**2))*self.g  # weight per unit length [N/m]
        
        return(L,m,w,d_outer)
    
    
    def reposition(self, r_center=None, heading=None, project=None, degrees=False, **kwargs):
        '''Adjusts dynamic cable position based on changed platform location or
        heading. It can call a custom "adjuster" function if one is
        provided. Otherwise it will just update the end positions.
        
        Parameters
        ----------
        r_center
            The x, y coordinates of the platform (undisplaced) [m].
        heading : float
            The absolute heading of the mooring line [deg or rad] depending on
            degrees parameter (True or False).
        project : FAModel Project, optional
            A Project-type object for site-specific information used in custom
            mooring line adjustment functions (mooring.adjuster).
        **kwargs : dict
            Additional arguments passed through to the designer function.
        '''
        
        # Adjust heading if provided
        if not heading == None:
            if degrees:
                self.heading = np.radians(heading)
            else:
                self.heading = heading
        
        phi = np.radians(90-self.heading) # heading in x-y radian convention [rad]
        
        # heading 2D unit vector
        u = np.array([np.cos(phi), np.sin(phi)])
        #print(u)
        r_center = np.array(r_center)[:2]
        # Set the updated fairlead location
        self.setEndPosition(np.hstack([r_center + self.rad_fair*u, self.z_fair]), 'b')
        
        
        # Run custom function to update the mooring design (and anchor position)
        # this would also szie the anchor maybe?
        if self.adjuster:
            self.adjuster(self, project, r_center, u, **kwargs)
        
        else: # otherwise just set the anchor position based on a set spacing
            self.setEndPosition(np.hstack([r_center + self.rad_anch*u, self.z_anch]), 'a', sink=True)
        
    
    
    def setEndPosition(self, r, end, sink=False):
        '''Set the position of an end of the dynamic cable.
        
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
        
        # >>> UPDATE FROM CABLEDESIGN2 <<<
        # go through dictionary and add up costs
        cost = 0
        # get cost of bare cable
        if 'cost' in self.dd['cable_type']:
            cost += self.dd['cable_type']['cost']*self.L 
        # get cost of buoyancy modules
        for bs in self.dd['buoyancy_sections']:
            if 'cost' in bs['module_props']:
                cost += bs['module_props']['cost']*bs['N_modules']
        # sum up the costs in the dictionary and return
        self.cost['total'] = cost
        return(cost)
        #return sum(self.cost.values()) 
    
    def updateTensions(self):
        ''' Gets tensions from subsystem and updates the max tensions dictionary if it is larger than a previous tension
        '''
        if not 'TAmax' in self.loads:
            self.loads['TAmax'] = 0
        if not 'TBmax' in self.loads:
            self.loads['TBmax'] = 0
        if not self.ss:
            if self.shared:
                self.createSubsystem(case=1)
            else:
                self.createSubsystem()
        # get anchor tensions
        if abs(self.ss.TA) > self.loads['TAmax']:
            self.loads['TAmax'] = deepcopy(self.ss.TA)
        # get TB tensions
        if abs(self.ss.TB) > self.loads['TBmax']:
            self.loads['TBmax'] = deepcopy(self.loads.TB)
            
        return(self.loads['TAmax'],self.loads['TBmax'])
    
    def createSubsystem(self, case=0,pristine=True,dd=None):
        ''' Creates a subsystem for cable and buoyancy section(s) configuration from the design dictionary
        
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
        '''
        
        # >>> USE SOME METHODS FROM CABLEDESIGN2 TO CONVERT BUOYANCY MODULES TO SECTIONS! <<<
        # set design dictionary as self.dd if none given, same with connectorList
        if not dd:
            dd = self.dd
            
        # check that cableType property is up to date
        cableType = self.makeCableType(dd['cable_type'])

        # Set up cable sections that will correspond to MoorPy Line objects
        #dd['sections'] = []  # initialize a list of these sections
        currentL = 0  # the length along the cable as we process the sections
        
        # lists used to make the subsystem
        lengths = []
        types = []
        
        # If no buoyancy sections, it's just one section of the bare cable
        if not 'buoyancy_sections' in dd:
            #self.dd['sections'].append({'type':self.cableType,'length':self.L})
            types.append(cableType)
            lengths.append(self.L) 
        
        if 'sections' in dd and dd['sections']: # this will be the case for marine growth or possibly other cases
            for i, sec in enumerate(dd['sections']):
                types.append(sec['type'])
                lengths.append(sec['length'])
        elif 'buoyancy_sections' in dd:       
            # Parse buoyancy sections to compute their properties and all lengths
            for i, bs in enumerate(dd['buoyancy_sections']):
                # check which end to start from ( need to flip locations of buoyancy modules if end A is at joint and end B is at platform)
                if isinstance(self.attached_to[0],Joint):
                    L_mid = self.L-bs['L_mid']
                else:
                    L_mid = bs['L_mid']
                # get buoyancy section information
                Ls,m,w,d_vol = self.calcEquivBuoyancy(bs) 
                
                # If this buoyancy section isn't at the very start of the cable
                if i > 0 or Ls/2 < L_mid:  
                    # Add a bare cable section before this buoyancy section
                    lengths.append(L_mid - Ls/2 - currentL)
                    types.append(cableType)
                    currentL = L_mid - Ls/2 # save the end location of the section
                    
                # create buoyancy section equivalent cable type dict
                buoyCableType = deepcopy(cableType)
                buoyCableType['d_vol'] = d_vol
                buoyCableType['m'] = m
                buoyCableType['w'] = w
                buoyCableType['name'] = cableType['name']+'_'+'buoy'+str(i)
                
                types.append(buoyCableType)
                lengths.append(Ls)
                
                #dd['sections'].append({'type':buoyCableType})
                #dd['sections'][-1]['length'] = Ls
                
                # update end location of the section 
                currentL += Ls 
               
                if i == len(dd['buoyancy_sections'])-1:
                    # this is the last section - add cable length at the end
                    #dd['sections'].append({'type':self.cableType})
                    #dd['sections'][-1]['length'] = self.L - bs['L_mid'] - Ls/2
                    types.append(cableType)
                    lengths.append(self.L - L_mid - Ls/2)
                    
                    currentL += self.L - L_mid - Ls/2
        '''
        currentL = 0
        # add buoyancy sections to design dictionary if it doesn't already exist
        if not 'sections' in dd:
            dd['sections'] = []
            if not self.dd['buoyancy_sections']:
                dd['sections'].append({'type':self.cableType,'length':self.L})
            else:
                for i,bs in enumerate(self.dd['buoyancy_sections']):
                    # get buoyancy section information
                    Ls,m,w,d_vol = self.calcEquivBuoyancy(bs) 
                    
                    if i == 0 and bs['L_mid']<Ls/2:
                        pass
                    else:            
                        # cable doesn't start with a buoyancy section of this is a mid section - need to add cable length section before buoyancy section
                        dd['sections'].append({'type':self.cableType})                
                        dd['sections'][-1]['length'] = bs['L_mid'] - Ls/2 - currentL
                        currentL = bs['L_mid'] - Ls/2 # save the end location of the section
                    # create buoyancy section equivalent cable type dict
                    buoyCableType = deepcopy(self.cableType)
                    buoyCableType['d_vol'] = d_vol
                    buoyCableType['m'] = m
                    buoyCableType['w'] = w
                    buoyCableType['name'] = self.cableType['name']+'_'+'buoy'+str(i)
                    dd['sections'].append({'type':buoyCableType})
                    dd['sections'][-1]['length'] = Ls
                    # update end location of the section 
                    currentL += Ls 
                   
                    if i == len(self.dd['buoyancy_sections'])-1:
                        # this is the last section - add cable length at the end
                        dd['sections'].append({'type':self.cableType})
                        dd['sections'][-1]['length'] = self.L - bs['L_mid'] - Ls/2
                        currentL += self.L - bs['L_mid'] - Ls/2
        '''        
        # # check if a subsystem already exists
        # if pristine:
        #     if self.ss:
        #         print('A subsystem for this Dynamic cable class instance already exists, this will be overwritten.')
        # else:
        #     #breakpoint()
        #     if self.ss_mod:
                # print('A modified subsystem for this Dynamic cable class instance already exists, this will be overwritten.')
        ss=Subsystem(depth=self.depth, rho=self.rho, g=self.g, 
                          span=self.span, rad_fair=self.rad_fair,
                          z_fair=self.z_fair)
        '''
        lengths = []
        types = []
        nsegs = []
        # run through each line section and collect the length and type
        for sec in dd['sections']:
            lengths.append(sec['length']) 
            types.append(sec['type']['name'])
            ss.lineTypes[types[-1]] = sec['type']  # points to existing type dict in self.dd for now
            nsegs.append(np.round(lengths[-1]/3))
            if nsegs[-1] > 25:
                nsegs[-1] = 25
            elif nsegs[-1] < 3:
                nsegs[-1] = 3
        '''
        # make the lines and set the points
        # breakpoint()
        ss.makeGeneric(lengths,types,suspended=case)
        ss.setEndPosition(self.rA,endB=0)
        ss.setEndPosition(self.rB,endB=1)
        #breakpoint()
        # note: next bit has similar code/function as Connector.makeMoorPyConnector <<<
        
        # # add in connector info to subsystem points
        # if case == 0: # has an anchor - need to ignore connection for first point
        #     startNum = 1
        # else: # no anchor - need to include all connections
        #     startNum = 0 

        # for i in range(startNum,len(self.ss.pointList)):                               
        #     point = self.ss.pointList[i]
        #     point.m = self.dd['connectors'][i]['m']
        #     point.v = self.dd['connectors'][i]['v']
        #     point.CdA = self.dd['connectors'][i]['CdA']
        # solve the system
        ss.initialize()
        try:
            ss.staticSolve(maxIter = 5000)
        except:
            breakpoint()
        #breakpoint()
        
        # save it in the object
        if pristine:
            self.ss = ss
            return(self.ss)
        else:
            self.ss_mod = ss
            return(self.ss_mod)
                  
    def addMarineGrowth(self, mgDict,updateDepths=False,tol=2):
        '''Creates a new design dictionary (does not overwrite old one) to account for marine 
        growth on the subystem, then calls createSubsystem() to recreate the cable

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
                buoy_th : list of thicknesses associated with buoyancy sections. This is used to simplify the marine growth modeling on buoys.
                        * In order of sections from end A to end B (also the order of sections listed in buoyancy_sections list in the dd)
                        example, if 3 buoyancy sections: - 0.1
                                                         - 0.05
                                                         - 0.1
                rho : list of densities for each thickness, or one density for all thicknesses, [kg/m^3] (optional - default is 1325 kg/m^3)
                }
        Returns
        -------
        changePoints : list
            List of point indices in the moorpy subsystem that are at the changeDepth
        changeDepths : list
            List of cutoff depths the changePoints should be located at

        '''
        def getMG(mgDict):
            # create a reference subsystem if it doesn't already exist
            if not self.ss:
                self.createSubsystem(pristine=1)
            # set location of reference subsystem
            oldLine = self.ss
            # set up variables
            LTypes = [] # list of line types for new lines (types listed are from reference object)
            Lengths = [] # lengths of each section for new line
            # Mats = [] # materials list for new line        
            sCount = [] # new list of connectors (need to add empty connector objects in between changeDepths)
            LThick = [] # list of mg thicknesses for new lines
            ln_raw = [] # list of line lengths from rA to current split in line (needed to determine length of new sections when there are multiple splits in one line section)
            # set up variables needed to check before/after of changeDepths
            changePoints = []
            changeDepths = [] # index of list that has the corresponding changeDepth
            
            sCount = 0
            buoyCount = -1
            # go through each line section
            for i in range(0,len(oldLine.lineList)):
                
                slthick = [] # mg thicknesses for the section (if rA is above rB, needs to be flipped before being added to full subsystem list LThick)
                slength = [] # line lengths for the section (if rA is above rB, needs to be flipped before being added to full subsystem list)
                schangeDepth = [] # changeDepths for the section (if rA is above rB, needs to be flipped before being added to full subsystem list)
                # set reference subsystem line section location
                ssLine = oldLine.lineList[i]
                
                # determine if it's a buoy section
                if ssLine.type['w'] < 0:
                    calcCDs = 0 # don't go through loop to find change depths of this section
                    # add to buoyCounter
                    buoyCount += 1
                else:
                    calcCDs = 1
                
                
                # add line material, type to list
                # Mats.append(ssLine.type['material'])
                LTypes.append(ssLine.type['name'])
                
                if calcCDs:
                           
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
                    while rAth==0: # and count1 < len(th):
                        if count1 == len(th):
                            breakpoint()
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
                            # Mats.append(ssLine.type['material'])
                            # add an empty connector object to list for split location
                            sCount += 1
                            changePoints.append(sCount)
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
                                        #print('\n\n',old_line,'\n\n')
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
                                
                else:
                    # this is a buoy section, add the prescribed buoy marine growth thickness
                    if isinstance(mgDict['buoy_th'],list):
                        # add the buoy thickness for this specific buoyancy section
                        LThick.append(mgDict['buoy_th'][buoyCount])
                    else:
                        # add the single buoy thickness
                        LThick.append(mgDict['buoy_th'])
                    
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
                
                sCount += 1
                # changePoints.append(sCount)
                                   
            # Set up list variables for pristine line info
            EA = []
            m = []
            d_ve_old = []
            cd = []
            cdAx = []
            d_nom_old = []
            ve_nom_adjust = []
                                                
            # create arrays
            # d_nom_old = np.zeros((len(LTypes),1))        
            # ve_nom_adjust = np.zeros((len(LTypes),1))
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
                        for j in range(0,len(LThick)):
                            thind = np.where(th[:][0]==LThick[j])
                            # assign rho_mg based on the rho_mg of the thickness
                            rho_mg[i] = mgDict['rho'][thind]             
                    
        
            nd = [] # list of dictionaries for new design dictionary sections part
            
            for j,ltyp in enumerate(LTypes):
                st =  deepcopy(oldLine.lineTypes)
                # find which old line linetypes have the same name as LTypes
                for k in st:
                    if st[k]['name'] == ltyp:
                        linekey = k
                # add in information for each line type without marine growth
                EA.append(st[linekey]['EA'])
                m.append(st[linekey]['m'])
                d_ve_old.append(st[linekey]['d_vol'])
                if not 'd_nom' in st[linekey]:
                    # calculate d_nom ***assume d_vol = d_nom
                    d_nom_old.append(st[linekey]['d_vol'])
                else:               
                    d_nom_old.append(st[linekey]['d_nom'])
                ve_nom_adjust.append(d_ve_old[-1]/d_nom_old[-1])
                # new dictionary for this line type
                nd.append({'type':{}, 'length':{}}) # new design dictionary
                ndt = nd[j]['type']

                if 'Cd' in st[linekey]:
                    cd.append(st[linekey]['Cd'])
                else:
                    cd.append(2)
                if 'CdAx' in st[linekey]:
                    cdAx.append(st[linekey]['CdAx'])
                else:
                    cdAx.append(0.5)
                if LThick[j] == 0:
                    nd[j]['type'] = deepcopy(st[linekey])
                    #nd[j]['type']['name'] = j
                else:
                    # calculate nominal diameter
                    d_nom_old[j] = d_ve_old[j]/ve_nom_adjust[j] # m
                    
                    # re-form dictionaries with marine growth values
                    
                    # calculate new line diameter that includes marine growth
                    ndt['d_nom'] = float(d_nom_old[j]+2*LThick[j]) #m
                    
                    
                    # check if this is a buoyant section or regular section
                    if st[linekey]['w'] < 0:
                        # this is a buoyant section - add marine growth on modules
                        
                        # step 1: get module info from design dictionary
                        bs = self.dd['buoyancy_sections'][buoyCount]
                        ls = bs['spacing'] # get length of segment between centers of 2 modules
                        Ls = bs['spacing']*(bs['N_modules']-1) # length of buoyancy segment
                        N_b = bs['N_modules']
                        d_b = bs['module_props']['d']
                        v_b = bs['module_props']['volume']
                        l_b = bs['module_props']['l']
                        m_b = bs['module_props']['m'] # total mass of buoy (not mass/m)
                        w_b = bs['module_props']['w']
                        # step 2: calc volume & weight of marine growth on one buoy (assume cylindrical buoy)
                        v_bmg = np.pi/4*((d_b+2*LThick[j])**2-d_b**2)*l_b 
                        v_bmgs = v_bmg + np.pi/2*((d_b+2*LThick[j])**2-self.d0**2)*LThick[j] # around cylinder + sides of cylinder
                        w_bmg = float(rho_mg[j]-self.rho)*self.g*v_bmgs
                        # step 3: calc volume & weight of marine growth on bare cable segment between 2 buoys (segment length - buoy length - 2*mg thickness from side wall marine growth of buoy)
                        v_bare_mg = np.pi/4*((self.d0+2*LThick[j])**2-self.d0**2)*(ls-l_b-2*LThick[j])*(N_b-1)
                        w_bare_mg = float(rho_mg[j]-self.rho)*self.g*v_bare_mg
                        # step 4: calc total weight of segment, including mg, buoy weight, and cable weight (buoy weight + cable weight combined in listed weight for the line section)
                        w_seg = w_bmg*N_b + w_bare_mg + st[linekey]['w']*Ls
                        # step 5: calc total weight of section/m
                        ndt['w'] = float(w_seg/Ls)
                        # get section mass in case it's needed
                        m_seg = v_bmg*rho_mg[j]*N_b + v_bare_mg*rho_mg[j] + st[linekey]['m']*Ls
                        ndt['m'] = float(m_seg/Ls)
                        
                        
                    else:                  
                        mu_mg[j] = 1
     
                        # calculate the new mass per meter including marine growth
                        growthMass = np.pi/4*(ndt['d_nom']**2-d_nom_old[j]**2)*rho_mg[j]*mu_mg[j] # marine growth mass
                        ndt['m'] =  float(growthMass + m[j]) # kg/m (total mass)
                        
                        # calculate the submerged weight per meter including marine growth
                        ndt['w'] = float(growthMass*(1-self.rho/rho_mg[j])*self.g + (m[j]-np.pi/4*d_ve_old[j]**2*self.rho)*self.g) # N/m
                        
                    # calculate new volume-equivalent diameter (cannot use regular chain/polyester conversion because of added marine growth)
                    ndt['d_vol'] = float(np.sqrt(4*((ndt['m']*self.g-ndt['w'])/self.rho/self.g)/np.pi))
                    
                    # calculate new increased drag coefficient from marine growth
                    # convert cd to cd for nominal diameter, then multiply by inverse of new ve_nom_adjust (ratio of d_nom with mg to d_ve with mg) to return to cd for volume equivalent diameter
                    ndt['Cd'] = float(cd[j]*ve_nom_adjust[j]*(ndt['d_nom']/ndt['d_vol']))
                    ndt['CdAx'] = float(cdAx[j]*ve_nom_adjust[j]*(ndt['d_nom']/ndt['d_vol']))
                    
                    # add line details to dictionary
                    ndt['name'] = oldLine.lineTypes[linekey]['name']
                    if 'MBL' in oldLine.lineTypes[linekey]:
                        ndt['MBL'] = oldLine.lineTypes[linekey]['MBL']
                    if 'MBR' in oldLine.lineTypes[linekey]:
                        ndt['MBR'] = oldLine.lineTypes[linekey]['MBR']
                    if 'cost' in oldLine.lineTypes[linekey]:
                        ndt['cost'] = oldLine.lineTypes[linekey]['cost']
                    ndt['EA'] = EA[j]
                    if 'EAd' in oldLine.lineTypes[linekey]:
                        ndt['EAd'] = oldLine.lineTypes[linekey]['EAd']
                    if 'EI' in oldLine.lineTypes[linekey]:
                        ndt['EI'] = oldLine.lineTypes[linekey]['EI']
                # add lengths                 
                nd[j]['length'] = Lengths[j]
            
            # fill out rest of new design dictionary
            nd1 = self.dd
            nd1['sections'] = nd
            
            # call createSubsystem() to make moorpy subsystem with marine growth
            if self.shared>0:
                self.createSubsystem(case=int(self.shared),dd=nd1,pristine=0)
            else:
                self.createSubsystem(dd=nd1,pristine=0)
                
            return(changeDepths,changePoints)
         
        ############################################################################################################################
        if updateDepths:
            mgDict1 = deepcopy(mgDict)
            cEq = tol + 1 
            ct = 0
            while(cEq>tol):
                changeDepths,changePoints = getMG(mgDict1)
                D = []
                for d in range(0,len(changeDepths)):
                    diff = mgDict['th'][changeDepths[d][0]][changeDepths[d][1]] - self.ss_mod.pointList[changePoints[d]].r[2]
                    D.append(diff)
                cEq = np.mean(D)
                # adjust depth to change based on difference between actual and desired change depth
                if cEq>tol:
                    mgDict1['th'][0][2] = mgDict1['th'][0][2] + cEq
                    for j in range(1,len(mgDict['th'])):
                        for k in range(1,3):
                            if ct < 4 and abs(cEq)<12:
                                mgDict1['th'][j][k] = mgDict1['th'][j][k] + cEq
                            elif (ct >= 4 and ct < 9) or abs(cEq)>=12:
                                # could be ping-ponging between two different things, try adding half
                                mgDict1['th'][j][k] = mgDict1['th'][j][k] + 0.5*cEq
                print('average difference between expected and actual change depth is: ',cEq)
                if ct > 10:
                    print('Could not meet tolerance in 10 attempts. Exiting loop')
                    break
                ct += 1
        else:
            changeDepths,changePoints = getMG(mgDict)
            return(changeDepths,changePoints)
            
            
        
    
    def symmetricalMirror(self,create_ss=False):
        '''
        Mirror a symmetrical cable design dictionary, allow to overwrite ss
        
        Parameters
        ----------
        create_ss : bool, optional
            Controls whether or not to call createSubsystem() function to make a new subsystem after updating dd

        Returns
        -------
        None.

        '''
        from copy import deepcopy
        if 'sections' in self.dd:
            sections1 = deepcopy(self.dd['sections'])
            sections1.reverse()
            newsections = []
            newsections.extend(sections1)
            newsections.extend(deepcopy(self.dd['sections']))
            self.dd['sections'] = newsections
        else:
            oldbs = self.dd['buoyancy_sections']
            # get half length
            halfL = self.dd['length']
            # double dd['length']
            self.dd['length'] = halfL*2 
            
            if oldbs[0]['L_mid'] == 0:
                # ends on a buoyancy section
                endB = 1 
            else:
                # ends on a bare cable section
                endB = 0
            
            
            bs = []
            if endB:
                # update buoyancy sections
                # add in all pre-mirrored buoyancy sections (start from end instead of middle)
                for k in range(len(oldbs)-1,-1,-1):
                    bs.append(deepcopy(oldbs[k]))
                    # update L_mid - change coordinates to x starting from end A
                    bs[-1]['L_mid'] = halfL - bs[-1]['L_mid']
                    if k == 0:
                        N0 = bs[-1]['N_modules']
                        # double middle buoyancy section length while keeping buoyancy the same - spacing needs to be changed by
                        # (2N0-2)/(2N0-1)*sp0 where sp0 is old spacing, N0 is old # of buoyancy sections 
                        bs[-1]['spacing'] = bs[-1]['spacing']*(2*N0-2)/(2*N0-1)
                        # double number of modules to double the buoyancy
                        bs[-1]['N_modules'] = N0*2
                        
                # add in new buoyancy sections (mirrored, in same order as initial list (middle to end))
                for k in range(1,len(oldbs)):
                    bs.append(deepcopy(oldbs[k]))
                    # update L_mid - change coordinates to x starting from end A
                    bs[-1]['L_mid'] = halfL + bs[-1]['L_mid']
            else:     
                # update buoyancy sections
                # add in all pre-mirrored buoyancy sections (start from end instead of middle)
                for k in range(len(oldbs)-1,-1,-1):     
                    bs.append(deepcopy(oldbs[k]))
                    # update L_mid - change coordinates to x starting from end A
                    bs[-1]['L_mid'] = halfL - bs[-1]['L_mid']
                # add in new buoyancy sections (mirrored, in same order as initial list (middle to end))
                for k in range(0,len(oldbs)):
                    bs.append(deepcopy(oldbs[k]))
                    # update L_mid - change coordinates to x starting from end A, and push L_mid up by bsEnd (length of bare cable middle section pre-mirroring)
                    bs[-1]['L_mid'] = halfL + bs[-1]['L_mid']
                
                    
                
            self.dd['buoyancy_sections'] = bs
        # reset length
        self.L = self.L*2
        # # reset rA to -span and assume same z fairlead as rB
        self.rA = [-self.span,0,self.rB[2]]
        self.rB = [0,0,self.rB[2]]
        # call createSubsystem function if asked to
        if create_ss:
            self.createSubsystem(case=1,pristine=1)
        
        '''
        # determine if given sections ends on a buoyancy section or a regular section
        Ls = []
        bEnd = 0
        currentL = 0
        if 'buoyancy_sections' in self.dd:
            for i, bs in enumerate(self.dd['buoyancy_sections']):
                lenseg = bs['N_modules']*bs['spacing']
                if lenseg/2 >= bs['L_mid']:
                    #buoyancy section
                    Ls.append(lenseg)
                else:
                    Ls.append(bs['L_mid']-lenseg/2) # regular section
                    Ls.append(lenseg) # buoyancy section
                    
                bsEnd = bs['L_mid'] + bs['N_modules']*bs['spacing']/2
                currentL = bsEnd
                if bsEnd  >= self.span/2:
                    # ends on a buoyancy section - find out where in the buoyancy section
                    if bsEnd == self.span/2:
                        # double the buoyancy section
                        bs['L_mid'] = bsEnd
                        bs['N_modules'] = bs['N_modules']*2 
                        Ls[-1] = Ls[-1]*2
                    else:
                        # keep buoyancy section the same
                        pass
                    #  halfL = sum(Ls)
                    bEnd = 1 # flag boolean for ending on buoyancy section
                    # currentL = bsEnd + bs['N_modules']*bs['spacing']/2
                    
            
            if bEnd == 0:
                # halfL = sum(Ls)
                # determine length of last regular section based on end of last buoyancy section
                regL = self.span/2 - bsEnd
                Ls.append(regL*2)
                # double this length by adding next buoyancy section regL*2 away from last one
                self.dd['buoyancy_sections'].append(deepcopy(bs))
                # adjust L_mid
                newbs = self.dd['buoyancy_sections'][-1]
                newbs['L_mid'] = bsEnd + regL*2 + bs['N_modules']*bs['spacing']/2
                # currentL = newbs['L_mid'] + bs['N_modules']*bs['spacing']/2
                Ls.append(bs['N_modules']*bs['spacing'])
            
            if len(self.dd['buoyancy_sections']) >= 2+bEnd:
                for ii in range(0,len(self.dd['buoyancy_sections'])-2):
                    j = len(self.dd['buoyancy_sections']) - ii - 2 - bEnd
                    # add new buoyancy section
                    self.dd['buoyancy_sections'].append(deepcopy(self.dd['buoyancy_sections'][j]))
                    # adjust L_mid
                    newbs = self.dd['buoyancy_sections'][-1]
                    if bEnd:
                        Ls.append(-2-ii*2) # regular section
                        newbs['L_mid'] = self.dd['buoyancy_sections'][-2]['L_mid'] +  self.dd['buoyancy_sections'][-2]['N_sections']*self.dd['buoyancy_sections'][-2]['spacing']/2 + Ls[-1] +newbs['N_modules']*newbs['spacing']/2
                        Ls.append(newbs['N_modules']*newbs['spacing']) # buoyancy section
                    else:
                        Ls.append(Ls[-4-ii*2]) # regular line
                        newbs['L_mid'] = self.dd['buoyancy_sections'][-2]['L_mid'] + Ls[-2-ii*2]/2 + Ls[-1] + newbs['N_modules']*newbs['spacing']/2
                        Ls.append(newbs['N_modules']*newbs['spacing']) # buoyancy section
                    #newbs['L_mid'] = self.dd['buoyancy_sections'][-2]['L_mid'] + Ls[]
        else:
            # no buoyancy sections
            pass
        
        '''
        
            
            
            
            
            
            
    """  Could have similar marine growth method as Mooring, but need to also consider buoyancy modules
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
        
        return(changeDepths,changePoints)
    """
