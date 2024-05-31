# class for a dynamic subsea power cable

import numpy as np
from copy import deepcopy
from moorpy.subsystem import Subsystem
from moorpy import helpers
from famodel.mooring.connector import Connector, Section
from famodel.famodel_base import Edge
from famodel.cables import cable_properties as cp

class DynamicCable(Edge):
    '''
    Class for a dynamic power cable. It inherits from Cable(Edge, dict)
    which describes the bare uniform cable before accessories are added.
    A DynamicCable object will likely be within a SubseaCable object, which could
    also include a StaticCable.
    End A of the dynamic cable could attach to a subsea Joint, or in the case
    of a suspended cable it would attach to another Platform.
    '''
    
    def __init__(self, id, dd=None, subsystem=None, anchor=None, rA=[0,0,0], rB=[0,0,0],
                 rad_anch=500, rad_fair=58, z_anch=-100, z_fair=-14, 
                 rho=1025, g=9.81,span=2000,length=2200,A=None,conductorSize=None, 
                 type='dynamic',zJTube=-30,voltage=66,powerRating=None,cable_type=None,
                 headingA=None,headingB=None,buoyancy_sections=None,zAnchor=-200):
        '''
        Parameters
        ----------
 
        '''
        Edge.__init__(self, id)  # initialize Edge base class
        # Design description dictionary for this dynamic cable
        self.dd = dd
        
        self.n_sec = 1
        
        self.span = span

        # Store the cable type properties dict here for easy access (temporary - may be an inconsistent coding choice)
        self.cableType = self.makeCableType(self.dd['cable_type'])  # Process/check it into a new dict

        # Save some constants for use when computing buoyancy module stuff
        self.d0 = self.cableType['d_vol']  # diameter of bare dynamic cable
        self.m0 = self.cableType['m']      # mass/m of bare dynamic cable
        self.w0 = self.cableType['w']      # weight/m of bare dynamic cable
        
        # Turn what's in dd into a list of buoyancy section info dicts
        self.buoyancySections = []
        if 'buoyancy_sections' in self.dd:
            for i, bs in enumerate(self.dd['buoyancy_sections']):
                for key in ['L_mid', 'module_props', 'N_modules', 'spacing']:
                    if not key in bs:  # make sure no entry is missing
                        raise Exception(f'Required entry {key} not found in buoyancy_sections entry')
                
                self.buoyancySections.append(bs)
        
        # MoorPy subsystem that corresponds to the dynamic cable
        self.ss = subsystem
        
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
            self.headingB = self.dd['headingB']
            self.headingA = 0

        else:
            self.headingA = 0
            self.headingB = 0
        
        
        self.L = length
        
        # relative positions (variables could be renamed)
        self.rad_anch = rad_anch
        self.rad_fair = rad_fair
        self.z_anch   = z_anch  
        self.z_fair   = z_fair
        
        self.adjuster = None  # custom function that can adjust the mooring
        
        self.shared = False # boolean for if the mooring line is a fully suspended cable
        self.symmetric = False # boolean for if the mooring line is a symmetric suspended cable
        
        # relevant site info
        self.rho = rho
        self.g = g
        
        # Dictionaries for addition information
        self.loads = {}
        self.reliability = {}
        self.cost = {}
    
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
    
    def updateSubsystem(self):
        '''Adjusts the subsystem properties when the buoyancy section info changes'''
        
        # make sure subsystem already exists
        if not self.ss:
            raise Exception('Subsystem does not yet exist in this DynamicCable')
        
        # make sure the number of buoyancy sections matches the subsystem
        if len(self.buoyancySections)*2 + 1 == self.ss.nLines:  
            case = 1  # typical case where buoyancy sections are in middle of cable
        elif len(self.buoyancySections)*2 == self.ss.nLines: 
            case = 0  # case where buoyancy section starts right at end A
        else:
            raise Exception("Number of buoyancy sections doesn't match subsystem")
        
        # Turn what's in dd into a list of buoyancy section info dicts
        for i, bs in self.buoyancySections:
        
            # bs contains 'L_mid', 'module_props', 'N_modules', 'spacing'
            
            # call function to calculate equivalent buoyancy 
            L,m,w,d_outer = self.calcEquivBuoyancy(bs)
            
            # update properties of the corresponding Subsystem Line
            self.ss.lineList[2*i+case].setL(L)
            self.ss.lineTypes[2*i+case]['m'] = m
            self.ss.lineTypes[2*i+case]['w'] = w
            self.ss.lineTypes[2*i+case]['d_vol'] = d_outer
            
            # still need to update lengths of Lines in between buoyancy sections!! <<<
    
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
        L = bs['N_modules'] * bs['spacing']  # length of section approximating buoyancy section
        
        # compute what diameter of buoyancy module is needed to achieve this buoyancy per unit length
        d_inner = self.d0  # inner diameter for buoyancy module [m]
        rho_buoy = bs['module_props']['density']  # constant density of buoyancy modules [kg/m^3]
        v_buoy = bs['module_props']['volume']  # displaced volume [m^3]
        
        # equivalent outer diameter of buoyancy section
        d_outer = np.sqrt(4*v_buoy/bs['spacing']/np.pi + d_inner**2) 
        
        # mass per meter of spread buoyancy module [kg/m]
        m_buoy = rho_buoy*v_buoy/bs['spacing']
        
        m = self.m0 + m_buoy  # mass per unit length of combined cable + spread buoyancy modules [kg/m]
        w = m*self.g - self.rho*(np.pi/4*(d_outer**2))*self.g  # weight per unit length [N/m]
        
        return(L,m,w,d_outer)
    
    def reposition(self, r_center=None, heading=None, project=None, degrees=False, **kwargs):
        '''Adjusts mooring position based on changed platform location or
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
            
        # heading 2D unit vector
        u = np.array([np.cos(self.heading), np.sin(self.heading)])
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
        
        # >>> UPDATE FROM CABLEDESIGN2 <<<
        
        # sum up the costs in the dictionary and return
        return sum(self.cost.values()) 
    
    def createSubsystem(self, case=0):
        ''' Create a subsystem for a line configuration from the design dictionary
        
        Parameters
        ----------
        case : int
            Selector shared/suspended cases:
                - 0 (default): end A is on the seabed
                - 1: assembly is suspended and end A is at another floating system
                - 2: the assembly is suspended and assumed symmetric, end A is the midpoint
        '''
        
        # >>> USE SOME METHODS FROM CABLEDESIGN2 TO CONVERT BUOYANCY MODULES TO SECTIONS! <<<
        self.dd['sections'] = []
        currentL = 0
        if not self.buoyancySections:
            self.dd['sections'].append({'type':self.cableType,'length':self.L})
        for i,bs in enumerate(self.buoyancySections):
            # get buoyancy section information
            Ls,m,w,d_vol = self.calcEquivBuoyancy(bs) 
            
            if i == 0 and bs['L_mid']<Ls/2:
                pass
            else:            
                # cable doesn't start with a buoyancy section of this is a mid section - need to add cable length section before buoyancy section
                self.dd['sections'].append({'type':self.cableType})                
                self.dd['sections'][-1]['length'] = bs['L_mid'] - Ls/2 - currentL
                currentL = bs['L_mid'] - Ls/2 # save the end location of the section
            # create buoyancy section equivalent cable type dict
            buoyCableType = deepcopy(self.cableType)
            buoyCableType['d'] = d_vol
            buoyCableType['m'] = m
            buoyCableType['w'] = w
            buoyCableType['name'] = self.cableType['name']+'_'+'buoy'+str(i)
            self.dd['sections'].append({'type':buoyCableType})
            self.dd['sections'][-1]['length'] = Ls
            # update end location of the section 
            currentL += Ls 
           
            if i == len(self.buoyancySections)-1:
                # this is the last section - add cable length at the end
                self.dd['sections'].append({'type':self.cableType})
                self.dd['sections'][-1]['length'] = self.L - bs['L_mid'] - Ls/2
                currentL += self.L - bs['L_mid'] - Ls/2

        
        # store in self.dd['sections'] ???
        # check if a subsystem already exists
        if self.ss:
            print('A subsystem for this Dynamic cable class instance already exists, this will be overwritten.')
            
        self.ss=Subsystem(depth=-self.dd['zAnchor'], rho=self.rho, g=self.g, span=self.dd['span'], rBfair=self.rB)
        lengths = []
        types = []
        nsegs = []
        # run through each line section and collect the length and type
        for sec in self.dd['sections']:
            lengths.append(sec['length'])
            types.append(sec['type']['name'])
            self.ss.lineTypes[types[-1]] = sec['type']  # points to existing type dict in self.dd for now
            nsegs.append(np.round(lengths[-1]/3))
            if nsegs[-1] > 25:
                nsegs[-1] = 25
            elif nsegs[-1] < 3:
                nsegs[-1] = 3
        # make the lines and set the points
        self.ss.makeGeneric(lengths,types,suspended=case,nsegs=nsegs)
        self.ss.setEndPosition(self.rA,endB=0)
        self.ss.setEndPosition(self.rB,endB=1)
        
        
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
        self.ss.staticSolve()
        
        return(self.ss)      
    

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
