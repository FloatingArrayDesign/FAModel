# New version of LineDesign that uses Subsystem

import moorpy as mp # type: ignore
#import moordesign.MoorSolve as msolve
from famodel.design.fadsolvers import dsolve2, dopt2, doptPlot
from moorpy.MoorProps import getAnchorProps # type: ignore
from moorpy.helpers import (loadLineProps, getLineProps,  # type: ignore
                            rotationMatrix, getFromDict) 

from famodel.mooring.mooring import Mooring

import numpy as np
import matplotlib.pyplot as plt
import yaml
import time



class LineDesign(Mooring):
    '''
    The LineDesign class inherits from Mooring, which includes a Subsystem. For some cases where offsets 
    need to be computed, this class will also add a System to look at N line case.
    
    - The dynamic component of the design process will utilize various dynamic amplification factors (DAFs)
      that will be used an input into this class to design the Moorings
    - Design variables are imported through the 'allVars' variable,designated by the 'Xindices' variable,
      and stored in the "X" variable
    - The objective is calculated using the method self.objectiveFun to evaluate the cost of the Mooring
    - Constraints are initialized by a global dictionary of all possible constraints: self.confundict
      - Each constraint in confundict has a corresponding function (method) to evaluate that constraint
      - The key is a string (e.g. "min_lay_length") that the user can pass in through a list, with a corresponding
        number to designate the Line in the Mooring for the constraint to apply to, and its quantitative limit
      - There is a member function (e.g. con_lay_length) that pertains to each constraint. It accepts the design
        vector X, updates the design to ensure the mooring system properties are updated, and evaluates the constraint.
        It returns a negative scalar if the constraint is not met by the quantity specified by the user
    
    Other notable capabilities
    - Shared vs anchored moorings designated by "shared" parameter
    - Anchor spacing can be a design variable. Turbine spacing for a shared line needs to be defined.
    - rBFair is the fairlead coordinates relative to the attached body's reference point for an anchored line in Quadrant 1.
      For example, rBFair can be something like [7.875,0,-21] or [5.57,5.57,-21] or [3.93,6.82,-21]
      For a shared line, the fairlead coordinates are assumed to be the same for both bodies but flipped
    - Design variables can be given one of four designations in Xindices: an integer, 'c' (constant), 
      's' (to be solved for), or 'r' (like a constant, but can be set as a ratio wrt another variable)
    
    Example allVars vector: X = [A or W0, L1, D1, <W1, L2, D2>...]    where <  > section repeats
        For anchor lines, the first entry is anchor spacing. For shared lines, the first entry is midpoint weight.
        Example Mooring: (anchor at spacing A)---- L1,D1-----(W1)------L2,D2------(W2)------L3,D3------(end B)
        Clump weights and buoyancy floats are not specified directly. They are both 'weights' and can have either a postitive or negative value
    
    '''
    
    def __init__(self, depth, lineProps=None, **kwargs):
        '''Creates a LineDesign Mooring object to be used for evaluating or optimizing a mooring line design.
        
        Parameters
        ----------
        depth : float
            Water depth
            
        Keyword Arguments
        -----------------
        solve_for : string
            Keyword indicating which built-in design algorithm to use, if any. Options are:
            'tension' - adjusts a line section length to achieve a target horizontal tension on the line.
            'offset' - adjusts a line section length to achieve a target mean offset considering all lines.
            'stiffness' - adjusts a line section length to achieve a target undisplaced surge stiffness considering all lines.
            'ghost' - adjusts anchor spacing to achieve a target minimum laid length - IN PROGRESS.
            'fancy' - adjusts a line section length to ensure mean offset is less than a target value - IN PROGRESS.
            'none' - makes no adjustment.
            All options except none require that one of the line section lengths be set as solved ('s') rather than fixed/variable.
        DAFs : float or float array, optional
            Dynamic amplification factors to use to scale up quasi-static predicted deviations from mean 
            values to approximate dynamic ones. Provide a scalar or an n+1 array where n is the number
            of line sections and the last entry is the DAF to be used for anchor loads. Default is 1.

        '''
        
        self.display = getFromDict(kwargs, 'display', default=0)
        
        # add the parameters set by the input settings dictionary
        self.name       = getFromDict(kwargs, 'name', dtype=str, default='no name provided')
        lineTypeNames   = getFromDict(kwargs, 'lineTypeNames' , dtype=str, shape=-1, default=[])
        
        
        # set up the mooring system object with the basics from the System class
        rho        = getFromDict(kwargs, 'rho', default=1025.0)
        g          = getFromDict(kwargs, 'g'  , default=9.81)
        self.depth = depth  # used?
        
        # ----- Set properties for Mooring object and its Subsystem -----
        # set model-specific parameters
        self.shared      = getFromDict(kwargs, 'shared', dtype=bool, default=False)
        self.span        = getFromDict(kwargs, 'span', default=0)  # [m] horizontal extent of mooring (formerly "spacing")
        
        # set remaining Mooring-specific parameters
        self.rBFair     = getFromDict(kwargs, 'rBFair', shape=-1, default=[0,0,0])  # [m] end coordinates relative to attached body's ref point
        self.nLines = len(lineTypeNames)                # number of sections in the mooring line
        
        
        
        # ============== set the design variable list ==============
        self.solve_for  = getFromDict(kwargs, 'solve_for', dtype=str, default='offset') # whether to solve for offsets assuming 3 lines, or solve for mean horizontal tension of this line (for use with shared array design tools)
        
        self.allVars    = getFromDict(kwargs, 'allVars' , shape=3*len(lineTypeNames))
        
        # set the design variable type list
        if 'Xindices' in kwargs:
            self.Xindices = list(kwargs['Xindices'])
            if not len(self.Xindices)==len(self.allVars):
                raise Exception("Xindices must be the same length as allVars")
        else:
            raise Exception("Xindices must be provided.")
            
        # find the largest integer to determine the number of desired design variables
        self.nX = 1 + max([ix for ix in self.Xindices if isinstance(ix, int)])
        
        # check for errors in Xindices
        for i in range(self.nX):
            if not i in self.Xindices:
                raise Exception(f"Design variable number {i} is missing from Xindices.")
        valid = list(range(self.nX))+['c','s','r','g']  # entries must be either design variable index or constant/solve/ratio flags
        for xi in self.Xindices:
            if not xi in valid:
                raise Exception(f"The entry '{xi}' in Xindices is not valid. Must be a d.v. index, 'c', 's', or 'r'.")
        
        # find the length solve index 's' and make sure it's valid
        sInds = [i for i,xi in enumerate(self.Xindices) if xi=='s']
        if len(sInds) == 1:
            if (sInds[0]-1)%3 == 0:
                self.iL = int((sInds[0]-1)/3)  # this is the line index whose length will be adjusted in the dsolve inner loop
            else:
                raise Exception("The 's' flag in Xindices must be at a line length (i.e. the 2nd, 5th, 8th...) position.")
        elif len(sInds) == 0:
            if self.solve_for in ['none', 'ghost']:
                self.iL = 0         # arbitrary line index. The index won't matter when solve_for = 'none'
            else:
                raise Exception("A single 's' flag for line length solving must be provided in Xindices")
        else:
            raise Exception("A single 's' flag for line length solving must be provided in Xindices")
        
        # check for 'r' variable option
        self.rInds = [i for i,xi in enumerate(self.Xindices) if xi=='r']
        for i in range(len(self.rInds)):
            if self.allVars[self.rInds[i]] >= 1.0 or self.allVars[self.rInds[i]] <= 0.0:
                raise Exception("The ratio variable needs to be between 1 and 0")
            
        
        # set up the mooring system for the specific configuration type
        '''
        Just makes the connections, sizing happens later.
        
        Example allVars vector: X = [A or W0, L1, D1, <W1, L2, D2>...]    where <  > section repeats
        For anchor lines, the first entry is anchor spacing. For shared lines, the first entry is midpoint weight.
        Example Mooring: (anchor at spacing A)---- L1,D1-----(W1)------L2,D2------(W2)------L3,D3------(end B)
        Clump weights and buoyancy floats are not specified directly. They are both 'weights' and can have either a postitive or negative value
        '''
        
        # first set the weight, length, and diameter lists based on the allVars inputs. Don't worry about design variables yet. Create the units list too.
        
        if self.shared==1:
            if self.span == 0:  raise Exception("For shared arrangements, a span must be provided to the Mooring object.")
            Ws = self.allVars[0::3].tolist()
        else:
            self.span = self.allVars[0]*10 - self.rBFair[0] # in tens of meters
            Ws = self.allVars[3::3].tolist()

        Ls = self.allVars[1::3].tolist()
        Ds = self.allVars[2::3].tolist()

        unitPattern = ['t', 'm', 'mm']
        self.allVarsUnits = [unitPattern[i % 3] for i in range(len(self.allVars))]
        if self.shared==0:
            self.allVarsUnits[0] = 'm'
        # if any of the input lengths are in ratio form, convert them to real value form
        # (this can currently only handle 1 ration variable per Mooring)
        if len(self.rInds) > 0:
            self.nsll_ratio = self.allVars[self.rInds[0]]
            self.allVars[self.rInds[0]] = self.nsll_ratio*self.span 
        Ls = self.allVars[1::3].tolist()    # reset the Ls variable

        # ====================================================================
        
        
        
        # ----- Initialize some objects -----
        if self.shared==1:
            shared=1
            shareCase=2  # assumed symmetric and we model half the shared line.
        elif self.shared==0:
            shared=0
            shareCase=0        
        # make a dummy design dictionary for Mooring to make a Subsystem with???
        dd = dict(subcomponents={})
        
        # Create subcomponents list: alternating Connectors and Sections
        # Pattern: [Connector, Section, Connector, Section, ..., Connector]
        # Total length: 2*nLines + 1 (nLines sections + nLines+1 connectors)
        dd['subcomponents'] = [{} for i in range(2*self.nLines + 1)]
        
        # the sizing function coefficients to use in the design
        self.lineProps = loadLineProps(lineProps)
        
        # Build alternating subcomponents list
        for i in range(self.nLines):
            # Connector at position 2*i (even indices: 0, 2, 4, ...)
            connector_idx = 2*i
            section_idx = 2*i + 1
            
            # Initialize connector properties (will be populated below)
            dd['subcomponents'][connector_idx] = {'m': 0, 'v': 0, 'CdA': 0}
            
            # Assign section properties
            dd['subcomponents'][section_idx]['type'] = getLineProps(Ds[i], 
                material=lineTypeNames[i], name=i, lineProps=self.lineProps)
            dd['subcomponents'][section_idx]['L'] = Ls[i]
        
        # Add final connector at end
        dd['subcomponents'][2*self.nLines] = {'m': 0, 'v': 0, 'CdA': 0}
        
        # Assign props of first connector if shared (midpoint weight)
        if self.shared==1:
            pointDict = self.getClumpMV(Ws[0])
            dd['subcomponents'][0]['m'] = pointDict['m']
            dd['subcomponents'][0]['v'] = pointDict['v']
        
        # Assign props for intermediate connectors
        for i in range(self.nLines-1):
            # Intermediate connectors are at positions 2, 4, 6, ... (2*(i+1))
            connector_idx = 2*(i+1)
            pointDict = self.getClumpMV(Ws[ i + 1*(self.shared==1)]) 
            
            dd['subcomponents'][connector_idx]['m'] = pointDict['m']
            dd['subcomponents'][connector_idx]['v'] = pointDict['v']
            # CdA could be added here if needed
        
        # General mooring dimension info
        dd['span'    ] = self.span
        dd['zAnchor' ] = -self.depth
        dd['rad_fair'] = np.abs(self.rBFair[0])
        dd['z_fair'  ] = self.rBFair[2]
        
        # super().__init__(depth=depth, rho=rho, g=g, lineProps=lineProps) # if we're a subsystem
        
        # Call Mooring init function (parent class)

        
        Mooring.__init__(self, dd=dd, rho=rho, g=g, shared=shared)        
        # The above will also create Mooring self parameters like self.rad_anch
        
        # Save a copy of the original anchoring radius to use with the 
        # solve_for=ghost option to adjust the chain length.
        self.rad_anch0 = float(self.rad_anch) 
        
        self.createSubsystem(case=int(shareCase))  
        if self.shared==1:
            self.ss.rA[2] = self.rBFair[2]
            
            # HARDCODING THIS FOR NOW (MIDPOINT WEIGHT MUST BE UPDATED)
            pointDict = self.getClumpMV(.5*Ws[0]) 
    
            self.dd['subcomponents'][0]['m'] = pointDict['m']
            self.dd['subcomponents'][0]['v'] = pointDict['v']
            
            self.ss.pointList[0].m = pointDict['m']
            self.ss.pointList[0].v = pointDict['v']
            
        self.ss.eqtol = getFromDict(kwargs, 'eqtol', default=0.002)  # position tolerance to use in equilibrium solves [m]

        # load a custom line props scaling dict if provided ??
        #self.ss.lineProps = lineProps 
        
        
        # identify number of line sections and initialize dynamic amplification factors
        self.DAFs = getFromDict(kwargs, 'DAFs', shape=self.nLines+2, default=1.0)   # dynamic amplication factor for each line section, and anchor forces (DAFS[-2] is for vertical load, DAFS[-1] is for horizontal load)
        self.Te0 = np.zeros([self.nLines,2])   # undisplaced tension [N] of each line section end [section #, end A/B]
        self.LayLen_adj = getFromDict(kwargs, 'LayLen_adj', shape=0, default=0.0) # adjustment on laylength... positive means that the dynamic lay length is greater than linedesign laylength 
        self.damage = getFromDict(kwargs, 'damage', shape = -1, default = 0.0) #Lifetime fatigue damage *(MBL/dT/dx)^m in list with same order as fatigue_headings
        self.fatigue_headings = getFromDict(kwargs, 'fatigue_headings', shape = -1, default = [0]) #loading directions for fatigue damage, same order as self.damage
        self.ms_fatigue_index = int(getFromDict(kwargs, 'ms_fatigue_index', shape = 0, default = 1)) #index of line in full moorpy system for fatigue damage evaluation. linelist follows the order in headings
        self.corrosion_mm = getFromDict(kwargs, 'corrosion_mm', default=0)             # [mm] the corrosion of line material over a 25 year lifetime
        
        # ----- Set solver and optimization settings -----
        
        self.x_target    = getFromDict(kwargs, 'x_target', default=0)                       # [m] target mean offset at rated load (e.g. from LinearSystem) - only used in solve_for offset or ghost
        self.x_mean_in   = getFromDict(kwargs, 'x_mean_in', default=0)  
        self.x_mean_out  = getFromDict(kwargs, 'x_mean_out', default=0)
        #self.x_mean_max  = getFromDict(kwargs, 'x_mean_max', default=self.x_mean)       # set the maximum tolerable mean offset to match the initial target mean offset << appears no longer really used
        self.x_ampl      = getFromDict(kwargs, 'x_ampl'  , default=10)                  # [m] expected wave-frequency motion amplitude about mean
        #self.x_extreme  = getFromDict(kwargs, 'xextreme'  , default=self.xmax)          # >>> same as below, but leaving for now for backward compatibility <<<
        #self.x_extr_pos = getFromDict(kwargs, 'x_extr_pos', default=self.xmax)          # [m] expected maximum extreme offset (mean + dynamics)
        #self.x_extr_neg = getFromDict(kwargs, 'x_extr_neg', default=-self.x_extr_pos)   # [m] expected maximum extreme negative offset (negative of xextreme unless provided separately)
        self.fx_target   = getFromDict(kwargs, 'fx_target')                # [N] the expected thrust force or target horizontal line tension
        self.kx_target   = getFromDict(kwargs, 'kx_target', default=0)                  # [N/m] the target horizontal line stiffness
        if self.solve_for == 'ghost':
            self.lay_length_target = getFromDict(kwargs, 'lay_target')  # [m] Target laid length - required when solve_for is ghost
        
        self.headings   = getFromDict(kwargs, 'headings' ,  shape=-1, default=[60, 180, 300])      # [deg] headings of the mooring lines (used only when solve_for is 'offset', 'stiffness', or 'fancy')
        
        # >>> TODO: add something that adjusts headings to give min/max offsets in -/+ x direction <<<
        
        self.noFail = False   # can be set to True for some optimizers to avoid failing on errors
        
        self.iter = -1  # iteration number of a given optimization run (incremented by updateDesign)
        self.log = dict(x=[], f=[], g=[], time=[], xe=[], a=[])  # initialize a log dict with empty values 
        
        # set the anchor type and initialize the horizontal and vertical capacities of the anchor
        self.anchorType = getFromDict(kwargs, 'anchorType', dtype=str, default='drag-embedment')
        self.anchorFx   = 0.0
        self.anchorFz   = 0.0            
        
        
        # ----- optimization stuff -----
        # get design variable bounds and last step size
        self.Xmin       = getFromDict(kwargs, 'Xmin'   , shape=self.nX)             # minimum bounds on each design variable
        self.Xmax       = getFromDict(kwargs, 'Xmax'   , shape=self.nX)             # maximum bounds on each design variable
        self.dX_last    = getFromDict(kwargs, 'dX_last', shape=self.nX, default=[]) # 'last' step size for each design variable
        if len(self.Xmin) != self.nX or len(self.Xmax) != self.nX or len(self.dX_last) != self.nX:
            raise Exception("The size of Xmin/Xmax/dX_last does not match the number of design variables")
        
        
        # initialize the vector of the last design variables, which each iteration will compare against
        self.Xlast = np.zeros(self.nX)
        
        # fill in the X0 value (initial design variable values) based on provided allVars and Xindices (uses first value if a DV has multiple in allVars)
        self.X0 = np.array([self.allVars[self.Xindices.index(i)] for i in range(self.nX)])
        self.X0Units = np.array(self.allVarsUnits)[[self.Xindices.index(i) for i in range(self.nX)]]  # corresponding units
        self.X_denorm = np.ones(self.nX)  # normalization factor for design variables
        self.obj_denorm = 1.0  # normalization factor for objective function


        # ----- Set up the constraint functions and lists -----
        
        if 'constraints' in kwargs:
            self.constraints = kwargs['constraints']
        else:
            self.constraints = []
                
        # a hard-coded dictionary that points to all of the possible constraint functions by name
        self.confundict = {"min_Kx"                : self.con_Kx,             # a minimum for the effective horizontal stiffness of the mooring
                           "max_offset"            : self.con_offset,         # a maximum for the horizontal offset in the extreme loaded position
                           "min_lay_length"        : self.con_lay_length,     # a minimum for the length of Line 1 on the seabed at x=x_extr_pos (replaces anchor_uplift)
                           "rope_contact"          : self.con_rope_contact,   # a margin off the seabed for Point 2 (bottom of Line 2) at x=x_extr_neg
                           "tension_safety_factor" : self.con_strength,       # a minimum ratio of MBL/tension for all lines in the Mooring at x=x_extr_pos
                           "min_sag"               : self.con_min_sag,        # a minimum for the lowest point's depth at x=x_extr_pos
                           "max_sag"               : self.con_max_sag,        # a maximum for the lowest point's depth at x=x_extr_neg
                           "max_total_length"      : self.con_total_length,   # a maximum line length    
                           "min_yaw_stiff"         : self.con_yaw_stiffness,  # a minimum yaw stiffness for the whole system about the extreme negative position
                           "max_damage"            : self.con_damage,         # a maximum fatigue damage for a specified mooring line (scales from provided damage from previous iteration)
                           "min_tension"           : self.con_min_tension,    # a minimum line tension
                           "min_angle"             : self.con_min_angle,      # a minimum inclination angle for the line
                           "max_angle"             : self.con_max_angle       # a maximum inclination angle for the line
                           }   
        
        # a hard-coded dictionary of the units associated with the constraints
        conUnitsDict = {"min_Kx"                : "N/m",
                        "max_offset"            : "m",
                        "min_lay_length"        : "m",
                        "rope_contact"          : "m",
                        "tension_safety_factor" : "-",
                        "min_sag"               : "m",
                        "max_sag"               : "m",
                        "max_total_length"      : "m",
                        "min_yaw_stiff"         : "Nm/deg",
                        "max_damage"            : "-",
                        "min_tension"           : "N",
                        "min_angle"             : "deg",
                        "max_angle"             : "deg"
                        }

        # add units to the constraints dictionary
        for con in self.constraints:
            if con['name'] in conUnitsDict:
                con['unit'] = conUnitsDict[con['name']]
            else:
                raise Exception(f"The constraint name '{con['name']}' is not recognized.")
        # set up list of active constraint functions
        self.conList = []
        self.convals = np.zeros(len(self.constraints))  # array to hold constraint values
        self.con_denorm = np.ones(len(self.constraints))  # array to hold constraint normalization constants
        self.con_denorm_default = np.ones(len(self.constraints))  # default constraint normalization constants
        
        for i, con in enumerate(self.constraints):    # for each list (each constraint) in the constraint dictionary

            # ensure each desired constraint name matches an included constraint function
            if con['name'] in self.confundict:

                # the constraint function for internal use (this would be called in UpdateDesign)
                def internalConFun(cc, ii):   # this is a closure so that Python doesn't update index and threshold
                    def conf_maker(X):
                        def func():
                            # compute the constraint value using the specified function
                            val = self.confundict[cc['name']](X, cc['index'], cc['threshold'])
                            
                            # record the constraint value in the list
                            self.convals[ii] = val / self.con_denorm[ii] # (normalized)
                            self.constraints[ii]['value'] = val  # save to dict (not normalized)
                            
                            return val
                        return func()
                    return conf_maker

                # make the internal function and save it in the constraints dictionary
                con['fun'] = internalConFun(con, i)

                # the externally usable constraint function maker
                def externalConFun(name, ii):   # this is a closure so that Python doesn't update index and threshold
                    def conf_maker(X):
                        def func():
                            # Call the updatedesign function (internally avoids redundancy)
                            try:
                                self.updateDesign(X)
                                
                                # get the constraint value from the internal list
                                val = self.convals[ii]
                            except:
                                val = -1000
                            
                            return val
                        return func()
                    return conf_maker

                # add the conf function to the conList
                self.conList.append(externalConFun(con['name'], i))
                
                # Save the default/recommended normalization constant
                
                if con['name'] in ['max_total_length']:
                    self.con_denorm_default[i] = con['threshold'] # sum([line.L for line in self.ss.lineList])
                
                elif con['name'] in ['min_Kx', 'tension_safety_factor', 'min_yaw_stiff', 'max_damage']:
                    self.con_denorm_default[i] = con['threshold']
                
                elif con['name'] in ['max_offset', 'min_lay_length', 'rope_contact', 'min_sag', 'max_sag', 'max_touchdown_range']:
                    self.con_denorm_default[i] = depth
                
            else:
                raise ValueError("Constraint parameter "+con['name']+" is not a supported constraint type.")

        
        
        # ensure each constraint is applicable for the type of mooring
        if self.shared==1:
            if any([con['name'] in ["min_lay_length", "rope_contact"] for con in self.constraints]):
                raise ValueError("You are using a constraint that will not work for a shared mooring line")
        else:
            if any([con['name'] in ["min_sag", "max_sag"] for con in self.constraints]):
                raise ValueError("You are using a sag constraint that will not work for an anchored mooring line")
        
        if not self.solve_for in ['none', 'ghost']:
            if any([con['name'] == "max_offset" for con in self.constraints]):
                raise ValueError("The offset constraints should only be used when solve_for is none or ghost")
        
        if self.solve_for == 'ghost':
            if any([con['name'] == "min_lay_length" for con in self.constraints]):
                print('Warning: having a min_lay_length cosntraint may conflict with lay_length_target in solve_for ghost.')
                #ind = [con['name'] for con in self.constraints].index('min_lay_length')
                #self.lay_length_target = self.constraints[ind]['threshold']
            if shared:
                raise Exception("solve_for ghost can't be used for shared lines")
            if not self.Xindices[0] == 'c':
                raise Exception("solve_for ghost requires the Xindices[0] to be 'c'.")
        
        
        # =============================================================
        
        # make the mooring system
        # self.makeGenericMooring( Ls, Ds, lineTypeNames, Ws, suspended=int(self.shared))
        
        if self.solve_for in ['none', 'offset'] and len(self.headings) == 0:
            raise Exception('When solve_for is none or offset, line headings must be provided.')
        
        # If needed, make a MoorPy System to use for determining offsets
        self.ms = None


        # These options require forces/stiffnesses of the whole mooring system
        if self.solve_for in ['none', 'offset', 'ghost']:
            
            self.ms = mp.System(depth=self.depth, rho=self.rho, g=self.g)
            #lineProps=lineProps)
            
            # Add a coupled body to represent the platform
            self.ms.addBody(-1, np.zeros(6), DOFs=[0,1])
            
            # Set up Subsystems at the headings
            for i, heading in enumerate(self.headings):

                rotMat = rotationMatrix(0, 0, np.radians(heading))
                
                # create end Points for the line
                self.ms.addPoint(1, np.matmul(rotMat, [self.rad_anch, 0, -self.depth]))
                self.ms.addPoint(1, np.matmul(rotMat, [self.rad_fair, 0, self.z_fair]), body=1)
                
                # Make subsystem and attach it
                ss = mp.Subsystem(mooringSys=self.ms, depth=self.depth, 
                        span=self.span, rBfair=[-self.rad_fair, 0, self.z_fair])

                # set up the Subsystem design, with references to the types in dd
                types = [self.dd['subcomponents'][i]['type'] for i in range(1, len(self.dd['subcomponents']), 2)]
                ss.makeGeneric(lengths=Ls, types=types)
                self.ms.lineList.append(ss)  # add the SubSystem to the System's lineList
                ss.number = i+1

                # attach it to the respective points
                self.ms.pointList[2*i+0].attachLine(i+1, 0)
                self.ms.pointList[2*i+1].attachLine(i+1, 1)
            
            self.ms.initialize()
        
        # initialize the created mooring system
        self.ss.initialize(daf_dict = {'DAFs': self.DAFs})
        self.ss.setOffset(0)
        self.updateDesign(self.X0, normalized=False)  # assuming X0/AllVars is not normalized                                                     
    
    
    def updateDesign(self, X, display=0, display2=0, normalized=True):
        '''updates the design with the current design variables using improved Fx/Kx solver methods
        '''
        start_time = time.time()
        
        # Design vector error checks
        if len(X)==0:                   # if any empty design vector is passed (useful for checking constraints quickly)
            return
        elif not len(X)==self.nX:
            raise ValueError(f"LineDesign.updateDesign passed design vector of length {len(X)} when expecting length {self.nX}")
        elif any(np.isnan(X)):
            raise ValueError("NaN value found in design vector")
        
        # If X is normalized, denormalize (scale) it up to the full values
        if normalized:
            X = X*self.X_denorm
        
        # If any design variable has changed, update the design and the metrics
        if not all(X == self.Xlast):
            
            self.Xlast = np.array(X)   # record the current design variables
        
            self.iter += 1
            
            if self.display > 2:
                print(f"Iteration {self.iter}")
            
            if self.display > 1:
                print("Updated design")
                print(X)
            
            
            # ----- Apply the design variables to update the design -----
            
            # update anchor spacing
            dvi = self.Xindices[0]              # design variable index - will either be an integer or a string
            if dvi in range(self.nX):           # only update if it's tied to a design variable (if it's an integer)
                if self.shared==1:                 # if it's a shared line, this would be the midpoint weight (we divide by two because we're simulating half the line)
                    
                    pointDict = self.getClumpMV(.5*X[dvi]) 
            
                    self.dd['subcomponents'][0]['m'] = pointDict['m']
                    self.dd['subcomponents'][0]['v'] = pointDict['v']
                    
                    self.ss.pointList[0].m = pointDict['m']
                    self.ss.pointList[0].v = pointDict['v']
                    
                    # arbitrary depth of self.depth/2. Will be equilibrated soon
                    #self.ss.pointList[0].setPosition([ -0.5*self.span, 0, -0.5*self.depth]) 
                    
                else:
                    # if it's an anchor line, this would be the anchor spacing
                    
                    dv = X[dvi]*10
                    self.ss.span = dv - np.abs(self.rBFair[0])                          # update the span of the ld.ss subsystem
                    self.dd['span'] = dv - np.abs(self.rBFair[0])                       # can update the ld.subsystem's design dictionary too
                    
                    self.ss.pointList[0].setPosition([-self.ss.span, 0, -self.depth])

                    self.setAnchoringRadius(dv)
            
            
            # Update section lengths or diameters 
            for i in range(self.nLines):
            
                # length
                dvi = self.Xindices[3*i+1]  # design variable index
                if dvi in range(self.nX):  # only update if it's tied to a design variable
                
                    # Modify section 1 length (if using ghost option)
                    if i==0 and self.solve_for=='ghost':
                        L_new = X[dvi] + (self.rad_anch - self.rad_anch0)
                    else:
                        L_new = X[dvi]
                
                    self.setSectionLength(L_new, i)
                    if self.ms:
                        for ss in self.ms.lineList:
                            ss.lineList[i].setL(L_new)
                    
                #elif dvi=='r':                  # if the line length is a ratio variable, update it to stay the same proportion of the updated anchor spacing
                #    self.ss.lineList[i].setL(self.nsll_ratio*self.span)
            
                # diameter
                dvi = self.Xindices[3*i+2]          # design variable index
                if dvi in range(self.nX):           # only update if it's tied to a design variable
                    lineType = getLineProps(X[dvi], 
                                  material=self.dd['subcomponents'][2*i+1]['type']['material'], 
                                  name=i, lineProps=self.lineProps)
                    # use the update method to preserve refs to the original dict - this 'points'/connects to the subsystem object too!
                    self.dd['subcomponents'][2*i+1]['type'].update(lineType)
            
            # apply corrosion to the mooring's MBL dictionary (which gets references in the getTenSF constraint in subsystem)
            self.addCorrosion(corrosion_mm=self.corrosion_mm)
            
            # update the intermediate points if they have any weight or buoyancy
            for i in range(self.nLines-1):
                dvi = self.Xindices[3*i+3]      # design variable index
                if dvi in range(self.nX):       # only update if it's tied to a design variable
             
                    pointDict = self.getClumpMV(X[dvi])

                    self.dd['subcomponents'][2*(i+1)]['m'] = pointDict['m']
                    self.dd['subcomponents'][2*(i+1)]['v'] = pointDict['v']
                    
                    self.ss.pointList[i+1].m = pointDict['m'] # update clump buoyancy
                    self.ss.pointList[i+1].v = pointDict['v'] # update clump mass
                    
                    if self.ms:  # also update things in the ms if there is one
                        for ss in self.ms.lineList:
                            ss.pointList[i+1].m = pointDict['m']
                            ss.pointList[i+1].v = pointDict['v']
            
            
            # ----- Screen design to make sure it's physically feasible -----

            # >>> TODO: check for negative line lengths that somehow get set <<<


            Lmax = 0.95*(self.ss.span + self.depth+self.rBFair[2])
            
            if sum([self.ss.lineList[i].L for i in range(self.nLines)]) > Lmax:     # check to make sure the total length of line is less than the maximum L shape (helpful for GA optimizations)
                
                if self.solve_for=='none':
                    self.x_mean_in = -1e3
                    self.x_mean_out = 1e3
                self.x_mean_eval = 1e3      # arbitrary high number to set the offset (and constraints)

                for i,con in enumerate(self.constraints):
                    val = -1e3
                    self.convals[i] = val / self.con_denorm[i] # (normalized)
                    self.constraints[i]['value'] = val  # save to dict (not normalized)

            
            else:            
            
                # ----- Length adjustment (seems sketchy) -----
                # print(self.ms.bodyList[0].r6[2])
                # set x0 as a 1D list of the line length to be solved for
                x0 = [self.ss.lineList[self.iL].L]
                
                # maximum length of the segment being sized to avoid fully slack
                Lmax = 0.99*(self.ss.span + self.depth+self.rBFair[2]) - sum([self.ss.lineList[i].L for i in range(self.nLines) if i != self.iL])
                
                # >>> may need a different Lmax for shared lines <<<
                
                if x0[0] >= Lmax:
                    x0[0] = 0.8*Lmax 
                
                
                # ----- Solver process -----
                
                # call dsolve2 to tune line length - eval function depends on solve_for
                # note: use a high stepfac so that dsolve2 detects a nonzero slope even when the slope is quite shallow
                if self.solve_for == "tension":
                    x, y, info = dsolve2(self.func_TH_L, x0, tol=[0.4*self.ss.eqtol], args=dict(direction='horizontal'),
                                        Xmin=[10], Xmax=[Lmax], dX_last=[10], maxIter=40,
                                        stepfac=100, display=self.display-1)
                elif self.solve_for == "offset":
                    
                    args = dict(xOffset=self.x_target, display=self.display-1)
                    
                    x, y, info = dsolve2(self.func_fx_L, x0, args=args, 
                                        tol=[0.4*self.ss.eqtol], Xmin=[10], Xmax=[Lmax], 
                                        dX_last=[10], stepfac=100, maxIter=40, 
                                        display=self.display-1)

                elif self.solve_for == "none":
                    pass
                    # >>> can remove this from if else block once solve_for error check is done in init func <<<
                
                elif self.solve_for == 'stiffness':
                    
                    x, y, info = dsolve2(self.func_kx_L, x0, args=dict(display=display2),
                                        tol=[0.4*self.ss.eqtol], Xmin=[10], Xmax=[Lmax], 
                                        dX_last=[10], stepfac=100, display=self.display-1)
                    
                    # this solves for the line length to meet a stiffness equality constraint
                    # which means that we can still have an offset constraint since the line 
                    # length isn't being solved for to meet a certain offset

                
                elif self.solve_for == 'fancy':    # a new option to allow lower mean offsets (need to rename!)
                    # Outer loop determines offset that gives target tension SF, inner loop adjusts line length to achieve said offset
                    def tuneLineLengthsForOffset(xCurrent, args):   # this function does the standard "offset"-mode solve, but now it can be done in the loop of another solve
                        
                        args = dict(xOffset=xCurrent[0], fx_target=self.fx_target)
                        
                        # tune line length until thrust force is balanced at this mean offset
                        x, y, info = dsolve2(self.func_fx_L, x0, args=args,
                                            tol=[0.4*self.ss.eqtol], Xmin=[10], Xmax=[Lmax],
                                            dX_last=[10], stepfac=100, display=0)
                    
                        stopFlag = False if info['success'] else True  # if the line length solve was unsuccessful, set the flat to stop the mean offset solve
                    
                        # check strength constraint at this offset + some dynamic additional offset
                        # (doing this manually here for now, and avoiding the strength constaint at higher levels >>> do not use tension_safety_factor! <<<)
                        '''This ensures the MBL of the line is always greater than the maximum tension the line feels times a safety factor'''
                        self.ss.lineList[self.iL].setL(x[0])                            # make sure the design is up to date (in terms of tuned line length)
                        self.ss.setOffset(xCurrent[0] + 10)  # offset the body the desired amount (current mean offset + wave offset)
                        cMin = self.ss.getMinSF(display=display) - 2.0          # compute the constraint value     
                        
                        print(f" xmax={xCurrent[0]:8.2f}  L={x[0]:8.3f}  dFx={y[0]:8.0f}  minSF={self.getMinSF():7.3f}")
                        #breakpoint()
                        
                        return np.array([cMin]), dict(status=1), stopFlag     # return the constraint value - we'll actually solve for this to be zero - finding the offset that just barely satisifes the SFs
                    
                    # solve for mean offset that will satisfy tension safety factor constraint (after dynamic wave offset is added)
                    x, y, info = dsolve2(tuneLineLengthsForOffset, [5], tol=[4*self.ss.eqtol], Xmin=[1], Xmax=[4*self.x_target], dX_last=[5], stepfac=10, display=1)
                
                
                elif self.solve_for=='ghost': 
                    '''Use very large anchor spacing and compute an imaginary
                    anchor spacing and line length based on the desired lay
                    length.'''
                    
                    # Compute the offset with the adjusted design variables
                    self.x_mean_out = self.getOffset(self.fx_target)
                    self.ms.bodyList[0].setPosition([0,0,0,0,0,0])  # ensure body is re-centered

                    # self.span and self.ss.span seems redundant. Does LD/Mooring need it?? <<<
                    
                    # figure out tension in least laid length scenario...
                    self.ss.setOffset(self.x_mean_out)  # apply max static offsets
                    self.ss.setDynamicOffset(self.x_mean_out + self.x_ampl)  # move to dynamic offset
                    max_anchor_tension = self.ss.TeD[0,0]  # save tension at anchor

                    # Set anchoring radius a bit larger than needed, and evaluate once (ss only)
                    length_to_add = 0.2 * self.rad_anch
                    new_length = self.dd['subcomponents'][1]['L'] + length_to_add/(1 + max_anchor_tension/self.ss.lineList[0].EA)
                    self.rad_anch = float(self.rad_anch + length_to_add)
                    self.ss.span  = self.rad_anch - self.rBFair[0]
                    self.ss.setEndPosition([-self.rad_anch, 0, -self.depth], endB=False)
                    Mooring.setSectionLength(self, new_length, 0)  # ss only, skip ms

                    # Figure out lay length
                    self.ss.setOffset(self.x_mean_out)  # apply max static offsets
                    self.ss.setDynamicOffset(self.x_mean_out + self.x_ampl)  # move to dynamic offset
                    max_anchor_tension = self.ss.TeD[0,0]  # save tension at anchor
                    min_lay_length = self.ss.getLayLength()  # record minimum lay length

                    # Adjust anchor positions to hit target
                    unused_length = min_lay_length - self.lay_length_target
                    new_length = self.dd['subcomponents'][1]['L'] - unused_length
                    new_spacing = self.rad_anch - unused_length*(1 + max_anchor_tension/self.ss.lineList[0].EA)
                    self.setAnchoringRadius(new_spacing)
                    self.setSectionLength(new_length, 0)

                    # Update the Subsystem solutions after the adjustments
                    self.ss.staticSolve()
                    for ss in self.ms.lineList:
                        ss.staticSolve()
                    
                    #print(f"{self.iter}  {self.ss.offset:6.2f}m offset,  {self.rad_anch:6.2f} rad_anch, {self.ss.lineList[0].L:6.2f} L")

                else:
                    raise Exception("solve_for must one of 'offset', 'tension', 'none', 'stiffness, 'fancy', or 'ghost'")
                
                
                if not self.solve_for in ['none', 'ghost']:
                    if info["success"] == False:
                        print("Warning: dsolve2 line length tuning solve did not converge.")
                        #breakpoint()    # <<<< handle non convergence <<<
                    else:
                        #>>>>> deal with nonzero y - penalize it somehow - for optimizer <<<<<
                        
                        # ensure system uses latest tuned line length
                        #self.ss.lineList[self.iL].setL(x[0])
                        self.setSectionLength(x[0], self.iL)
                
                
                # ----- Compute (or set) high and low mean offsets -----
                # (solve for the offsets at which the horizontal mooring reactions balance with fx_target)
                if self.solve_for in ['none', 'ghost']:
                    self.x_mean_out = self.getOffset(self.fx_target) 
                    self.x_mean_in = -self.getOffset(-self.fx_target)
                    if self.display > 1: print(f" Found offsets  x_mean_out: {self.x_mean_out:.2f},  x_mean_in: {self.x_mean_in:.2f}")
                    self.ms.bodyList[0].setPosition([0,0,0,0,0,0])  # ensure body is re-centered
                
                # x_mean_in is the offset when the input headings are flipped, representing the opposite loading direction.
                # This will only be worst-case/best-case offsets when one of the input headings is either completely upwind or completely downwind.
                
                
                # ----- Evaluate system state and constraint values at offsets -----
                
                # Evaluate any constraints in the list, at the appropriate displacements.
                # The following function calls will fill in the self.convals array.
                
                
                # ZERO OFFSET:
                self.ss.setOffset(0)
                
                # get undisplaced tensions of each line section and anchors
                for i, line in enumerate(self.ss.lineList):
                    self.Te0[i,0] = np.linalg.norm(line.fA)
                    self.Te0[i,1] = np.linalg.norm(line.fB)
                    

                
                # Call any constraints that evaluate at the undisplaced position
                #self.calcCurvature()
                for con in self.constraints:
                    if con['offset'] == 'zero':
                        con['fun'](X)
                
                
                # MAX OFFSET: 
                self.ss.setOffset(self.x_mean_out)  # apply static offsets
                self.ss.setDynamicOffset(self.x_mean_out + self.x_ampl)  # move to dynamic offset
                # save maximum anchor tensions for use in cost calculations (includes DAF)
                self.anchorFx = self.ss.anchorFxD
                self.anchorFz = self.ss.anchorFzD
                
                self.min_lay_length = self.ss.getLayLength()  # record minimum lay length
                #print(f"{self.iter}  {self.ss.offset:6.2f}m offset,  {self.rad_anch:6.2f} rad_anch, {self.ss.lineList[0].L:6.2f} L")
                #print(f"Min lay length is {self.min_lay_length}")
                self.x_mean_eval = float(self.x_mean_out)    # the x_mean value to evaluate if there's an offset constraint
                
                # Call any constraints needing a positive displacement
                for con in self.constraints:
                    if con['offset'] == 'max':
                        con['fun'](X)

                
                # MIN OFFSET: 
                self.ss.setOffset(-self.x_mean_in)  # apply static offset
                self.ss.setDynamicOffset(-self.x_mean_in + -self.x_ampl)  # peak offset

                self.max_lay_length = self.ss.getLayLength()  # record maximum lay length

                self.x_mean_eval = float(self.x_mean_in)    # the x_mean value to evaluate if there's an offset constraint
                
                # Call any constraints needing a negative displacement
                for con in self.constraints:
                    if con['offset'] == 'min':
                        con['fun'](X)
                
                
                # OTHER: 
                self.ss.setOffset(0)  # restore to zero offset and static EA            
                # or at least set back to static states
                
                # Call any constraints that depend on results across offsets
                for con in self.constraints:
                    if con['offset'] == 'other' or con['offset'] == 'zero':
                        con['fun'](X)
                
                ############################################################  
            
            # ----- Evaluate objective function -----
            
            # Calculate the cost from all components in the Mooring
            self.lineCost = 0.0
            for line in self.ss.lineList:
                self.lineCost += line.L*line.type['cost']
            
            # Adjust cost for active length in case of ghost option
            if self.solve_for == 'ghost':
            
                # the length beyond the minimum lay length is not used
                unused_length = self.min_lay_length - self.lay_length_target
                
                # just adjust costs from first section
                self.lineCost -= unused_length * self.ss.lineList[0].type['cost']
                
                # also adjust and record anchor position ???  (would be nice to show on plots)

            
            # calculate anchor cost (using anchor forces calculated when the mooring's constraints were analyzed)
            if self.shared==1:
                self.anchorCost = 0.0
                self.anchorMatCost = 0.0
                self.anchorInstCost = 0.0
                self.anchorDecomCost = 0.0
            else:
                self.anchorMatCost, self.anchorInstCost, self.anchorDecomCost = getAnchorProps(self.anchorFx, self.anchorFz, type=self.anchorType)
                self.anchorCost = self.anchorMatCost + self.anchorInstCost + self.anchorDecomCost
            
            # calculate weight/float cost
            self.wCost = 0.0
            self.WF = 1.0   # weight factor: a multiplier for the weight cost per unit mass (kg)
            for point in self.ss.pointList:
                if point.number > 1 and point.number < self.nLines+1:
                    self.wCost += abs(point.m + point.v*self.rho)*self.WF
            
            # if it's shared, we need to double the line costs since it's mirrored
            if self.shared==1:
                self.lineCost = self.lineCost*2
                self.wCost = self.wCost*2     
                
            # total cost for all 3 moorings
            self.cost = self.lineCost + self.anchorCost + self.wCost
            
            if self.display > 1:
                print(' Cost is ',self.cost)


            # >>> dynamic_L = self.ss.lineList[0].L - self.min_lay_length   #(for line [0] only...)
            
            self.obj_val = self.cost / self.obj_denorm  # normalize objective function value
            
            
            # ----- write to log -----
            
            # log the iteration number, design variables, objective, and constraints
            self.log['x'].append(list(X))
            self.log['f'].append(list([self.obj_val*self.obj_denorm]))
            self.log['g'].append(list(self.convals*self.con_denorm))
            self.log['time'].append(time.time() - start_time)
            self.log['xe'].append(self.ss.lineList[self.iL].L)
            self.log['a'].append((self.ss.span + self.rBFair[0])/10)
        
            # TODO: log relevant tuned values (line length, lay length, etc.) for each solve_for option <<<
        
        # provide some output
        if self.display > 5:
            f = self.objectiveFun(X, display=1)
            
            Fx = self.fB_L[0] # get horizontal force from mooring on body

            print(f"Fx: {Fx/1e3:7.1f} vs target of {self.fx_target/1e3:7.1f}")
            
            print("Line lengths are ")
            for line in self.ss.lineList:
                print(line.L)
                
            print("Line input diameters are ")
            for lineType in self.lineTypes.values():
                print(lineType['input_d'])
                
            print(f"Cost is {f}")    
                
            self.evaluateConstraints(X, normalized=False, display=1)
            
            self.plotProfile()
            plt.show()
    
    
    def objectiveFun(self, X, display=0, normalized=True):
        '''objective of the optimization. Set to minimize cost'''
        
        self.updateDesign(X, display=display, normalized=normalized)

        if display > 1:
            print(f"Cost is {self.cost:.1f} and objective value is {self.obj_val:.3f}.")

        return float(self.obj_val)  # return a copy
    
    
    def evaluateConstraints(self, X, display=0, normalized=True):
        '''Update the design (if necessary) and display the constraint
        values.'''

        self.updateDesign(X, display=display, normalized=normalized)
        
        if display > 1:
            for i, con in enumerate(self.constraints):
                print(f" Constraint {i:2d} value of {con['value']:8.2f} "
                      +f"for {con['name']}: {con['threshold']} of {con['index']} at {con['offset']} displacement.")
        
        return np.array(self.convals)  # return a copy
        
        
    def setNormalization(self):
        '''Set normalization factors for optimization 
        (based on initial design state).'''
        
        # design variables
        self.X_denorm = np.array(self.Xlast)
        # objective
        self.obj_denorm = self.cost
        # constraints
        self.con_denorm = self.con_denorm_default
    
    
    def clearNormalization(self):
        '''Clear any normalization constants to unity so no scaling is done.'''
        self.X_denorm = np.ones(self.nX)
        self.obj_denorm = 1.0
        self.con_denorm = np.ones(len(self.constraints))
    
    
    def optimize(self, gtol=0.03, maxIter=40, nRetry=0, plot=False, display=0, stepfac=4, method='dopt'):
        '''Optimize the design variables according to objectve, constraints, bounds, etc.
        '''
        
        # set the display value to use over the entire process
        self.display = display
        
        # reset iteration counter
        self.iter = -1
        
        # clear optimization progress tracking lists
        self.log['x'] = []
        self.log['f'] = []
        self.log['g'] = []
        self.log['time'] = []
        self.log['xe'] = []
        self.log['a'] = []
        
        def eval_func(X, args):
            '''Mooring object evaluation function condusive with MoorSolve.dopt2'''
            
            f = self.objectiveFun(X, display=display)
            g = self.evaluateConstraints(X, display=display)
            oths = dict(status=1)
            Fx = self.ss.fB_L[0]
            
            return f, g, [self.ss.lineList[self.iL].L], Fx, oths, False
        
        
        # set noFail for GAs in case they come up with crazy designs (avoids exceptions)
        if 'GA' in method:
            self.noFail = True
        else:
            self.noFail = False
        
        # Set starting point to normalized value
        X0 = self.X0 / self.X_denorm
        dX_last = self.dX_last / self.X_denorm
        Xmax = self.Xmax / self.X_denorm
        Xmin = self.Xmin / self.X_denorm
        
        
        # call optimizer to perform optimization
        # --------- dopt method: Newton Iteration -----------
        if method=='dopt':
        
            if display > 0:  print("\n --- Beginning LineDesign2 optimize iterations using DOPT2 ---")

            X, min_cost, infodict = dopt2(eval_func, X0, tol=4*self.ss.eqtol, a_max=1.4, maxIter=maxIter, stepfac=stepfac,
                                                Xmin=Xmin, Xmax=Xmax, dX_last=dX_last, display=4) #self.display)
            
        
            # Retry procedure if things don't work
            for i in range(nRetry):
                print(f" Mooring optimization attempt {i} was UNSUCCESFUL")
                print(f" Message from dopt on attempt {i}: {infodict['message']}")
            
                self.updateDesign(X)                                # update the mooring using the optimized design variables
                G = self.evaluateConstraints(X)              # evaluate the constraints of the mooring
                
                # check how far the constraints are off
                c_rel = G / np.array([con[2] for con in self.constraints], dtype=float) # get relative value of constraints (denominator converts dict values to a np.array)
                i_kx = [i for i,con in enumerate(self.constraints) if con=='min_Kx'][0] # index of Kx constraint

                if display > 1: print(f' stiffness is {c_rel[i_kx]*100+100.:5.1f}% of target')
                c_rel[i_kx] = 0.0            
                # zero the kx constraint since it's okay to break it (that's why we iterate with LinearSystem)

      
            
                if np.min(c_rel) < -gtol:
                
                    # try to catch some known problem cases
                    if stepfac==10:
                        print(' retrying optimization with step size (stepfac) boosted from 10 to 100')
                        stepfac = 100
                            
                    else:
                        self.updateDesign(X)                             # make sure it's left at the optimized state
                        break  # out of ideas, so that's the best we can do with this design problem
                        
                    # rerun the optimizer with modified settings
                    X, min_cost, infodict = dopt2(eval_func, X0, tol=0.001, a_max=1.4, maxIter=maxIter,
                                                        Xmin=Xmin, Xmax=Xmax, dX_last=dX_last)
    
    
                else:                                                # if successful
                    self.updateDesign(X)                             # make sure it's left at the optimized state
                    break                                            # exit the retry loop
            
                if i==nRetry-1:    # re-check the design if we did all retries, since otherwise it won't be done by the loop above
                    self.updateDesign(X)                                  # update the mooring using the optimized design variables
                    G = self.evaluateConstraints(X, display=0)     # evaluate the constraints of the mooring
                    
                    # check how far the constraints are off
                    #c_rel = G / np.fromiter(self.constraints.values(), dtype=float)         # get relative value of constraints (denominator converts dict values to a np.array)
                    c_rel = G / np.array([con[2] for con in self.constraints], dtype=float) # get relative value of constraints (denominator converts dict values to a np.array)
                    i_kx = [i for i,con in enumerate(self.constraints) if con=='min_Kx'][0] # index of Kx constraint
                    if display > 1: print(f' stiffness is {c_rel[i_kx]*100+100.:5.1f}% of target')
                    c_rel[i_kx] = 0.0
        
        
        # --------- COBYLA method -----------
        elif method in ['COBYLA', 'SLSQP']:
            
            from scipy.optimize import minimize    
            
            if self.display > 0:  print("\n --- Beginning LineDesign2 optimize iterations using COBYLA ---")

            condict = [dict(type="ineq", fun=con) for con in self.conList]
            cons_tuple = tuple(condict)
            
            if method=='COBYLA':
                result = minimize(self.objectiveFun, X0, constraints=cons_tuple,
                                  method="COBYLA", options={'maxiter':maxIter,
                                  'disp':True, 'rhobeg':0.1, 'catol':0.001}) # 'rhobeg':10.0
            
            elif method=='SLSQP':
                result = minimize(self.objectiveFun, X0, constraints=cons_tuple,
                                  method='SLSQP', bounds = list(zip(Xmin, Xmax)),
                                  options={'maxiter':maxIter, 'eps':0.02,
                                  'ftol':1e-6, 'disp': True, 'iprint': 99})
            
            X = result.x
        
        
        # --------- Bayesian method -----------
        elif method == 'bayesian':
            
            from bayes_opt import BayesianOptimization
            from scipy.optimize import NonlinearConstraint 
            
            if self.display > 0:  print("\n --- Beginning LineDesign2 optimize iterations using Bayesian Optimization ---")
            
            # --- make list of decision variable names ---
            # design parameter names [A or W0, L1, D1, <W1, L2, D2>...]
            '''
            param_names = []
            for i in range( (2+self.nX)//3):
                param_names = param_names + [f'W{i}', f'L{i+1}', f'D{i+1}']
            if not self.shared:
                param_names[0] = 'A'
            '''
            dvnames = [str(i) for i in range(self.nX)]
            
            # --- set up constraints ---
            def constraint_function(**kwargs):
                
                # Reconstruct decision variable vector
                X = np.zeros(self.nX)  
                for i in range(self.nX):
                    X[i] = kwargs[dvnames[i]]
                
                # Call and evaluate each constraint?
                G = self.evaluateConstraints(X, display=0)
                
                return G
            
            # Make constraint objects (valid values are from zero to infinity)
            zeros = np.zeros(len(self.conList))
            constraint = NonlinearConstraint(constraint_function, zeros, zeros+np.inf)
            
            # Make objective function to maximize
            def negativeObjectiveFun(**kwargs):
                
                # Reconstruct decision variable vector
                X = np.zeros(self.nX)  
                for i in range(self.nX):
                    X[i] = kwargs[dvnames[i]]
                    
                # Negative version of objective function
                return -1*self.objectiveFun(X, display=0)
            
            # Bounded region of parameter space
            pbounds = {}
            for i in range(self.nX):
                pbounds[dvnames[i]] = (Xmin[i], Xmax[i])
            
            # Set up optimizer
            optimizer = BayesianOptimization(
                f=negativeObjectiveFun, constraint=constraint,
                pbounds=pbounds, verbose=2, random_state=1)
                
            # Find some valid starting points using a random search
            pts = 0  # how many valid starting points found so far
            for i in range(1000):
                x = optimizer.space.random_sample()
                print(x)
                if optimizer.constraint.allowed(optimizer.space.probe(x)[1]):
                    print('Registering start point')
                    print(x)
                    optimizer.space.probe(x)
                    pts += 1
                    
                    if pts > 3:  # <<< How many valid starting points to ensure
                        break
                        
            # Do the optimization
            optimizer.maximize(
                init_points=4,   # <<< Total number of start points before iterating
                n_iter=maxIter)        # <<< Number of points to iterate through
        
            print(optimizer.max)
            
            X = np.array(list(optimizer.max['params'].values()))
        
        
        # ----- CMNGA -----
        elif method=='CMNGA':
            from cmnga import cmnga  # type: ignore
            
            bounds = np.array([[Xmin[i], Xmax[i]] for i in range(len(Xmin))])
            
            X, min_cost, infoDict = cmnga(self.objectiveFun, bounds, self.conList, 
                dc=0.03, nIndivs=14, nRetry=500, maxGens=20, maxNindivs=600 )
        
        
        # --------- Genetic Algorithm ----------
        elif method=='GA':
            
            # import the GA from scipy to save from importing if other optimizers are selected
            from scipy.optimize import differential_evolution, NonlinearConstraint

            # initialize storage variables for the GA, including an iterator variable to track through LineDesign
            n = 100000
            self.XsGA = np.zeros([self.nX, n])
            self.CsGA = np.zeros([len(self.constraints), n])
            self.FsGA = np.zeros([3,n])

            # initialize some GA parameters
            self.popsize = 2
            self.maxiter = 40

            self.popsize = 15
            self.maxiter = 1000

            # bounds
            bounds = [(Xmin[i], Xmax[i]) for i in range(len(Xmin))]
            # constraints
            constraints = tuple([ NonlinearConstraint(self.conList[i], 0, np.inf) for i in range(len(self.conList)) ])
            
            # run the GA
            result = differential_evolution(self.objectiveFun, bounds, maxiter=self.maxiter,
                                            constraints=constraints, popsize=self.popsize, tol=0.1, disp=True, polish=False)
            
            # this doesn't require the initial design variable vector, it searches the whole design space initially

            # set the number of individuals in the population (for some reason, it means NP = popsize*N, where N is # of DV's)
            if self.popsize==1:
                self.NP = self.popsize*self.nX + 1
            else:
                self.NP = self.popsize*self.nX

            # organize the stored variables better (trim the excess zeros)
            XsGA = np.zeros([len(self.XsGA), len(np.trim_zeros(self.XsGA[0,:]))])
            for i in range(len(self.XsGA)):
                XsGA[i,:] = np.trim_zeros(self.XsGA[i,:])
            self.XsGA = np.array(XsGA)
            maxCsGA = 0
            for i in range(len(self.CsGA)):
                if len(np.trim_zeros(self.CsGA[i,:])) > maxCsGA:
                    maxCsGA = len(np.trim_zeros(self.CsGA[i,:]))
            CsGA = np.zeros([len(self.CsGA), maxCsGA])
            for i in range(len(self.CsGA)):
                CsGA[i,:len(np.trim_zeros(self.CsGA[i,:]))] = np.trim_zeros(self.CsGA[i,:])
            self.CsGA = np.array(CsGA)
            FsGA = np.zeros([len(self.FsGA),len(self.CsGA[0])])
            for i in range(len(FsGA)):
                FsGA[i,:] = np.array(self.FsGA[i,0:len(self.CsGA[0])])
            self.FsGA = FsGA
            
            X = result.x

    
        # --------- Particle Swarm Algorithm ----------
        elif method=='PSO':

            from pyswarm import pso

            xopt, fopt = pso(self.objectiveFun, Xmin, Xmax, f_ieqcons=self.getCons4PSO, 
                      maxiter=50, swarmsize=1000, debug=True)

            # TODO: Either implement a change in the pyswarm.pso function to rerun if a feasible design isn't found in the first generation OR
            # implement a try/except statement in updateDesign for when solveEquilibrium errors occur (which will happen if a feasible design isn't found)

            X = xopt

        
        else:
            raise Exception('Specified optimization method not recognized.')
        
        # make sure it's left at the optimized state
        self.updateDesign(X)

        # save a couple extra metrics
        #self.infodict['weight'] = -self.fB_L[2]
        
        # check whether optimization passed or failed based on constraint values
        self.evaluateConstraints(X, display = 5)
        if min(self.convals) > -0.01:
            self.success = True
        else:
            self.success = False
      
        if plot:
            self.plotOptimization()

        return X, self.cost #, infodict
    

    def getCons4PSO(self, X):
        conList = []
        for con in self.constraints:
            conList.append(con['value'])
        return conList
    

    def plotOptimization(self, layout="tall", return_fig=False):
        '''Plot the optimization trajectory, including design variables, constraints and cost.
        
        Parameters
        ----------
        layout : str
            "tall" (default) or "grid" layout for subplots. The grid will place all d.v.s in the first column, all constraints in the second column, and cost in the third column.
        '''
        if len(self.log['x']) == 0:
            print("No optimization trajectory saved (log is empty). Nothing to plot.")
            return
        
        Xs = np.array(self.log['x'])
        Fs = np.array(self.log['f'])
        Gs = np.array(self.log['g'])
        
        n_dv  = len(self.X0)
        n_con = len(self.constraints)
        
        if layout=="tall":
            n_rows = n_dv + n_con + 1
            fig, axes = plt.subplots(n_rows, 1, sharex=True, figsize=[6, 1.5*n_rows])
            axes = axes.reshape(-1, 1)
        elif layout=="grid":
            n_rows = max(n_dv, n_con, 1)
            fig, axes = plt.subplots(n_rows, 3, sharex=True, figsize=[12, 1.5*n_rows])
            if n_rows == 1:
                axes = axes[np.newaxis, :]                

        # --- Column 1: Design variables ---
        for i in range(n_dv):
            ax = axes[i, 0]
            ax.plot(Xs[:, i], color='blue')
            ax.set_ylabel(f"d.v.{i+1} ({self.X0Units[i]})", rotation=90, fontsize=10, va='center')

        # --- Column 2 / stacked: Constraints ---
        for i, con in enumerate(self.constraints):
            idx = i if layout == "grid" else n_dv + i
            ax = axes[idx, 1 if layout == "grid" else 0]
            ax.axhline(0, color=[0.5,0.5,0.5])
            tol = 0.005 * (max(Gs[:, i])-min(Gs[:, i]))
            color = 'green' if Gs[-1, i] >= -tol else 'red'
            ax.plot(Gs[:, i], color=color)
            if con['name']=='tension_safety_factor':
                con_name = 'SF'
            else:
                con_name = con['name']
            ax.set_ylabel(f"{con_name} ({con['unit']})",
                        rotation=90, 
                        va='center', fontsize=10)
            # Show threshold value inside plot
            ax.text(0.98, 0.90, f"{con['threshold']}",
                    transform=ax.transAxes,
                    va='top', ha='right', fontsize=8, color='black')
            
        # --- Column 3 / stacked: Cost ---
        if layout == "grid":
            ax_cost = axes[0, 2]
        else:
            ax_cost = axes[-1, 0]
        ax_cost.plot(Fs/1e6, color='black')
        ax_cost.set_ylabel('cost (M$)', rotation=90, va='center', fontsize=10)
    
        # remove unused axes if layout='grid'
        if layout=="grid":
            for i in range(n_dv, n_rows):
                axes[i, 0].axis('off')
            
            for i in range(n_con, n_rows):
                axes[i, 1].axis('off')  
            
            for i in range(2, n_rows):
                axes[i, 2].axis('off')

        # --- X labels only on bottom subplots ---
        for i in range(n_rows):
            for j in range(axes.shape[1]):
                if axes[i, j].has_data():
                    if layout == "tall":
                        if i == n_rows-1:  # only bottom row
                            axes[i, j].set_xlabel("function evaluations")
                        else:
                            axes[i, j].set_xlabel("")
                    elif layout == "grid":
                        if (i == n_dv-1 and j == 0) or (i == n_con-1 and j == 1) or (i == 0 and j == 2):
                            axes[i, j].set_xlabel("function evaluations")
                        else:
                            axes[i, j].set_xlabel("")


        plt.tight_layout()

        if return_fig:
            return fig, axes

    def plotGA(self):
        '''A function dedicated to plotting relevant GA outputs'''

        # determine how many "generations" the GA went through
        gens = [0]; m=3
        gens.append(self.NP*m)
        feasible=False
        while len(gens) < self.maxiter+1:
            while not feasible and len(gens) < self.maxiter+1:
                m=2
                nextgen = gens[-1] + (self.NP*m)
                gens.append(nextgen)
                for f in self.FsGA[0,gens[-2]:nextgen]:
                    if f>0:
                        feasible=True
                        m=1
            if len(gens) < self.maxiter+1:
                nextgen = gens[-1] + (self.NP*m)
                gens.append(nextgen)
        if len(gens) != self.maxiter+1: raise ValueError('Something is not right')


        #Ls = self.allVars[1::3].tolist()
        #Ds = self.allVars[2::3].tolist()
        #Ws = self.allVars[3::3].tolist()

        # set the x-axis vector of each individual that was evaluated
        iters = np.arange(1, self.iter+1 + 1, 1)

        # plot the change in design variables across each individual
        chainL = [self.XsGA[0,i] for i in range(len(self.XsGA[0]))]
        chainD = [self.XsGA[1,i] for i in range(len(self.XsGA[1]))]
        #polyL = [self.XsGA[2,i] for i in range(len(self.XsGA[2]))]
        #polyD = [self.XsGA[3,i] for i in range(len(self.XsGA[3]))]

        fig, ax = plt.subplots(2,1, sharex=True)
        ax[0].plot(iters, chainL, label='chain')
        #ax[0].plot(iters, polyL, label='polyester')
        ax[1].plot(iters, chainD, label='chain')
        #ax[1].plot(iters, polyD, label='polyester')
        ax[1].set_xlabel('individual evaluated')
        ax[1].set_ylabel('line diameter (mm)')
        ax[0].set_ylabel('line length (m)')
        ax[0].legend()
        ax[1].legend()
        #for i in range(len(gens)):
            #ax[0].axvline(x=gens[i], color='k')

        # plot the change in each constraint of each individual across the optimization
        Cnames = ['lay_length','rope_contact','offset','strength0','strength1']
        Cline = np.zeros_like(self.CsGA)
        fig, ax = plt.subplots(len(self.CsGA), 1, sharex=True)
        for i in range(len(self.CsGA)):
            for j in range(len(self.CsGA[i])):
                if self.CsGA[i,j] < -9000:
                    ax[i].plot(iters[j], 0, 'rx')
                    Cline[i,j] = 0
                else:
                    ax[i].plot(iters[j], self.CsGA[i,j], 'bo')
                    Cline[i,j] = self.CsGA[i,j]
            ax[i].set_ylabel(f'{Cnames[i]}')
            ax[i].plot(iters, Cline[i,:], 'g')
            ax[i].plot(iters, np.zeros(len(iters)), 'r')
        ax[-1].set_xlabel('individual evaluated')
        
        # plot the change in objective (cost) of each individual across the optimization
        Fnames = ['Line Cost', 'Anchor Cost', 'Total Cost']
        fig, ax = plt.subplots(1, 1, sharex=True)
        ax.plot(iters, self.FsGA[0,:], label='Line Cost')
        ax.plot(iters, self.FsGA[1,:], label='Anchor Cost')
        ax.plot(iters, self.FsGA[2,:], label='Total Cost')
        ax.set_ylabel('Cost ($)')
        ax.set_xlabel('individual evaluated')
        ax.legend()
        '''
        # to calculate all the iterations (individuals) that had all nonzero constraints
        a=[]
        for j in range(len(ld.CsGA[0])):
            if np.all(ld.CsGA[:,j]>0):
                a.append(j)
        
        # attempting to plot only the nonzero points on the plot
        for i in range(len(ld.FsGA)):
            for j in range(len(ld.FsGA[i])):
                if ld.FsGA[i,j]==np.nan:
                    ld.FsGA[i,j] = None
        '''



    def storeGA(self, val, i, type='X', name='', index=0):
        '''function to store the design variable vector, constraint values, and objective results for each iteration, based on self.iter,
        where self.iter is updated every time updateDesign is called'''

        #if method=='GA':
        if type=='X':
            self.XsGA[:,i] = val
        
        elif type=='C':
            confunnames = [self.confundict[con[0]].__name__ for con in self.constraints]
            for c in range(len(np.unique(confunnames))):
                if name==confunnames[c]:
                    self.CsGA[c+index,i] = val
        
        elif type=='F':
            self.FsGA[:,i] = val


    """
    def checkGA(self, type='normal'):
        '''function to check the feasibility of a design, mostly used in a GA, to ensure that LineDesign can even run it.
        More specifically, if the GA comes up with a design with sum of line lengths longer than span+depth, it will return False'''

        total_linelength = sum([self.ss.lineList[i].L for i in range(self.nLines)])
        
        if type=='normal':
            Lmax0 = self.span-self.rBFair[0] + self.depth+self.rBFair[2]                             # maximum possible line length allowable in equilibrium position
            if total_linelength > Lmax0:
                return False
            else:
                return True

        elif type=='offset':
            Lmax1 = self.span-self.rBFair[0]-self.x_mean_out-self.x_ampl + self.depth+self.rBFair[2]     # maximum possible line length allowable in offset position
            if total_linelength > Lmax1:
                return False
            else:
                return True
    """


    # :::::::::: solver functions ::::::::::
    
    # the original function from LineDesign, for tuning the line's horizontal tension
    def func_TH_L(self, Xl, args):
        '''Apply specified section L, return the horizontal pretension error.'''
        self.setSectionLength(Xl[0], self.iL)
        # option to setOffset?
        self.ss.staticSolve()
        if args['direction']=='horizontal':
            Fx = abs(self.ss.fB_L[0])  # horizontal fairlead tension
        elif args['direction']=='norm':
            Fx = np.linalg.norm(self.ss.fB_L)
        
        return np.array([Fx - self.fx_target]), dict(status=1) , False
    
    
    def func_kH_L(self, Xl, args):
        '''Apply specified section L, return the horizontal stiffness error.'''
        self.ss.lineList[self.iL].setL(Xl[0])
        # option to setOffset?
        self.staticSolve()
        Kx = self.KB_L[0,0]  # horizontal inline stiffness
        
        return np.array([Kx - self.kx_target]), dict(status=1) , False
    
    
    def func_fx_L(self, Xl, args):
        '''Apply specified section L, return the Fx error when system is offset.'''
        '''Function for solving line length that achieves equilibrium at a specified offset. 
        Expects xOffset, fx_target, and angles as keys in args dictionary.
        Receives line length and returns net force at xOffset.'''
        
        if self.ms:
            for ss in self.ms.lineList:
                ss.lineList[self.iL].setL(Xl[0])
            self.ms.bodyList[0].setPosition([args['xOffset'], 0,0,0,0,0])
            self.ms.solveEquilibrium()
            Fx = -self.ms.bodyList[0].getForces()[0]
        else:
            self.ss.lineList[self.iL].setL(Xl[0])
            self.ss.setOffset(args['xOffset'])
            Fx = np.abs(self.ss.fB_L[0])  # horizontal fairlead tension.
        
        if 'display' in args:
            if args['display'] > 2: 
                print(f" Xl is {Xl[0]:6.3f}  and Fx is {Fx/1e3:10.0f} kN  so error is  {(Fx+self.fx_target)/1e3:8.0f} kN")
        
        return np.array([Fx - self.fx_target]), dict(status=1), False
    
    
    def func_kx_L(self, Xl, args):    # evaluate how close the system horizontal stiffness is compared to the kx_target
            
        for ss in self.ms.lineList:  # go through each Subsystem
            ss.lineList[self.iL].setL(Xl[0])  # update the section length
            
        # option to setOffset?
        self.ms.bodyList[0].setPosition([0, 0,0,0,0,0])  # apply offset
        self.ms.solveEquilibrium()
        Kx = self.ms.getCoupledStiffness()[0,0]  # mooring system stiffness in x
        
        if 'display' in args:
            if args['display'] > 1: 
                print(f" Xl is {Xl[0]:6.3f}  and Kx is {Kx/1e3:10.0f} kN/m  so error is  {(Kx+self.kx_target)/1e3:8.0f} kN/m")
                
        return np.array([Kx - self.kx_target]), dict(status=1), False
    
    
                
    def func_fx_x(self, X, args):
    
        self.ms.bodyList[0].setPosition([X[0], 0,0,0,0,0])  # apply offset
        self.ms.solveEquilibrium()
        FxMoorings = self.ms.bodyList[0].getForces()[0]  # net mooring force in x
        FxApplied = args['FxApplied']
        
        return np.array([FxApplied + FxMoorings]), dict(status=1), False
        
        
    
    def step_fx_x(self, X, args, Y, oths, Ytarget, err, tols, iter, maxIter):
        ''' this now assumes tols passed in is a vector'''
        
        FxMoorings = self.ms.bodyList[0].getForces()[0]  # net mooring force in x
        FxApplied = args['FxApplied']
        
        dY = FxApplied + FxMoorings
        
        Kx = self.ms.bodyList[0].getStiffnessA(lines_only=True)[0,0]
        
        if Kx > 0:
            dX = dY/Kx
            
        else:  # backup case, just move 10 m
            
            dX = np.sign(dY)*10
        
        return np.array([dX])


    
    
    def setAnchoringRadius(self, a):
        '''Sets the anchoring radius, including of any LineDesign MoorPy 
        System. Input is the anchoring radius from platform centerline [m].
        '''
        
        if a < 0:
            raise Exception("The value passed to setAnchoringRadius must be positive.")
        
        self.rad_anch = float(a)
        
        self.dd['span'] = self.rad_anch - self.rBFair[0]
        self.ss.span  = float(self.dd['span'])
        
        self.ss.setEndPosition([-self.rad_anch, 0, -self.depth], endB=False)
        
        # Now handle the MoorPy system, if there is one, moving the anchor points
        if self.ms:
            for i, heading in enumerate(self.headings):
                rotMat = rotationMatrix(0, 0, np.radians(heading))
                self.ms.pointList[2*i].setPosition(np.matmul(rotMat, [self.rad_anch, 0, -self.depth]))
                
                # set subsystem span if needed... <<<
                self.ms.lineList[i].span = float(self.dd['span'])
    
    def setSectionLength(self, L, i):
        '''Sets the length of a section, including in the MoorPy System if there 
        is one. Overrides Mooring.setSectionLength'''
        
        # First call the Mooring version of this method, which handles the subsystem
        Mooring.setSectionLength(self, L, i)
        
        # Now handle the MoorPy system, if there is one
        if self.ms:
            for ss in self.ms.lineList:
                ss.lineList[i].setL(L)
    

    # ::::::::::::::::::::::::::::::: constraint functions :::::::::::::::::::::::::::::::
    
    # Each should return a scalar C where C >= 0 is valid and C < 0 is violated.
    
    def con_Kx(self, X, index, value, display=0):
        '''This ensures Kx, the effective horizontal stiffness, is greater than a given value. 
        Note: this constraint doesn't use the index input.'''
        
        Kx = self.ss.KB_L[0,0]     # get effective horizontal stiffness at current/undisplaced position
        c = Kx - value
                
        return c
    
    
    def con_total_length(self, X, index, value):
        '''This ensures that the total length of the Mooring does not result in a fully slack Mooring
        (ProfileType=4) in its negative extreme mean position'''
        # ['max_line_length', index, value] # index and value are completely arbitrary right now
        
        Lmax = 0.95*(self.ss.span + self.depth+self.rBFair[2])

        total_linelength =  sum([self.ss.lineList[i].L for i in range(self.nLines)])
        c = Lmax-total_linelength

        return c
    
    # ----- offset constraints -----
    
    def getOffset(self, FxApplied, headings=[]):
        '''Computes the horizontal offset of the body in response to an 
        applied horizontal force, considering all mooring lines, by solving
        for offset at which mooring reaction force equals FxApplied.'''
        
        # Ensure everything is switched back to status stiffnesses
        self.ms.revertToStaticStiffness()
        
        # Solve for the surge offset that matches the applied force
        '''
        x, y, info = dsolve2(self.func_fx_x, [0], step_func=self.step_fx_x,
                             args=dict(FxApplied=FxApplied, 
                             heading=headings), tol=[0.01], Xmin=[-1e5], 
                             Xmax=[1e5], dX_last=[10], stepfac=4, display=0)
                             
        return x[0]
        '''
        self.ms.bodyList[0].f6Ext = np.array([FxApplied, 0,0, 0,0,0])
        self.ms.solveEquilibrium(DOFtype='both')
        #offset = self.ms.bodyList[0].r6[0]
        #self.ms.bodyList[0].f6Ext = [0,0,0,0,0,0]
        #self.ms.bodyList[0].setPosition([0,0,0,0,0,0])
        #self.ms.solveEquilibrium()
        #return offset
        return self.ms.bodyList[0].r6[0]
        
    
    def con_offset0(self, X, index, value):
        '''This ensures that the system does not offset by a certain amount in its unloaded position'''
        
        # placeholder, this method may not make sense as-is
        return value - self.getOffset(0)
    
    
    def con_offset(self, X, index, value):
        '''This ensures that the system does not offset by a certain amount in its fully loaded position'''

        return value - abs(self.x_mean_eval)
    
    # ----- lay length constraints -----
    
    def con_lay_length(self, X, index, threshold, display=0):
        '''This ensures there is a minimum amount of line on the seabed at the +extreme displaced position.'''
        return self.ss.getLayLength(iLine=index) - threshold  + self.ss.LayLen_adj

    def con_max_td_range(self, X, index, threshold, display=0):
        '''Ensures the range of motion of the touchdown point betweeen the
        range of offsets is less then a certain distance.
        This constraint is for the system as a whole (index is ignored) and 
        must have offset='other' so that it's evaluated at the end.'''
        return threshold - (self.max_lay_length - self.min_lay_length)
    
    
    # ----- rope contact constraints -----
    
    def con_rope_contact(self, X, index, threshold, display=0):
        '''Ensures the first line node doesn't touch the seabed by some 
        minimum clearance.'''

        return self.ss.getPointHeight(index) - threshold    # compute the constraint value     
    
    
    # ----- strength constraints -----
    
    def con_strength(self, X, index, threshold, display=0):
        '''This ensures the MBL of the line is always greater than the maximum 
        tension the line feels times a safety factor.'''
        return self.ss.getTenSF(index) - threshold
    
    def con_min_tension(self, X, index, threshold, display = 0):
        '''Ensures that the minimum line tension is above a threshold'''
        return self.ss.getMinTen(index) - threshold
    
    def con_curvature(self, X, index, threshold, display=0):
        '''Ensure that the MBR of the cable is always greater than the maximum 
        actual curvature times a safety factor.'''
        return self.ss.getCurvSF(index) - threshold    
    

    def getDamage(self, index, display=0):    
        ''' method to predict fatigue damage based on previous iteration'''
        
        damage = self.damage
        
        if sum(damage) == 0:
            raise ValueError("Fatigue damage from previous iteration was not provided")
        

        sumdamage = 0

        
        #fatigue_headings are loading direction for fatigue dynamic factor calculation. must match order of damage in self.damage
        for i, ang in enumerate(self.fatigue_headings):
            
            #apply fx_target at direction in fatigue_headings
            self.ms.bodyList[0].f6Ext = np.array([self.fx_target*np.cos(np.radians(ang)), self.fx_target*np.sin(np.radians(ang)),0, 0,0,0])
            self.ms.solveEquilibrium(DOFtype='both')
            
            #store offset
            offsetx = self.ms.bodyList[0].r6[0]
            offsety = self.ms.bodyList[0].r6[1]
            
            #tension 1
            Ten1 = max(np.linalg.norm(self.ms.lineList[self.ms_fatigue_index].lineList[index].fA),np.linalg.norm(self.ms.lineList[self.ms_fatigue_index].lineList[index].fB))
            
            #set force back to zero  
            self.ms.bodyList[0].f6Ext = [0,0,0,0,0,0]
            
            #add dx to previous offset to get dtdx (slope of tension-displacement curve)
            dx = 0.5
            self.ms.bodyList[0].setPosition(np.array([offsetx + dx*np.cos(np.radians(ang)),offsety+dx*np.sin(np.radians(ang)),0,0,0,0]))     # move the body by the change in distance
            self.ms.solveEquilibrium()
            
            #tension 1
            Ten2 = max(np.linalg.norm(self.ms.lineList[self.ms_fatigue_index].lineList[index].fA),np.linalg.norm(self.ms.lineList[self.ms_fatigue_index].lineList[index].fB))
            
            #slope of tension-displacement curve at fx_target applied at ang
            dTdx = (Ten2 - Ten1)/dx
            
            #ratio is based on fatigue damage equation (Tension/MBL)^m, where m = 3 for chain
            MBL_corroded = self.ms.lineList[self.ms_fatigue_index].lineList[index].type['MBL'] * ( (self.ms.lineList[self.ms_fatigue_index].lineList[index].type['d_nom'] - (self.corrosion_mm/1000)) / self.ms.lineList[self.ms_fatigue_index].lineList[index].type['d_nom'] )**2
            ratio = (dTdx/ MBL_corroded)**3
            
            #ratio is multipled by the inputted previous iteration damage*MBL1/dTdx1
            sumdamage = sumdamage + ratio * damage[i]

            
        return sumdamage 


    def con_damage(self, X, index, threshold, display=0):
        '''constraint method to ensure the scaled fatigue damage meets required fatigue damage'''
        
        return threshold - self.getDamage(index, display=display)  
    
    
    def getYawStiffness(self, x_offset, display=0):
        '''method to calculate the yaw stiffness of the whole mooring system using an analytical equation'''

        yawstiff = 0
        # calculate stiffness in different situations
        for i, ang in enumerate(self.headings):
            spacing_x = self.span*np.cos(np.radians(ang)) - x_offset         # x distance from offset fairlead to anchor point
            spacing_y = self.span*np.sin(np.radians(ang))                    # y distance from offset fairlead to anchor point
            spacing_xy= np.linalg.norm([spacing_x, spacing_y])                  # radial distance from offset fairlead to anchor point
            self.setPosition(spacing_xy-self.span)
            tau0 = self.ss.fB_L[0]       # calculate the horizontal tension on the body from the 1 line
            
            # analytic equation for yaw stiffness for each mooring line heading
            yawstiff += (-tau0/spacing_xy)*self.ss.rBFair[0]**2 + -tau0*self.ss.rBFair[0]

        self.ss.setOffset(0)  # restore to zero offset and static EA
        
        return yawstiff

    
    def con_yaw_stiffness0(self, X, value, display=0):
        '''constraint method to ensure the yaw stiffness of the mooring system represented by this line design meets a certain yaw stiffness requirement,
        quasi-statically, and in the undisplaced position'''

        c = self.getYawStiffness(x_offset=0, display=display) - value              # compute the constraint value   
        
        return c

    def con_yaw_stiffness(self, X, index, value, display=0):
        '''constraint method to ensure the yaw stiffness of the mooring system represented by this line design meets a certain yaw stiffness requirement,
        quasi-statically, and in the extreme displaced position'''
        
        try:
            bodyPosition = np.array([-self.x_mean_in-self.x_ampl, 0,0,0,0,0])
            c = self.getYawStiffness(x_offset=bodyPosition[0], display=display) - value  # compute the constraint value    

        except Exception as e:  
            if self.noFail:
                c = -60000
            else:
                raise(e)

        return c
    

    # ----- shared line sag constraints -----
    
    def con_min_sag(self, X, index, threshold, display=0):
        '''Ensure the lowest point of a line section is below 
        a minimum depth.'''
        return threshold - self.ss.getSag(index)
    
    def con_max_sag(self, X, index, threshold, display=0):
        '''Ensures the lowest point of a line section is above
        a certain maximum depth.'''
        return self.ss.getSag(index) - threshold
    
    # ----- angle constraints -----
    def con_min_angle(self, X, index, threshold, display=0):
        '''Ensure the angle of a line section is above a minimum value.'''
        return self.ss.getAng(index) - threshold

    def con_max_angle(self, X, index, threshold, display=0):
        '''Ensure the angle of a line section is below a maximum value.'''
        return threshold - self.ss.getAng(index)

    # ----- utility functions -----
            
    def plotProfile(self, Xuvec=[1,0,0], Yuvec=[0,0,1], ax=None, color=None, title="", slack=False, displaced=True, figsize=(6,4), label=None):
        '''Plot the mooring profile in undisplaced and extreme displaced positions

        Parameters
        ----------
        Xuvec : list, optional
            plane at which the x-axis is desired. The default is [1,0,0].
        Yuvec : lsit, optional
            plane at which the y-axis is desired. The default is [0,0,1].            
        ax : axes, optional
            Plot on an existing set of axes
        color : string, optional
            Some way to control the color of the plot ... TBD <<<
        title : string, optional
            A title of the plot. The default is "".
        slack : bool, optional
            If false, equal axis aspect ratios are not enforced to allow compatibility in subplots with axis constraints.
        dispalced : bool, optional
            If true (default), displaced line profiles are also plotted.

        Returns
        -------
        fig : figure object
            To hold the axes of the plot
        ax: axis object
            To hold the points and drawing of the plot

        '''
      
        # if axes not passed in, make a new figure
        if ax == None:
            fig, ax = plt.subplots(1,1, figsize=figsize)
            ax.set_xlabel('Horizontal distance (m)')
            ax.set_ylabel('Depth (m)')
        else:
            fig = plt.gcf()   # will this work like this? <<<


        if displaced:
            offsets = [0, self.x_mean_out+self.x_ampl, -self.x_mean_in-self.x_ampl]
        else:
            offsets = [0]
        
        for x in offsets:
            
            alph = 1 if x==0 else 0.5  # make semi-transparent for offset profiles
            
            self.ss.setOffset(x)

            #ax.plot(self.rB[0], self.rB[2],'ko',markersize = 2)  # fairlead location
            ax.plot(x, 0,'ko',markersize = 2)  # platform ref point location
            # self.ss.drawLine2d(0,ax)
            for i, line in enumerate(self.ss.lineList):
                if i != 0:
                    label = None
                if color==None:     # alternate colors so the segments are visible
                    if line.type['material'][0]=='c':
                        line.drawLine2d(0, ax, color=[.1, 0, 0], alpha=alph, Xuvec=Xuvec, Yuvec=Yuvec, Xoff=-self.rad_fair, label=label)
                        if self.shared==1:  # plot other half too if it's a shared line where only half is modeled <<<
                            line.drawLine2d(0, ax, color=[.1, 0, 0], alpha=alph, Xuvec=-np.array(Xuvec), Yuvec=Yuvec, Xoff=-self.span-self.rad_fair, label=label) 
                    elif 'nylon' in line.type['material']:
                        line.drawLine2d(0, ax, color=[.8,.8,.2], alpha=alph, Xuvec=Xuvec, Yuvec=Yuvec,Xoff=-self.rad_fair, label=label)
                    else:
                        line.drawLine2d(0, ax, color=[.3,.5,.5], alpha=alph, Xuvec=Xuvec, Yuvec=Yuvec,Xoff=-self.rad_fair, label=label)
                        if self.shared==1:  # plot other half too if it's a shared line where only half is modeled <<<
                            line.drawLine2d(0, ax, color=[.3,.5,.5], alpha=alph, Xuvec=-np.array(Xuvec), Yuvec=Yuvec, Xoff=-self.span-self.rad_fair, label=label)
                else:
                    line.drawLine2d(0, ax, color=color, alpha=alph, Xuvec=Xuvec, Yuvec=Yuvec,Xoff=-self.rad_fair, label=label)
                    if self.shared==1:  # plot other half too if it's a shared line where only half is modeled <<<
                        line.drawLine2d(0, ax, color=color, alpha=alph, Xuvec=-np.array(Xuvec), Yuvec=Yuvec, Xoff=-self.span-self.rad_fair, label=label)

            '''
            # plot points/weights/floats along the line >>> needs to be updated to account for Xuvec and Yuvec <<<
            for point in self.pointList:
                if point.number > 1 and point.number < self.nLines+1:
                    if point.v > 0:
                        ax.plot(point.r[0],point.r[2],'yo',markersize=5)
                    elif point.m > 0:
                        ax.plot(point.r[0],point.r[2],'ko',markersize=5)
                    else:
                        ax.plot(point.r[0],point.r[2],'bo',markersize=1)
            '''
            
        # make legend entries available
        if displaced:
            if not color==None:
                ax.plot(np.nan, np.nan, color=color, alpha=1, label="undisplaced") 
                ax.plot(np.nan, np.nan, color=color, alpha=0.5, label="displaced") 
            
        #ax.plot([self.ss.lineList[0].rA[0], 0], [-self.depth, -self.depth], color='k')
        # only force equal aspect ratio if "slack" keyword isn't specified (so that sharex=True, sharey-True plots are possible)
        if not slack:
            ax.axis("equal")
            
        ax.set_title(title)
        #ax.set_ylim(-1,1)
        
        self.ss.setOffset(0)  # return to undisplaced position
        self.ss.solveEquilibrium(tol=self.ss.eqtol)
        
        return fig, ax  # return the figure and axis object in case it will be used later to update the plot

    
    def plotCurves(self, ax=[], color="k", title=""):
        '''Plot key performance curves for the mooring as a function of offset

        Parameters
        ----------
        ax : axes, optional
            Plot on an existing set of axes
        title : string, optional
            A title of the plot. The default is "".

        Returns
        -------
        fig : figure object
            To hold the axes of the plot
        ax: axis object
            To hold the points and drawing of the plot

        '''
      
        # if axes not passed in, make a new figure
        if len(ax) == 0:
            fig, ax = plt.subplots(2,1, sharex=True)
            newFig=True
        else:
            if not len(ax) == 2:
                raise Exception("ax provided to plotCurves must be a list of 2 axes.")
            fig = plt.gcf()
            newFig = False
        
        x = np.linspace(-self.x_mean_in-self.x_ampl, self.x_mean_out+self.x_ampl, 50)
        
        Fx = np.zeros(len(x))
        Ts = np.zeros([len(x), len(self.ss.lineList)])
        
        # calculate values at each offset point
        for i in range(len(x)):                                      # go through each offset point
            
            self.ss.setOffset(x[i])  # offset the desired amount
            
            Fx[i] = self.ss.fB_L[0]  # get horizontal mooring force
            
            for j in range(len(self.ss.lineList)):                 # get upper end tension of each line segment
                Ts[i,j] = self.ss.lineList[j].TB
        
        # plots        
        ax[0].plot(x, -Fx/1e3, c=color)
        
        for j in range(len(self.ss.lineList)):    
            ax[1].plot(x, Ts[:,j]/1e3, c=color, dashes=[5-0.5*j, 0.5*j], label=f"segment {j+1}")            
            
        ax[0].set_ylabel("Fx (kN)") 
        ax[1].set_ylabel("Tension (kN)") 
        if newFig:  ax[1].legend()
        ax[1].set_xlabel("Offset (m)")    
        #fig.set_title(title)
        
        self.ss.setOffset(0)  # restore to undisplaced position
        
        return fig, ax  # return the figure and axis object in case it will be used later to update the plot

        
    def dump(self):
        '''Puts info about the mooring into a dictionary and returns it.'''
        
        self.objectiveFun([])  # ensure things are calculated
        
        info = dict(arrangement={}, design={}, performance={}, cost={})    # the dictionary and its top-level entries
        
        info['arrangement']['name']  = self.name
        '''
        info['design']['X'         ] = self.Xlast                                                       # design variables
        info['design']['Gdict'     ] = self.evaluateConstraints([])[1]                                  # dict of constraint names and values of evaluated constraint functions
        info['design']['Ls'        ] = [line.L                            for line in self.ss.lineList]    # length of each segment
        info['design']['Ds'        ] = [line.type['input_d']              for line in self.ss.lineList]    # *input* diameter  of each segment 
        info['design']['lineTypes' ] = [line.type['name']                 for line in self.ss.lineList]    # line type of each segment (may be redundant with what's in arrangement)
        info['design']['anchorType'] = self.anchorType                                                  # (may be redundant with what's in arrangement)
        info['design']['span'   ] = self.span                                                     # platform-platform of platfom-anchor horizontal span just in case it's changed
        info['design']['Ltot'      ] = sum([line.L for line in self.ss.lineList])                          # total mooring length
    
        info['performance']['Fx'] = self.fB_L[0]
        info['performance']['Kx'] = self.bodyList[0].getStiffness(tol=self.eqtol)[0,0]
    
        info['cost']['total'  ] = self.cost
        info['cost']['line'   ] = self.lineCost
        if not self.shared:
            info['cost']['anchor' ] = self.anchorMatCost
            info['cost']['install'] = self.anchorInstCost   # eventually should sort out if this represents the total installation cost
            info['cost']['decom'  ] = self.anchorDecomCost
        '''
        
        # this version converts out of numpy format for yaml export (should make a better system for this)
        info['design']['X'         ] = self.Xlast.tolist()                                              # design variables
        #info['design']['Gdict'     ] = self.evaluateConstraints([])[1]                                  # dict of constraint names and values of evaluated constraint functions
        info['design']['Ls'        ] = [float(line.L              )       for line in self.ss.lineList]    # length of each segment
        info['design']['Ds'        ] = [float(line.type['d_nom'])       for line in self.ss.lineList]    # *input* diameter  of each segment 
        info['design']['lineTypes' ] = [str(line.type['material'])        for line in self.ss.lineList]    # line type of each segment (may be redundant with what's in arrangement)
        info['design']['anchorType'] = self.anchorType                                                  # (may be redundant with what's in arrangement)
        info['design']['span'   ] = float(self.span)                                              # platform-platform of platfom-anchor horizontal span just in case it's changed
        info['design']['Ltot'      ] = float(sum([line.L for line in self.ss.lineList]))                   # total mooring length
    
        info['performance']['Fx'] = float(self.fB_L[0]    )
        info['performance']['Kx'] = float(self.KB_L[0,0])
    
        info['cost']['total'  ] = float(self.cost    )
        info['cost']['line'   ] = float(self.lineCost)
        if not self.shared==1:
            info['cost']['anchor' ] = float(self.anchorMatCost  )
            info['cost']['install'] = float(self.anchorInstCost )  # eventually should sort out if this represents the total installation cost
            info['cost']['decom'  ] = float(self.anchorDecomCost)
        
        
        return info
    

    def adjustConstraint(self, key, value):
        '''Modifies the value of an existing constraint'''
        for con in self.constraints:
            if con[0] == key:
                con[2] = value
                
    @staticmethod
    def getClumpMV(weight, rho=1025.0, g=9.81, **kwargs):
        
        '''A function to provide a consistent scheme for converting a clump weight/float magnitude to the 
        mass and volume to use in a MoorPy Point.'''
        
        if weight >= 0:                          # if the top point of the intermediate line has a clump weight
            pointvol = 0.0
            pointmass = weight*1000.0           # input variables are in units of tons (1000 kg), convert to kg
        else:
            pointvol = -weight*1200.0/rho  # input variables are still in tons. Assume additional 20% of BM mass
            pointmass = -weight*200.0
    
        return dict(m=pointmass, v=pointvol)



if __name__ == '__main__':
    
    # Example case from Stein
    '''
    settings = {}
    settings['rBFair'] = [58,0,-14]
    settings['x_ampl'] = 10     # xmax value is designed to be the "target" offset, used for solve_for = 'tension'
    settings['fx_target'] = 1.95e6
    settings['solve_for'] = 'none'
    settings['headings'] = [60, 180, 300]
     
    settings['name'] = 'chain-poly-chain'
    settings['lineTypeNames'] = ['chain','polyester','chain']
    settings['anchorType'] = 'suction'
    settings['allVars'] = [1000/10, 100, 120, 0, 800, 200, 0, 100, 120]
    settings['Xindices'] = ['c', 0, 'c', 'c', 1,    2, 'c', 'c', 'c']
    settings['Xmin'] = [10, 10, 10]
    settings['Xmax'] = [500, 10000, 500]
    settings['dX_last'] = [10, 10, 10]
    
    settings['constraints'] = [dict(name='rope_contact'  , index=1, threshold=5 , offset='min'),
                               dict(name='max_offset'    , index=0, threshold=60, offset='max')] 
    
    for j in range(len(settings['lineTypeNames'])):
        settings['constraints'].append(dict(name='tension_safety_factor', index=j, threshold=2.0, offset='max'))

    depth = 766.765
    ld = LineDesign(depth, **settings)

    ld.setNormalization()  # turn on normalization (important for COBYLA etc)

    start_time = time.time()
    #X, min_cost = ld.optimize(maxIter=12, plot=False, display=3, stepfac=4, method='dopt')
    X, min_cost = ld.optimize(maxIter=10, plot=True, display=3, stepfac=4, method='COBYLA')
    print("optimize time: {:8.2f} seconds".format((time.time() - start_time)))
    ld.objectiveFun(X, display=2)
    ld.evaluateConstraints(X, display=0)
    ld.updateDesign(X, display=0)
    ld.plotProfile()
    plt.show()
    '''


    depth = 200

    settings = {}
    settings['rBFair'] = [58,0,-14]
    settings['x_ampl'] = 10
    settings['fx_target'] = 1.95e6
    settings['headings'] = [60, 180, 300]

    settings['solve_for'] = 'none'
    #settings['solve_for'] = 'ghost'
    
    settings['name'] = 'DEA-chain-poly'  # <<< semitaut option
    settings['lineTypeNames'] = ['chain','polyester']
    settings['anchorType'] = 'drag-embedment'
    '''
    settings['allVars'] = [800/10, 400, 120, 0, 400, 200,]   
    settings['Xindices'] = ['c', 0, 1, 'c', 2, 3]
    settings['Xmin'] = [10, 10, 10, 10]
    settings['Xmax'] = [10000, 500, 800, 500]
    settings['dX_last'] = [10, 10, 10, 10]
    '''
    settings['allVars'] = [1000/10, 800, 120, 0, 80, 200,]   
    settings['Xindices'] = ['c', 0, 1, 'c', 'c', 2]
    settings['Xmin'] = [400, 10, 10]
    settings['Xmax'] = [2000, 500, 500]
    settings['dX_last'] = [10, 10, 10]
    
    '''
    settings['name'] = 'DEA-chain'  # <<< catenary option
    settings['lineTypeNames'] = ['chain']
    settings['anchorType'] = 'drag-embedment'
    settings['lay_target'] = 200
    settings['allVars'] = [1000/10, 1000, 120]    
    settings['Xindices'] = ['c', 0, 1]
    settings['Xmin'] = [500, 50]
    settings['Xmax'] = [1500, 300]
    settings['dX_last'] = [10, 10]
    
    
    #settings['solve_for'] = 'offset'
    settings['solve_for'] = 'tension'
    settings['Xindices'] = ['c', 's', 0]
    settings['Xmin'] = [10]
    settings['Xmax'] = [500]
    settings['dX_last'] = [10]
    settings['x_target'] = 34.560922734165835
    settings['x_mean_out'] = 34.560922734165835
    settings['x_mean_in'] = 60
    '''
    settings['constraints'] = [dict(name='min_lay_length', index=0, threshold=20, offset='max'),
                               dict(name='max_offset'    , index=0, threshold=25, offset='max'),
                               dict(name='rope_contact'  , index=1, threshold=5 , offset='min')]
    
    for j in range(len(settings['lineTypeNames'])):
        settings['constraints'].append(dict(name='tension_safety_factor', index=j, threshold=2.0, offset='max'))

    


    ld = LineDesign(depth, **settings)

    ld.setNormalization()  # turn on normalization (important for COBYLA etc)

    start_time = time.time()
    #X, min_cost = ld.optimize(maxIter=20, plot=False, display=3, stepfac=4, method='dopt')
    #X, min_cost = ld.optimize(maxIter=40, plot=True, display=3, stepfac=4, method='COBYLA')
    #X, min_cost = ld.optimize(maxIter=40, plot=True, display=3, stepfac=4, method='CMNGA')
    #X, min_cost = ld.optimize(maxIter=40, plot=True, display=3, stepfac=4, method='PSO')
    X, min_cost = ld.optimize(maxIter=40, plot=True, display=0, stepfac=4, method='bayesian')
    
    print('')
    print('Analyzing Results:')
    print( "    optimize time:                 {:8.2f} seconds".format((time.time() - start_time)))
    print( '    design variables (normalized):   ', [f"{x:8.3f}" for x in X])
    print( '    design variables (denormalized): ', [f"{x:8.2f}" for x in X*ld.X_denorm])
    print(f'    solved line length:               {ld.ss.lineList[ld.iL].L:8.2f} m')
    print('')

    ld.objectiveFun(X, display=2)
    ld.evaluateConstraints(X, display=2)
    ld.updateDesign(X, display=0)
    ld.plotProfile()
    plt.show()

    a = 2
