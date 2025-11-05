# Draft Cable Design class adapted from LineDesign - R

import moorpy as mp
from moorpy.subsystem import Subsystem

from famodel.cables.dynamic_cable import DynamicCable
import famodel.cables.cable_properties as cprops

from fadesign.fadsolvers import dopt2, doptPlot
from moorpy.helpers import getLineProps, getFromDict
from copy import deepcopy

import numpy as np
import matplotlib.pyplot as plt
import yaml
import time


class CableDesign(DynamicCable):
    '''
    The Dynamic Cable class inherits the properties of MoorPy's Subsystem class
    (i.e. solving for equilibrium) for the purpose of quasi-static design and
    analysis. Eventually the DynamicCable class will live in FAModel, and this 
    code will be streamlined to inherit it and just add design methods.
    
    Example allVars vector: X = [span, L, <B, Lmid, Ls>...]    
        where <  > section repeats and is composed of
        B - net buoyancy provided by all modules on this section [kN]
        Lmid - the buoyancy section midpoint along the cable arc length [m]
        Ls - the length of this buoyancy section (centered about the midpoint) [m]
    Xindices
        specify the design variable number, or optional key characters:
        c - constant, will not be changed
        r - the AllVars value will be interpreted as a ratio to the total length
        In other words, the actual value will be the specified value times L.
    '''

    def __init__(self, depth, cableType, buoyType, n=3, i_buoy=[1], mgdict = None, **kwargs):
        '''Creates a DynamicCable object to be used for evaluating or
        optimizing a dynamic cable design.

        Parameters
        ----------
        depth : float
            Water depth
        span : float
            Horizontal distance of dynamic cable [m].
        n : int
            Number of sections (typically alternating: cable, cable+buoyancy, ...)
        i_buoy : list
            List of section indices that can have buoyancy modules.
        cableType : dict
            Dictionary of bare cable properties.
        buoyType : dict
            Dictionary of buoyancy module properties.
        name : string
            Name of dictionary entry in cableProps yaml file to get data from.
        X0 : array
            Initial design variable values (length n).
        offset : float
            Maximum mean/steady offset in surge [m].
        '''
        
        self.display = getFromDict(kwargs, 'display', default=0)
        
        # add the parameters set by the input settings dictionary
        self.name       = getFromDict(kwargs, 'name', dtype=str, default='no name provided')
        
        # set up the mooring system object with the basics from the System class
        rho        = getFromDict(kwargs, 'rho', default=1025.0)
        g          = getFromDict(kwargs, 'g'  , default=9.81)
        
        # ----- set model-specific parameters -----

        self.shared  = getFromDict(kwargs, 'shared',  default=0)      # flag to indicate shared line
        self.rBFair     = getFromDict(kwargs, 'rBFair', shape=-1, default=[0,0,0])  # [m] end coordinates relative to attached body's ref point
        self.nLines = n       # number of sections
        self.i_buoy = i_buoy  # index of any sections that have buoyancy modules
        self.bs = [0 for i in range(0,len(i_buoy))]
        
        #-------set marine growth parameters---------------------
        self.MG = getFromDict(kwargs, 'MG', default = False)
        if self.MG:
            if mgdict == None:
                raise Exception('mgdict must be provied if MG == True')
            else:
                self.mgdict = mgdict
     
        # ============== set the design variable list ==============
        self.ignore_static = getFromDict(kwargs, 'ignore_static', default=False)
        
        self.allVars    = getFromDict(kwargs, 'allVars' , shape=2 + 3*len(self.i_buoy))
        
        # set the design variable type list
        if 'Xindices' in kwargs:
            self.Xindices = list(kwargs['Xindices'])
            if not len(self.Xindices)==len(self.allVars):
                raise Exception("Xindices must be the same length as allVars")
        else:
            raise Exception("Xindices must be provided.")
        

        # number of design variables (the design vector is the length of each
        # find the largest integer to determine the number of desired design variables
        self.nX = 1 + max([ix for ix in self.Xindices if isinstance(ix, int)])


        # check for errors in Xindices
        for i in range(self.nX):
            if not i in self.Xindices:
                raise Exception(f"Design variable number {i} is missing from Xindices.")
        # entries must be either design variable index or constant/solve/ratio flags
        valid = list(range(self.nX))+['c','r']
        for xi in self.Xindices:
            if not xi in valid:
                raise Exception(f"The entry '{xi}' in Xindices is not valid. Must be a d.v. index, 'c', or 'r'.")
        
        # check for 'r' variable option
        self.rInds = [i for i,xi in enumerate(self.Xindices) if xi=='r']
        for i in range(len(self.rInds)):
            if self.allVars[self.rInds[i]] >= 1.0 or self.allVars[self.rInds[i]] <= 0.0:
                raise Exception("The ratio variable needs to be between 1 and 0")
        
        
        # ----- Initialize some objects -----
        
        self.span = self.allVars[0]
        self.L = self.allVars[1]
        
        # Store the bare cable type by itself for easy access (TODO: reduce redundancy)
        self.cableType = cableType
        self.buoyType = buoyType
        
        # make a dummy design dictionary for Mooring to make a Subsystem with???
        dd = {}
        
        # The bare cable properties dict
        dd['cable_type'] = cableType
        
        #length properties
        dd['length'] = self.L
        
        #span
        dd['span'] = self.span
        
        # Buoyancy section properties

        for i in range(len(i_buoy)):
            
            # Net buoyancy per buoyancy module [N]
            F_buoy = (rho - buoyType['density'])*g*buoyType['volume']
            
            # Buoyancy
            if self.shared == 2 and i ==0:
                N_modules = 1000*self.allVars[3*i+2] / (F_buoy) / 2  # split buoyancy force across the full length
            else:
                N_modules = 1000*self.allVars[3*i+2] / F_buoy  # my not be an integer, that's okay


            # L_mid (position along cable)
            if self.Xindices[3*i + 3] == 'r':
                
                #set equal to ratio * cable length
                L_mid = self.allVars[3*i+3] * self.L         
            else:
                L_mid = self.allVars[3*i+3]
            
            # Spacing
            spacing = self.allVars[3*i+4] / (N_modules - 1)
            
            if N_modules > 0:
                if not 'buoyancy_sections' in dd:
                    dd['buoyancy_sections'] = []
                dd['buoyancy_sections'].append(dict(L_mid=L_mid, 
                                                    module_props=buoyType, 
                                                    N_modules = N_modules,
                                                    spacing = spacing))
        
        # Call Mooring init function (parent class)
        if self.shared == 1:
            
            DynamicCable.__init__(self, 'designed cable', dd=dd, 
                     rA=[self.span,0,self.rBFair[2]], rB=self.rBFair,
                     rad_anch=self.span, rad_fair=self.rBFair[0], z_anch=-depth, 
                     z_fair=self.rBFair[2], rho=rho, g=g, span=self.span, length=self.L, shared = self.shared)  # arbitrary initial length
            
        elif self.shared == 2:
            DynamicCable.__init__(self, 'designed cable', dd=dd, 
                     rA=[-0.5*self.span-self.rBFair[0], 0, -1], rB=self.rBFair,
                     rad_anch=self.span, rad_fair=self.rBFair[0], z_anch=-depth, 
                     z_fair=self.rBFair[2], rho=rho, g=g, span=self.span, length=self.L, shared = self.shared)  # arbitrary initial length

        else:
            DynamicCable.__init__(self, 'designed cable', dd=dd, 
                     rA=[self.span,0,-depth], rB=self.rBFair,
                     rad_anch=self.span, rad_fair=self.rBFair[0], z_anch=-depth, 
                     z_fair=self.rBFair[2], rho=rho, g=g, span=self.span, length=self.L, shared = self.shared)  # arbitrary initial length
                     
        
        
        # now make Subsystem, self.ss
        self.createSubsystem(case=int(self.shared))
        self.ss.eqtol= 0.05  # position tolerance to use in equilibrium solves [m]

        # amplification factors etc.
        self.DAFs = getFromDict(kwargs, 'DAFs', shape=self.nLines+2, default=1.0)   # dynamic amplication factor for each line section, and anchor forces (DAFS[-2] is for vertical load, DAFS[-1] is for horizontal load)
        self.Te0 = np.zeros([self.nLines,2])   # undisplaced tension [N] of each line section end [section #, end A/B]
        self.LayLen_adj = getFromDict(kwargs, 'LayLen_adj', shape=0, default=0.0) # adjustment on laylength... positive means that the dynamic lay length is greater than linedesign laylength
        self.damage = getFromDict(kwargs, 'damage', shape = -1, default = 0.0) #Lifetime fatigue damage from previous iteration for wind/wave headings in self.headings and the 180 degree reverse of self.headings
        #self.unload(f'{configuration}.dat') 


        # ----- Set solver and optimization settings -----
        self.x_mean  = getFromDict(kwargs, 'offset', default=0)
        self.x_ampl  = getFromDict(kwargs, 'x_ampl', default=10)  # [m] expected wave-frequency motion amplitude about mean

        self.eqtol = 0.002    # position tolerance to use in equilibrium solves [m]
        self.noFail = False   # can be set to True for some optimizers to avoid failing on errors
        
        self.iter = -1  # iteration number of a given optimization run (incremented by updateDesign)
        self.log = dict(x=[], f=[], g=[])  # initialize a log dict with empty values 
        
        
        # ----- optimization stuff -----
        # get design variable bounds and last step size
        self.Xmin       = getFromDict(kwargs, 'Xmin'   , shape=self.nX, default=np.zeros(self.nX))      # minimum bounds on each design variable
        self.Xmax       = getFromDict(kwargs, 'Xmax'   , shape=self.nX, default=np.zeros(self.nX)+1000) # maximum bounds on each design variable
        self.dX_last    = getFromDict(kwargs, 'dX_last', shape=self.nX, default=[]) # 'last' step size for each design variable

        if len(self.Xmin) != self.nX or len(self.Xmax) != self.nX or len(self.dX_last) != self.nX:
            raise Exception("The size of Xmin/Xmax/dX_last does not match the number of design variables")
        
        #set up initial design variable values from allVars input 
        self.X0 = np.array([self.allVars[self.Xindices.index(i)] for i in range(self.nX)])

        # initialize the vector of the last design variables, which each iteration will compare against
        self.Xlast = np.zeros(self.nX)
        
        self.X_denorm = np.ones(self.nX)  # normalization factor for design variables
        self.obj_denorm = 1.0  # normalization factor for objective function


        # ----- set up the constraint functions and lists -----

        if 'constraints' in kwargs:
            self.constraints = kwargs['constraints']
        else:
            self.constraints = {}
            #raise ValueError('A constraints dictionary must be passed when initializing a new Mooring')

        # a hard-coded dictionary that points to all of the possible constraint functions by name
        self.confundict = {"max_total_length"      : self.con_total_length,   # maximum total length of combined line sections
                           "min_lay_length"        : self.con_lay_length,     # minimum length of a line section on the seabed
                           "max_lay_length"        : self.con_max_lay_length,     # minimum length of a line section on the seabed
                           "tension_safety_factor" : self.con_strength,       # minimum ratio of MBL/tension for a section
                           "overall_tension_safety_factor" : self.con_overall_strength,       # minimum ratio of MBL/tension for all sections
                           "curvature_safety_factor":self.con_curvature,      # minimum ratio of curvature_limit/curvature for a section
                           "max_curvature"          :self.con_max_curvature,      # minimum ratio of curvature_limit/curvature for a section
                           "min_sag"               : self.con_min_sag,        # minimum for the lowest point of a section
                           "max_sag"               : self.con_max_sag,        # maximum for the lowest point of a section
                           "max_hog"               : self.con_max_hog,        # maximum for the highest point of a section
                           "max_touchdown_range"   : self.con_max_td_range,   # maximum for the lowest point of a section
                           }

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
                            self.updateDesign(X)
                            
                            # get the constraint value from the internal list
                            val = self.convals[ii]
                            
                            return val
                        return func()
                    return conf_maker

                # add the conf function to the conList
                self.conList.append(externalConFun(con['name'], i))
                
                # Save the default/recommended normalization constant
                
                if con['name'] in ['max_total_length']:
                    self.con_denorm_default[i] = con['threshold'] # sum([line.L for line in self.ss.lineList])
                
                elif con['name'] in ['tension_safety_factor', 'curvature_safety_factor']:
                    self.con_denorm_default[i] = 4*con['threshold']
                
                elif con['name'] in ['max_curvature']:
                    self.con_denorm_default[i] = 4*con['threshold']
                
                elif con['name'] in ['min_lay_length', 'min_sag', 'max_sag', 'max_hog', 'max_touchdown_range']:
                    self.con_denorm_default[i] = depth
                
            else:
                raise ValueError("Constraint parameter "+con['name']+" is not a supported constraint type.")



        # ----- Set up the cable properties -----
        '''
        # For now, this will load the cable properties YAML and manually add
        # the selected cable type to the MoorPy system.
        with open(cableProps) as file:
            source = yaml.load(file, Loader=yaml.FullLoader)
        
        # Get dictionary of the specified cable type from the yaml
        di = source['cable_types'][name]  
        
        cableType = self.makeCableType(di, name)  # Process/check it into a new dict
        # ^^^ I forget why this is done
        
        # Save some constants for use when computing buoyancy module stuff
        
        self.d0 = cableType['d_vol']  # diameter of bare dynamic cable
        self.m0 = cableType['m']      # mass/m of bare dynamic cable
        self.w0 = cableType['w']      # weight/m of bare dynamic cable
        
        #self.rho_buoy = cableType['rho_buoy']  # aggregate density of buoyancy modules [kg/m^3]
        '''

        '''
        # ----- set up the dynamic cable in MoorPy -----
        
        lengths = self.allVars[:self.nLines] # Length of each section [m] (first n entries of allVars)
        types = [] 
        
        # Set up line types list
        for i in range(self.nLines):
            # Give the buoyancy sections their own type so they can be adjusted independently
            if i in self.i_buoy:
                types.append(deepcopy(cableType))
            
            else:  # All bare cable sections can reference the same type
                types.append(cableType)  
        
        # call to the Subsystem method to put it all together
        if self.shared == True:
            
            # set second platform connection at the same coordinates as the first platform connection
            self.rAFair = self.rBFair
            self.makeGeneric(lengths, types, suspended = 1)
        else:
            self.makeGeneric(lengths, types)
        '''
        # initialize and equilibrate this initial cable
        self.ss.initialize()
        self.ss.maxIter = 5000
        self.ss.setOffset(0)
        self.updateDesign(self.X0, normalized=False)  # assuming X0/AllVars is not normalized



    def updateDesign(self, X, display=0, display2=0, normalized=True):
        '''updates the design with the current design variables (X).
        
        Example allVars vector: X = [span, L, <B, Lmid, Ls>...]    
            where <  > section repeats and is composed of
            B - net buoyancy provided by all modules on this section [N]
            Lmid - the buoyancy section midpoint along the cable arc length
            Ls - the length of this buoyancy section (centered about the midpoint)
        Xindices
            specify the design variable number, or optional key characters:
            c - constant, will not be changed
            r - the AllVars value will be interpreted as a ratio to the total length
            In other words, the actual value will be the specified value times L.
        
        If self.shared==2, then the buoyancy sections are measured from the
        center point of the cable, and are assumed to be mirrored on both sides.
        '''
        
        # Design vector error checks
        if len(X)==0:                   # if any empty design vector is passed (useful for checking constraints quickly)
            return
        elif not len(X)==self.nX:
            raise ValueError(f"DynamicCable.updateDesign passed design vector of length {len(X)} when expecting length {self.nX}")
        elif any(np.isnan(X)):
            raise ValueError("NaN value found in design vector")
        
        # If X is normalized, denormalize (scale) it up to the full values
        if normalized:
            X = X*self.X_denorm
        
        
        # If any design variable has changed, update the design and the metrics
        if not all(X == self.Xlast):
        
            self.Xlast = np.array(X)   # record the current design variables

            if self.display > 1:
                print("Updated design")
                print(X)
            
            self.iter += 1
            
            # ----- Apply the design variables to update the design -----
            
            # Update span
            dvi = self.Xindices[0]              # design variable index - will either be an integer or a string
            if dvi in range(self.nX):           # only update if it's tied to a design variable (if it's an integer)
                self.span = X[dvi]
                self.dd['span'] = self.span
            
            # Update total cable length
            dvi = self.Xindices[1]
            if dvi in range(self.nX):
                self.L = X[dvi]
                self.dd['length'] = X[dvi]  # some redundancy - need to streamline DynamicCable
            
            
            # Update each buoyancy section
            for i in range(len(self.i_buoy)):
                
                bs = self.dd['buoyancy_sections'][i]  # shorthand for the buoyancy section dict
                
                # Net buoyancy per buoyancy module [N]
                F_buoy = (self.rho - bs['module_props']['density'])*self.g*bs['module_props']['volume']
                
                
                # Buoyancy
                dvi = self.Xindices[3*i+2]  # buoyancy design variable index
                if dvi in range(self.nX):  # only update if it's tied to a design variable
                    if self.shared == 2 and i == 0:
                        bs['N_modules'] = 1000*X[dvi] / F_buoy / 2  # my not be an integer, that's okay
                    else:
                        bs['N_modules'] = 1000*X[dvi] / F_buoy  # my not be an integer, that's okay

                
                # L_mid (position along cable)
                dvi = self.Xindices[3*i+3]  # buoyancy design variable index
                if dvi in range(self.nX):  # only update if it's tied to a design variable
                    bs['L_mid'] = X[dvi]
                elif dvi == 'r':
                    bs['L_mid'] = self.allVars[3*i+3] * self.L 
                    
                # Spacing
                dvi = self.Xindices[3*i+4]  # buoyancy design variable index
                if dvi in range(self.nX):  # only update if it's tied to a design variable
                    length = X[dvi]
                else:
                    length = self.allVars[3*i+4]
                    
                bs['spacing'] = length / (bs['N_modules'] - 1)
                
                #store the buoyancy module spacing
                self.bs[i] = bs['spacing'] 
            
            
            
            # get these design dictionary changes applied in DynamicCable
            if len(self.i_buoy) > 0:
                self.updateSubsystem()
            else:
                self.ss.lineList[0].setL(self.L)
            
            '''
            for i in range(self.nLines):     # go through each section

                # update the section length from the design variables
                dvi = self.Xindices[i] #design variable index
                
                # only update if design variable is in list (not constant)
                if dvi in range(self.nX):
                    L=X[dvi]
                    self.ss.lineList[i].setL(L)

                # if the line has buoyancy, apply the buoyancy design variable
                if i in self.i_buoy:
                    
                    # check if design variable is in list (not constant)
                    dvi = self.Xindices[self.nLines + self.i_buoy.index(i)]
                    if dvi in range(self.nX):
                        B = X[dvi]                   # take the buoyancy per unit length design variable [N/m]
            
                    #handle cases where buoyancy is fixed
                    else:
                        B = self.allVars[self.nLines + i - 1]
    
                    # compute what diameter of buoyancy module is needed to achieve this buoyancy per unit length
                    d_inner = self.d0  # inner diameter for buoyancy module [m]
                    rho_buoy = self.rho_buoy                                    # constant density of buoyancy modules
                    
                    d_outer = np.sqrt(((4*B)/((self.rho-self.rho_buoy)*np.pi*self.g))+d_inner**2) # required outer diameter of buoyancy modules (assuming spread over section length) [m]
                    m_buoy = rho_buoy*(np.pi/4*(d_outer**2 - d_inner**2))       # mass per meter of spread buoyancy module [kg/m]
                    m = self.m0 + m_buoy                                        # mass per unit length of combined cable + spread buoyancy modules [kg/m]
                    w = m*self.g - self.rho*(np.pi/4*(d_outer**2))*self.g       # weight per unit length [N/m]
                    
                    # update line properties
                    self.ss.lineTypes[i]['m'] = m
                    self.ss.lineTypes[i]['w'] = w
                    self.ss.lineTypes[i]['d_vol'] = d_outer
            '''
            
            
            # ----- evaluate constraints -----
            # Evaluate any constraints in the list, at the appropriate displacements.
            # The following function calls will fill in the self.convals array.
            
            #increase solveEquilibrium tolerance
            self.ss.eqtol = 0.05
            
            
            if self.MG:
                self.addMarineGrowth(self.mgdict)
                self.ss_mod.eqtol = 0.05
                self.ss_mod.maxIter = 5000
                
         
            # ZERO OFFSET: 
            self.ss.setOffset(0)
            self.ss.calcCurvature()
            
            if self.MG: 
                self.ss_mod.setOffset(0)
                self.ss_mod.calcCurvature()
            
            # Save tensions # these aren't used anywhere... not saving the MG tensions
            for i, line in enumerate(self.ss.lineList):
                self.Te0[i,0] = np.linalg.norm(line.fA)
                self.Te0[i,1] = np.linalg.norm(line.fB)
            
            # Call any constraints that evaluate at the undisplaced position
            for con in self.constraints:
                if con['offset'] == 'zero':
                    con['fun'](X)
            
            # MAX OFFSET:
            self.ss.setOffset(self.x_mean+self.x_ampl)  # apply static + dynamic offsets
            self.ss.calcCurvature()
            
            if self.MG: 
                self.ss_mod.setOffset(self.x_mean+self.x_ampl)
                self.ss_mod.calcCurvature()
            
            # Call any constraints needing a positive displacement
            for con in self.constraints:
                if con['offset'] == 'max':
                    con['fun'](X)
                    
            self.min_lay_length = self.ss.getLayLength()  # record minimum lay length
            
            if self.MG:
                self.min_lay_length = min([self.ss.getLayLength(), self.ss_mod.getLayLength()])
            
            # MIN OFSET: 
            self.ss.setOffset(-self.x_mean-self.x_ampl)  # apply static + dynamic offsets
            self.ss.calcCurvature()
            
            if self.MG: 
                self.ss_mod.setOffset(-self.x_mean-self.x_ampl)
                self.ss_mod.calcCurvature()
            
            # Call any constraints needing a negative displacement
            for con in self.constraints:
                if con['offset'] == 'min':
                    con['fun'](X)
                    
            self.max_lay_length = self.ss.getLayLength()  # record maximum lay length
            if self.MG:
                self.max_lay_length = max([self.ss.getLayLength(), self.ss_mod.getLayLength()])
            
            # OTHER:
            self.ss.setOffset(0)  # restore to zero offset and static EA
            if self.MG: 
                self.ss_mod.setOffset(0)
                
            # or at least set back to static states
            
            # Call any constraints that depend on results across offsets
            for con in self.constraints:
                if con['offset'] == 'other':
                    con['fun'](X)
            
            
            # --- evaluate objective function ---
            
            # calculate the cost of each section
            self.cost = {}
            
            if self.ignore_static:  # option to ignore static portion of cable in cost calcs
                L = self.L - self.min_lay_length
            else:
                L = self.L            
            
            self.cost['cable'] = L*self.cableType['cost']
            
            self.cost['buoyancy'] = 0
            
            if 'buoyancy_sections' in self.dd:
                for bs in self.dd['buoyancy_sections']:
                    self.cost['buoyancy'] += bs['N_modules']*self.buoyType['cost']
                
            self.cost['total'] = self.cost['cable'] + self.cost['buoyancy']
            
            self.obj_val = self.cost['total'] / self.obj_denorm  # normalize objective function value
            
            # could also add a cost for touchdown protection sleeve based on 
            # the touchdown point range of motion
            # e.g.  c_touchdwown_protection = (self.max_lay_length - self.min_lay_length) * cost_factor
            
            
            # ----- write to log -----
            
            # log the iteration number, design variables, objective, and constraints
            self.log['x'].append(list(X))
            self.log['f'].append(list([self.obj_val]))
            self.log['g'].append(list(self.convals))
            
            
            # provide some output?
            if display > 5:
                f = self.objectiveFun(X, display=1)
                print("Line lengths are ")
                for line in self.ss.lineList:
                    print(line.L)

                print(f"Cost is {f}")
                self.evaluateConstraints(X, display=1)
                self.ss.plotProfile()
                plt.show()

    
    def objectiveFun(self, X, display=0, normalized=True):
        '''Update the design (if necessary) and return the objective function
        (cost) value.'''

        self.updateDesign(X, display=display, normalized=normalized)

        if display > 1:
            print(f"Cost is {self.cost['total']:.1f} and objective value is {self.obj_val:.3f}.")

        return self.obj_val


    def evaluateConstraints(self, X, display=0, normalized=True):
        '''Update the design (if necessary) and display the constraint
        values.'''

        self.updateDesign(X, display=display, normalized=normalized)
        
        if display > 1:
            for i, con in enumerate(self.constraints):
                print(f" Constraint {i:2d} value of {con['value']:8.2f} "
                      +f"for {con['name']}: {con['threshold']} of {con['index']} at {con['offset']} displacement.")
        
        return self.convals


    def setNormalization(self):
        '''Set normalization factors for optimization 
        (based on initial design state).'''
        
        # design variables
        self.X_denorm = np.array(self.Xlast)
        # objective
        self.obj_denorm = self.cost['total']
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
        
        # reset iteration counter
        self.iter = -1
        
        # clear optimization progress tracking lists
        self.log['x'] = []
        self.log['f'] = []
        self.log['g'] = []
        
        # set combined objective+constraints function for dopt
        def eval_func(X):
            '''DynamicCable object evaluation function'''

            self.updateDesign(X)
            f = self.obj
            g = np.array(self.convals)  # needs to be a copy to not pass by ref
            oths = dict(status=1)

            return f, g, [], [], oths, False

        # set the display value to use over the entire process
        self.display = display
        self.method = method
        
        # Set starting point to normalized value
        X0 = self.X0 / self.X_denorm
        dX_last = self.dX_last / self.X_denorm
        Xmax = self.Xmax / self.X_denorm
        Xmin = self.Xmin / self.X_denorm

        # call optimizer to perform optimization
        if method=='dopt':

            if display > 0:  print("\n --- Beginning CableDesign2 optimize iterations using DOPT2 ---")

            X, min_cost, infodict = dopt2(eval_func, X0, tol=0.001, a_max=1.4, maxIter=maxIter, stepfac=stepfac,
                                          Xmin=Xmin, Xmax=Xmax, dX_last=dX_last, display=self.display)
        
        elif method in ['COBYLA', 'SLSQP']:

            from scipy.optimize import minimize

            if self.display > 0:  print("\n --- Beginning CableDesign2 optimize iterations using COBYLA ---")

            condict = [dict(type="ineq", fun=con) for con in self.conList]
            cons_tuple = tuple(condict)

            if method=='COBYLA':
                result = minimize(self.objectiveFun, X0, constraints=cons_tuple, method="COBYLA",
                                  options={'maxiter':maxIter, 'disp':True, 'rhobeg':0.1})
                                  #options={'maxiter':maxIter, 'disp':True, 'rhobeg':10.0})
            
            elif method=='SLSQP':
                result = minimize(self.objectiveFun, X0, constraints=cons_tuple, method='SLSQP',
                                  bounds = list(zip(Xmin, Xmax)), 
                                  options={'maxiter':maxIter, 'eps':0.02,'ftol':1e-6, 'disp': True, 'iprint': 99})
            
            X = result.x
        
        #elif method=='CMNGA':
        #    from cmnga import cmnga

        #    bounds = np.array([[self.Xmin[i], self.Xmax[i]] for i in range(len(self.Xmin))])

        #    X, min_cost, infoDict = cmnga(self.objectiveFun, bounds, self.conList, dc=0.2, nIndivs=12, nRetry=100, maxGens=20, maxNindivs=500 )

            #, maxIter=maxIter, stepfac=stepfac, Xmin=self.Xmin, Xmax=self.Xmax, dX_last=self.dX_last, display=self.display)

        else:
            raise Exception('Optimization method unsupported.')
        
        # make sure it's left at the optimized state
        self.updateDesign(X)
            
        # plot
        if plot:
            self.plotOptimization()

        return X, self.cost['total'] # , infodict
    
    
    def plotOptimization(self):
    
        if len(self.log['x']) == 0:
            print("No optimization trajectory saved (log is empty). Nothing to plot.")
            return
            
        fig, ax = plt.subplots(len(self.X0)+1+len(self.constraints),1, sharex=True, figsize=[6,8])
        fig.subplots_adjust(left=0.4)
        Xs = np.array(self.log['x'])
        Fs = np.array(self.log['f'])
        Gs = np.array(self.log['g'])
        
        for i in range(len(self.X0)):
            ax[i].plot(Xs[:,i])
            #ax[i].axhline(self.Xmin[i], color=[0.5,0.5,0.5], dashes=[1,1])
            #ax[i].axhline(self.Xmax[i], color=[0.5,0.5,0.5], dashes=[1,1])

        ax[len(self.X0)].plot(Fs)
        ax[len(self.X0)].set_ylabel("cost", rotation='horizontal')

        for i, con in enumerate(self.constraints):
            j = i+1+len(self.X0)
            ax[j].axhline(0, color=[0.5,0.5,0.5])
            ax[j].plot(Gs[:,i])
            ax[j].set_ylabel(f"{con['name']}({con['threshold']})", 
                           rotation='horizontal', labelpad=80)

        ax[j].set_xlabel("function evaluations")

    # ::::::::::::::::::::::::::::::: constraint functions :::::::::::::::::::::::::::::::

    # Each should return a scalar C where C >= 0 means valid and C < 0 means violated.


    def con_total_length(self, X, index, threshold):
        '''This ensures that the total length of the Mooring does not result in a fully slack Mooring
        (ProfileType=4) in its negative extreme mean position'''
        # ['max_line_length', index, threshold] # index and threshold are completely arbitrary right now

        Lmax = (self.span-self.rBFair[0]-self.x_mean + self.depth+self.ss.rBFair[2]) # (3-14-23) this method might now be deprecated with more recent updates to ensure the combined line lengths aren't too large
        total_linelength =  sum([self.ss.lineList[i].L for i in range(self.nLines)])
        c = Lmax-total_linelength

        return c
        
        
    def con_lay_length(self, X, index, threshold, display=0):
        '''This ensures there is a minimum amount of line on the seabed at the +extreme displaced position.'''
        
        if self.MG:
            minlaylength = min([self.ss.getLayLength(iLine=index),self.ss_mod.getLayLength(iLine=index)])
        else:
            minlaylength = self.ss.getLayLength(iLine=index)
            
        return  minlaylength - threshold  + self.LayLen_adj
    
    def con_max_lay_length(self, X, index, threshold, display=0):
        '''This ensures there is a minimum amount of line on the seabed at the +extreme displaced position.'''
        
        if self.MG:
            minlaylength = min([self.ss.getLayLength(iLine=index),self.ss_mod.getLayLength(iLine=index)])
        else:
            minlaylength = self.ss.getLayLength(iLine=index)
            
        return  threshold - minlaylength 

    def con_max_td_range(self, X, index, threshold, display=0):
        '''Ensures the range of motion of the touchdown point betweeen the
        range of offsets is less then a certain distance.
        This constraint is for the system as a whole (index is ignored) and 
        must have offset='other' so that it's evaluated at the end.'''
        return threshold - (self.max_lay_length - self.min_lay_length)

    """
    def con_buoy_contact(self, X, index, threshold, display=0):
        '''This ensures the first line node doesn't touch the seabed by some minimum clearance in the +extreme displaced position.'''
        return self.getPointHeight(index) - threshold
        <<<< seems funny <<<
    """

    def con_strength(self, X, index, threshold, display=0):
        '''This ensures the MBL of the line is always greater than the maximum 
        tension the line feels times a safety factor.'''
        if self.MG:
            minsf = min([self.ss.getTenSF(index),self.ss_mod.getTenSF(index)])
        else:
            minsf = self.ss.getTenSF(index)
        return minsf - threshold
    
    def con_overall_strength(self, X, index, threshold, display=0):
        '''This ensures the MBL of the line is always greater than the maximum 
        tension the line feels times a safety factor. *** checks all line sections ***'''
        
        sfs = []

        #check both ss_mod and ss if there's marine growth
        if self.MG:
            
            #iterate through linelist and append safety factors
            for index in range(0, len(self.ss_mod.lineList)):    
                minsf = self.ss_mod.getTenSF(index)
                sfs.append(minsf - threshold)
        
        for index in range(0, len(self.ss.lineList)):    
            minsf = self.ss.getTenSF(index)
            sfs.append(minsf - threshold)
                
        return min(sfs)
    
    
    
    def con_curvature(self, X, index, threshold, display=0):
        '''Ensure that the MBR of the cable is always greater than the maximum 
        actual curvature times a safety factor.'''
        if self.MG:
            mincsf = min([ self.ss.getCurvSF(index),  self.ss_mod.getCurvSF(index)])
        else: 
            mincsf = self.ss.getCurvSF(index)
        return mincsf - threshold
    
    def con_max_curvature(self, x, index, threshold, display=0):
        '''Ensures that the MBR divided by the maximum curvature over the 
        entire cable is greater than a threshold safety factor.
        
        >>> make a single set of cable props for the line overall
        >>> then there will be more for the buoyancy sections '''
        if self.MG:
            maxks = max([max(self.ss.Ks), max(self.ss_mod.Ks)])
        else:
            maxks = max(self.ss.Ks)
        return 1 /( self.cableType['MBR'] * maxks ) - threshold
    

    def con_min_sag(self, X, index, threshold, display=0):
        '''Ensure the lowest point of a line section is below 
        a minimum depth.'''
        if self.MG:
            minsag = min([self.ss.getSag(index), self.ss_mod.getSag(index)])
        else:
            minsag = self.ss.getSag(index)
        return threshold - minsag
    
    def con_max_sag(self, X, index, threshold, display=0):
        '''Ensures the lowest point of a line section is above
        a certain maximum depth.'''
        if self.MG:
            maxsag = max([self.ss.getSag(index), self.ss_mod.getSag(index)])
        else:
            maxsag = self.ss.getSag(index)
        return maxsag - threshold
    
    def con_max_hog(self, X, index, threshold, display=0):
        '''Ensures the highest point of a line section is below 
        a certain maximum depth '''
        if self.MG:
            maxhog = max([self.ss.getHog(index), self.ss_mod.getHog(index)])
        else:
            maxhog = self.ss.getHog(index)
        return threshold - maxhog
    
    # ----- utility functions -----

    def plotProfile(self, iPoint=1, Xuvec=[1,0,0], Yuvec=[0,0,1], ax=None, color=None, title="", slack=False, displaced=True, figsize=(6,4)):
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
        displaced : bool, optional
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
            
            if self.MG:
                fig1, ax1 = plt.subplots(1,1, figsize=figsize)
                ax1.set_xlabel('Horizontal distance (m)')
                ax1.set_ylabel('Depth (m)')
        else:
            fig = plt.gcf()   # will this work like this? <<<

        
        if displaced:
            offsets = [0, self.x_mean+self.x_ampl, -self.x_mean-self.x_ampl]
        else:
            offsets = [0]

        for x in offsets:
            alph = 1 if x==0 else 0.5  # make semi-transparent for offset profiles
            
            self.ss.setOffset(x)

            ax.plot(x, 0,'ko',markersize = 2)  # draw platform reference point

            if self.shared == 2:  # plot other half too if it's a shared line where only half is modeled <<<
                for i, line in enumerate(self.ss.lineList):
                    if i in self.i_buoy:
                        self.ss.lineList[i].color = [.6,.6,.0]
                    else:
                        self.ss.lineList[i].color = [.3,.5,.5]
                self.ss.drawLine2d(0, ax, color = 'self', Xoff = -self.ss.span/2)
                
                #store ss cos_th before plotting the flipped half cable
                self.ss.cos_th = -self.ss.cos_th
                self.ss.drawLine2d(0, ax, color = 'self', Xoff = self.ss.span/2)
                self.ss.cos_th = -self.ss.cos_th
             
            else:
                for i, line in enumerate(self.ss.lineList):
                    if color==None:     # alternate colors so the segments are visible
                        if i in self.i_buoy:
                            line.drawLine2d(0, ax, color=[.6,.6,.0], alpha=alph, Xuvec=Xuvec, Yuvec=Yuvec)
                        else:
                            line.drawLine2d(0, ax, color=[.3,.5,.5], alpha=alph, Xuvec=Xuvec, Yuvec=Yuvec)
                    else:
                            line.drawLine2d(0, ax, color=color, alpha=alph, Xuvec=Xuvec, Yuvec=Yuvec)
              
            if self.MG:
                alph = 1 if x==0 else 0.5  # make semi-transparent for offset profiles
                
                self.ss_mod.setOffset(x)

                ax1.plot(x, 0,'ko',markersize = 2)  # draw platform reference point

                if self.shared == 2:  # plot other half too if it's a shared line where only half is modeled <<<
                    for i, line in enumerate(self.ss_mod.lineList):
                        
                        #check if linetype name has buoy in it (**** this is highly dependent on naming convention)
                        if line.type['name'].split("_")[-1][:4] == 'buoy':
                            self.ss_mod.lineList[i].color = [.6,.6,.0]
                        else:
                            self.ss_mod.lineList[i].color = [.3,.5,.5]
                    self.ss_mod.drawLine2d(0, ax1, color = 'self', Xoff = -self.ss.span/2)
                    
                    #store ss cos_th before plotting the flipped half cable
                    self.ss_mod.cos_th = -self.ss_mod.cos_th
                    self.ss_mod.drawLine2d(0, ax1, color = 'self', Xoff = self.ss_mod.span/2)
                    self.ss_mod.cos_th = -self.ss_mod.cos_th
                 
                else:
                    for i, line in enumerate(self.ss_mod.lineList):
                        if color==None:     # alternate colors so the segments are visible
                        
                            #check if linetype name has buoy in it (**** this is highly dependent on naming convention)
                            if line.type['name'].split("_")[-1][:4] == 'buoy':
                                line.drawLine2d(0, ax1, color=[.6,.6,.0], alpha=alph, Xuvec=Xuvec, Yuvec=Yuvec)
                            else:
                                line.drawLine2d(0, ax1, color=[.3,.5,.5], alpha=alph, Xuvec=Xuvec, Yuvec=Yuvec)
                        else:
                                line.drawLine2d(0, ax1, color=color, alpha=alph, Xuvec=Xuvec, Yuvec=Yuvec)
            
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


        self.ss.setOffset(0)  # set back to its neutral position
        
        if self.MG:
            
            if not slack:
                ax1.axis("equal")

            ax1.set_title(title + " Marine Growth")
            
            self.ss_mod.setOffset(0)

            return fig, ax, fig1, ax1  # return the figure and axis object in case it will be used later to update the plot
        
        else:
            return fig, ax

    def plotCurves(self, ax=[], color="k", title=""):
        '''Plot key performance curves for the cable as a function of offset

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
        
        x = np.linspace(-self.x_mean_high-self.x_ampl, self.x_mean_low+self.x_ampl, 50)
        
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

    """
    def makeCableType(self, di, name):
        '''sets up a cableType dictinoary by reading in from a dictionary a
        a specified name entry.'''
        
        # a few calculations
        d = float(di['d'])  # [m]
        m = float(di['m'])  # [kg/m]
        w = (m - np.pi/4*d**2 *self.rho)*self.g

        # make and fill in a cableType dictionary, which will go in MoorPy's lineTypes dictionary
        cableType = dict(name=name)
        cableType['d_vol'] = float(di['d'])  # [m]
        cableType['m']     = m  # [kg/m
        cableType['w']     = w  # [N/m] wet weight per unit length]
        cableType['EA']    = getFromDict(di, 'EA')  # [N] axial stiffness
        cableType['EI']    = getFromDict(di, 'EI' , default=0)  # [N-m^2] bending stiffness
        cableType['MBL']   = getFromDict(di, 'MBL', default=0)  # [N] minimum breaking load
        cableType['MBR']   = getFromDict(di, 'MBR', default=0)  # [m] minimum bend radius
        cableType['A_con'] = getFromDict(di, 'A'  , default=0)  # [mm^2] conductor area
        cableType['dynamic'] = getFromDict(di, 'dynamic', dtype=bool, default=True)
        cableType['DC']      = getFromDict(di, 'DC'     , dtype=bool, default=False)
        cableType['cable_cost'] = getFromDict(di, 'cable_cost', default=0)   # $/m dynamic cable cost 
        cableType['buoy_cost']  =  getFromDict(di, 'buoy_cost', default=0)   # cost of each module
        cableType['buoy_length'] = getFromDict(di, 'buoy_length', default=0) # meters for each buoyancy module
        cableType['L_BM'] = getFromDict(di, 'L_BM', default=0)               # [m] center to center spacing between two buoyancy modules
        cableType['D_BM'] = getFromDict(di, 'D_BM', default=0)               # [m] Diameter of buoyancy module
        cableType['V_BM'] = getFromDict(di, 'V_BM', default=0)               # [m] volume of buoyancy module
        cableType['rho_buoy'] = getFromDict(di, 'rho_buoy', default=500)     # [kg/m^3] aggregate density of buoyancy module
        if cableType['V_BM'] <= 0:
            raise Exception("Volume of buoyancy module must be greater than zero")
        
        return cableType
    """
    def updateHyroCoeffs(self, C_dnc = 1.2, C_dnb = 1.2, C_dab1 = 1, C_dab2 = 0, C_dac = 0, C_anb = 1, C_anc = 1, C_aab = 0.5 , C_aac = 0):
        '''
        

        Parameters
        ----------
        C_dnc : Normal drag coeff for the cable. The default is 1.2.
        C_dnb : Normal drag coeff for the buoyancy module. The default is 1.2.
        C_dab1 : Drag coefficient for exposed ends of buoyancy module. The default is 1.
        C_dab2 : Axial drag coefficient for buoyancy module (skin friction). The default is 0.
        C_dac : Axial drag coefficient for cable (skin friction). The default is 0.
        C_anb : Normal added mass coefficient for buoyancy module. The default is 1.
        C_anc : Normal added mass coefficient for cable. The default is 1.
        C_aab : Axial added mass coefficient for buoyancy module. The default is 0.5 (assumed sphere added mass coeff).
        C_aac : Axial added mass coefficient for cable. The default is 0.

        Returns
        -------
        None.

        '''
        #iterate through list of line properties
        buoycount = -1
        for i in (range(0, len(self.ss.lineTypes))):
            linetype = self.ss.lineTypes[i]
            if linetype['name'].split("_")[-1][:4] == 'buoy':
                buoycount += 1 
                deq = linetype['d_vol'] # volume equiv diameter for buoy section
                dc = self.cableType['d_vol'] # diameter of cable
                db = self.buoyType['d'] # diameter of buoy
                Lbs = self.bs[buoycount]
                if Lbs == 0:
                    ValueError('Buoyancy module spacing is zero')
                Lb = self.buoyType['l']
                
                self.ss.lineTypes[i]['Cd'] = 1/(Lbs * deq)*(C_dnc * dc * (Lbs - Lb) + C_dnb * db * Lb)
                self.ss.lineTypes[i]['CdAx'] = 1 / (Lbs * deq) *(C_dab1 * (db**2 - dc**2)/4 + C_dab2 * db * Lb + C_dac * dc * (Lbs - Lb))
                self.ss.lineTypes[i]['Ca'] = 1 / (Lbs * deq**2) * (C_anb * db**2 * Lb + C_anc * dc**2 *(Lbs - Lb))
                self.ss.lineTypes[i]['CaAx'] =  1 / (Lbs * deq**2) * (C_aab * db**2 * Lb + C_aac * dc**2 *(Lbs - Lb))
            else:
                self.ss.lineTypes[i]['CdAx'] = 0.0
                self.ss.lineTypes[i]['Ca'] = 1.0
        if self.MG:
            for i in (range(0, len(self.ss_mod.lineTypes))):
                linetype = self.ss_mod.lineTypes[i]
                if linetype['name'].split("_")[-1][:4] == 'buoy':
                    buoycount += 1 
                    deq = linetype['d_vol'] # volume equiv diameter for buoy section
                    dc = self.cableType['d_vol'] # diameter of cable
                    db = self.buoyType['d'] # diameter of buoy
                    Lbs = self.bs[buoycount]
                    if Lbs == 0:
                        ValueError('Buoyancy module spacing is zero')
                    Lb = self.buoyType['l']
                    
                    self.ss_mod.lineTypes[i]['Cd'] = 1/(Lbs * deq)*(C_dnc * dc * (Lbs - Lb) + C_dnb * db * Lb)
                    self.ss_mod.lineTypes[i]['CdAx'] = 1 / (Lbs * deq) *(C_dab1 * (db**2 - dc**2)/4 + C_dab2 * db * Lb + C_dac * dc * (Lbs - Lb))
                    self.ss_mod.lineTypes[i]['Ca'] = 1 / (Lbs * deq**2) * (C_anb * db**2 * Lb + C_anc * dc**2 *(Lbs - Lb))
                    self.ss_mod.lineTypes[i]['CaAx'] =  1 / (Lbs * deq**2) * (C_aab * db**2 * Lb + C_aac * dc**2 *(Lbs - Lb))
                else:
                    self.ss_mod.lineTypes[i]['CdAx'] = 0.0
                    self.ss_mod.lineTypes[i]['Ca'] = 1.0
         

# ----- Main Script -----
if __name__ == '__main__':

    # EXAMPLE

    depth = 800
    configuration = 'Humboldt'

    settings = {}
    settings['rBFair'] = [0,0,-14]  # relative attachment coordinate on FOWT [m]
    settings['span'] = 950  # relative attachment coordinate on FOWT [m]

    settings['offset'] = 80  # mean surge offsets in either direction [m]
    settings['x_ampl'] = 5   # additional dynamic surge amplitude about the mean [m]

    
    # design variables: initial values, min and max bounds
    settings['Xindices'] = ['c', 0,      1,     2,  'c'] # order of design variables. multiple line segments can have the same design variable. 'c' flag means that it stays constant
    #                      span  L     B1[kN]   Lmid1  Spread
    settings['allVars'] = [950, 1100, 100, 613, 300]  # must be the same length as Xindices
    settings['Xmin'] =    [100, 100, 100] # must be same length as # of design variables
    settings['Xmax'] =    [1200, 800, 1000] # must be same length as # of design variables
    settings['dX_last'] = [10, 10, 10] # must be same length as # of design variables

    # set up constraints
    settings['constraints'] = [dict(name='min_lay_length', index=0, threshold=     80, offset='max'),  # ensure there is at least 20 m of cable along the seabed
                               dict(name='max_sag',        index=1, threshold=5-depth, offset='min')]  # ensure the start of the buoyancy section stays 5 m off the seabed

    # also add a tension safety factor constraint for each section

    for i in range(3):
        settings['constraints'].append(dict(name='tension_safety_factor', index=i, threshold=2.0, offset='max'))

    # add a curvature safety factor constraint for each offset of the cable or section of the cable
    for i in range(3):
        settings['constraints'].append(dict(name='curvature_safety_factor', index=i, threshold=2.0, offset='min'))

    # add a maximum touchdown point range of motion constraint
    settings['constraints'].append(dict(name='max_touchdown_range', index=0, threshold=50.0, offset='other'))
    
    # load property coefficients
    cable_props = cprops.loadCableProps(None) # load default property scaling coefficients
    cableType = cprops.getCableProps(400, 'dynamic_cable_66', cableProps=cable_props)
    
    buoy_props = cprops.loadBuoyProps(None) # load default property scaling coefficients
    buoyType = cprops.getBuoyProps(1, 'Buoyancy_750m', buoyProps=buoy_props)
    
    
    #set up the object
    dc = CableDesign(depth, cableType, buoyType, n=3, i_buoy=[1], **settings)

    #plot initial design
    dc.plotProfile(title='initial (X0)')
    dc.setNormalization()
    
    X, min_cost = dc.optimize(maxIter=3, plot=False, display=2, stepfac=4, method='COBYLA')
    #X, min_cost = dc.optimize(maxIter=8, plot=False, display=1, stepfac=4, method='SLSQP')
    #X, min_cost = dc.optimize(maxIter=2, plot=False, display=1, stepfac=4, method='dopt')

    dc.objectiveFun(X, display=2)
    dc.evaluateConstraints(X, display=2)
    dc.updateDesign(X, display=0)
    dc.plotProfile(title= 'dopt')
    dc.plotOptimization()
    #dc.unload('Humboldt.dat')
    
    
    plt.show()
