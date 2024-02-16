# class for a mooring line

import numpy as np


class Mooring():
    '''
    Class for a floating array mooring line (anchored or shared).
    The idea is to have a common object that can be used for 2D
    layout optimization as well as (optionally) 3D design.
    Work in progress. Eventually will inherit from Edge.
    '''
    
    def __init__(self, dd=None, subsystem=None, anchor=None, rA=[0,0,0], rB=[0,0,0],
                 rad_anch=500, rad_fair=58, z_anch=-100, z_fair=-14):
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
                                  length
                                 }
                         }
                     connectors: {}
                     rAnchor
                     zAnchor
                     rFair
                     zFair
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
        
        # Design description dictionary for this Mooring
        self.dd = dd
        
        # MoorPy subsystem that corresponds to the mooring line
        self.subsystem = subsystem
        
        # end point absolute coordinates, to be set later
        self.rA = rA
        self.rB = rB
        self.heading = 0
        
        # relative positions (variables could be renamed)
        self.rad_anch = rad_anch
        self.rad_fair = rad_fair
        self.z_anch   = z_anch  
        self.z_fair   = z_fair
        
        self.adjuster = None  # custom function that can adjust the mooring
        
        # Dictionaries for addition information
        self.loads = {}
        self.reliability = {}
        self.cost = {}
    
    
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
            
            if self.subsystem:
                self.subsystem.setEndPosition(self.rA, False, sink=sink)
            
        elif end in ['b', 'B', 1]:
            self.rB = np.array(r)
            
            if self.subsystem:
                self.subsystem.setEndPosition(self.rB, True, sink=sink)
                
        else:
            raise Exception('End A or B must be specified with either the letter, 0/1, or False/True.')
    
    
    def getCost(self):
        
        # mooring line cost
        line_cost = 0
        if self.subsystem:
            for line in self.subsystem.lineList:
                line_cost += line.getCost()
        else:
            pass # could have some proxy cost calculation
        
        self.cost['line'] = line_cost
        
        # anchor cost  (it should already be computed)
        
        # sum up the costs in the dictionary and return
        return sum(self.cost.values()) 
        
    def createSubsystem(self, case=0):
        from moorpy.subsystem import Subsystem
        ''' Create a subsystem for a line configuration from the design dictionary
        
        Parameters
        ----------
        case : int
            Selector shared/suspended cases:
                - 0 (default): end A is on the seabed
                - 1: assembly is suspended and end A is at another floating system
                - 2: the assembly is suspended and assumed symmetric, end A is the midpoint
        '''
        # check if a subsystem already exists
        if self.subsystem:
            print('A subsystem for this Mooring class instance already exists, this will be overwritten.')
        # create subsystem
        self.subsystem=Subsystem(depth=-self.dd['zAnchor'], spacing=self.dd['rAnchor'], rBfair=self.rB)
        lengths = []
        types = []
        # run through each line section and collect the length and type
        for i in range(0,len(self.dd['sections'])):
            lengths.append(self.dd['sections'][i]['length'])
            types.append(self.dd['sections'][i]['type']['name'])
            self.subsystem.lineTypes[types[-1]] = self.dd['sections'][i]['type']

        
        # make the lines and set the points 
        self.subsystem.makeGeneric(lengths,types,suspended=case)
        self.subsystem.setEndPosition(self.rA,endB=0)
        self.subsystem.setEndPosition(self.rB,endB=1)
        return(self.subsystem)
        # solve the system
        self.subsystem.staticSolve()
