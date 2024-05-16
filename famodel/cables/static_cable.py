# class for a static subsea power cable

import numpy as np
from copy import deepcopy
from moorpy.subsystem import Subsystem
from moorpy import helpers
from famodel.mooring.connector import Connector, Section
from famodel.famodel_base import Edge

class StaticCable(Edge):
    '''
    Class for a static power cable. It inherits from Cable(Edge, dict)
    which describes the bare uniform cable before accessories are added.
    A StaticCable object will likely be within a SubseaCable object, which could
    also include a DynamicCable.
    Both ends will be attached to subsea joints.
    The main thing done by this class is describing the cable routing on the
    seabed.
    '''
    
    def __init__(self, dd=None, subsystem=None, anchor=None, rA=[0,0,0], rB=[0,0,0],
                 rad_anch=500, rad_fair=58, z_anch=-100, z_fair=-14, 
                 rho=1025, g=9.81, id=None):
        '''
        Parameters
        ----------
 
        '''
        Edge.__init__(self, id)  # initialize Edge base class
        # Design description dictionary for this dynamic cable
        self.dd = dd
        
        self.n_sec = 1
        
        # Store the cable type properties dict here for easy access (temporary - may be an inconsistent coding choice)
        self.cableType = self.makeCableType(self.dd['cable_type'])  # Process/check it into a new dict

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
        
        self.shared = False # boolean for if the mooring line is a fully suspended cable
        self.symmetric = False # boolean for if the mooring line is a symmetric suspended cable
        
        # relevant site info
        self.rho = rho
        self.g = g
        
        # Dictionaries for addition information
        self.loads = {}
        self.reliability = {}
        self.cost = {}
    
    
    def setSectionLength(self, L, i):
        '''Sets length of section, including in the subdsystem if there is
        one.'''
        
        # >>> PROBABLY NEED TO REVISE HOW THIS WORKS TO CONSIDER BUOYANCY MODULES <<<
        
        self.dd['sections'][i]['L'] = L  # set length in dd (which is also Section/subcomponent)
        
        if self.ss:  # is Subsystem exists, adjust length there too
            self.ss.lineList[i].setL(L)
    
    
    
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
        elif end in ['b', 'B', 1]:
            self.rB = np.array(r)
        else:
            raise Exception('End A or B must be specified with either the letter, 0/1, or False/True.')
    
    
    def getCost(self):
        
        # >>> UPDATE <<<
        
        # sum up the costs in the dictionary and return
        return sum(self.cost.values()) 
    
    # maybe cable routing methods would go here (currently drafted in cable.py)
