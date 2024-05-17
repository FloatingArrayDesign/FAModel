# class for a power cable in a floating array

import numpy as np
from copy import deepcopy
from moorpy.subsystem import Subsystem
from moorpy import helpers

from famodel.cables.dynamic_cable import DynamicCable
from famodel.cables.static_cable import StaticCable
from famodel.cables.components import Joint
from famodel.famodel_base import Edge


class SubseaCable(Edge):
    '''Class for an entire subsea cable that transmits power between two 
    turbines or between a turbine and a substation. This can include a
    dynamic cable section (at each end), a static cable section, and
    subsea joints that connect them. Possible scenarios:
    DynamicCable - Joint - StaticCable - Joint - DynamicCable
    DynamicCable - Joint - StaticCable - Joint
    DynamicCable - Joint
                   Joint - StaticCable - Joint - DynamicCable
                                         Joint - DynamicCable
    DynamicCable (fully suspended)
    '''
    
    def __init__(self, id, dd=None):
        ''' '
        dd : dict
            Dictionary of cable information (see ontology).
        '''
        
        Edge.__init__(self, id)  # initialize Edge base class
        
        self.dd = dd  # save design dictionary
        
        self.n_sec = len(self.dd['cables'])
        
        # Turn what's in dd and turn it into Sections and Connectors
        for i, joi in enumerate(self.dd['joints']):
            if joi:
                Jid = id+'_'+joi+str(i)
            else:
                Jid = id+'_'+i
            self.dd['joints'][i] = Joint(Jid,**self.dd['joints'][i])
        
        for i, sec in enumerate(self.dd['cables']):
            Cid = id+'_'+sec+str(i)
            if sec['type'] == 'static':
                self.dd['cables'][i] = StaticCable(Cid, **self.dd['cables'][i])
            else:
                self.dd['cables'][i] = DynamicCable(Cid, **self.dd['cables'][i])
        
        # Connect them and store them in self(Edge).subcomponents!
        subcons = []  # list of node-edge-node... to pass to the function
        for i in range(self.n_sec):
            subcons.append(self.dd['joints'][i])
            subcons.append(self.dd['cables'][i])
        subcons.append(self.dd['joints'][-1])
        self.addSubcomponents(subcons)  # Edge method to connect and store em
        
        # Indices of connectors and sections in self.subcomponents list
        self.i_con = list(range(0, 2*self.n_sec+1, 2))
        self.i_sec = list(range(1, 2*self.n_sec+1, 2))
        
        '''
        self.system = system
        
        self.AttachA   =d['AttachA']
        self.AttachB   =d['AttachB']
        self.DynCableA =self.system.dynamicCableConfigs[d['DynCableA']]
        self.DynCableB =self.system.dynamicCableConfigs[d['DynCableB']]
        self.headingA  =np.radians(d['headingA'])
        self.headingB  =np.radians(d['headingB'])
        self.cableType =d['cableType']
        '''
        
        # static cable end points
        self.rA = np.zeros(3)
        self.rB = np.zeros(3)
        
        # cable routing info
        self.x = []
        self.y = []
        self.r = []
        
        self.L = 0  # total length (to be computed) [m]
        
    