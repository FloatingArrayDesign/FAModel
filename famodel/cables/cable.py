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
    
    def __init__(self, id, d=None):
        ''' '
        d : dict
            Dictionary of cable information (see ontology).
        '''
        
        Edge.__init__(self, id)  # initialize Edge base class
        
        self.dd = {'joints':[],'cables':[]}  # save design dictionary
        self.n_sec = len(d['cables'])
        
        # Turn what's in dd and turn it into Sections and Connectors (if there is more than one section)
        if len(d['cables'])>1:
            for i, joi in enumerate(d['joints']):
                if joi:
                    Jid = id+'_'+d['joints'][i]['type']+str(i)
                else:
                    Jid = id+'_'+str(i)
                self.dd['joints'].append(Joint(Jid,**d['joints'][i]))
            
            for i, sec in enumerate(d['cables']):
                Cid = id+'_'+sec['cable_type']['name']+str(i)
                if sec['type'] == 'static':
                    self.dd['cables'].append(StaticCable(Cid, dd=d['cables'][i], **d['cables'][i]))
                else:
                    self.dd['cables'].append(DynamicCable(Cid, dd=d['cables'][i], **d['cables'][i]))
            
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
        
        else:
            # just create the singular cable object as a dynamic cable
            self.dd = d
            Cid = id+'_'+self.dd['cables'][0]['cable_type']['name']+str(0)
            cabObj = DynamicCable(Cid, dd=self.dd['cables'][0],**self.dd['cables'][0])
            self.dd['cables'][0]['object'] = cabObj
        
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
        
    