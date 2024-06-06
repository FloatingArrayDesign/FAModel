# class for cable components

import numpy as np
from famodel.famodel_base import Node, Edge


class Joint(Node, dict):
    '''Subsea joint for power cables, such as might join dynamic and static cables together.
    '''
    def __init__(self,id, r=None,**kwargs):
        '''
        Connectors inherit from dict, so properties can be passed in as arguments
        and they will be assigned like dictionary entries. 
        
        Parameters
        ----------
        dd : dictionary, optional
            Design dictionary The default is None.
        r : list, optional
            x,y,z location of the connector. The default is [0,0,0].
        kwargs
            Additional optional parameters, such as m
        '''
        from famodel.project import getFromDict
        dict.__init__(self, **kwargs)  # initialize dict base class (will put kwargs into self dict)
        Node.__init__(self, id)  # initialize Node base class
        # Joint position and orientation

        # set defaults if they weren't provided
        self['m'  ] = getFromDict(self, 'm'  , default=0)
        
        # MoorPy Point Object for Joint
        self.mpConn = None
        
        # Dictionary for failure probability
        self.failure_probability = {}
        
    #this might be useful for the ends of dynamic cables
    def makeMoorPyConnector(self, ms):
        '''Create a MoorPy connector object in a MoorPy system
        Parameters
        ----------
        ms : class instance
            MoorPy system
        
        Returns
        -------
        ms : class instance
            MoorPy system 

        '''
        from famodel.project import getFromDict
        # create connector as a point in MoorPy system
        ms.addPoint(1,self['r'])
        # assign this point as mpConn in the anchor class instance
        self.mpConn = ms.pointList[-1]

        self.mpConn.m = getFromDict(self,'m',default=10000)
        self.mpConn.v = getFromDict(self,'v',default=0)
        self.mpConn.CdA = getFromDict(self,'CdA',default=0)

        return(ms)

"""
class Cable(Edge, dict):
    '''A length of a subsea power cable product (i.e. same cross section of
    the actual cable. This is just a very simple cable description.
    Buoyancy modules and other appendages might exist along its length, described
    in other classes. (Note this is different from how Mooring Sections are
    modeled). The cable sectional specs should be stored in Cable['type'].
    '''
    def __init__(self, **kwargs):
        '''
        '''

        dict.__init__(self, **kwargs)  # initialize dict base class
        Edge.__init__(self, 'no name')  # initialize Edge base class
        
        # <<< have a ['type'] entry that stores a type, which could be used for sizing...
        
        # set defaults if they weren't provided
        self['L'   ] = getFromDict(self, 'L'  , default=0)
        
        # if the type dict wasn't provided, set as none to start with
        if not 'type' in self:
            self['type'] = None
"""