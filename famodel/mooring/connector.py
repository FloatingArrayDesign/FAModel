# class for a connector

import numpy as np
from famodel.famodel_base import Node, Edge
from moorpy.helpers import getFromDict

class Connector(Node, dict):
    '''
    '''
    def __init__(self,id, r=[0,0,0], **kwargs):
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
            Additional optional parameters, such as m, v, CdA...
        '''

        dict.__init__(self, **kwargs)  # initialize dict base class (will put kwargs into self dict)
        Node.__init__(self, id)  # initialize Node base class
        
        # <<< have a ['type'] entry that stores a type, which could be used for sizing...
        
        # Connector position and orientation
        self.r = np.array(r)  # x, y, z coordinates of connector [m]
        
        # set defaults if they weren't provided
        self['m'  ] = getFromDict(self, 'm'  , default=0)
        self['v'  ] = getFromDict(self, 'v'  , default=0)
        self['CdA'] = getFromDict(self, 'CdA', default=0)
        
        # MoorPy Point Object for Connector
        self.mpConn = None
        
        # dictionary of failure probabilities
        self.failure_probability = {}
        
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
        # create connector as a point in MoorPy system
        ms.addPoint(1,self.r)
        # assign this point as mpConn in the anchor class instance
        self.mpConn = ms.pointList[-1]

        self.mpConn.m = self['m']
        self.mpConn.v = self['v']
        self.mpConn.CdA = self['CdA']
        

        return(ms)


class Section(Edge, dict):
    '''
    '''
    def __init__(self,id, **kwargs):
        '''
        '''

        dict.__init__(self, **kwargs)  # initialize dict base class
        Edge.__init__(self, id)  # initialize Edge base class
        
        
        # <<< have a ['type'] entry that stores a type, which could be used for sizing...


        # set defaults if they weren't provided
        self['L'   ] = getFromDict(self, 'L'  , default=0)
        
        # if the type dict wasn't provided, set as none to start with
        if not 'type' in self:
            self['type'] = None