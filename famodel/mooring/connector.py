# class for a connector

import numpy as np

class Connector():
    '''
    '''
    def __init__(self, dd=None, r=[0,0,0]):
        '''

        Parameters
        ----------
        dd : dictionary, optional
            Design dictionary The default is None.
        r : list, optional
            x,y,z location of the connector. The default is [0,0,0].

        Returns
        -------
        None.

        '''

        # Design description dictionary for this Platform
        self.dd = dd
        
        # Connector position and orientation
        self.r = np.array(r)  # x, y, z coordinates of connector [m]
        
        if self.dd:
            self.m = self.dd['mass']
            self.v = self.dd['volume']
            if 'CdA' in self.dd:
                self.CdA = self.dd['CdA']
            else:
                self.CdA = None
        else:
            self.CdA = None
            self.m = None
            self.v = None
        
        # MoorPy Point Object for Connector
        self.mpConn = None
        
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

        # add mass if available
        if self.dd['design']['mass']:
            self.mpConn.m = self.dd['design']['mass']
        # set connector volume and axial Cd if available
        if self.dd['design']['volume']:
            self.mpConn.v = self.dd['design']['volume']
        if self.dd['design']['CdA']:
            self.mpConn.CdA = self.dd['design']['CdA']

        return(ms)
        