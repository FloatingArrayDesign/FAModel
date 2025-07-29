# class for a connector

import numpy as np
from famodel.famodel_base import Node, Edge
from moorpy.helpers import getFromDict, getPointProps, loadPointProps

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
        
        # cost dictionary
        self.cost = {}
        
        self.getProps()
    
    
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
        ms.addPoint(0,self.r)
        
        # assign this point as mpConn in the anchor class instance
        self.mpConn = ms.pointList[-1]

        self.mpConn.m = self['m']
        self.mpConn.v = self['v']
        self.mpConn.CdA = self['CdA']
        
        # set point type in ms
        self.getProps()
        
        return(ms)
    
    
    def getProps(self):
        '''
        Wrapper function to get moorpy point props dictionary
        and set the point type in the moorpy system
        (if it exists)'''
        # get point type information if possible
        pt = {}
        if 'type' in self:
            try:
                design = {f"num_c_{self['type']}":1}
                details = {}
                if self['m']>0:
                    details['m'] = self['m']
                if self['v']>0:
                    details['v'] = self['v']
                if self['CdA']>0:
                    details['CdA'] = self['CdA']
                if self.mpConn:
                    pt = self.mpConn.sys.setPointType(design,**details)
                else:
                    props = loadPointProps(None)
                    pt = getPointProps(design, Props=props, **details)
                self.required_safety_factor = pt['FOS']
            except:
                pass
            
        return(pt)
    
    
    def getCost(self,update=True, fx=0.0, fz=0.0, peak_tension=None, MBL=None):
        '''Get cost of the connector from MoorPy pointProps.
        Wrapper for moorpy's getCost_and_MBL helper function'''
        if update:
            if self.mpConn:
                # use pointProps to update cost
                try:
                    self.getProps()
                    self.cost['materials'], MBL, info = self.mpConn.getCost_and_MBL(fx=fx, fz=fz, peak_tension=peak_tension)
                except:
                    print('Warning: unable to find cost from MoorPy pointProps, cost dictionary not updated')
        # if update == False, just return existing costs
                
        return sum(self.cost.values())


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
        
        # MoorPy Line object for the section
        self.mpLine = None

    
    def makeMoorPyLine(self, ms):
        '''Create a MoorPy Line object in a MoorPy system.
        If this section is attached to connectors that already have associated
        MoorPy point objects, then those attachments will also be made in 
        MoorPy.
        
        Parameters
        ----------
        ms : MoorPy System object
            The MoorPy system to create the Line object in.
        '''
        
        # See if this section is attached to any already-created MoorPy Points
        pointA = 0
        if self.attached_to[0]:  # if an end A attachment
            if self.attached_to[0].mpConn:  # if it has a MoorPy point object
                pointA = self.attached_to[0].mpConn.number  # get its number
        pointB = 0
        if self.attached_to[1]:
            if self.attached_to[1].mpConn:
                pointB = self.attached_to[1].mpConn.number
        
        # Create a Line for the section in MoorPy system
        ms.addLine(self['L'], self['type'], pointA=pointA, pointB=pointB)
        
