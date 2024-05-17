"""Anchor class for FAModel, containing information and key methods anchors of mooring lines
    Work in progress
"""
import moorpy as mp
import numpy as np
from famodel.famodel_base import Node

class Anchor(Node):
    
    def __init__(self, dd=None, ms=None, r=[0,0,0], aNum=None,id=None):
        '''
        Parameters
        ----------
        dd: dictionary
            Design dictionary that contains all information on an anchor for a mooring line/shared mooring
            {
                type:
                design:
                    d: # diameter
                    L: # length
                    A: # Area
                    t: # Thickness
                    h: # Embedment depth
                    m: # Mass
                    UHC: # Ultimate Holding Capacity
                soil_type:
                angle: # seabed angle 
                cost:
                    totCost: #total cost
                    matCost: # material cost
                    instCost: # installation cost
                    decomCost: # decomissioning cost                          
                }
        ms: system object
            MoorPy system object the anchor is in
            
        r: list
            Location of anchor in x,y,z
        '''
        # Initialize as a node
        Node.__init__(self,id)
        
        # Design description dictionary for this Anchor
        self.dd = dd
        
        # MoorPy system this anchor is in
        self.ms = ms
        
        # x,y,z location of anchor
        self.r = r
        
        # anchor index in array mooring list (only used for shared moorings/anchors)
        self.aNum = aNum
        
        # list of mooring lines connected to the anchor
        self.mooringList = {}
        
        # MoorPy anchor object
        self.mpAnchor = None
        
        # anchor capacity
        self.anchorCapacity = None
        
        # Dictionaries for additional information
        self.loads = {}
        '''
        {
           ff: # horizontal maximum anchor loads [N]
           fz: # vertical maximum anchor loads [N]
           method: # dynamic or static
            }
        '''
        self.soilProps = {}
        # self.cost = {}
        
    def makeMoorPyAnchor(self, ms):
        '''Create a MoorPy anchor object in a moorpy system
        Parameters
        ----------
        ms : class instance
            MoorPy system
        
        Returns
        -------
        ms : class instance
            MoorPy system 

        '''
        # create anchor as a fixed point in MoorPy system
        ms.addPoint(1,self.r)
        # assign this point as mpAnchor in the anchor class instance
        self.mpAnchor = ms.pointList[-1]

        # add mass if available
        if 'm' in self.dd['design'] and self.dd['design']['m']:
            self.mpAnchor.m = self.dd['design']['m']
        # set anchor diameter
        if 'd' in self.dd['design'] and self.dd['design']['d']:
            self.mpAnchor.d = self.dd['design']['d']
        # set the point as an anchor entity
        self.mpAnchor.entity= {'type': 'anchor', 'anchor_type':self.dd['type']}

        return(ms)
            
    def getMPForces(self, lines_only=False, seabed=True, xyz=False):   
        '''Find forces on anchor using MoorPy Point.getForces method and stores in loads dictionary
        Parameters
        ----------
        lines_only : boolean
            Calculate forces from just mooring lines (True) or not (False). Default is false
        seabed : boolean
            Include effect of seabed pushing up the anchor (True) or not (False) Default is false
        xyz : boolean
            Return forces in x,y,z DOFs (True) or only the enabled DOFs (False)
        '''
        # call getForces method from moorpy point object
        loads = self.mpAnchor.getForces(lines_only=lines_only, seabed=seabed, xyz=xyz)
        self.loads['ff'] = np.sqrt(loads[0]**2+loads[1]**2)
        self.loads['fz'] = loads[2]
        # loads determined from moorpy are static (?)
        self.loads['method'] = 'static'
        
    def getCost(self):
        '''find costs of anchor from MoorProps and store in design dictionary

        '''
        # check if anchor loads are available
        if not self.loads:
            # if not, check if theres a moorpy anchor object and calculate loads from that
            if self.mpAnchor:
                print("Need anchor loads to obtain cost, using getMPForces to determine loads in MoorPy")
                self.getMPForces()
            elif self.ms:
                print('Need anchor loads to obtain cost, creating a MoorPy anchor object and using getMPForces to determine loads in MoorPy')
                self.makeMoorPyAnchor(self.ms)
                self.getMPForces()
            else:
                raise Exception("Need anchor loads to obtain cost")
        # check again if there are loads
        if self.loads:
            c = self.dd['cost'] # set location for clarity
            # calculate individual costs and total cost for the anchor
            c['matCost'], c['instCost'], c['decomCost'] = mp.Point.getcost(self.mpAnchor)
            c['totCost'] = c['matCost'] + c['instCost'] + c['decomCost']

    
    def getMass(self,uhc_mode):
        '''find mass and/or UHC of anchor from MoorProps and store in design dictionary
        Parameters
        ----------
        uhc_mode : boolean
            True : obtain UHC from mass
            False : obtain Masss and UHC from loads
        '''
        if uhc_mode: # if looking for UHC given mass
            if self.dd['design']['m']: # check anchor mass is given
                self.dd['design']['UHC'], self.dd['design']['m'], info = mp.MoorProps.getAnchorMass(uhc_mode=1, mass_int=self.dd['design']['m'], anchor=self.dd['type'], soil_type=self.anchorProps['soil_type'])
            else:
                raise Exception("Need anchor mass to calculate UHC when uhc_mode = 1")
        else: # if looking for mass and UHC given loads
            if self.loads: # check the loads section exists
                self.dd['design']['UHC'], self.dd['design']['m'], info = mp.MoorProps.getAnchorMass(uhc_mode=0, fx=self.loads['ff'], fz=self.loads['fz'], anchor=self.dd['type'],soil_type=self.dd['soil_type'],method=self.loads['method'])
            elif self.mpAnchor:
                print("Need anchor loads to obtain mass, using getMPForces to determine loads in MoorPy")
                self.getMPForces()
                self.dd['design']['UHC'], self.dd['design']['m'], info = mp.MoorProps.getAnchorMass(uhc_mode=0, fx=self.loads['ff'], fz=self.loads['fz'], anchor=self.dd['type'],soil_type=self.dd['soil_type'],method=self.loads['method'])
            elif self.ms:
                print('Need anchor loads to obtain mass, creating a MoorPy anchor object and using getMPForces to determine loads in MoorPy')
                self.makeMoorPyAnchor(self.ms)
                self.getMPForces()
                self.dd['design']['UHC'], self.dd['design']['m'], info = mp.MoorProps.getAnchorMass(uhc_mode=0, fx=self.loads['ff'], fz=self.loads['fz'], anchor=self.dd['type'],soil_type=self.dd['soil_type'],method=self.loads['method'])
            else:
                raise Exception("Need anchor loads to obtain mass")