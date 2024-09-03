"""Anchor class for FAModel, containing information and key methods for anchors of mooring lines
    Work in progress
"""
import moorpy as mp
import numpy as np
from famodel.famodel_base import Node
from famodel.mooring.mooring import Mooring

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
                    # all geometric info from yaml file
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
               
        # MoorPy anchor object
        self.mpAnchor = None
        
        # anchor capacity
        self.anchorCapacity = None
        
        # Dictionaries for additional information
        self.loads = {}
        '''
        {
           ff:      # horizontal maximum anchor loads [N]
           fz:      # vertical maximum anchor loads [N]
           Tm:      # load at the mudline [N]
           theta_m: # angle of load at the mudline [rad]
           method:  # dynamic or static method of calculation
            }
        '''
        self.soilProps = {}
        self.failure_probability = {}
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
    
    
    def getAnchorCapacity(self,ground_conds=None,installAdj=1,profile=None,loads=None):
        '''
        Calls anchor capacity functions developed by Felipe Moreno for the correct anchor type 

        Parameters
        ----------
        ground_conds : dict, optional
            Ground conditions ex: UCS,Em,phi,gamma,effective stress,etc. The default is None.
            If no dict provided, the ground conds will be pulled from the dd['soil_properties']
        installAdj : float, optional
            Adjustment to the capacity based on installation (dummy variable for now, but future installation functions
                                                              will dictate this value)
        profile : 2D array,  optional
            2d array of depths (m) and corresponding undrained shear strength (Pa). Su must not be zero 
            (change to small value such as .001), but z must always start at 0. Ex: array([z1,Su1],[z2,Su2],...)
            Used only for driven pile and drilled and grouted pile anchors.

        Returns
        -------
        capacity : float
            Capacity of the anchor
        info : dict
            Dictionary of other information on the anchor such as mass, cost, etc. May vary by type of anchor.

        '''
        anchType = self.dd['type'] 
        geom = self.dd['design']# geometric anchor information

        if not ground_conds: 
            ground_conds = self.dd['soil_properties']
         
        for key,prop in ground_conds:
            if len(prop)>1:
                print('Warning: Only homogeneous soils are supported at this time. Only the first item in a property list will be used.')
                break               
        soil = self.dd['soil_type'] # soil type
        
        if not loads:
            loads = self.loads
        
       
        # logic to determine what functions to call based on anchor type and model type...
        if anchType == 'SEPLA':
            from .capacity_plate import getCapacityPlate
            if 'clay' in soil or 'mud' in soil:
                if 'Su0' in ground_conds and 'k' in ground_conds and 'gamma' in ground_conds:
                    # # get line type of connected mooring 
                    # for att in self.attached:
                    #     if isinstance(att,Mooring):
                    #         mtype = att.dd['sections'][0]['type']['material']
                    #         md = att.dd['sections'][0]['type']['d_nom']               
                    results = getCapacityPlate(geom['A'], geom['Hs'], geom['beta'], soil_type='clay', gamma=ground_conds['gamma'],
                                               Su0=ground_conds['Su0'], k=ground_conds['k'],H=loads['Ha'],V=loads['Va'])
                else:
                    raise Exception('Ground conditions dictionary needs Su0, k, gamma, information for clay SEPLA')
            else:
                print(f'Warning: Soil type {soil} is not compatible with plate anchors (SEPLA)')
                                
        elif anchType == 'DEA':
            from .capacity_plate import getCapacityPlate
            if 'clay' in soil or 'mud' in soil:
                if 'Su0' in ground_conds and 'k' in ground_conds:
                    results = getCapacityPlate(geom['A'],geom['Hs'],soil_type='clay',Su0=ground_conds['Su0'],k=ground_conds['k'])
                else:
                    raise Exception('Ground conditions dictionary needs Su0, k information for clay DEA')
            else:
                print(f'Warning: Soil type {soil} is not compatible with drag-embedment anchors')
        elif anchType == 'SCA':
            from .capacity_suction import getCapacitySuction
            if 'sand' in soil:
                if 'phi' in ground_conds and 'beta' in ground_conds:
                    results = getCapacitySuction(geom['D'],geom['L'],soil_type='sand',phi=ground_conds['phi'],beta=ground_conds['beta'],rhows=geom['rhows'])
                else:
                    raise Exception('Ground conditions dictionary needs phi and beta information for sand suction pile anchor')
            elif 'clay' in soil or 'mud' in soil:
                if 'Su0' in ground_conds and 'k' in ground_conds and 'alpha' in ground_conds:# and 'gamma_sub' in ground_conds:
                    results = getCapacitySuction(geom['D'],geom['L'],Su0=ground_conds['Su0'],
                                       k=ground_conds['k'],alpha=ground_conds['alpha'],rhows=geom['rhows'])
                else:
                    raise Exception('Ground conditions dictionary needs Su0, k, and alpha information for clay suction pile anchor')
            else:
                print(f'Warning: Soil type {soil} is not compatible with suction pile anchor')
        elif anchType == 'helical':
            from .capacity_helical import getCapacityHelical
            if 'sand' in soil:
                if 'phi' in ground_conds and 'gamma' in ground_conds and 'alpha_star' in ground_conds:
                    results = getCapacityHelical(geom['D'],geom['L'],geom['d'],soil_type='sand',gamma=ground_conds['gamma'],
                                                 alpha_star=ground_conds['alpha_star'],phi=ground_conds['phi'])
                else:
                    raise Exception('Ground conditions dictionary needs phi, gamma and beta information for clay helical pile anchor')
            elif 'clay' in soil or 'mud' in soil:
                if 'Su0' in ground_conds and 'k' in ground_conds and 'gamma' in ground_conds and 'alpha_star' in ground_conds:
                    getCapacityHelical(geom['D'],geom['L'],geom['d'],soil_type='clay',Su0=ground_conds['Su0'],k=ground_conds['k'],
                                       gamma=ground_conds['gamma'],alpha_star=ground_conds['alpha_star']) 
                else:
                    raise Exception('Ground conditions dictionary needs Su0, k, gamma, and alpha_star information for clay helical pile anchor')
            else:
                print(f'Warning: Soil type {soil} is not compatible with helical pile anchor')
        elif anchType == 'torpedo':
            from .capacity_torpedo import getCapacityTorpedo
            if 'clay' in soil or 'mud' in soil:
                if 'Su0' in ground_conds and 'k' in ground_conds and 'alpha' in ground_conds:
                    # # get line type of connected mooring 
                    # for att in self.attached:
                    #     if isinstance(att,Mooring):
                    #         mtype = att.dd['sections'][0]['type']['material']
                    #         md = att.dd['sections'][0]['type']['d_nom']
                    results = getCapacityTorpedo(geom['D1'],geom['D2'],geom['L1'],geom['L2'],geom['Hp'],ground_conds['Su0'],
                                                 ground_conds['k'],ground_conds['alpha'],)
                else:
                    raise Exception('Ground conditions dictionary needs Su0, k, and alpha information')
            else:
                print('Warning: Soil type {soil} is not compatible with torpedo pile anchor')
        elif anchType == 'driven': # driven pile anchor
            # check soil
            if 'weak_rock' in soil:
                from capacity_drivenrock import getCapacityDrivenRock
                if profile or ('UCS' in ground_conds and 'Em' in ground_conds and 'depth' in ground_conds):
                    if not profile:
                        profile = [list(x) for x in list(zip(ground_conds['depth'],ground_conds['UCS'],ground_conds['Em']))]
                        # for x in profile:
                        #     x.append('Reese')
                    results = getCapacityDrivenRock(profile,geom['L'],geom['D'],geom['t'],geom['E'])
                else:
                    raise Exception('Ground conditions dictionary needs UCS, Em, and depth information for weak rock driven pile anchor')
            elif 'sand' in soil:
                from capacity_drivensoil import getCapacityDrivenSoil
                if profile or ('gamma' in ground_conds and 'phi' in ground_conds and 'depth' in ground_conds):
                    if not profile:
                        profile = [list(x) for x in list(zip(ground_conds['depth'],ground_conds['phi'],ground_conds['gamma']))]
                        # for x in profile:
                        #     x.extend('Matlock',0.02)
                    results = getCapacityDrivenSoil(profile,geom['L'],geom['D'],geom['t'],geom['E'],soil_type='sand')
                else:
                    raise Exception('Ground conditions dictionary needs phi, gamma, S_eff, k, and depth information for sand driven pile anchor')
            elif 'clay' in soil or 'mud' in soil:
                from capacity_drivensoil import getCapacityDrivenSoil
                if profile or ('Su' in ground_conds and 'gamma' in ground_conds and 'depth' in ground_conds):
                    if not profile:
                        profile = [list(x) for x in list(zip(ground_conds['depth'],ground_conds['Su'],ground_conds['gamma']))]
                        # for x in profile:
                        #     x.extend('Matlock',0.02)
                    results = getCapacityDrivenSoil(profile,geom['L'],geom['D'],geom['t'],geom['E'],soil_type='clay')
                else:
                    raise Exception('Ground conditions dictionary needs Su, gamma,S_eff and depth information for clay driven pile anchor')
            else:
                print(f'Warning: Soil type {soil} is not compatible with driven pile anchors')
                        
        elif anchType == 'dandg_pile':  # drill and grout pile
            from .capacity_dandg import getCapacityDandG
            # check for correct soil
            if 'rock' in soil:
                # check for correct ground properties
                if profile or ('UCS' in ground_conds and 'Em' in ground_conds and 'depth' in ground_conds):
                    if not profile:
                        profile = [list(x) for x in list(zip(ground_conds['depth'],ground_conds['UCS'],ground_conds['Em']))]
                        # for x in profile:
                        #     x.append('Reese')
                    results = getCapacityDandG(profile,geom['L'],geom['D'],geom['t'],geom['E'])
                else:
                    raise Exception('Ground conditions dictionary need UCS, Em, and depth information for drill and grout pile')
            else:
                print(f'Warning: soil type {soil} is not compatible with drill and grout pile')
        else:
            raise Exception(f'Anchor type {anchType} is not supported at this time')
            
        
        # capacity = cap*installAdj ??? OR is installAdj an input to the capacity functions?
        # save capacity 
        # return(capacity,info)
            
    def getMPForces(self, lines_only=False, seabed=True, xyz=False,capacity_loads=False):   
        '''Find forces on anchor using MoorPy Point.getForces method and stores in loads dictionary
        Parameters
        ----------
        lines_only : boolean
            Calculate forces from just mooring lines (True) or not (False). Default is false
        seabed : boolean
            Include effect of seabed pushing up the anchor (True) or not (False). Default is true
        xyz : boolean
            Return forces in x,y,z DOFs (True) or only the enabled DOFs (False). Default is false
        '''
        # call getForces method from moorpy point object
        loads = self.mpAnchor.getForces(lines_only=lines_only, seabed=seabed, xyz=xyz)
        self.loads['Hm'] = np.sqrt(loads[0]**2+loads[1]**2)
        self.loads['Vm'] = loads[2]
        self.loads['theta_m'] = np.arctan(self.loads['Hm']/self.loads['Vm'])
        Tm =  np.sqrt(loads[0]**2+loads[1]**2+loads[2]**2)
    
        if capacity_loads:
            from .capacity_load import getAnchorLoad
            if 'zlug' in self.dd['design']:
                if self.dd['design']['zlug'] != 0:
                    # get line type
                    for att in self.attached:
                        if isinstance(att,Mooring):
                            mtype = att.dd['sections'][0]['type']['material']
                            md = att.dd['sections'][0]['type']['d_nom']
                    soil = self.dd['soil_type']
                    ground_conds = self.dd['soil_properties']
                    if 'clay' in soil or 'mud' in soil:
                        loadresults = getAnchorLoad(Tm,self.loads['theta_m'],self.dd['design']['zlug'],mtype,md,soil,
                                                    Su0=ground_conds['Su0'],k=ground_conds['k']) # output Ha and Va       
                    else:
                        loadresults = getAnchorLoad(Tm,self.loads['theta_m'],self.dd['design']['zlug'],mtype,md,soil) # output Ha and Va       
                    
                    #
        
        
        # loads determined from moorpy are static
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