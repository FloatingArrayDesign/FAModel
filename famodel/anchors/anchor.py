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
        self.anchorCapacity = {}
        
        # Dictionaries for additional information
        self.loads = {}
        '''
        {
           Hm:      # horizontal maximum anchor loads at mudline [N]
           Vm:      # vertical maximum anchor loads at mudline [N]
           thetam: # angle of load at the mudline [rad]
           Ha:      # horizontal maximum loads at lug
           Va:      # vertical maximum loads at lug
           thetaa:  # angle of load at lug
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
        loads : dict, optional
            Dictionary of loads on the anchor at the lug point. If not provided, will use the loads dictionary property 
            of the anchor. If this is empty and it is needed for the capacity function (i.e. driven piles) then 
            the anchor.getMPForces() function will be called.

        Returns
        -------
        results : dict
            Dictionary of capacity of the anchor (generally a max force [kN] in H and V, but can be a max displacement (driven, dandg piles))

        '''
        anchType = self.dd['type'] 
        geom = self.dd['design']# geometric anchor information

        if not ground_conds: 
            ground_conds = self.dd['soil_properties']
         
        for key,prop in ground_conds.items():
            if len(prop)>1:
                print('Warning: Only homogeneous soils are supported at this time. Only the first item in a property list will be used.')
                break               
        soil = self.dd['soil_type'] # soil type
        
        # loads only used for driven and drilled and grouted piles
        if not loads:
            loads = self.loads
        
       
        # logic to determine what functions to call based on anchor type and soil type...
        if anchType == 'SEPLA' or anchType == 'DEA' or anchType == 'DEPLA' or anchType == 'VLA' or anchType == 'plate':
            from .anchors_famodel.capacity_plate import getCapacityPlate
            if 'clay' in soil or 'mud' in soil:
                # write or overwrite beta in geom dictionary from loads function
                if anchType != 'DEA':
                    if not 'beta' in geom:
                        if not 'thetaa' in loads:
                            loads = self.getMPForces(lug_loads=True)
                        geom['beta'] = 90 - loads['thetaa']
                else:
                    geom['beta'] = 0
                if 'Su0' in ground_conds and 'k' in ground_conds and 'gamma' in ground_conds:             
                    results = getCapacityPlate(geom['A'], geom['beta'], geom['zlug'], 'clay', ground_conds['gamma'][0],
                                               Su0=ground_conds['Su0'][0], k=ground_conds['k'][0])
                else:
                    raise Exception('Ground conditions dictionary needs Su0, k, gamma information for clay plate anchors')
            else:
                print(f'Warning: Soil type {soil} is not compatible with plate anchors (SEPLA/DEPLA/DEA/VLA)')
                                
        elif anchType == 'suction' or anchType == 'suction_pile':
            from .anchors_famodel.capacity_suction import getCapacitySuction
            print('suction pile')
            # check loads have been calculated (needed for capacity function in this case)
            if not 'Ha' in loads:
                # call getMPForces function 
                loads = self.getMPForces(lug_loads=True)
            if 'sand' in soil:
                if 'phi' in ground_conds and 'beta' in ground_conds:
                    results = getCapacitySuction(geom['D'],geom['L'],geom['zlug'],loads['Ha']/1000,loads['Va']/1000,'sand',ground_conds['gamma'][0],
                                                 phi=ground_conds['phi'][0],beta=ground_conds['beta'][0])
                else:
                    raise Exception('Ground conditions dictionary needs phi and beta information for sand suction pile anchor')
            elif 'clay' in soil or 'mud' in soil:
                if 'Su0' in ground_conds and 'k' in ground_conds and 'alpha' in ground_conds:# and 'gamma_sub' in ground_conds:
                    results = getCapacitySuction(geom['D'],geom['L'],geom['zlug'],loads['Ha']/1000,loads['Va']/1000,'clay',ground_conds['gamma'][0],
                                                 Su0=ground_conds['Su0'][0],k=ground_conds['k'][0],alpha=ground_conds['alpha'][0])
                    results['Horizontal max.'] = results['Horizontal max.']
                    results['Vertical max.'] = results['Vertical max.']
                    
                else:
                    raise Exception('Ground conditions dictionary needs Su0, k, and alpha information for clay suction pile anchor')
            else:
                print(f'Warning: Soil type {soil} is not compatible with suction pile anchor')
        elif anchType == 'helical' or anchType == 'helical_pile':
            from .anchors_famodel.capacity_helical import getCapacityHelical
            if 'sand' in soil:
                if not 'alpha_star' in ground_conds:
                    ground_conds['alpha_star'] = ground_conds['alpha']
                if 'phi' in ground_conds and 'gamma' in ground_conds:
                    results = getCapacityHelical(geom['D'],geom['L'],geom['d'],'sand',ground_conds['gamma'][0],
                                                 ground_conds['alpha_star'][0],phi=ground_conds['phi'][0])
                    results['Vertical max.'] = results['Capacity']
                else:
                    raise Exception('Ground conditions dictionary needs phi, gamma and beta information for clay helical pile anchor')
            elif 'clay' in soil or 'mud' in soil:
                if not 'alpha_star' in ground_conds:
                    ground_conds['alpha_star'] = ground_conds['alpha']
                if 'Su0' in ground_conds and 'k' in ground_conds and 'gamma' in ground_conds and 'alpha_star' in ground_conds:
                    results = getCapacityHelical(geom['D'],geom['L'],geom['d'],'clay',ground_conds['gamma'][0],ground_conds['alpha_star'][0],
                                       Su0=ground_conds['Su0'][0],k=ground_conds['k'][0])
                    results['Vertical max.'] = results['Capacity']
                else:
                    raise Exception('Ground conditions dictionary needs Su0, k, gamma, and alpha_star information for clay helical pile anchor')
            else:
                print(f'Warning: Soil type {soil} is not compatible with helical pile anchor')
        elif anchType == 'torpedo' or anchType == 'torpedo_pile':
            from .anchors_famodel.capacity_torpedo import getCapacityTorpedo
            if 'clay' in soil or 'mud' in soil:
                if 'Su0' in ground_conds and 'k' in ground_conds and 'alpha' in ground_conds:
                    results = getCapacityTorpedo(geom['D1'],geom['D2'],geom['L1'],geom['L2'],geom['zlug'],'clay',ground_conds['Su0'][0],
                                                 ground_conds['k'][0],ground_conds['alpha'][0])
                    results['Horizontal max.'] = results['Horizontal max.']
                    results['Vertical max.'] = results['Vertical max.']
                else:
                    raise Exception('Ground conditions dictionary needs Su0, k, and alpha information')
            else:
                print('Warning: Soil type {soil} is not compatible with torpedo pile anchor')
        elif anchType == 'driven' or anchType == 'driven_pile': # driven pile anchor
            # check loads have been calculated (needed for capacity function in this case)
            if not 'Ha' in loads:
                # call getMPForces function 
                loads = self.getMPForces(lug_loads=True)
            # check soil
            if 'weak_rock' in soil:
                from .anchors_famodel.capacity_drivenrock import getCapacityDrivenRock
                
                if not profile:
                    if 'UCS' in ground_conds and 'Em' in ground_conds:
                        profile = [[0,ground_conds['UCS'][0],ground_conds['Em'][0]],[75,ground_conds['UCS'][0],ground_conds['Em'][0]]] #profile = [list(x) for x in list(zip(ground_conds['depth'],ground_conds['UCS'],ground_conds['Em']))]
                    else:
                        raise Exception('Ground conditions dictionary needs UCS, Em, and depth information for weak rock driven pile anchor')
                results = getCapacityDrivenRock(profile,geom['L'],geom['D'],geom['zlug'],loads['Va'],loads['Ha'])
                   
            elif 'sand' in soil:
                from .anchors_famodel.capacity_drivensoil import getCapacityDrivenSoil
                if profile or ('gamma' in ground_conds and 'phi' in ground_conds):
                    if not profile:
                        profile = [[0,ground_conds['phi'][0],ground_conds['gamma'][0]],[75,ground_conds['phi'][0],ground_conds['gamma'][0]]] #profile = [list(x) for x in list(zip(ground_conds['depth'],ground_conds['phi'],ground_conds['gamma']))]

                    results = getCapacityDrivenSoil(profile,'sand',geom['L'],geom['D'],geom['zlug'],loads['Va'],loads['Ha'])
                else:
                    raise Exception('Ground conditions dictionary needs phi, gamma, and depth information for sand driven pile anchor')
            elif 'clay' in soil or 'mud' in soil:
                from .anchors_famodel.capacity_drivensoil import getCapacityDrivenSoil
                #if profile or ('Su' in ground_conds and 'gamma' in ground_conds and 'depth' in ground_conds) or ('Su0' in ground_conds and 'k' in ground_conds):
                if not profile:
                    if 'Su' in ground_conds and 'depth' in ground_conds and 'gamma' in ground_conds:
                        profile = [list(x) for x in list(zip(ground_conds['depth'],ground_conds['Su'],ground_conds['gamma']))]
                    elif 'Su0' in ground_conds and 'k' in ground_conds and 'gamma' in ground_conds:
                        Su = ground_conds['Su0'][0]+ground_conds['k'][0]*75
                        profile = [[0,ground_conds['Su0'][0],ground_conds['gamma'][0]],[75,Su,ground_conds['gamma'][0]]]
                    else:
                        raise Exception('Ground conditions dictionary needs information for clay driven pile anchor')

                results = getCapacityDrivenSoil(profile,'clay',geom['L'],geom['D'],geom['zlug'],loads['Va'],loads['Ha'])
                    
            else:
                print(f'Warning: Soil type {soil} is not compatible with driven pile anchors')
                        
        elif anchType == 'dandg_pile' or anchType == 'dandg':  # drill and grout pile
            from .anchors_famodel.capacity_dandg import getCapacityDandG
            # check for correct soil
            if 'rock' in soil:
                # check loads have been calculated (needed for capacity function in this case)
                if not loads:
                    # call getMPForces function 
                    loads = self.getMPForces(lug_loads=True)
                # check for correct ground properties
                if profile or ('UCS' in ground_conds and 'Em' in ground_conds):
                    if not profile:
                        profile = [[0,ground_conds['UCS'][0],ground_conds['Em'][0]],[75,ground_conds['UCS'][0],ground_conds['Em'][0]]] #[list(x) for x in list(zip(ground_conds['depth'],ground_conds['UCS'],ground_conds['Em']))]
                        # for x in profile:
                        #     x.append('Reese')
                    results = getCapacityDandG(profile,geom['L'],geom['D'],geom['zlug'],loads['Va'],loads['Ha'])
                else:
                    raise Exception('Ground conditions dictionary need UCS and Em information for drill and grout pile')
            else:
                print(f'Warning: soil type {soil} is not compatible with drill and grout pile')
        else:
            raise Exception(f'Anchor type {anchType} is not supported at this time')
            
        
        # capacity = cap*installAdj ??? OR is installAdj an input to the capacity functions?
        # save capacity 
        if 'dandg' in anchType or 'driven' in anchType: # will take in dandg, dandg_pile, driven, driven_pile
            self.anchorCapacity['Lat_max'] = results[2]['Lateral displacement']
            self.anchorCapacity['Rot_max'] = results[2]['Rotational displacement']
        else:
            if 'Horizontal max.' in results:
                self.anchorCapacity['Ha_max'] = results['Horizontal max.']*1000 # [N]
            self.anchorCapacity['Va_max'] = results['Vertical max.']*1000 # [N]
        return(results)
            
    def getMPForces(self, max_force=False,lines_only=False, seabed=True, xyz=False,mud_loads=None,lug_loads=False):   
        '''Find forces on anchor at mudline using the platform.getWatchCircle method or MoorPy Point.getForces method. 
        Optionally, get forces at anchor lug location with getTransferLoad function in capacity_loads.py.
        Stores in loads dictionary
        Parameters
        ----------
        max_force : boolean, optional
            Find and save the maximum force on the anchor (True) or just get force at the current MoorPy system state (False)           
        lines_only : boolean, optional
            Calculate forces from just mooring lines (True) or not (False). Default is false
        seabed : boolean, optional
            Include effect of seabed pushing up the anchor (True) or not (False). Default is true
        xyz : boolean, optional
            Return forces in x,y,z DOFs (True) or only the enabled DOFs (False). Default is false
        lug_loads : boolean, optional
            Return loads from getTransferLoad function (forces at anchor lug)
        mud_loads : dictionary, optional
            Dictionary of mudline loads, which will be used as inputs to the getTransferLoad()
            
        '''
        if max_force:
            pass
            # # find platform associated with this anchor
            # for att in self.attachments.items():
            #     if isinstance(att['obj'],Mooring):
            #         for attM in att['obj'].attached_to:
            #             if isinstance(attM,Platform):
            #                 locx,locy,maxVals = attM.getWatchCircle()
            #                 load = maxVals['minTenSF']*att['obj'].dd['']
        # call getForces method from moorpy point object
        else:
            loads = self.mpAnchor.getForces(lines_only=lines_only, seabed=seabed, xyz=xyz)
            self.loads['Hm'] = np.sqrt(loads[0]**2+loads[1]**2) # mudline forces in [N]
            self.loads['Vm'] = loads[2] # [N]
            self.loads['thetam'] = np.degrees(np.arctan(self.loads['Vm']/self.loads['Hm'])) # [deg]
            Tm =  np.sqrt(self.loads['Hm']**2+self.loads['Vm']**2) # [N]
        
        
        if lug_loads:
            from .anchors_famodel.capacity_load import getTransferLoad
            if 'zlug' in self.dd['design']:
                if self.dd['design']['zlug'] > 0:
                    # get line type
                    for att in self.attachments.values():
                        if isinstance(att['obj'],Mooring):
                            mtype = att['obj'].dd['sections'][0]['type']['material']
                            md = att['obj'].dd['sections'][0]['type']['d_nom']
                            mw = att['obj'].dd['sections'][0]['type']['w']
                    soil = self.dd['soil_type']
                    ground_conds = self.dd['soil_properties']
                    if 'clay' in soil or 'mud' in soil:
                        # Tm, thetam, zlug, line_type, d, soil_type, Su0=None, k=None, w=None
                        loadresults = getTransferLoad(Tm/1000,self.loads['thetam'],self.dd['design']['zlug'],mtype,md,
                                                       'clay',Su0=ground_conds['Su0'][0],k=ground_conds['k'][0],w=mw/1000) # output Ha and Va    (convert weight to kN/m)   
                    elif 'sand' in soil:
                            soil = 'sand'
                            #loadresults = getAnchorLoad(Tm/1000,self.loads['thetam'],self.dd['design']['zlug'],md,
                                                           #soil,gamma=ground_conds['gamma']) # output Ha and Va  (convert weight to kN/m)
                    if 'rock' in soil:
                        raise ValueError('zlug should be <= 0 for rock.')
                    else:
                        self.loads['Ha'] = loadresults['H']*1000 # [N]
                        self.loads['Va'] = loadresults['V']*1000 # [N]
                        self.loads['thetaa'] = loadresults['angle'] # [deg]
                else:
                    # Ha = Hm because zlug is at mudline or above
                    self.loads['Ha'] = self.loads['Hm'] # [N]
                    self.loads['Va'] = self.loads['Vm'] # [N]
                    self.loads['thetaa'] = self.loads['thetam'] # [deg]

        
        # loads determined from moorpy are static
        self.loads['method'] = 'static'
        
        return(self.loads)
    
    def getFS(self):
        '''
        Compute safety factor for loads on the anchor

        Returns
        -------
        FS : dict
            Dictionary of safety factors (often horizontal and vertical load SFs, but could be displacement SFs (drilled and grouted/driven piles))

        '''
        if not 'Ha' in self.loads:
            self.getMPForces(lug_loads=True)
        if not self.anchorCapacity:
            self.getAnchorCapacity()
         
        # look for load dictionary key in capacity dictionary
        FS = {}
        for Lkey,Lval in self.loads.items():
            for Ckey,Cval in self.anchorCapacity.items():
                if Lkey in Ckey:
                    FS[Lkey] = Cval/Lval
                    
                    
        return(FS)
            
        
               
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