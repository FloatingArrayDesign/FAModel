"""Anchor class for FAModel, containing information and key methods for anchors of mooring lines
    Work in progress
"""
import moorpy as mp
import numpy as np
from famodel.famodel_base import Node
from famodel.mooring.mooring import Mooring
import famodel.platform.platform 


class Anchor(Node):
    
    def __init__(self, dd=None, ms=None, r=[0,0,0], aNum=None, id=None, 
                 g=9.81, rho=1025):
        '''
        Parameters
        ----------
        dd: dictionary
            Design dictionary that contains all information on an anchor for a mooring line/shared mooring
            {
                type:   # anchor type (plate,suction_pile,torpedo_pile,helical_pile,driven_pile,dandg_pile)
                design: # all geometric info from yaml file, only need to include info relevant to your anchor type
                     A    plate anchor area
                     D    anchor diameter (or helix diameter for helical piles)
                     D1   torpedo anchor wing diameter
                     D2   torpedo anchor shaft diameter
                     d    helical pile shaft diameter 
                     L    pile anchor length
                     L1   torpedo anchor wing length 
                     L2   torpedo anchor shaft length
                     zlug padeye z elevation (+ down into the soil)
                     beta angle of plate anchor after keying (optional)
                soil_type:    <<< not design information
                angle: # seabed angle   <<< not design information
                cost:
                    matCost: # material cost
                    instCost: # installation cost
                    decomCost: # decomissioning cost                          
            }
        ms: system object
            MoorPy system object the anchor is in
            
        r: list
            Location of anchor in x,y,z
            
        aNum: int
            entry number in project.anchorList dictionary (may remove...)
        id: str/int
            unique id of this object
        g: float
            acceleration due to gravity in m/s^2
        rho: float
            density of water in kg/m^3
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
        
        # get environment info
        self.g = g # acceleration due to gravity (m/s^2)
        self.rho = rho # density of fluid (kg/m^3)
        
        # anchor mass
        if 'mass' in self.dd['design']:
            self.mass = self.dd['design']['mass']
        else:
            self.mass = None
        
        # Dictionaries for additional information
        # anchor capacity
        self.anchorCapacity = {}
        
        # anchor costs
        self.cost = {}
        
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
        self.mpAnchor.entity= {'type': 'anchor'}
        if 'type' in self.dd:
            self.mpAnchor.entity['anchor_type']=self.dd['type']
        
        return(ms)
    
    
    def getAnchorCapacity(self,ground_conds=None,installAdj=1,profile=None,loads=None,plot=True):
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
            Dictionary of loads on the anchor at the lug point in [N]. If not provided, will use the loads dictionary property 
            of the anchor. If this is empty and it is needed for the capacity function (i.e. driven piles) then 
            the anchor.getLugForces() function will be called.

        Returns
        -------
        results : dict
            Dictionary of capacity of the anchor (generally a max force [N] in H and V, but can be a max displacement (driven, dandg piles))

        '''
        # - - - - set details - - - -
        anchType = self.dd['type'] 
        geom = self.dd['design']# geometric anchor information

        if not ground_conds: 
            ground_conds = self.dd['soil_properties']
         
        for key,prop in ground_conds.items():
            if isinstance(prop,list) or isinstance(prop,np.ndarray):
                if len(prop)>1:
                    print('Warning: Only homogeneous soils are supported at this time. Only the first item in a property list will be used.')
                    break 
                else:
                    ground_conds[key] = prop[0]
        soil = self.dd['soil_type'] # soil type
        
        if loads:
            # find out if mudline loads or anchor loads
            if not 'Ha' in loads:
                # get loads at lug
                loads = self.getLugForces(mudloads=loads,plot=plot)
        else:
            loads = self.loads
        
        
       
        # logic to determine what functions to call based on anchor type and soil type...
        
        # - - - - plate anchors - - - - 
        if anchType == 'SEPLA' or anchType == 'DEA' or anchType == 'DEPLA' or anchType == 'VLA' or anchType == 'plate':
            from .anchors_famodel.capacity_plate import getCapacityPlate
            if 'clay' in soil or 'mud' in soil:
                # write or overwrite beta in geom dictionary from loads function
                if anchType != 'DEA':
                    if not 'beta' in geom:
                        if not 'thetaa' in loads:
                            loads = self.getLugForces(plot=plot)
                        geom['beta'] = 90 - loads['thetaa']
                else:
                    geom['beta'] = 0
                if 'Su0' in ground_conds and 'k' in ground_conds and 'gamma' in ground_conds:             
                    results = getCapacityPlate(geom['A'], geom['beta'], geom['zlug'], 'clay', ground_conds['gamma'],
                                               Su0=ground_conds['Su0'], k=ground_conds['k'])
                else:
                    raise Exception('Ground conditions dictionary needs Su0, k, gamma information for clay plate anchors')
            else:
                print(f'Warning: Soil type {soil} is not compatible with plate anchors (SEPLA/DEPLA/DEA/VLA)')
         
        # - - - - suction buckets - - - -                       
        elif 'suction' in anchType:
            from .anchors_famodel.capacity_suction import getCapacitySuction
            # check loads have been calculated (needed for capacity function in this case)
            if not 'Ha' in loads:
                # call getMPForces function 
                loads = self.getLugForces(plot=plot)
            if 'sand' in soil:
                if 'phi' in ground_conds and 'Dr' in ground_conds:
                    results = getCapacitySuction(geom['D'], geom['L'], geom['zlug'],
                                                 loads['Ha']/1000, loads['Va']/1000,
                                                 'sand', ground_conds['gamma'],
                                                 phi=ground_conds['phi'],
                                                 Dr=ground_conds['Dr'], plot=plot)
                else:
                    raise Exception('Ground conditions dictionary needs phi and relative density information for sand suction pile anchor')
            elif 'clay' in soil or 'mud' in soil:
                if 'Su0' in ground_conds and 'k' in ground_conds and 'alpha' in ground_conds:# and 'gamma_sub' in ground_conds:
                    results = getCapacitySuction(geom['D'],geom['L'], geom['zlug'],
                                                 loads['Ha']/1000, loads['Va']/1000,
                                                 'clay', ground_conds['gamma'],
                                                 Su0=ground_conds['Su0'],
                                                 k=ground_conds['k'], plot=plot)
                    results['Horizontal max.'] = results['Horizontal max.']
                    results['Vertical max.'] = results['Vertical max.']
                    
                else:
                    raise Exception('Ground conditions dictionary needs Su0, k, and alpha information for clay suction pile anchor')
            else:
                print(f'Warning: Soil type {soil} is not compatible with suction pile anchor')
                
        # - - - - helical piles - - - -
        elif 'helical' in anchType:
            from .anchors_famodel.capacity_helical import getCapacityHelical
            if 'sand' in soil:
                if 'phi' in ground_conds and 'gamma' in ground_conds:
                    results = getCapacityHelical(geom['D'], geom['L'], geom['d'],
                                                 geom['zlug'], 'sand', 
                                                 ground_conds['gamma'],
                                                 phi=ground_conds['phi'], 
                                                 Dr=ground_conds['Dr'])
                    results['Vertical max.'] = results['Capacity']
                else:
                    raise Exception('Ground conditions dictionary needs phi, gamma and relative density information for clay helical pile anchor')
            elif 'clay' in soil or 'mud' in soil:
                if not 'alpha_star' in ground_conds:
                    ground_conds['alpha_star'] = ground_conds['alpha']
                if 'Su0' in ground_conds and 'k' in ground_conds and 'gamma' in ground_conds:
                    results = getCapacityHelical(geom['D'], geom['L'], geom['d'],
                                                 geom['zlug'], 'clay', 
                                                 ground_conds['gamma'],
                                                 Su0=ground_conds['Su0'], 
                                                 k=ground_conds['k'])
                    results['Vertical max.'] = results['Capacity']
                else:
                    raise Exception('Ground conditions dictionary needs Su0, k, gamma, and alpha_star information for clay helical pile anchor')
            else:
                print(f'Warning: Soil type {soil} is not compatible with helical pile anchor')
        
        # - - - - torpedo piles - - - - 
        elif 'torpedo' in anchType:
            from .anchors_famodel.capacity_torpedo import getCapacityTorpedo
            if 'clay' in soil or 'mud' in soil:
                if 'Su0' in ground_conds and 'k' in ground_conds and 'alpha' in ground_conds:
                    results = getCapacityTorpedo(geom['D1'], geom['D2'], 
                                                 geom['L1'], geom['L2'], 
                                                 geom['zlug'], 'clay', 
                                                 ground_conds['Su0'],
                                                 ground_conds['k'], 
                                                 ground_conds['alpha'])
                    results['Horizontal max.'] = results['Horizontal max.']
                    results['Vertical max.'] = results['Vertical max.']
                else:
                    raise Exception('Ground conditions dictionary needs Su0, k, and alpha information')
            else:
                print('Warning: Soil type {soil} is not compatible with torpedo pile anchor')
        
        # - - - - driven piles - - - - 
        elif 'driven' in anchType: # driven pile anchor
            # check loads have been calculated (needed for capacity function in this case)
            if not 'Ha' in loads:
                # call getLugForces function 
                loads = self.getLugForces(plot=plot)
            H_inc = loads['Ha']*0.1 # increment H by 10% of Ha load in the while loops to back-calc max H from displacements
            H = 0
            # check soil
            if 'weak_rock' in soil:
                from .anchors_famodel.capacity_drivenrock import getCapacityDrivenRock
                
                if not profile:
                    if 'UCS' in ground_conds and 'Em' in ground_conds:
                        profile = [[0,ground_conds['UCS'],ground_conds['Em']], 
                                   [75,ground_conds['UCS'],ground_conds['Em']]] #profile = [list(x) for x in list(zip(ground_conds['depth'],ground_conds['UCS'],ground_conds['Em']))]
                    else:
                        raise Exception('Ground conditions dictionary needs UCS, Em, and depth information for weak rock driven pile anchor')
                
                y, z, results = getCapacityDrivenRock(profile, geom['L'], geom['D'], 
                                                      geom['zlug'], loads['Va'], 
                                                      loads['Ha'], plot=plot)
                
                # loop through, calling capacity with larger H values until a displacement value goes above limit
                while results['Lateral displacement']< 0.05*geom['D'] and results['Rotational displacement'] < 0.25:
                    # increment H 
                    H += H_inc
                    # call capacity function       
                    y, z, results = getCapacityDrivenRock(profile, geom['L'], 
                                                          geom['D'], geom['zlug'], 
                                                          loads['Va'], H=H, plot=plot)
                    
                   
            elif 'sand' in soil:
                from .anchors_famodel.capacity_drivensoil import getCapacityDrivenSoil
                if profile or ('gamma' in ground_conds and 'Dr' in ground_conds and 'phi' in ground_conds):
                    if not profile:
                        profile = [[0,ground_conds['phi'],ground_conds['gamma'],ground_conds['Dr']], 
                                   [75,ground_conds['phi'],ground_conds['gamma'],ground_conds['Dr']]] #profile = [list(x) for x in list(zip(ground_conds['depth'],ground_conds['phi'],ground_conds['gamma']))]

                    y, z, results = getCapacityDrivenSoil(profile, 'sand', 
                                                          geom['L'], geom['D'], 
                                                          geom['zlug'], loads['Va'], 
                                                          loads['Ha'], plot=plot)
                    if geom['zlug'] > 0:
                        # need to check bending moment if lug is below mudline (+ zlug)
                        # loop through, calling capacity with larger H values until a displacement value goes above limit
                        while results['Lateral displacement']<= 0.05*geom['D'] and results['Bending moment'] <= results['Plastic moment']:
                            # increment H by 10% of load
                            H += H_inc
                            # call capacity function     
                            y, z, results = getCapacityDrivenSoil(profile,'clay', 
                                                                  geom['L'], geom['D'],
                                                                  geom['zlug'], loads['Va'], 
                                                                  H=H, plot=plot)
                            
                    else:
                        while results['Lateral displacement']<= 0.05*geom['D'] and results['Rotational displacement'] <= 0.25:
                            # increment H by 10% of load
                            H += H_inc  
                            # call capacity function
                            y, z, results = getCapacityDrivenSoil(profile, 'clay', 
                                                                  geom['L'], geom['D'], 
                                                                  geom['zlug'], loads['Va'],
                                                                  H=H, plot=plot)
                else:
                    raise Exception('Ground conditions dictionary needs phi, gamma, and depth information for sand driven pile anchor')
            elif 'clay' in soil or 'mud' in soil:
                from .anchors_famodel.capacity_drivensoil import getCapacityDrivenSoil
                #if profile or ('Su' in ground_conds and 'gamma' in ground_conds and 'depth' in ground_conds) or ('Su0' in ground_conds and 'k' in ground_conds):
                if not profile:
                    if 'Su' in ground_conds and 'depth' in ground_conds and 'gamma' in ground_conds:
                        profile = [list(x) for x in list(zip(ground_conds['depth'],ground_conds['Su'],ground_conds['gamma']))]
                    elif 'Su0' in ground_conds and 'k' in ground_conds and 'gamma' in ground_conds:
                        Su = ground_conds['Su0']+ground_conds['k']*75
                        profile = [[0,ground_conds['Su0'],ground_conds['gamma']],[75,Su,ground_conds['gamma']]]
                    else:
                        raise Exception('Ground conditions dictionary needs information for clay driven pile anchor')

                y, z, results = getCapacityDrivenSoil(profile,'clay',geom['L'],geom['D'],geom['zlug'],loads['Va'],loads['Ha'], plot=plot)
                
                if geom['zlug'] > 0:
                    # need to check bending moment if lug is below mudline (+ zlug)
                    # loop through, calling capacity with larger H values until a displacement value goes above limit
                    while results['Lateral displacement']<= 0.05*geom['D'] and results['Bending moment'] <= results['Plastic moment']:
                        # increment H by 10% of load
                        H += H_inc
                        # call capacity function     
                        y, z, results = getCapacityDrivenSoil(profile,'clay',geom['L'],geom['D'],geom['zlug'],loads['Va'], H=H, plot=plot)
                
                else:
                    while results['Lateral displacement']<= 0.05*geom['D'] and results['Rotational displacement'] <= 0.25:
                        # increment H by 10% of load
                        H += H_inc  
                        # call capacity function
                        y, z, results = getCapacityDrivenSoil(profile,'clay',geom['L'],geom['D'],geom['zlug'],loads['Va'], H=H, plot=plot)
                        
                    
            else:
                print(f'Warning: Soil type {soil} is not compatible with driven pile anchors')
                
        # - - - - drilled and grouted piles - - - -                
        elif 'dandg' in anchType:  # drill and grout pile
            from .anchors_famodel.capacity_dandg import getCapacityDandG
            # check for correct soil
            if 'rock' in soil:
                # check loads have been calculated (needed for capacity function in this case)
                if not 'Ha' in loads:
                    # call getMPForces function 
                    loads = self.getLugForces(plot=plot)
                # check for correct ground properties
                if profile or ('UCS' in ground_conds and 'Em' in ground_conds):
                    if not profile:
                        profile = [[0,ground_conds['UCS'],ground_conds['Em']],[75,ground_conds['UCS'],ground_conds['Em']]] #[list(x) for x in list(zip(ground_conds['depth'],ground_conds['UCS'],ground_conds['Em']))]

                    # call capacity function once to get displacement values
                    y, z, results = getCapacityDandG(profile,geom['L'],geom['D'], 
                                                     geom['zlug'], loads['Va'], 
                                                     loads['Ha'], plot=plot)
                    H_inc = loads['Ha']*0.1 # increment H by 10% of Ha load
                    H = H_inc # start H at 10% of Ha load
                    # loop through, calling capacity with larger H values until a displacement value goes above limit
                    while results['Lateral displacement']< 0.05*geom['D'] and results['Rotational displacement'] < 0.25:
                        # call capacity function       
                        y, z, results = getCapacityDandG(profile, geom['L'], geom['D'], 
                                                         geom['zlug'], loads['Va'], 
                                                         H=H, plot=plot)
                        # increment H 
                        H += H_inc
                else:
                    raise Exception('Ground conditions dictionary need UCS and Em information for drill and grout pile')
            else:
                print(f'Warning: soil type {soil} is not compatible with drill and grout pile')
                
        # - - - - anchor type not recognized or supported - - - - 
        else:
            raise Exception(f'Anchor type {anchType} is not supported at this time')
            
        # - - - - save relevant results in dictionary using common terms - - - - 
        # capacity = cap*installAdj ??? OR is installAdj an input to the capacity functions?
        # save capacity 
        if 'dandg' in anchType or 'driven' in anchType: # will take in dandg, dandg_pile, driven, driven_pile
            self.anchorCapacity['Lat_max'] = results['Lateral displacement'] # [deg]
            if 'Rotational displacement' in results:
                self.anchorCapacity['Rot_max'] = results['Rotational displacement'] # [deg]
            elif 'Bending moment' in results:
                self.anchorCapacity['Mbend_max'] = results['Bending moment']
            self.anchorCapacity['Va_max'] = results['Axial capacity'] # [N]
            self.anchorCapacity['Ha_max'] = H

        else:
            if 'Horizontal max.' in results:
                self.anchorCapacity['Ha_max'] = results['Horizontal max.']*1000 # [N]
            self.anchorCapacity['Va_max'] = results['Vertical max.']*1000 # [N]
            
        if 'dandg' in anchType or 'driven' in anchType or 'helical' in anchType or 'torpedo' in anchType: 
            self.mass = results['Pile weight']*1000/self.g # mass in [kg]
        elif 'plate' in anchType or 'DEA' in anchType or 'SEPLA' in anchType or 'DEPLA' in anchType or 'VLA' in anchType:
            self.mass = results['Plate weight']*1000/self.g # mass in [kg]
        elif 'suction' in anchType:
            self.mass = results['Weight Pile']*1000/self.g # mass in [kg]
            
        # add on extra for drag-embedment anchors (flukes)
        if 'DEA' in anchType:
            self.mass *= 1.75
        
            
        return(results)
            
    def getMudlineForces(self, max_force=False,lines_only=False, seabed=True, xyz=False,project=None):   
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
            
        '''
        Platform = famodel.platform.platform.Platform
        if max_force:
            if project:
                # get watch circle of platform(s)
                project.arrayWatchCircle()
            else:
                # find platform associated with this anchor
                for att in self.attachments.values():
                    if isinstance(att['obj'],Mooring):
                        for attM in att['obj'].attached_to:
                            if isinstance(attM,Platform):
                                locx,locy,maxVals = attM.getWatchCircle()
        # call getForces method from moorpy point object
        else:
            loads = self.mpAnchor.getForces(lines_only=lines_only, seabed=seabed, xyz=xyz)
            self.loads['Hm'] = np.sqrt(loads[0]**2+loads[1]**2) # mudline forces in [N]
            self.loads['Vm'] = loads[2] # [N]
            self.loads['thetam'] = np.degrees(np.arctan(self.loads['Vm']/self.loads['Hm'])) # [deg]
            self.loads['mudline_load_type'] = 'current_state'
        
        # loads determined from moorpy are static
        self.loads['method'] = 'static'
        
        return(self.loads)
    
    def getLugForces(self,mudloads=None,max_force=True,plot=True):
        '''
        Find forces on an anchor at the lug point based on the mudline forces and angles. Calls getTransferFunction script

        Parameters
        ----------
        mudloads : dict, optional
            Dictionary of max mudline forces. The default is None.

        Returns
        -------
        loads: dict
            Dictionary of loads at the lug point [N]

        '''
        from .anchors_famodel.capacity_load import getTransferLoad
        
        nolugload = False
        
        if not mudloads:        
            if not self.loads:
                # get max mudline forces first
                self.getMudlineForces(max_force=max_force)
            elif not 'mudline_load_type' in self.loads:
                raise KeyError("Loads dictionary must specify 'mudline_load_type'='current_state' or 'mudline_load_type'='max', where 'max' indicates the loads are maximum loads.")
            elif max_force and self.loads['mudline_load_type'] != 'max':
                # need max forces, not current state
                self.getMudlineForces(max_force=True)
            mudloads = self.loads
        else:
            # check syntax
            if not 'Hm' in mudloads or not 'Vm' in mudloads:
                raise KeyError('Mudline load dictionary must have Hm and Vm for horizontal load and vertical load (in [N]) at the mudline')
            if not 'thetam' in mudloads:
                mudloads['thetam'] = np.degrees(np.arctan(mudloads['Vm']/mudloads['Hm']))
                
        def makeEqual_TaTm(mudloads):
            mudloads['Ha'] = mudloads['Hm'] # [N]
            mudloads['Va'] = mudloads['Vm'] # [N]
            mudloads['thetaa'] = mudloads['thetam'] # [deg]
        
        if 'zlug' in self.dd['design']:
            if self.dd['design']['zlug'] > 0:
                # get line type
                for att in self.attachments.values():
                    if isinstance(att['obj'],Mooring):
                        mtype = att['obj'].dd['sections'][0]['type']['material']
                        if not 'chain' in mtype:
                            print('No chain on seafloor, setting Ta=Tm')
                            nolugload = True
                            break
                        else:
                            md = att['obj'].dd['sections'][0]['type']['d_nom']
                            mw = att['obj'].dd['sections'][0]['type']['w']
                soil = self.dd['soil_type']
                ground_conds = self.dd['soil_properties']
                # update soil conds as needed to be homogeneous
                for key,prop in ground_conds.items():
                    if isinstance(prop,list) or isinstance(prop,np.ndarray):
                        if len(prop)>1:
                            print('Warning: Only homogeneous soils are supported at this time. Only the first item in a property list will be used.')
                            break 
                        else:
                            ground_conds[key] = prop[0]
                            
                Tm =  np.sqrt(mudloads['Hm']**2+mudloads['Vm']**2) # [N]
                if 'clay' in soil or 'mud' in soil and not nolugload:
                    # Tm, thetam, zlug, line_type, d, soil_type, Su0=None, k=None, w=None
                    try:
                        loadresults = getTransferLoad(Tm/1000,mudloads['thetam'],
                                                      self.dd['design']['zlug'],mtype,md,
                                                      'clay',Su0=ground_conds['Su0'],
                                                      k=ground_conds['k'],w=mw/1000,
                                                      plot=plot) # output Ha and Va    (convert weight to kN/m)   
                    except Exception as e:
                        print(e)
                        print('Unable to get loads at anchor lug location. Setting Ta = Tm')
                        nolugload = True
                elif 'sand' in soil and not nolugload:
                        soil = 'sand'
                        try:
                            loadresults = getTransferLoad(Tm/1000, self.loads['thetam'], 
                                                        self.dd['design']['zlug'],
                                                        mtype, md, soil,
                                                        gamma=ground_conds['gamma'], 
                                                        phi=ground_conds['phi'],
                                                        delta=ground_conds['delta'], 
                                                        w=mw/1000,plot=plot) # output Ha and Va  (convert weight to kN/m)
                        except Exception as e:
                            print(e)
                            print('Unable to get loads at anchor lug location. Setting Ta = Tm')
                            nolugload = True
                elif 'rock' in soil and not nolugload:
                    raise ValueError('zlug should be <= 0 for rock.')
                    
                # if loadresults['V']<0:
                #     # results are invalid
                #     print('Warning: invalid results for the combination of anchor ',self.dd['type'],' soil ',soil,' and loads ',mudloads,'. Setting Ha=Hm, Va=Vm, thetaa=thetam')
                #     makeEqual_TaTm(mudloads)
                if nolugload:
                    makeEqual_TaTm(mudloads)
                else:
                    mudloads['Ha'] = loadresults['H']*1000 # [N]
                    mudloads['Va'] = loadresults['V']*1000 # [N]
                    mudloads['thetaa'] = loadresults['angle'] # [deg]
            else:
                # Ha = Hm because zlug is at mudline or above
                makeEqual_TaTm(mudloads)
        else:
            print('No zlug given, assuming loads at mudline = loads at anchor lug')
            makeEqual_TaTm(mudloads)
            
        if not 'method' in mudloads:
            # assume mudloads are static unless told otherwise
            # loads determined from moorpy are static
            mudloads['method'] = 'static'
        else:
            mudloads['method'] = mudloads['method']
        
        return(self.loads)
    
    def getFS(self, loads=None, acceptance_crit=None):
        '''
        Compute safety factor for loads on the anchor
        
        Parameters
        ----------
        loads : dict, optional
            Dictionary of loads on the anchor. 
        acceptance_crit : dict, optional
            Dictionary of acceptable factors of safety for each load type.
            Key is the load type, and value is the minimum acceptable safety factor.
            Default is None (in which case no comparison between FS and acceptance criteria is calculated)

        Returns
        -------
        FS : dict
            Dictionary of safety factors (often horizontal and vertical load SFs, but could be displacement SFs (drilled and grouted/driven piles))
        acceptance : dict
            Dictionary of bools that state whether the FS>=acceptance_crit for each load
        

        '''
        if not loads:
            if not 'Ha' in self.loads:
                self.getLugForces()
            loads = self.loads
        if not self.anchorCapacity:
            self.getAnchorCapacity()
         
        # look for load dictionary key in capacity dictionary
        FS = {}
        acceptance = {}
        for Lkey,Lval in loads.items():
            for Ckey,Cval in self.anchorCapacity.items():
                if Lkey in Ckey:
                    if Lval == 0:
                        FS[Lkey] = float('inf')
                    else:
                        FS[Lkey] = Cval/Lval
                    if acceptance_crit and Lkey in acceptance_crit:
                        acceptance[Lkey] = acceptance_crit[Lkey]<=FS[Lkey]
                    
        if acceptance_crit:
            return(FS,acceptance[Lkey])
        else:
            return(FS)
            
        
               
    def getCost(self,costDict='default'):
        '''find costs of anchor and store in design dictionary
        
        Parameters
        ----------
        costDict : dictionary or yaml, optional
            Dictionary of various costs for anchors. Sub costs that can be included are: 
            material : material costs
            
        '''
        if isinstance(costDict,str) and costDict != 'default':
            import yaml
            costDict = yaml.load(costDict, Loader=yaml.FullLoader)
        anchType = self.dd['type']
        if costDict == 'default':
            matCostDict = {'DEA':5.705,'suction_pile':4.435,'gravity':1.905} # mean values from Task 49 Design Basis ranges
            instCostDict = {}
            decomCostDict = {}
        else:
            matCostDict = costDict['material']
            if 'install' in costDict:
                instCostDict = costDict['install']
            if 'decom' in costDict:
                decomCostDict = costDict['decom']
        keyFail = True
        # check if mass info is available
        if not self.mass:
            if 'soil_properties' in self.dd:               
                # need mass - call capacity functions
                self.getAnchorCapacity(plot=False)
            else:
                print('Soil properties needed to calculate anchor mass for cost. Setting cost to 0.')
                self.mass = 0
            
        # sort by type of anchor
        for Ckey,Cval in matCostDict.items():
            if anchType in Ckey:
                self.cost['materials'] = matCostDict[Ckey]*self.mass
                # self.cost['install'] = instCostDict[Ckey]
                # self.cost['decom'] = decomCostDict[Ckey]
                keyFail = False
        # raise error if anchType not found in cost dictionary
        if keyFail:
            raise KeyError(f'anchor type {anchType} not found in material cost dictionary')
            
        return(sum(self.cost.values()))

           
    def getSize(self,startGeom,loads=None,minfs={'Ha':1.8,'Va':1.8},inc_pct=10):
        '''
        

        Parameters
        ----------
        startGeom : dict
            Dictionary of required geometric values to start with. These will be increased by
            inc_pct for each run of the while loop
        loads : TYPE, optional
            DESCRIPTION. The default is None.
        inc_pct : float
            Percent to increase the geometric properties by in each iteration of the while loop

        Returns
        -------
        None.

        '''
        def fs_check(fs,minfs):
            '''
            Checks if all fs > minfs

            Parameters
            ----------
            fs : dict
                factors of safety
            minfs : dict
                minimum factors of safety

            Returns
            -------
            result: bool
                True if all fs > minfs, else False

            '''
            if fs['Ha'] > minfs['Ha'] and fs['Va'] > minfs['Va']:
                return(True)
            else:
                return(False)
            
        
        anchType = self.dd['type']
        if not loads:
            loads = self.loads
            
        if not 'Ha' in loads:
            self.getLugForces(mudloads=loads)
            
        if not 'zlug' in self.dd['design']:
            self.dd['design']['zlug'] = 0
        
        geom = self.dd['design']
        geom.update(startGeom)

        try:
            self.getAnchorCapacity(loads=loads,plot=False)
        except:
            for cKey in startGeom.keys():
                geom[cKey]*=2
            self.getAnchorCapacity(loads=loads,plot=False)
        fs = self.getFS(loads=loads)
        # check if it's already over the minimum factor of safety
        last_fs = fs_check(fs,minfs)
        signswitch = False
        if last_fs:
            # reduce sizing to below fs so we can increase slowly later
            while last_fs:
                for cKey in startGeom.keys():
                    geom[cKey] *= .75
                    self.getAnchorCapacity(loads=loads,plot=False)
        
            
            
        while not last_fs:
            
            for cKey in startGeom.keys():
                geom[cKey] *= (1+inc_pct/100)
                
            self.getAnchorCapacity(loads=loads,plot=False)
            fs = self.getFS(loads=loads)
            last_fs = fs_check(fs,minfs)
            

            
                
                
        return(geom,fs)
                
    def getSuctionSize(self,D,L,loads=None,minfs={'Ha':1.6,'Va':2},LD_con=5):
        '''
        

        Parameters
        ----------
        D : float
            Diameter of suction bucket
        L : float
            Length of suction bucket
        loads : TYPE, optional
            DESCRIPTION. The default is None.
        minfs : dict,optoinal
            Minimum factors of safety in horizontal and vertical directions
        LD_con : float
            Constraint for L/D parameter

        Returns
        -------
        None.

        '''
        from scipy.optimize import minimize
        anchType = self.dd['type']
        if not loads:
            loads = self.loads
            
        if not 'Ha' in loads:
            loads = self.getLugForces(mudloads=loads)
           
        loads['Ha'] = minfs['Ha']*loads['Ha']
        loads['Va'] = minfs['Va']*loads['Va']
                   
        if not 'zlug' in self.dd['design']:
            self.dd['design']['zlug'] = (2/3)*L
            
        # Define the objective function: Minimize |UC - 1| (aim for UC to be 1)
        def objective(vars):
            D, L = vars
            self.dd['design']['D'] = D
            self.dd['design']['L'] = L
            self.dd['design']['zlug'] = (2/3)*L
            results = self.getAnchorCapacity(plot=False)
            return abs(results['UC'] - 1)  
        
        def conFun(vars,LD_con):
            D, L = vars
            if L/D >= LD_con:
                conval = -1 
            else:
                conval = 1 
            
            return(conval)
        
        # Initial guess for D and L
        initial_guess = [D, L]       # Input values for D and L
        
        # Bounds for D and L (adjust as needed)
        bounds = [(1, 5), (5, 50)]   # Bounds for D and L
        
        # constraints
        constraints = [{'type':'ineq','fun':conFun,'args':(LD_con,)}]
        
        # Run the optimization to find D and L that satisfy UC close to 1
        solution = minimize(objective, initial_guess, bounds=bounds,method="COBYLA",
                            constraints=constraints)
        
        # Extract the optimized values of D and L
        self.dd['design']['D'], self.dd['design']['L'] = solution.x
        self.dd['design']['zlug'] = (2/3)*self.dd['design']['L'] 
        results = self.getAnchorCapacity(plot=False)    
          
       
            
        # # check if anchor loads are available
        # if not self.loads:
        #     # if not, check if theres a moorpy anchor object and calculate loads from that
        #     if self.mpAnchor:
        #         print("Need anchor loads to obtain cost, using getMPForces to determine loads in MoorPy")
        #         self.getLugForces()
        #     elif self.ms:
        #         print('Need anchor loads to obtain cost, creating a MoorPy anchor object and using getMPForces to determine loads in MoorPy')
        #         self.makeMoorPyAnchor(self.ms)
        #         self.getLugForces()
        #     else:
        #         raise Exception("Need anchor loads to obtain cost")
        # # check again if there are loads
        # if self.loads:
        #     c = self.dd['cost'] # set location for clarity
        #     # calculate individual costs and total cost for the anchor
        #     c['matCost'], c['instCost'], c['decomCost'] = mp.Point.getcost(self.mpAnchor)
        #     c['totCost'] = c['matCost'] + c['instCost'] + c['decomCost']

    
    # def getMass(self,uhc_mode):
    #     '''find mass and/or UHC of anchor from MoorProps and store in design dictionary
    #     Parameters
    #     ----------
    #     uhc_mode : boolean
    #         True : obtain UHC from mass
    #         False : obtain Masss and UHC from loads
    #     '''
    #     if uhc_mode: # if looking for UHC given mass
    #         if self.dd['design']['m']: # check anchor mass is given
    #             self.dd['design']['UHC'], self.dd['design']['m'], info = mp.MoorProps.getAnchorMass(uhc_mode=1, mass_int=self.dd['design']['m'], anchor=self.dd['type'], soil_type=self.anchorProps['soil_type'])
    #         else:
    #             raise Exception("Need anchor mass to calculate UHC when uhc_mode = 1")
    #     else: # if looking for mass and UHC given loads
    #         if self.loads: # check the loads section exists
    #             self.dd['design']['UHC'], self.dd['design']['m'], info = mp.MoorProps.getAnchorMass(uhc_mode=0, fx=self.loads['ff'], fz=self.loads['fz'], anchor=self.dd['type'],soil_type=self.dd['soil_type'],method=self.loads['method'])
    #         elif self.mpAnchor:
    #             print("Need anchor loads to obtain mass, using getMPForces to determine loads in MoorPy")
    #             self.getLugForces()
    #             self.dd['design']['UHC'], self.dd['design']['m'], info = mp.MoorProps.getAnchorMass(uhc_mode=0, fx=self.loads['ff'], fz=self.loads['fz'], anchor=self.dd['type'],soil_type=self.dd['soil_type'],method=self.loads['method'])
    #         elif self.ms:
    #             print('Need anchor loads to obtain mass, creating a MoorPy anchor object and using getMPForces to determine loads in MoorPy')
    #             self.makeMoorPyAnchor(self.ms)
    #             self.getLugForces()
    #             self.dd['design']['UHC'], self.dd['design']['m'], info = mp.MoorProps.getAnchorMass(uhc_mode=0, fx=self.loads['ff'], fz=self.loads['fz'], anchor=self.dd['type'],soil_type=self.dd['soil_type'],method=self.loads['method'])
    #         else:
    #             raise Exception("Need anchor loads to obtain mass")
