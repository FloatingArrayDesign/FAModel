"""Anchor class for FAModel, containing information and key methods for anchors of mooring lines
    Work in progress
"""
import moorpy as mp
import numpy as np
from scipy.optimize import minimize
from famodel.famodel_base import Node
from famodel.mooring.mooring import Mooring
import matplotlib.pyplot as plt
from collections import defaultdict
import famodel.platform.platform 
import shapely as sh

class Anchor(Node):
    
    def __init__(self, dd=None, ms=None, r=[0,0,0], aNum=None, id=None,
                 g=9.81, rho=1025, profile_map=None):
        '''
        Initialize an Anchor object.

        Parameters
        ----------
        dd : dict
            Design dictionary containing all information on the anchor.
        ms : MoorPy system object
            MoorPy system instance.
        r : list of float
            Anchor position coordinates (x, y, z) (m)
        aNum : int, optional
            Index in anchor list.
        id : str or int, optional
            Unique anchor identifier.
        g : float, optional
            Gravity.
        rho : float, optional
            Water density.
        profile_map : list of dict, optional
            Full soil profile map for selecting local soil layers.
        '''

        from famodel.famodel_base import Node
        Node.__init__(self, id)

        self.dd = dd
        self.ms = ms
        self.r = r
        self.aNum = aNum
        self.g = g
        self.rho = rho

        if dd and 'type' in dd:
            self.anchType = dd['type']
        else:
            self.anchType = 'suction'
            print(f"[Anchor] No type provided. Defaulting to 'suction'.")

        self.soil_type = None
        self.soil_profile = None
        self.profile_name = None
        self.soil_type_list = []

        self.mpAnchor = None
        self.capacity_format = None
        self.mass = dd.get('design', {}).get('mass', None) if dd else None
        self.point_num = 0  # initialize point number
        
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
        self.safety_factors = {} # calculated safety factor
        self.safety_factors_required = {} # minimum allowable safety factor
        
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
        self.loads = {}
        self.anchorCapacity = {}
        self.cost = {}
        self.failure_probability = {}
        self.env_impact = {}

        # Assign soil profile if map is provided
        if profile_map is not None:
            if len(profile_map) == 1:
                self.setSoilProfile(profile_map)
            elif len(profile_map) >= 4:
                self.interpolateSoilProfile(profile_map)
            else:
                raise ValueError("profile_map must contain either 1 or ≥4 CPTs for soil assignment.")

    def setSoilProfile(self, profile_map):
        '''
        Assign a soil profile directly from a single CPT.
        Assumes profile_map is a list with only one entry.
        '''
        if len(profile_map) != 1:
            raise ValueError("setSoilProfile expects a profile_map with exactly one CPT.")

        cpt = profile_map[0]
        self.soil_profile = cpt['layers']
        self.profile_name = cpt.get('name', 'CPT_Assigned')

        # Extract soil types from layers
        soil_types = [layer['soil_type'] for layer in self.soil_profile]
        self.soil_type_list = list(set(soil_types))
        self.soil_type = soil_types[0] if len(self.soil_type_list) == 1 else 'mixed'

        # Group layers by soil type
        soilProps = defaultdict(list)
        for layer in self.soil_profile:
            layer_copy = layer.copy()
            soil_type = layer_copy.pop('soil_type')
            soilProps[soil_type].append(layer_copy)
        self.soilProps = dict(soilProps)

        print(f"[Anchor] Assigned soil profile from {self.profile_name} with soil types {self.soil_type_list}")


    def interpolateSoilProfile(self, profile_map):
        '''
        Interpolate a soil profile from the 4 nearest CPTs in profile_map.
        '''
        if len(profile_map) < 4:
            raise ValueError("interpolateSoilProfile requires at least 4 CPTs.")

        x_anchor, y_anchor = self.r[0], self.r[1]

        # Sort CPTs by distance
        distances = [np.hypot(p['x'] - x_anchor, p['y'] - y_anchor) for p in profile_map]
        idx_sorted = np.argsort(distances)
        CPTs = [profile_map[i] for i in idx_sorted[:4]]

        # Inverse distance weighting
        x = np.array([cpt['x'] for cpt in CPTs])
        y = np.array([cpt['y'] for cpt in CPTs])
        d = np.hypot(x - x_anchor, y - y_anchor)
        w = 1 / np.maximum(d, 1e-3)**2
        w /= np.sum(w)

        # Interpolate layer-by-layer (assumes same layer structure)
        layers_list = [cpt['layers'] for cpt in CPTs]
        n_layers = len(layers_list[0])
        interpolated_layers = []

        for i in range(n_layers):
            base_layer = layers_list[0][i]
            layer = {'soil_type': base_layer['soil_type']}

            for key in base_layer:
                if key == 'soil_type':
                    continue
                if all(key in l[i] for l in layers_list):
                    vals = [l[i][key] for l in layers_list]
                    layer[key] = np.dot(w, vals)

            interpolated_layers.append(layer)

        self.soil_profile = interpolated_layers
        self.profile_name = "Interpolated_2D"

        # Extract soil types
        soil_types = [layer['soil_type'] for layer in self.soil_profile]
        self.soil_type_list = list(set(soil_types))
        self.soil_type = soil_types[0] if len(self.soil_type_list) == 1 else 'mixed'

        # Group interpolated layers by soil type
        soilProps = defaultdict(list)
        for layer in self.soil_profile:
            layer_copy = layer.copy()
            soil_type = layer_copy.pop('soil_type')
            soilProps[soil_type].append(layer_copy)
        self.soilProps = dict(soilProps)

        print(f"[Anchor] Interpolated soil profile: {self.profile_name} with soil types {self.soil_type_list}")

    def makeMoorPyAnchor(self, ms):
        '''
        Create a MoorPy anchor object in a MoorPy system.

        Parameters
        ----------
        ms : MoorPy system instance
            The MoorPy system to add the anchor to.

        Returns
        -------
        ms : MoorPy system instance
            The updated MoorPy system with the anchor added.
        '''       
        anchType = self.anchType or 'suction'

        # Create anchor as a fixed point in MoorPy system
        ms.addPoint(1, self.r)

        # Assign this point as mpAnchor in the anchor class instance
        self.mpAnchor = ms.pointList[-1]

        # Set mass if available
        if 'mass' in self.dd.get('design', {}):
            self.mpAnchor.m = self.dd['design']['mass']

        # Set diameter if available
        if 'd' in self.dd.get('design', {}):
            self.mpAnchor.d = self.dd['design']['d']
            
        # Set dummy design to get PointType from MoorPy
        design = {f"num_a_{anchType}": 1}
        pointType = ms.setPointType(design, source=None)
        self.mpAnchor.entity = pointType

        return ms
       
    def getLineProperties(self):
        '''
        Retrieve line_type, diameter and unit weight from attached mooring.

        Returns
        -------
        line_type : str
            Type of mooring line ('chain' or 'wire')
        d : float
            Nominal diameter (m)
        w : float
            Unit weight (N/m)
        '''
        for att in self.attachments.values():
            if isinstance(att['obj'], Mooring):
                mtype = att['obj'].dd['sections'][0]['type']['material'].lower()
                if 'chain' not in mtype:
                    print('No chain below seafloor, setting Ta=Tm (no load transfer).')
                    return mtype, None, None, True
                else:
                    d_nom = att['obj'].dd['sections'][0]['type']['d_nom']
                    w_nom = att['obj'].dd['sections'][0]['type']['w']
                    return 'chain', d_nom, w_nom, False
        raise ValueError('No mooring line attachment found for anchor.')

    def getMudlineForces(self, max_force=False, lines_only=False, seabed=True, xyz=False, project=None):
        '''
        Find forces on anchor at mudline using the platform.getWatchCircle method
        or the MoorPy Point.getForces method. Optionally computes the maximum force 
        based on platform excursion using the project's arrayWatchCircle method or
        the attached platform's getWatchCircle method.

        Parameters
        ----------
        max_force : bool, optional
            If True, computes the maximum expected force on the anchor 
            using platform excursion. Default is False.
        lines_only : bool, optional
            Calculate forces from just mooring lines (True) or not (False). Default is False.
        seabed : bool, optional
            Include effect of seabed pushing up the anchor (True) or not (False). Default is True.
        xyz : bool, optional
            Return forces in x, y, z DOFs (True) or only the enabled DOFs (False). Default is False.
        project : object, optional
            Project object that can run arrayWatchCircle(). Used only if max_force is True.

        Returns
        -------
        dict
            Dictionary containing mudline forces.
        '''
        Platform = famodel.platform.platform.Platform

        if max_force:
            if project:
                project.arrayWatchCircle()
            else:
                for att in self.attachments.values():
                    if isinstance(att['obj'], Mooring):
                        for attM in att['obj'].attached_to:
                            if isinstance(attM, Platform):
                                locx, locy, maxVals = attM.getWatchCircle()
                                Hm = np.sqrt(maxVals[0]**2 + maxVals[1]**2)
                                Vm = maxVals[2]
                                thetam = np.degrees(np.arctan2(Vm, Hm))
                                self.loads['Hm'] = Hm
                                self.loads['Vm'] = Vm
                                self.loads['thetam'] = thetam
                                self.loads['mudline_load_type'] = 'max_force'
                                break
        else:
            loads = self.mpAnchor.getForces(lines_only=lines_only, seabed=seabed, xyz=xyz)
            Hm = np.sqrt(loads[0]**2 + loads[1]**2)
            Vm = loads[2]
            thetam = np.degrees(np.arctan2(Vm, Hm))
            self.loads['Hm'] = Hm
            self.loads['Vm'] = Vm
            self.loads['thetam'] = thetam
            self.loads['mudline_load_type'] = 'current_state'

        self.loads['method'] = 'static'
        return self.loads
       
    def getLugForces(self, Hm, Vm, zlug, line_type=None, d=None, w=None, plot=True):
        '''
        Calculate the lug forces Ha and Va based on mudline loads using local soil profile.

        Parameters
        ----------
        Hm : float
            Horizontal mudline load (N)
        Vm : float
            Vertical mudline load (N)
        zlug : float
            Padeye embedment depth (m)
        line_type : str, optional
            Type of mooring line ('chain' or 'wire')
        d : float, optional
            Mooring line diameter (m)
        w : float, optional
            Mooring line unit weight (N/m)
        plot : bool, optional
            Whether to plot the load transfer profile

        Returns
        -------
        Ha : float
            Horizontal load at lug (N).
        Va : float
            Vertical load at lug (N).
        '''
        from .anchors_famodel.capacity_load import getTransferLoad
        from .anchors_famodel.support_plots import plot_load

        # Ensure soil profile is available
        if self.soil_profile is None or self.soil_type is None:
            raise ValueError("Anchor soil profile or soil type is not assigned. Use setSoilProfile first.")

        soil_profile = self.soil_profile
        soil_type = self.soil_type

        # Determine mudline depth
        z0 = soil_profile[0]['top']

        # Load transfer if padeye is embedded
        if zlug > z0:
            if line_type is None or d is None or w is None:
                try:
                    line_type, d, w = self.getLineProperties()
                except ValueError:
                    print('[Warning] No mooring attachment found. Trying anchor-level line properties...')
                    line_type = getattr(self, 'line_type', None)
                    d = getattr(self, 'd', None)
                    w = getattr(self, 'w', None)

                    if any(v is None for v in [line_type, d, w]):
                        print('[Fallback] Using default chain properties.')
                        line_type = 'chain'
                        d = 0.16
                        w = 5500.0

            layers, loads = getTransferLoad(
                profile_map=[{'layers': self.soil_profile}],
                Tm=np.sqrt(Hm**2 + Vm**2),
                thetam=np.degrees(np.arctan2(Vm, Hm)),
                zlug=zlug,
                line_type=line_type,
                d=d,
                w=w,
                plot=plot
            )

            Ta = loads['Ta']
            thetaa = loads['thetaa']
            Ha = Ta*np.cos(np.deg2rad(thetaa))
            Va = Ta*np.sin(np.deg2rad(thetaa))

        else:
            Ha = Hm
            Va = Vm
            
        if plot == True:
            plot_load(layers, loads['drag_values'], loads['depth_values'], 
                      loads['Tm'], loads['thetam'], loads['Ta'], 
                      loads['thetaa'], zlug=zlug)

        return layers, Ha, Va

    def getCapacityAnchor(self, Hm, Vm, zlug, line_type=None, d=None, w=None, mass_update=False, plot=False):
        '''
        Calculate anchor capacity based on anchor type and local soil profile.
    
        Parameters
        ----------
        Hm : float
            Horizontal mudline load (N)
        Vm : float
            Vertical mudline load (N)
        zlug : float
            Padeye embedment depth (m)
        line_type : str, optional
            Type of mooring line ('chain' or 'wire')
        d : float, optional
            Mooring line diameter (m)
        w : float, optional
            Mooring line unit weight (N/m)
        mass_update : bool, optional
            Whether to update the mass when is not assigned
        plot : bool, optional
            Whether to plot the load transfer and pile geometry
    
        Returns
        -------
        results : dict
            Capacity results dictionary from the selected capacity function.
        '''

        # - - - - set details - - - -
        anchType = self.dd['type'] 
        geom = self.dd['design']# geometric anchor information

        if not ground_cons:
            soil = next(iter(self.soilProps.keys()), None) # soil type
            ground_conds = self.soilProps[soil]
        else:
            soil = next(iter(ground_cons.keys()))
            ground_conds = ground_cons[soil]
         
        for key,prop in ground_conds.items():
            if isinstance(prop,list) or isinstance(prop,np.ndarray):
                if len(prop)>1:
                    print('Warning: Only homogeneous soils are supported at this time. Only the first item in a property list will be used.')
                    break 
                else:
                    ground_conds[key] = prop[0]
        
        
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
                            # calculate thetaa from Ha and Va
                            loads['thetaa'] = np.arctan2(loads['Va'],loads['Ha'])
                            # loads = self.getLugForces(plot=plot)
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
        self.mass = results['Weight']*1000/self.g # mass in [kg]
            
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
    
    def getLugForces(self, mudloads=None, max_force=True, plot=False):
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
                soil = next(iter(self.soilProps.keys()), None)
                ground_conds = self.soilProps[soil]
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
        
        return mudloads
    
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
        acceptance_margin : dict
            Dictionary of difference between FS and acceptance criteria for each load type
        

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
        acceptance_margin = {}
        for Lkey,Lval in loads.items():
            for Ckey,Cval in self.anchorCapacity.items():
                if Lkey in Ckey:
                    if Lval == 0:
                        FS[Lkey] = float('inf')
                    else:
                        FS[Lkey] = Cval/Lval
                    if acceptance_crit and Lkey in acceptance_crit:
                        if Lval == 0 or acceptance_crit[Lkey] == 0:
                            acceptance[Lkey] = True
                        else:
                            acceptance[Lkey] = acceptance_crit[Lkey]<=FS[Lkey]
                            acceptance_margin[Lkey] = FS[Lkey] - acceptance_crit[Lkey]
                    
        if acceptance_crit:
            return(FS,acceptance,acceptance_margin)
        else:
            return(FS)
            
    def makeBuffer(self, buff_rad=50):
        point = sh.Point(self.r[:2])
        buff = point.buffer(buff_rad)
        return buff
               
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
           
        capacity_dispatch = {
            'suction': getCapacitySuction,
            'sepla': getCapacityPlate,
            'dea': getCapacityPlate,
            'depla': getCapacityPlate,
            'vla': getCapacityPlate,
            'plate': getCapacityPlate,
            'torpedo': getCapacityTorpedo,
            'helical': getCapacityHelical,
            'driven': getCapacityDriven,
            'dandg': getCapacityDandG
        }
    
        print(f'[Debug] mass_update = {mass_update}')
        anchType_clean = self.anchType.lower().replace(' ', '')
        capacity_func = capacity_dispatch.get(anchType_clean)
        if capacity_func is None:
            raise ValueError(f"Unknown anchor type '{self.anchType}' for anchor capacity calculation.")
    
        if self.soil_profile is None or self.soil_type is None:
            raise ValueError("Soil profile or soil type not set for this anchor.")
    
        soil_profile = self.soil_profile
        soil_type = self.soil_type
        z0 = soil_profile[0]['top']
    
        # Load transfer if padeye is embedded below mudline
        
        if line_type is None or d is None or w is None:
            try:
                line_type, d, w = self.getLineProperties()
            except ValueError:
                print('[Warning] No mooring attachment found. Trying anchor-level line properties...')
                line_type = getattr(self, 'line_type', None)
                d = getattr(self, 'd', None)
                w = getattr(self, 'w', None)

                if any(v is None for v in [line_type, d, w]):
                    print('[Fallback] Using default chain properties.')
                    line_type = 'chain'
                    d = 0.16
                    w = 5500.0
    
        if zlug > z0:
            layers, Ha, Va = self.getLugForces(
                Hm, Vm,
                zlug=zlug,
                line_type=line_type,
                d=d,
                w=w,
                plot=False
            )

            Ta = np.sqrt(Ha**2 + Va**2)
            thetaa = np.degrees(np.arctan2(Va, Ha))
            
            print(f'Input Hm = {Hm}, Vm = {Vm}, zlug = {zlug}')
            print(f'Output Ha = {Ha}, Va = {Va}, zlug = {zlug}')
            print(f'Output Ta = {Ta}, thetaa = {(thetaa)}')
            print(f"[Branch Check] Entered {'zlug>z0' if zlug>z0 else 'else'} for anchor {self.anchType}")

        else:
            Ha = Hm
            Va = Vm
            Ta = np.sqrt(Ha**2 + Va**2)
            thetaa = np.degrees(np.arctan2(Va, Ha))
            print(f'[Direct assign] Ha = {Ha}, Va = {Va}, Ta = {Ta}, thetaa = {thetaa}')
            print(f"[Branch Check] Entered {'zlug>z0' if zlug>z0 else 'else'} for anchor {self.anchType}")


    
        # --- Call the appropriate capacity function ---
        if anchType_clean in ['sepla', 'dea', 'depla', 'vla', 'plate']:
            self.capacity_format = 'plate'
            B = self.dd['design']['B']
            L = self.dd['design']['L']
            print(f"[Final Check] Ha = {Ha}, Va = {Va}, anchor = {self.anchType}")
            beta = 90.0 - np.degrees(np.arctan2(Va, Ha))
            self.dd['design']['beta'] = beta 
            layers, results = capacity_func(
                profile_map=[{'name': self.profile_name, 'layers': self.soil_profile}],
                location_name=self.profile_name,
                B=B, L=L, zlug=zlug,
                beta=beta,
                Ha=Ha, Va=Va,
                plot=plot
            )
            
        elif anchType_clean == 'suction':
            self.capacity_format = 'envelope'
            D = self.dd['design']['D']
            L = self.dd['design']['L']
            zlug = self.dd['design']['zlug']
            layers, results = capacity_func(
                profile_map=[{'name': self.profile_name, 'layers': self.soil_profile}],
                location_name=self.profile_name,
                D=D, L=L, zlug=zlug,
                Ha=Ha, Va=Va,
                thetalug=5, psilug=7.5,
                plot=plot
            )
            
        elif anchType_clean == 'torpedo':
            self.capacity_format = 'envelope'
            D1 = self.dd['design']['D1']
            D2 = self.dd['design']['D2']
            L1 = self.dd['design']['L1']
            L2 = self.dd['design']['L2']
            ballast = self.dd['design'].get('ballast', 0.0)
            layers, results = capacity_func(
                profile_map=[{'name': self.profile_name, 'layers': self.soil_profile}],
                location_name=self.profile_name,
                D1=D1, D2=D2, L1=L1, L2=L2,
                zlug=zlug,
                ballast=ballast,
                Ha=Ha, Va=Va,
                plot=plot
            )

        elif anchType_clean == 'helical':
            self.capacity_format = 'component'
            D = self.dd['design']['D']     
            L = self.dd['design']['L']     
            d = self.dd['design']['d']     
            zlug = self.dd['design']['zlug']
            layers, results = capacity_func(
                profile_map=[{'name': self.profile_name, 'layers': self.soil_profile}],
                location_name=self.profile_name,
                D=D, L=L, d=d,
                zlug=zlug,
                Ha=Ha, Va=Va,
                plot=plot
            )

        elif anchType_clean in ['driven', 'dandg']:
            self.capacity_format = 'component'
            L = self.dd['design']['L']
            D = self.dd['design']['D']
            zlug = self.dd['design']['zlug']
            layers, y, z, results = capacity_func(
                profile_map=[{'name': self.profile_name, 'layers': self.soil_profile}],
                location_name=self.profile_name,
                L=L, D=D, zlug=zlug,
                Ha=Ha, Va=Va,
                plot=plot
            )
    
        else:
            raise ValueError(f"Anchor type '{self.anchType}' not supported.")
    
        # --- Store results ---
        self.anchorCapacity = {
            'Hmax': results.get('Horizontal max.', np.nan),
            'Vmax': results.get('Vertical max.', np.nan),
            'Ha': Ha,
            'Va': Va,
            'zlug': zlug,
            'z0': z0}
        
        # Correct UC format
        if anchType_clean in ['suction', 'torpedo', 'plate', 'sepla', 'dea', 'depla', 'vla']:
            self.anchorCapacity['UC'] = results.get('Unity check', np.nan)
        
        elif anchType_clean in ['helical', 'driven', 'dandg']:
            self.anchorCapacity['Unity check (horizontal)'] = results.get('Unity check (horizontal)', np.nan)
            self.anchorCapacity['Unity check (vertical)'] = results.get('Unity check (vertical)', np.nan)
        
        # Copy over lateral and rotational displacements
        if 'Lateral displacement' in results:
            self.anchorCapacity['Lateral displacement'] = results['Lateral displacement']
        if 'Rotational displacement' in results:
            self.anchorCapacity['Rotational displacement'] = results['Rotational displacement']
        
        # Weight calculated via dimensions
        if mass_update == False:
            if 'Weight pile' in results:
                self.anchorCapacity['Weight pile'] = results['Weight pile']
            if 'Weight plate' in results:
                self.anchorCapacity['Weight plate'] = results['Weight plate']
        else:
            if 'Weight pile' in results:
                if self.mass is None:
                    self.mass = results['Weight pile']/self.g
                self.anchorCapacity['Weight pile'] = self.mass*self.g
            if 'Weight plate' in results:
                if self.mass is None:
                    self.mass = results['Weight plate']/self.g
                self.anchorCapacity['Weight plate'] = self.mass*self.g
                
        # print(f"[DEBUG] Stored Lateral displacement in anchorCapacity: {self.anchorCapacity['Lateral displacement']:.6f}")
         
    def getSizeAnchor(self, geom, geomKeys, geomBounds=None, loads=None,
                      lambdap_con=[4, 8], zlug_fix=True, safety_factor={'SF_combined': 1.0}, plot=False):
        '''
        Generalized optimization method for all anchor types, using dictionary-based safety factors.
        '''
    
        anchType_clean = self.dd['type'].lower().replace('', '')
    
        if loads is None:
            loads = self.loads
    
        Hm = loads['Hm']
        Vm = loads['Vm']
    
        line_type = getattr(self, 'line_type', 'chain')
        d = getattr(self, 'd', 0.16)
        w = getattr(self, 'w', 5000.0)
    
        def update_zlug():
            if anchType_clean == 'suction' and not zlug_fix and 'zlug' not in geomKeys:
                self.dd['design']['zlug'] = (2/3)*self.dd['design']['L']
            elif anchType_clean in ['driven', 'helical'] and not zlug_fix:    
                ratio = self.dd['design'].get('zlug_ratio', self.dd['design']['zlug']/self.dd['design']['L'])
                self.dd['design']['zlug_ratio'] = ratio
                self.dd['design']['zlug'] = ratio*self.dd['design']['L']
    
        def get_lambda():
            if anchType_clean == 'torpedo':
                L = self.dd['design']['L1'] + self.dd['design']['L2']
                A_wing = (self.dd['design']['D1'] - self.dd['design']['D2']) * self.dd['design']['L1']
                A_shaft = self.dd['design']['D2'] * L
                D = (A_wing + A_shaft) / L
            elif anchType_clean in ['driven', 'dandg', 'helical', 'suction']:
                L = self.dd['design']['L']
                D = self.dd['design']['D']
            elif anchType_clean in ['plate', 'sepla', 'dea', 'depla', 'vla']:
                L = self.dd['design']['L']
                D = self.dd['design']['B']
            else:
                raise ValueError(f'lambda not defined for anchor type: {anchType_clean}')
            return L/D
    
        def constraint_lambda_min(vars):
            return get_lambda() - lambdap_con[0]
    
        def constraint_lambda_max(vars):
            return lambdap_con[1] - get_lambda()
    
        if anchType_clean in ['suction', 'torpedo', 'plate', 'sepla', 'dea', 'depla', 'vla']:
            target_UC = 1.0/safety_factor.get('SF_combined', 1.0)
    
            def objective_uc(vars):
                for i, key in enumerate(geomKeys):
                    self.dd['design'][key] = vars[i]
                update_zlug()
                self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                                       line_type=line_type, d=d, w=w, mass_update=True, plot=False)
                UC = self.anchorCapacity.get('UC', 2.0)
                return (UC - target_UC)**2
    
            def constraint_uc_envelope(vars):
                return self.anchorCapacity.get('UC', 0.0) - target_UC
    
            constraints_uc = [
                {'type': 'ineq', 'fun': constraint_lambda_min},
                {'type': 'ineq', 'fun': constraint_lambda_max},
                {'type': 'ineq', 'fun': constraint_uc_envelope},
            ]
    
            result_uc = minimize(
                objective_uc,
                geom,
                method='COBYLA',
                bounds=geomBounds if geomBounds else None,
                constraints=constraints_uc,
                options={'rhobeg': 0.1, 'catol': 0.01, 'maxiter': 500}
            )
    
            endGeom = dict(zip(geomKeys, result_uc.x))
            self.dd['design'].update(endGeom)
            update_zlug()
            self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                                   line_type=line_type, d=d, w=w, mass_update=True, plot=plot)
    
            print('\nFinal Optimized Anchor (UC-based):')
            print('Design:', self.dd['design'])
            print('Capacity Results:', self.anchorCapacity)
            return
    

        def termination_condition():
            UC_h = self.anchorCapacity['Ha'] / self.anchorCapacity['Hmax']
            UC_v = self.anchorCapacity['Va'] / self.anchorCapacity['Vmax']
            disp_lat = abs(self.anchorCapacity.get('Lateral displacement', 0.0))
            disp_rot = abs(self.anchorCapacity.get('Rotational displacement', 0.0))
            limit_lat = 0.05 * self.dd['design']['D']
            limit_rot = 5.0
    
            if UC_h <= 1.0 and UC_v <= 1.0 and disp_lat <= limit_lat and disp_rot <= limit_rot:
                print('[Termination Condition Met] All four limits satisfied.')
                return 'terminate'
    
            return 'continue'
    
        def is_valid(value):
            return np.isfinite(value) and not np.isnan(value) and abs(value) < 1e6
        
        if anchType_clean in ['helical', 'driven', 'dandg']:
            L0, D0 = geom if len(geom) == 2 else [5.0, 1.0]
            self.dd['design']['L'] = L0
            self.dd['design']['D'] = D0
            #self.dd['design']['t'] = max(0.05, 0.1 * D0)
            update_zlug()
            self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                                   line_type=line_type, d=d, w=w, mass_update=True, plot=False)
        
            UC_h = self.anchorCapacity['Ha']/self.anchorCapacity['Hmax']
            UC_v = self.anchorCapacity['Va']/self.anchorCapacity['Vmax']
            disp_lat = abs(self.anchorCapacity.get('Lateral displacement', 0.0))
            disp_rot = abs(self.anchorCapacity.get('Rotational displacement', 0.0))
            limit_disp = 0.05*D0
            limit_rot = 5.0
            direction = 'shrink' if (UC_h <= 1.0 and UC_v <= 1.0 and disp_lat <= limit_disp and disp_rot <= limit_rot) else 'grow'
    
            max_iter = 200
            iter_count = 0
    
            if direction == 'shrink':
                for D in np.arange(D0, 0.49, -0.05):
                    self.dd['design']['D'] = D
                    #self.dd['design']['t'] = max(0.02, 0.1*D)
                    for L in np.arange(L0, 1.95, -0.25):
                        self.dd['design']['L'] = L
                        update_zlug()
                        self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                                               line_type=line_type, d=d, w=w, mass_update=True, plot=False)
                        UC_h = self.anchorCapacity['Ha']/self.anchorCapacity['Hmax']
                        UC_v = self.anchorCapacity['Va']/self.anchorCapacity['Vmax']
                        disp_lat = abs(self.anchorCapacity.get('Lateral displacement', 0.0))
                        disp_rot = abs(self.anchorCapacity.get('Rotational displacement', 0.0))
                        print(f'[Iter {iter_count}] L={L:.2f}, D={D:.2f}, UC_h={UC_h:.3f}, UC_v={UC_v:.3f}, lat={disp_lat:.3f} m, rot={disp_rot:.3f} deg')
                        iter_count += 1
                        if not all(is_valid(v) for v in [UC_h, UC_v, disp_lat, disp_rot]):
                            continue
                        if termination_condition():
                            print(f'\nTermination criteria met.')
                            print('Design:', self.dd['design'])
                            print('Capacity Results:', self.anchorCapacity)
                            return
            else:
                for D in np.arange(D0, 3.05, 0.05):
                    self.dd['design']['D'] = D
                    #self.dd['design']['t'] = max(0.02, 0.1*D)
                    for L in np.arange(L0, 50.25, 0.25):
                        self.dd['design']['L'] = L
                        update_zlug()
                        self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                                               line_type=line_type, d=d, w=w, mass_update=True, plot=False)
                        UC_h = self.anchorCapacity['Ha']/self.anchorCapacity['Hmax']
                        UC_v = self.anchorCapacity['Va']/self.anchorCapacity['Vmax']
                        disp_lat = abs(self.anchorCapacity.get('Lateral displacement', 0.0))
                        disp_rot = abs(self.anchorCapacity.get('Rotational displacement', 0.0))
                        print(f'[Iter {iter_count}] L={L:.2f}, D={D:.2f}, UC_h={UC_h:.3f}, UC_v={UC_v:.3f}, lat={disp_lat:.3f} m, rot={disp_rot:.3f} deg')
                        iter_count += 1
                        status = termination_condition()
                        if status == 'terminate':
                            print(f'Termination criteria met.')
                            print('Design:', self.dd['design'])
                            print('Capacity Results:', self.anchorCapacity)
                            return
                        elif status == 'continue':
                            continue
                    status = termination_condition()
                    if status == 'terminate':
                        print(f'\nTermination criteria met.')
                        print('Design:', self.dd['design'])
                        print('Capacity Results:', self.anchorCapacity)
                        return
    
            print('[Warning] While-loop search reached bounds without meeting criteria.')
    
        else:
            raise ValueError(f"Anchor type '{anchType_clean}' not supported for safety factor input.")
                   
    def getSizeAnchor2(self, geom, geomBounds=None, loads=None,
                       lambdap_con=[3, 6], zlug_fix=True, safety_factor={'SF_combined': 1.0}, plot=False):
        '''
        Grid-based optimization method for envelope anchors (suction, torpedo, plate).
        Evaluates UC over a grid of L and D, and selects the point closest to target UC.
        '''
        import matplotlib.pyplot as plt
        from matplotlib import cm
        import matplotlib.colors as mcolor
        import numpy as np

        anchType_clean = self.dd['type'].lower().replace('', '')

        if loads is None:
            loads = self.loads

        Hm = loads['Hm']
        Vm = loads['Vm']

        line_type = getattr(self, 'line_type', 'chain')
        d = getattr(self, 'd', 0.16)
        w = getattr(self, 'w', 5000.0)

        if anchType_clean not in ['suction', 'torpedo', 'plate']:
            raise ValueError(f"Grid-based getSizeAnchor only supports envelope anchors, not '{anchType_clean}'")

        UC_target = 1.0/safety_factor.get('SF_combined', 1.0)

        # Unpack bounds and generate grid
        L_vals = np.linspace(geomBounds[0][0], geomBounds[0][1], 10)
        D_vals = np.linspace(geomBounds[1][0], geomBounds[1][1], 10)

        L_grid, D_grid = np.meshgrid(L_vals, D_vals)
        UC_grid = np.full_like(L_grid, np.nan, dtype=float)
        mask = np.full_like(L_grid, False, dtype=bool)

        best_UC, best_L, best_D = None, None, None
        results = []

        for i in range(D_grid.shape[0]):  # loop over D
            for j in range(D_grid.shape[1]):  # loop over L
                D = D_grid[i, j]
                L = L_grid[i, j]
                lambdap = L/D

                if not (lambdap_con[0] <= lambdap <= lambdap_con[1]):
                    continue

                mask[i, j] = True
                self.dd['design']['L'] = L
                self.dd['design']['D'] = D

                if anchType_clean == 'suction' and not zlug_fix:
                    self.dd['design']['zlug'] = (2/3)*L

                try:
                    self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                                           line_type=line_type, d=d, w=w, 
                                           mass_update=True, plot=False)
                    UC = self.anchorCapacity.get('UC', np.nan)
                    results.append({
                        'L': L,
                        'D': D,
                        'UC': UC})

                    if UC > 1e-2 and UC < 10.0:
                        UC_grid[i, j] = UC
                        # Find UC closest to target
                        if best_UC is None or abs(UC - UC_target) < abs(best_UC - UC_target):
                            best_UC = UC
                            best_L = L
                            best_D = D

                except:
                    continue

        # Update best result
        # if best_L is not None and best_D is not None:
        self.dd['design']['L'] = best_L
        self.dd['design']['D'] = best_D
        if anchType_clean == 'suction' and not zlug_fix:
            self.dd['design']['zlug'] = (2/3)*best_L

        self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                               line_type=line_type, d=d, w=w, 
                               mass_update=True, plot=plot)

        print('\nFinal Optimized Anchor (Grid-based):')
        print('Design:', self.dd['design'])
        print('Capacity Results:', self.anchorCapacity)

        # else:
        #     print('[Warning] No valid combination found in the grid.')

        # Optional plot

        if plot:
            fig, ax = plt.subplots(figsize=(6, 8))
            vmin, vmax = 0.01, 10
            levels = np.logspace(np.log10(vmin), np.log10(vmax), 21)
            cp = ax.contourf(D_grid, L_grid, UC_grid, levels=levels, cmap='coolwarm', norm=mcolor.LogNorm(vmin=vmin, vmax=vmax))
            fig.colorbar(cp, ax=ax, label='Unity check (UC)')
            ax.contour(D_grid, L_grid, UC_grid, levels=levels, colors='k', linewidths=0.3, alpha=0.3)
            ax.contour(D_grid, L_grid, UC_grid, levels=[1.0], colors='red', linewidths=2, linestyles='--')
            ax.set_xlabel('Diameter (m)')
            ax.set_ylabel('Length (m)')
            ax.set_title('Unity Check (UC')
            ax.plot(best_D, best_L, 'ro', label='Best match')
            ax.annotate('Best match', (best_D, best_L), textcoords="offset points", xytext=(10,10), ha='center', color='red')
            ax.legend()
            plt.grid(True)
            plt.tight_layout()
            plt.show()
            
        #UC_target = 1.0
        closest = min(results, key=lambda x: abs(x['UC'] - UC_target))
        print("Closest to UC_target:")
        print(closest)
            
        return results
    
    def getSizeAnchor_BO(self, 
                         geom=[10.0, 2.0],
                         geomKeys=['L', 'D'],
                         geomBounds=[(5.0, 15.0), (1.0, 4.0)],
                         loads=None,
                         lambdap_con=[3, 6],
                         zlug_fix=False,
                         safety_factor={'SF_combined': 1.0},
                         n_calls=25,
                         plot=False,
                         verbose=True):
        '''
        Bayesian optimization to find (D, L) for UC closest to UC_target.
        Uses scikit-optimize for surrogate model and efficient sampling.
        '''
        from skopt import gp_minimize
        from skopt.space import Real
        from skopt.utils import use_named_args
        import numpy as np

        if loads is None:
            loads = self.loads

        Hm = loads['Hm']
        Vm = loads['Vm']

        line_type = getattr(self, 'line_type', 'chain')
        d = getattr(self, 'd', 0.16)
        w = getattr(self, 'w', 5000.0)

        UC_target = 1.0 / safety_factor.get('SF_combined', 1.0)

        # Define the search space
        space  = [
            Real(geomBounds[1][0], geomBounds[1][1], name='D'),
            Real(geomBounds[0][0], geomBounds[0][1], name='L')
        ]

        @use_named_args(space)
        def objective(**params):
            D = params['D']
            L = params['L']

            # Apply lambda constraint
            lambdap = L/D
            if not (lambdap_con[0] <= lambdap <= lambdap_con[1]):
                return 100.0

            self.dd['design']['D'] = D
            self.dd['design']['L'] = L
            if not zlug_fix:
                self.dd['design']['zlug'] = (2/3)*L

            try:
                self.getCapacityAnchor(
                    Hm=Hm,
                    Vm=Vm,
                    zlug=self.dd['design']['zlug'],
                    line_type=line_type,
                    d=d,
                    w=w,
                    mass_update=True,
                    plot=False
                )
                UC = self.anchorCapacity.get('UC', np.nan)
            except:
                UC = np.nan

            if verbose:
                print(f"Evaluated D={D:.3f}, L={L:.3f} -> UC={UC:.3f}")

            if not np.isfinite(UC):
                return 100.0

            if UC < UC_target:
                return (UC_target - UC)**2 * 0.5  # less penalty for overdesign
            else:
                return (UC - UC_target)**2 * 10   # higher penalty for failure

        # Run Bayesian optimization
        res = gp_minimize(
            objective,
            space,
            x0=[geom[1], geom[0]],
            n_calls=n_calls,
            random_state=42,
            verbose=verbose
        )

        # Best result
        best_D, best_L = res.x
        self.dd['design']['D'] = best_D
        self.dd['design']['L'] = best_L
        if not zlug_fix:
            self.dd['design']['zlug'] = (2/3)*best_L

        self.getCapacityAnchor(
            Hm=Hm,
            Vm=Vm,
            zlug=self.dd['design']['zlug'],
            line_type=line_type,
            d=d,
            w=w,
            mass_update=True,
            plot=plot
        )
        UC = self.anchorCapacity.get('UC', np.nan)

        print('\nBayesian Optimized Anchor:')
        print('Design:', self.dd['design'])
        print('Capacity Results:', self.anchorCapacity)
        print(f'Best UC: {UC:.4f} (target: {UC_target})')

        results = {'D': best_D, 'L': best_L, 'UC': UC, 'result': res}

        return results
    # PATCH for GRADIENT method: wrap getCapacityAnchor in safe evaluator
    def safe_get_uc(self, Hm, Vm, zlug, line_type, d, w, verbose=False):
        try:
            self.getCapacityAnchor(Hm, Vm, zlug, line_type, d, w, True, False)
            return self.anchorCapacity.get('UC', np.nan)
        except Exception as e:
            if verbose:
                print(f"[Safe Error] {str(e)}")
            return np.nan

    def getSizeAnchor_gradient(self, 
                               geom=[10.0, 2.0],
                               geomKeys=['L', 'D'],
                               geomBounds=[(5.0, 15.0), (1.0, 4.0)],
                               loads=None,
                               lambdap_con=[3, 6],
                               zlug_fix=False,
                               safety_factor={'SF_combined': 1.0},
                               step_size=0.2,
                               tol=0.05,
                               max_iter=30,
                               verbose=True):
        '''
        Gradient-based optimization with early stopping to match UC_target.
        '''
        import numpy as np
    
        if loads is None:
            loads = self.loads
    
        Hm = loads['Hm']
        Vm = loads['Vm']
    
        line_type = getattr(self, 'line_type', 'chain')
        d = getattr(self, 'd', 0.16)
        w = getattr(self, 'w', 5000.0)
    
        UC_target = 1.0 / safety_factor.get('SF_combined', 1.0)
    
        L, D = geom
    
        for iter in range(max_iter):
            lambdap = L / D
            if not (lambdap_con[0] <= lambdap <= lambdap_con[1]):
                if verbose:
                    print(f"[Iter {iter}] λ = {lambdap:.2f} out of bounds. Terminating.")
                break
    
            self.dd['design']['L'] = L
            self.dd['design']['D'] = D
            if not zlug_fix:
                self.dd['design']['zlug'] = (2/3)*L
    
            UC0 = self.safe_get_uc(Hm, Vm, self.dd['design']['zlug'], line_type, d, w, verbose=verbose)
    
            if not np.isfinite(UC0):
                break
    
            if verbose:
                print(f"[Iter {iter}] L={L:.2f}, D={D:.2f}, UC={UC0:.3f}")
    
            if abs(UC0 - UC_target) < tol:
                print("Early stopping: UC within tolerance.")
                break
    
            # Gradient estimate
            delta = 0.1
            UC_L = self.safe_get_uc(Hm, Vm, (2/3)*(L + delta), line_type, d, w, verbose=verbose)
            UC_D = self.safe_get_uc(Hm, Vm, (2/3)*L, line_type, d, w, verbose=verbose)
    
            grad_L = (UC_L - UC0)/delta if np.isfinite(UC_L) else 0.0
            grad_D = (UC_D - UC0)/delta if np.isfinite(UC_D) else 0.0
    
            # Update
            L -= step_size * grad_L
            D -= step_size * grad_D
            L = np.clip(L, geomBounds[0][0], geomBounds[0][1])
            D = np.clip(D, geomBounds[1][0], geomBounds[1][1])
    
            if not (lambdap_con[0] <= L/D <= lambdap_con[1]):
                if verbose:
                    print("Terminated: lambda constraint violated after update.")
                break
    
        self.dd['design']['L'] = L
        self.dd['design']['D'] = D
        self.dd['design']['zlug'] = (2/3)*L
        self.getCapacityAnchor(Hm, Vm, self.dd['design']['zlug'], line_type, d, w, True, True)
    
        print('\nGradient Optimized Anchor:')
        print('Design:', self.dd['design'])
        print('Capacity Results:', self.anchorCapacity)
    
        return {'D': D, 'L': L, 'UC': self.anchorCapacity.get('UC', np.nan)}
   
    def getSafetyFactor(self):
        '''
        Calculate the safety factor based on the unity checks stored in capacity results.

        Returns
        -------
        dict
            Dictionary containing safety factors.
        '''

        # - - - - Objective and Constraint Functions 
        
        # Define the objective function: Minimize weight of anchor (cost is dependent on weight)
        def objective(vars, args):

            geomKeys = args['geomKeys']
            input_loads = args['input_loads']
            fix_zlug = args['fix_zlug']

            newGeom = dict(zip(geomKeys,vars))
            self.dd['design'].update(newGeom)
            if 'suction' in self.dd['type'] and not fix_zlug:
                self.dd['design']['zlug'] = (2/3)*newGeom['L']
            
            if 'Hm' in input_loads or 'Vm' in input_loads:
                anchor_loads = self.getLugForces(mudloads=input_loads)
                input_loads = dict(Ha=anchor_loads['Ha'], Va=anchor_loads['Va'])    # overwrite the input_loads dictionary
            # get results
            results = self.getAnchorCapacity(loads=input_loads, plot=False)
                
            return(results['Weight'])
        
        # constraint for suction bucket sizing only. May add more constraints for other anchors in the future...
        def conFun_LD(vars, geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs):
            newGeom = dict(zip(geomKeys, vars))
            self.dd['design'].update(newGeom)

            if 'Hm' in input_loads or 'Vm' in input_loads:
                anchor_loads = self.getLugForces(mudloads=input_loads)
                input_loads = dict(Ha=anchor_loads['Ha'], Va=anchor_loads['Va'])    # overwrite the input_loads dictionary

            results = self.getAnchorCapacity(loads=input_loads, plot=False)
            
            convalA = newGeom['L']/newGeom['D'] - LD_con[0]
            convalB = LD_con[1] - newGeom['L']/newGeom['D']
            conval = min([convalA,convalB])
            # if newGeom['L']/newGeom['D'] >= LD_con[0] and newGeom['L']/newGeom['D'] <= LD_con[1]:
            #     conval = 1 
            # else:
            #     conval = -1 
            
            return(conval)
        # constraint to ensure unity check > 1 for suction buckets
        def conFun_Suction(vars, geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs):
            if 'Hm' in input_loads or 'Vm' in input_loads:
                anchor_loads = self.getLugForces(mudloads=input_loads)
                input_loads = dict(Ha=anchor_loads['Ha'], Va=anchor_loads['Va'])    # overwrite the input_loads dictionary
            results = self.getAnchorCapacity(loads=input_loads, plot=False)
            #conval = results['UC'] - 1
            conval = 1 - results['UC']
            # convalB = 1 - results['UC'] 
            return(conval)
        
        def conFun_DandG(vars, geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs):

            newGeom = dict(zip(geomKeys, vars))
            self.dd['design'].update(newGeom)
            if 'Hm' in input_loads or 'Vm' in input_loads:
                anchor_loads = self.getLugForces(mudloads=input_loads)
                input_loads = dict(Ha=anchor_loads['Ha'], Va=anchor_loads['Va'])    # overwrite the input_loads dictionary
            results = self.getAnchorCapacity(loads=input_loads, plot=False)

            return np.array([0.05*newGeom['D'] - results['Lateral displacement'] , 0.25 - results['Rotational displacement']])
            
        def conFunH(vars, geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs):
            # if 'suction' in self.dd['type']:
            #     results = self.getAnchorCapacity(plot=False)
            #     conval = results['UC'] - 1
            #     # if results['UC'] < 1:
            #     #     conval = -1*(results['UC'])
            # else:
            if 'Hm' in input_loads or 'Vm' in input_loads:
                anchor_loads = self.getLugForces(mudloads=input_loads)
                input_loads = dict(Ha=anchor_loads['Ha'], Va=anchor_loads['Va'])    # overwrite the input_loads dictionary
                minfs = dict(Ha=minfs['Hm'], Va=minfs['Vm'])
            FS, _, _ = self.getFS(loads=input_loads, acceptance_crit=minfs)
            conval = FS['Ha'] - 1
                # for key,val in FS.items():
                    
                #     if val/minfs[key]<1:
                #         if -1*(1-val/minfs[key]) < conval:
                #             conval = -1*(1-val/minfs[key])
            return(conval)
        
        def conFunV(vars, geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs):
            if 'Hm' in input_loads or 'Vm' in input_loads:
                anchor_loads = self.getLugForces(mudloads=input_loads)
                input_loads = dict(Ha=anchor_loads['Ha'], Va=anchor_loads['Va'])    # overwrite the input_loads dictionary
                minfs = dict(Ha=minfs['Hm'], Va=minfs['Vm'])
            FS, _, _ = self.getFS(loads=input_loads, acceptance_crit=minfs)
            # special case for DEAs
            if minfs['Va'] == 0:
                 conval = 1
            else:
                conval = FS['Va'] - 1
                
            # print('FS_V',FS['Va'])
            return(conval)
        
        def conBounds(vars, geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs):

            newGeom = dict(zip(geomKeys, vars))
            self.dd['design'].update(newGeom)

            if 'Hm' in input_loads or 'Vm' in input_loads:
                anchor_loads = self.getLugForces(mudloads=input_loads)
                input_loads = dict(Ha=anchor_loads['Ha'], Va=anchor_loads['Va'])    # overwrite the input_loads dictionary
            results = self.getAnchorCapacity(loads=input_loads, plot=False)

            bound_L_lower = newGeom['L'] - geomBounds[0][0]
            bound_L_upper = geomBounds[0][1] - newGeom['L']
            bound_D_lower = newGeom['D'] - geomBounds[1][0]
            bound_D_upper = geomBounds[1][1] - newGeom['D']

            return np.array([bound_L_lower, bound_L_upper, bound_D_lower, bound_D_upper])
        
        # - - - - - Setup & Optimization
        from scipy.optimize import minimize
        from copy import deepcopy

        anchType = self.dd['type']

        # loads['Ha'] = minfs['Ha']*loads['Ha']
        # loads['Va'] = minfs['Va']*loads['Va']
        startGeom = dict(zip(geomKeys,geom))
        print('start geometry: ',startGeom)
        # apply initial guess geometry
        self.dd['design'].update(startGeom)
                   
        if not 'zlug' in self.dd['design']:
            if 'suction' in anchType and not fix_zlug:
                self.dd['design']['zlug'] = (2/3)*startGeom['L']
            else:
                self.dd['design']['zlug'] = 0
        
        # if zlug is fixed, remove it from design variables
        if fix_zlug and 'zlug' in geomKeys:
            zlug_loc = geomKeys.index('zlug')
            startGeom.pop('zlug')
            geomKeys.remove('zlug')
            geom.pop(zlug_loc)
            if geomBounds:
                geomBounds.pop(zlug_loc)

        if not loads:
            loads = self.loads
            
        if not 'Ha' in loads:
            loads = self.getLugForces(mudloads=loads)
            
        # suction bucket needs to be loads*FS because of capacity envelope calculations in capacity function
        if ('Hm' in loads and 'Vm' in loads) and ('Hm' in minfs and 'Vm' in minfs):
            input_loads = {'Hm':loads['Hm']*minfs['Hm'], 'Vm':loads['Vm']*minfs['Vm']}
        else:
            input_loads = {'Ha':loads['Ha']*minfs['Ha'],'Va':loads['Va']*minfs['Va']}

           

        
        # Initial guess for geometry
        initial_guess = geom  # [val for val in startGeom.values()]       # Input values for geometry
        # geomKeys = [key for key in startGeom.keys()]
        
        # Bounds and constraints
        if 'suction' in anchType:
            # bounds = [(1, 7), (5, 50),()]   # Bounds for D and L
            # constraints
            
            constraints = [{'type':'ineq','fun':conFun_LD,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)},
                           {'type':'ineq','fun':conFun_Suction,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)},
                           {'type':'ineq','fun':conFunH,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)},
                           {'type':'ineq','fun':conFunV,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)},
                           {'type':'ineq','fun':conBounds,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)}]
        
        elif 'dandg' in anchType:
            constraints = [{'type':'ineq','fun':conFun_LD,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)},
                           {'type':'ineq','fun':conFun_DandG,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)},
                           {'type':'ineq','fun':conFunH,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)},
                           {'type':'ineq','fun':conFunV,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)},
                           {'type':'ineq','fun':conBounds,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)}]
        
        else:
            constraints = [{'type':'ineq','fun':conFunH,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)},
                           {'type':'ineq','fun':conFunV,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)}]
        
        # Run the optimization to find sizing that satisfy UC close to 1
        print('optimizing anchor size')
        
        if 'suction' in anchType or 'dandg' in anchType:
            solution = minimize(objective, initial_guess, args=dict(geomKeys=geomKeys, input_loads=input_loads, fix_zlug=fix_zlug, LD_con=LD_con, geomBounds=geomBounds, minfs=minfs),
                                method="COBYLA", constraints=constraints, options={'rhobeg':0.1, 'catol':0.001})
        else:
            solution = minimize(objective, initial_guess, args=dict(geomKeys=geomKeys, input_loads=input_loads, fix_zlug=fix_zlug, LD_con=LD_con, geomBounds=geomBounds, minfs=minfs),
                                method="COBYLA", constraints=constraints, options={'rhobeg':0.1, 'catol':0.001})
        
        FS, acceptance, FSdiff = self.getFS(loads=input_loads, acceptance_crit=minfs)
        
        # adjust starting value if you're far off from the acceptance criteria (in either direction)
        if FSdiff_max:
            count = 0
            while count<10 and (np.any([abs(FSdiff[key])>FSdiff_max[key] for key in FSdiff.keys()]) or np.any([diff<0 for diff in FSdiff.values()])):
                if np.any([diff<.02 for key,diff in FSdiff.items() if minfs[key]>0]) and np.all([diff>=0 for diff in FSdiff.values()]):
                    # exit loop if you're as close as can be on one of the FS even if other is above diff requirements UNLESS an FS is below minimum reqiured FS
                    break
                print('Factor of Safety not close enough to minimum factor of safety, trying again with adjusted initial guess.')
                print(FS)
                # calculate new percent difference of FS from min fs
                diffPCT = [FSdiff[key]/FS[key] for key in FSdiff]
                # create adjustment coefficient based on this or .25, whichever is lower
                adjust_coeff = np.min([np.min(diffPCT),0.25])
                # adjust initial guess values by adjustment coefficient
                for i,val in enumerate(initial_guess):
                    initial_guess[i] = val - val*adjust_coeff
                # update zlug for suction buckets as needed to be 2/3L
                if 'suction' in anchType and not fix_zlug:
                    zlug_loc = geomKeys.index('zlug')
                    L_loc = geomKeys.index('L')
                    initial_guess[zlug_loc] = (2/3)*initial_guess[L_loc]

                print('new initial guess',initial_guess)
                # re-run optimization
                if 'suction' in anchType or 'dandg' in anchType:
                    solution = minimize(objective, initial_guess, args=dict(geomKeys=geomKeys, input_loads=input_loads, fix_zlug=fix_zlug, LD_con=LD_con, geomBounds=geomBounds, minfs=minfs),
                                        method="COBYLA", constraints=constraints, options={'rhobeg':0.1, 'catol':0.001})
                else:
                    solution = minimize(objective, initial_guess, args=dict(geomKeys=geomKeys, input_loads=input_loads, fix_zlug=fix_zlug, LD_con=LD_con, geomBounds=geomBounds, minfs=minfs),
                                        method="COBYLA", constraints=constraints, options={'rhobeg':0.1, 'catol':0.001})
                # re-determine FS and diff from minFS
                FS, acceptance, FSdiff = self.getFS(loads=input_loads, acceptance_crit=minfs)  
                count += 1
        
        # Extract the optimized values of geometry
        endGeom = dict(zip(geomKeys,solution.x))
        print('End geometry: ',endGeom)
        self.dd['design'].update(endGeom)
        if 'suction' in anchType and not fix_zlug:
            self.dd['design']['zlug'] = (2/3)*self.dd['design']['L'] 
        results = self.getAnchorCapacity(loads=input_loads,plot=plot)            
            
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

            SF = 1.0/UC if UC != 0 else np.inf

            return {'SF_combined': SF} 
                            
    def getCostAnchor(self, ms=None):
        '''
        Assign material cost using a Point object and getCost_and_MBL().
        '''

        # Create or use existing MoorPy system
        if ms is None:
            ms = mp.System()

            # Create MoorPy Point using makeMoorPyAnchor
            self.makeMoorPyAnchor(ms)

        # Check if mass is assigned
        if self.mass is None:
            if 'Weight pile' in self.anchorCapacity:
                self.mass = self.anchorCapacity['Weight pile'] / self.g
            elif 'Weight plate' in self.anchorCapacity:
                self.mass = self.anchorCapacity['Weight plate'] / self.g
            else:
                raise KeyError("Missing 'Weight pile' or 'Weight plate' in anchorCapacity. \
                Run getCapacityAnchor() before getCostAnchor(), or define self.mass explicitly.")
        
        # Assign mass to MoorPy point
        self.mpAnchor.m = self.mass

        cost, MBL, info = self.mpAnchor.getCost_and_MBL()

        # Store results
        self.cost = {
            'Material cost': cost,
            'MBL': MBL,
            'unit_cost': cost/self.mpAnchor.m }

        return self.cost

    def getCombinedPlot(self):
        '''
        Create a plot showing the suction pile and the inverse catenary overlay in the same coordinate system.
        '''
        from anchors_famodel.capacity_load  import getTransferLoad
        from anchors_famodel.capacity_plots import plot_suction

        if self.anchType.lower() != 'suction':
            raise NotImplementedError("getCombinedPlot only supports suction piles.")

        # Extract design inputs
        design = self.dd['design']
        D = design['D']
        L = design['L']
        zlug = design['zlug']

        if self.soil_profile is None or self.soil_type is None:
            raise ValueError("Soil profile or type not assigned. Use setSoilProfile first.")

        soil_profile = self.soil_profile
        soil_type = self.soil_type
        z0 = soil_profile[0]['top']

        Hm = self.loads['Hm']
        Vm = self.loads['Vm']
        thetam = self.loads.get('thetam', np.degrees(np.arctan2(Vm, Hm)))

        line_type = getattr(self, 'line_type', 'chain')
        d = getattr(self, 'd', 0.16)
        w = getattr(self, 'w', 5000.0)

        # Get inverse catenary path
        layers, result = getTransferLoad(
            profile_map=[{'layers': self.soil_profile}],
            Tm=np.sqrt(Hm**2 + Vm**2),
            thetam=thetam,
            zlug=zlug,
            line_type=line_type,
            d=d,
            w=w,
            plot=False
        )

        drag_values  = np.array(result['drag_values'])
        depth_values = -np.array(result['depth_values'])[::-1]

        x_start = D/2 + drag_values[0]
        z_start = zlug
        drag_transformed  = x_start - drag_values
        depth_transformed = z_start + (depth_values- depth_values[0])

        # Plot suction pile
        plot_suction(soil_profile, L, D, z0=z0, zlug=zlug, title='Suction Pile and Mooring Line Load Path')


        # Overlay inverse catenary path
        plt.plot(drag_transformed, depth_transformed, color='b', lw=2.0, label='Inverse catenary')
        plt.plot(drag_transformed[-1], depth_transformed[-1], 'ro', label='Mudline end')
        plt.plot( drag_transformed[0],  depth_transformed[0], 'go', label='Embedded end')

        n = 2e6
        Tm = result['Tm']
        Ta = result['Ta']
        thetaa = result['thetaa']

        plt.arrow(drag_transformed[-1], depth_transformed[-1],
                  Tm*np.cos(np.deg2rad(thetam))/n, -Tm*np.sin(np.deg2rad(thetam))/n,
                  head_width=0.25, head_length=0.5, color='r', label='Mudline load')
    
        plt.arrow(drag_transformed[0], depth_transformed[0],
                  Ta*np.cos(np.deg2rad(thetaa))/n, -Ta*np.sin(np.deg2rad(thetaa))/n,
                  head_width=0.25, head_length=0.5, color='g', label='Padeye load')

        xmax = max(drag_transformed[-1] + D, 2*D)
        plt.xlim(-D, xmax)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()
