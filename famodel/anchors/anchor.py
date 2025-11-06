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
                 g=9.81, rho=1025, profile_map=None, display=0):
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

        # Initialize as a node
        Node.__init__(self, id)

        # Design description dictionary for this Anchor
        self.dd = dd
        # MoorPy system this anchor is in
        self.ms = ms
        # x,y,z location of anchor
        self.r = r
        # anchor index in array mooring list 
        self.aNum = aNum
        
        self.g = g
        self.rho = rho

        if dd and 'type' in dd:
            self.anchType = dd['type']
        else:
            self.anchType = 'suction'
            print(f"[Anchor] No type provided. Defaulting to 'suction'.")
        
        # raise errors/warnings if the anchor type is not what it needs to be
        anchor_type_options = ['suction', 'sepla', 'dea', 'depla', 'vla', 'plate', 'torpedo', 'helical', 'driven', 'drilled']
        if self.anchType.lower() not in anchor_type_options:
            raise ValueError(f"The anchor 'type' {self.anchType} needs to explicitly be one of {anchor_type_options} (Case not sensitive)")
        # if self.anchType not in ['drag-embedment', 'gravity', 'suction', 'SEPLA', 'VLA', 'driven']:
        #     print('Warning: The anchor type provided does not have any cost coefficients. This will default to a suction pile')

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
        
        self.display = display

    def setSoilProfile(self, profile_map):
        '''
        Assign a soil profile directly from a single CPT.
        Assumes profile_map is a list with only one entry.
        '''
        if len(profile_map) != 1:
            raise ValueError("setSoilProfile expects a profile_map with exactly one CPT.")

        cpt = profile_map[0]
        self.soil_profile = profile_map               
        self.profile_name = cpt.get('name', 'CPT_Assigned')
    
        # Extract soil types from layers
        layers = cpt['layers']
        soil_types = [layer['soil_type'] for layer in layers]
        self.soil_type_list = list(set(soil_types))
        self.soil_type = soil_types[0] if len(self.soil_type_list) == 1 else 'mixed'
    
        # Group layers by soil type
        soilProps = defaultdict(list)
        for layer in layers:
            layer_copy = layer.copy()
            soil_type = layer_copy.pop('soil_type')
            soilProps[soil_type].append(layer_copy)
        self.soilProps = dict(soilProps)

        if self.display > 0: print(f"[Anchor] Assigned soil profile from {self.profile_name} with soil types {self.soil_type_list}")

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

        self.soil_profile = [{
            'name': 'Interpolated_2D',
            'x': x_anchor,
            'y': y_anchor,
            'layers': interpolated_layers}]
        self.profile_name = "Interpolated_2D"

        # Extract soil types
        layers = self.soil_profile[0]['layers']
        soil_types = [layer['soil_type'] for layer in layers]
        self.soil_type_list = list(set(soil_types))
        self.soil_type = soil_types[0] if len(self.soil_type_list) == 1 else 'mixed'

        # Group interpolated layers by soil type
        soilProps = defaultdict(list)
        for layer in self.soil_profile[0]['layers']:
            layer_copy = layer.copy()
            soil_type = layer_copy.pop('soil_type')
            soilProps[soil_type].append(layer_copy)
        self.soilProps = dict(soilProps)

        if self.display > 0: print(f"[Anchor] Interpolated soil profile: {self.profile_name} with soil types {self.soil_type_list}")
        
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

        '''
        # Create anchor as a fixed point in MoorPy system
        ms.addPoint(1, self.r)

        # Assign this point as mpAnchor in the anchor class instance
        self.mpAnchor = ms.pointList[-1]
        '''

        # create anchor as a fixed body in MoorPy system and assign to mpAnchor property
        # r6 = [self.r[0],self.r[1],self.r[2],0,0,0]
        # self.mpAnchor = ms.addBody(1,r6)
        self.mpAnchor = ms.addPoint(1,self.r)

        # Set mass if available
        if 'mass' in self.dd.get('design', {}):
            self.mpAnchor.m = self.dd['design']['mass']

        # Set diameter if available
        if 'd' in self.dd.get('design', {}):
            self.mpAnchor.d = self.dd['design']['d']
            
        # Set dummy design to get PointType from MoorPy
        if anchType not in list(ms.pointProps['AnchorProps'].keys()):
            anchType = 'suction'    # default to a suction pile just to get the MoorPy pointProps working
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
                mtype = att['obj'].sections()[0]['type']['material'].lower()
                if 'chain' not in mtype:
                    print('No chain below seafloor, setting Ta=Tm (no load transfer).')
                    return mtype, None, None, True
                else:
                    d_nom = att['obj'].sections()[0]['type']['d_nom']
                    w_nom = att['obj'].sections()[0]['type']['w']
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
            # loads = self.mpAnchor.getForces(lines_only=lines_only, seabed=seabed, xyz=xyz)
            # MoorPy Body.getForces does not accept seabed/lines_only in current API.
            # Get forces (total), optionally post-process seabed if needed.
            loads = self.mpAnchor.getForces()
            Hm = np.sqrt(loads[0]**2 + loads[1]**2)
            Vm = loads[2]
            thetam = np.degrees(np.arctan2(Vm, Hm))
            self.loads['Hm'] = Hm
            self.loads['Vm'] = Vm
            self.loads['thetam'] = thetam
            self.loads['mudline_load_type'] = 'current_state'

        self.loads['method'] = 'static'
        return self.loads
       
    def getLugForces(self, Hm, Vm, zlug, line_type=None, d=None, w=None, plot=False, display=0):
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
        z0 = soil_profile[0]['layers'][0]['top']
        
        # Load transfer if padeye is embedded
        if zlug > z0:
            # Check if padeye is embedded in rock
            if any(layer.get('soil_type') == 'rock' for layer in self.soil_profile[0]['layers']):
                raise ValueError('[Warning] Padeye depth is embedded in rock. Embedded line in rock is not possible.')
             

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
                profile_map=self.soil_profile,
                Tm=np.sqrt(Hm**2 + Vm**2),
                thetam=np.degrees(np.arctan2(Vm, Hm)),
                zlug=zlug,
                line_type=line_type,
                d=d,
                w=w,
                plot=plot, display=display)

            Ta = loads['Ta']
            thetaa = loads['thetaa']
            Ha = Ta*np.cos(np.deg2rad(thetaa))
            Va = Ta*np.sin(np.deg2rad(thetaa))

        else:
            Ha = Hm
            Va = Vm
            layers = self.soil_profile[0]['layers']  

            
        if plot == True:
            plot_load(layers, loads['drag_values'], loads['depth_values'], 
                      loads['Tm'], loads['thetam'], loads['Ta'], 
                      loads['thetaa'], zlug=zlug)

        return layers, Ha, Va

    def getCapacityAnchor(self, Hm, Vm, zlug, line_type=None, d=None, w=None, mass_update=False, plot=False, display=0):
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

        from .anchors_famodel.capacity_plate import getCapacityPlate
        from .anchors_famodel.capacity_suction import getCapacitySuction
        from .anchors_famodel.capacity_torpedo import getCapacityTorpedo
        from .anchors_famodel.capacity_helical import getCapacityHelical
        from .anchors_famodel.capacity_driven import getCapacityDriven
        from .anchors_famodel.capacity_drilled import getCapacityDrilled
           
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
            'drilled': getCapacityDrilled}
        
        if self.display > 0: print('[DEBUG] profile_name:', self.profile_name)
        if self.display > 0: print('[DEBUG] soil_profile passed as profile_map:')
        for entry in self.soil_profile:
            if self.display > 0: print(entry.get('name'), list(entry.keys()))


        if self.display > 1: print(f'[Debug] mass_update = {mass_update}')
        anchType_clean = self.dd['type'].lower().replace(' ', '')
        capacity_func = capacity_dispatch.get(anchType_clean)
        if capacity_func is None:
            raise ValueError(f"Unknown anchor type '{self.anchType}' for anchor capacity calculation.")
    
        if self.soil_profile is None or self.soil_type is None:
            raise ValueError("Soil profile or soil type not set for this anchor.")
    
        soil_profile = self.soil_profile
        soil_type = self.soil_type
        z0 = soil_profile[0]['layers'][0]['top']  
    
        if line_type is None or d is None or w is None:
            try:
                line_type, d, w = self.getLineProperties()
            except ValueError:
                if self.display > 0: print('[Warning] No mooring attachment found. Trying anchor-level line properties...')
                line_type = getattr(self, 'line_type', None)
                d = getattr(self, 'd', None)
                w = getattr(self, 'w', None)

                if any(v is None for v in [line_type, d, w]):
                    if self.display > 0: print('[Fallback] Using default chain properties.')
                    line_type = 'chain'
                    d = 0.16
                    w = 5500.0
    
        # Load transfer if padeye is embedded below mudline
        if zlug > z0:
            layers, Ha, Va = self.getLugForces(
                Hm, Vm,
                zlug=zlug,
                line_type=line_type,
                d=d,
                w=w,
                plot=False, display=display)

            Ta = np.sqrt(Ha**2 + Va**2)
            thetaa = np.degrees(np.arctan2(Va, Ha))
            
            if self.display > 0: print(f"[Branch Check] Entered {'zlug>z0' if zlug>z0 else 'else'} for anchor {self.anchType}")

        else:
            Ha = Hm
            Va = Vm
            Ta = np.sqrt(Ha**2 + Va**2)
            thetaa = np.degrees(np.arctan2(Va, Ha))

            if self.display > 0: print(f"[Branch Check] Entered {'zlug>z0' if zlug>z0 else 'else'} for anchor {self.anchType}")

        # --- Call the appropriate capacity function ---
        if anchType_clean in ['sepla', 'dea', 'depla', 'vla', 'plate']:
            self.capacity_format = 'plate'
            B = self.dd['design']['B']
            L = self.dd['design']['L']
            if self.display > 1: print(f"[Final Check] Ha = {Ha}, Va = {Va}, anchor = {self.anchType}")
            beta = 90.0 - np.degrees(np.arctan2(Va, Ha))
            self.dd['design']['beta'] = beta 
            layers, results = capacity_func(
                profile_map=self.soil_profile,
                location_name=self.profile_name,
                B=B, L=L, zlug=zlug,
                beta=beta,
                Ha=Ha, Va=Va,
                plot=plot, display=display)
            
        elif anchType_clean == 'suction':
            self.capacity_format = 'envelope'
            D = self.dd['design']['D']
            L = self.dd['design']['L']
            zlug = self.dd['design']['zlug']
            layers, results = capacity_func(
                profile_map=self.soil_profile,
                location_name=self.profile_name,
                D=D, L=L, zlug=zlug,
                Ha=Ha, Va=Va,
                thetalug=5, psilug=7.5,
                plot=plot, display=display)
            
        elif anchType_clean == 'torpedo':
            self.capacity_format = 'envelope'
            D1 = self.dd['design']['D1']
            D2 = self.dd['design']['D2']
            L1 = self.dd['design']['L1']
            L2 = self.dd['design']['L2']
            ballast = self.dd['design'].get('ballast', 0.0)
            layers, results = capacity_func(
                profile_map=self.soil_profile,
                location_name=self.profile_name,
                D1=D1, D2=D2, L1=L1, L2=L2,
                zlug=zlug,
                ballast=ballast,
                Ha=Ha, Va=Va,
                plot=plot, display=display)

        elif anchType_clean == 'helical':
            self.capacity_format = 'component'
            D = self.dd['design']['D']     
            L = self.dd['design']['L']     
            d = self.dd['design']['d']     
            zlug = self.dd['design']['zlug']
            layers, results = capacity_func(
                profile_map=self.soil_profile,
                location_name=self.profile_name,
                D=D, L=L, d=d,
                zlug=zlug,
                Ha=Ha, Va=Va,
                plot=plot, display=display)

        elif anchType_clean == 'driven':
            self.capacity_format = 'component'
            L = self.dd['design']['L']
            D = self.dd['design']['D']
            zlug = self.dd['design']['zlug']
            layers, y, z, results = capacity_func(
                profile_map=self.soil_profile,
                location_name=self.profile_name,
                L=L, D=D, zlug=zlug,
                Ha=Ha, Va=Va,
                plot=plot, display=display)
        
        elif anchType_clean == 'drilled':
            self.capacity_format = 'component'           
            L = self.dd['design']['L']
            D = self.dd['design']['D']
            zlug = self.dd['design']['zlug']
            layers, y, z, results = capacity_func(
                profile_map=self.soil_profile,
                location_name=self.profile_name,
                L=L, D=D, zlug=zlug,
                Ha=Ha, Va=Va,
                plot=plot, display=display)

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
        
        elif anchType_clean in ['helical', 'driven', 'drilled']:
            self.anchorCapacity['Unity check (horizontal)'] = results.get('Unity check (horizontal)', np.nan)
            self.anchorCapacity['Unity check (vertical)'] = results.get('Unity check (vertical)', np.nan)
        
        # Copy over lateral and rotational displacements
        if 'Lateral displacement' in results:
            self.anchorCapacity['Lateral displacement'] = results['Lateral displacement']
        if 'Rotational displacement' in results:
            self.anchorCapacity['Rotational displacement'] = results['Rotational displacement']
        
        # Weight calculated via dimensions
        if not mass_update:
            if 'Weight pile' in results:
                self.anchorCapacity['Weight pile'] = results['Weight pile']
            if 'Weight plate' in results:
                self.anchorCapacity['Weight plate'] = results['Weight plate']
        else:
            if 'Weight pile' in results:
                self.mass = results['Weight pile']/self.g
                self.anchorCapacity['Weight pile'] = results['Weight pile']
            if 'Weight plate' in results:
                self.mass = results['Weight plate']/self.g
                self.anchorCapacity['Weight plate'] = results['Weight plate']
                
        # print(f"[DEBUG] Stored Lateral displacement in anchorCapacity: {self.anchorCapacity['Lateral displacement']:.6f}")
         
    def getSizeAnchor(self, geom, geomKeys, geomBounds=None, loads=None, lambdap_con=[4, 8],
                       zlug_fix=True, safety_factor={}, plot=False, display=0):
        '''
        Generalized optimization method for all anchor types, using dictionary-based safety factors.
        '''
        self.display = display
    
        anchType_clean = self.dd['type'].strip().lower()
        print(f"[Debug] Anchor type parsed: '{anchType_clean}'")
    
        if loads is None:
            loads = self.loads
        
        sf_Hm = safety_factor.get('Hm', safety_factor.get('SF_horizontal', 1.0))
        sf_Vm = safety_factor.get('Vm', safety_factor.get('SF_vertical', 1.0))
        sf_uc = safety_factor.get('SF_combined', max(sf_Hm, sf_Vm))  # conservative by default

        Hm = loads['Hm']*sf_Hm
        Vm = loads['Vm']*sf_Vm
    
        line_type = getattr(self, 'line_type', 'chain')
        d = getattr(self, 'd', 0.16)
        w = getattr(self, 'w', 5000.0)
    
        def update_zlug():
            if 'suction' in anchType_clean and not zlug_fix and 'zlug' not in geomKeys:
                self.dd['design']['zlug'] = (2/3)*self.dd['design']['L']
            elif np.any([name in anchType_clean for name in ['driven', 'helical']]) and not zlug_fix:    
                ratio = self.dd['design'].get('zlug_ratio', self.dd['design']['zlug']/self.dd['design']['L'])
                self.dd['design']['zlug_ratio'] = ratio
                self.dd['design']['zlug'] = ratio*self.dd['design']['L']
            elif 'drilled' in anchType_clean:
                self.dd['design']['zlug'] = 0
    
        def get_lambda():
            if 'torpedo' in anchType_clean:
                L = self.dd['design']['L1'] + self.dd['design']['L2']
                A_wing = (self.dd['design']['D1'] - self.dd['design']['D2']) * self.dd['design']['L1']
                A_shaft = self.dd['design']['D2'] * L
                D = (A_wing + A_shaft) / L
            elif np.any([name in anchType_clean for name in ['driven', 'drilled', 'helical', 'suction']]):
                L = self.dd['design']['L']
                D = self.dd['design']['D']
            elif np.any([name in anchType_clean for name in ['plate', 'sepla', 'dea', 'depla', 'vla']]):
                L = self.dd['design']['L']
                D = self.dd['design']['B']
            else:
                raise ValueError(f'lambda not defined for anchor type: {anchType_clean}')
            return L/D
    
        def constraint_lambda_min(vars):
            for i, key in enumerate(geomKeys):
                self.dd['design'][key] = vars[i]
            update_zlug()
            self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                                    line_type=line_type, d=d, w=w, mass_update=True, plot=False, display=display)
            return get_lambda() - lambdap_con[0]
    
        def constraint_lambda_max(vars):
            return lambdap_con[1] - get_lambda()
        
        def constraint_bounds(vars):
            con_bound_return = np.zeros(len(geomKeys)*2)
            for i,var in enumerate(geomKeys):
                con_bound_return[2*i] = self.dd['design'][var] - geomBounds[i][0]
                con_bound_return[2*i+1] = geomBounds[i][1] - self.dd['design'][var]
            return con_bound_return
    
        if np.any([name in anchType_clean for name in ['suction', 'torpedo', 'plate', 'sepla', 'dea', 'depla', 'vla']]):
            target_UC = 1.0/sf_uc
    
            def objective_uc(vars):
                '''
                for i, key in enumerate(geomKeys):
                    self.dd['design'][key] = vars[i]
                update_zlug()
                self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                                       line_type=line_type, d=d, w=w, mass_update=True, plot=False, display=display)
                '''
                #UC = self.anchorCapacity.get('UC', 2.0)
                #return (UC - target_UC)**2
                #return self.anchorCapacity.get('Weight pile')
                if any(name in anchType_clean for name in ['plate', 'sepla', 'dea', 'depla', 'vla']):
                    return self.anchorCapacity.get('Weight plate')
                else:
                    return self.anchorCapacity.get('Weight pile')
    
            def constraint_uc_envelope(vars):
                return self.anchorCapacity.get('UC', 0.0) - target_UC
    
            constraints_uc = [
                {'type': 'ineq', 'fun': constraint_lambda_min},
                {'type': 'ineq', 'fun': constraint_lambda_max},
                {'type': 'ineq', 'fun': constraint_uc_envelope},
                {'type': 'ineq', 'fun': constraint_bounds},
            ]
    
            result_uc = minimize(
                objective_uc,
                geom,
                method='COBYLA',
                constraints=constraints_uc,
                options={'rhobeg': 0.1, 'catol': 0.01, 'maxiter': 500}
            )
    
            endGeom = dict(zip(geomKeys, result_uc.x))
            self.dd['design'].update(endGeom)
            update_zlug()
            self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                                   line_type=line_type, d=d, w=w, mass_update=True, plot=plot, display=display)
    
            print('\nFinal Optimized Anchor (UC-based):')
            print('Design:', self.dd['design'])
            print('Capacity Results:', self.anchorCapacity)
            return
        
        def near_border():
            UC_h = self.anchorCapacity['Ha']/self.anchorCapacity['Hmax']
            UC_v = self.anchorCapacity['Va']/self.anchorCapacity['Vmax']
            disp_lat = abs(self.anchorCapacity.get('Lateral displacement', 0.0))
            disp_rot = abs(self.anchorCapacity.get('Rotational displacement', 0.0))
            limit_lat = 0.05*self.dd['design']['D']  # 10% of the pile diameter
            limit_rot = 10.0                         # 10 deg
        
            near_UC_h = 0.95 <= UC_h <= 1.0
            near_UC_v = 0.95 <= UC_v <= 1.0
            near_disp_lat = 0.95*limit_lat <= disp_lat <= limit_lat
            near_disp_rot = 4.75 <= disp_rot <= limit_rot
        
            return near_UC_h or near_UC_v or near_disp_lat or near_disp_rot

        def termination_condition():
            UC_h = self.anchorCapacity['Ha']/self.anchorCapacity['Hmax']
            UC_v = self.anchorCapacity['Va']/self.anchorCapacity['Vmax']
            disp_lat = abs(self.anchorCapacity.get('Lateral displacement', 0.0))
            disp_rot = abs(self.anchorCapacity.get('Rotational displacement', 0.0))
            limit_lat = 0.05*self.dd['design']['D']  # 10% of the pile diameter
            limit_rot = 10.0                         # 10 deg
        
            all_satisfied = (UC_h <= 1.0 and UC_v <= 1.0 and disp_lat <= limit_lat and disp_rot <= limit_rot)
        
            if all_satisfied:
                if near_border():
                    if self.display > 0: print('[Termination] All criteria satisfied and near border.')
                    return 'terminate'
                else:
                    if self.display > 0: print('[Safe but not near border] Continue shrinking...')
                    return 'continue'
            return 'continue'
        
        def termination_condition_drilled():
            UC_v = self.anchorCapacity['Va']/self.anchorCapacity['Vmax']
            disp_lat = abs(self.anchorCapacity.get('Lateral displacement', 0.0))
            disp_rot = abs(self.anchorCapacity.get('Rotational displacement', 0.0))
            limit_lat = 0.05*self.dd['design']['D']  # 10% of the pile diameter
            limit_rot = 10.0                         # 10 deg
        
            all_satisfied = (UC_v <= 1.0 and disp_lat <= limit_lat and disp_rot <= limit_rot)
        
            if all_satisfied:
                if near_border():
                    if self.display > 0: print('[Termination] All criteria satisfied and near border.')
                    return 'terminate'
                else:
                    if self.display > 0: print('[Safe but not near border] Continue shrinking...')
                    return 'continue'
            return 'continue'
           
        def is_valid(value):
            return np.isfinite(value) and not np.isnan(value) and abs(value) < 1e6
        
        if anchType_clean in ['helical', 'driven']:
            L0, D0 = geom if len(geom) == 2 else [5.0, 1.0]
            self.dd['design']['L'] = L0
            self.dd['design']['D'] = D0
            Lmin, Lmax = geomBounds[0]
            Dmin, Dmax = geomBounds[1]
            update_zlug()
            self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                                   line_type=line_type, d=d, w=w, mass_update=True, plot=False, display=display)
        
            UC_h = self.anchorCapacity['Ha']/self.anchorCapacity['Hmax']
            UC_v = self.anchorCapacity['Va']/self.anchorCapacity['Vmax']
            disp_lat = abs(self.anchorCapacity.get('Lateral displacement', 0.0))
            disp_rot = abs(self.anchorCapacity.get('Rotational displacement', 0.0))
            limit_disp = 0.10*D0    # 10% of the pile diameter
            limit_rot = 10.0        # 10 deg 
            direction = 'shrink' if (UC_h <= 1.0 and UC_v <= 1.0 and disp_lat <= limit_disp and disp_rot <= limit_rot) else 'grow'
    
            max_iter = 200
            iter_count = 0
    
            if direction == 'shrink':
                for L in np.arange(L0, Lmin - 1e-6, -0.25):
                    self.dd['design']['L'] = L
                    for D in np.arange(Dmax, Dmin - 1e-6, -0.05):
                        if L/D > lambdap_con[1] or L/D < lambdap_con[0]:
                            continue
                        self.dd['design']['L'] = L
                        update_zlug()
                        self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                                               line_type=line_type, d=d, w=w, mass_update=True, plot=False, display=display)
                        UC_h = self.anchorCapacity['Ha']/self.anchorCapacity['Hmax']
                        UC_v = self.anchorCapacity['Va']/self.anchorCapacity['Vmax']
                        disp_lat = abs(self.anchorCapacity.get('Lateral displacement', 0.0))
                        disp_rot = abs(self.anchorCapacity.get('Rotational displacement', 0.0))
                        if self.display > 0: print(f'[Iter {iter_count}] L={L:.2f}, D={D:.2f}, UC_h={UC_h:.3f}, UC_v={UC_v:.3f}, lat={disp_lat:.3f} m, rot={disp_rot:.3f} deg')
                        iter_count += 1
                        if not all(is_valid(v) for v in [UC_h, UC_v, disp_lat, disp_rot]):
                            continue
                        if termination_condition():
                            print(f'\nTermination criteria met.')
                            print('Design:', self.dd['design'])
                            print('Capacity Results:', self.anchorCapacity)
                            return
            elif direction == 'grow':
                for L in np.arange(L0, Lmax + 1e-6, 0.25):
                    self.dd['design']['L'] = L
                    for D in np.arange(Dmin, Dmax + 1e-6, 0.05):   
                        if L/D > lambdap_con[1] or L/D < lambdap_con[0]:
                            continue
                        self.dd['design']['L'] = L
                        update_zlug()
                        self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                                               line_type=line_type, d=d, w=w, mass_update=True, plot=False, display=display)
                        UC_h = self.anchorCapacity['Ha']/self.anchorCapacity['Hmax']
                        UC_v = self.anchorCapacity['Va']/self.anchorCapacity['Vmax']
                        disp_lat = abs(self.anchorCapacity.get('Lateral displacement', 0.0))
                        disp_rot = abs(self.anchorCapacity.get('Rotational displacement', 0.0))
                        if self.display > 0: print(f'[Iter {iter_count}] L={L:.2f}, D={D:.2f}, UC_h={UC_h:.3f}, UC_v={UC_v:.3f}, lat={disp_lat:.3f} m, rot={disp_rot:.3f} deg')
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
            else:
                raise ValueError(f"Unknown optimization direction: {direction}")
                
            if self.display > 0: print('[Warning] While-loop search reached bounds without meeting criteria.')
                    
        if 'drilled' in anchType_clean:
            L0, D0 = geom if len(geom) == 2 else [5.0, 1.0]
            self.dd['design']['L'] = L0
            self.dd['design']['D'] = D0
            Lmin, Lmax = geomBounds[0]
            Dmin, Dmax = geomBounds[1]
            update_zlug()
            self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                                   line_type=line_type, d=d, w=w, mass_update=True, plot=False, display=display)
        
            UC_v = self.anchorCapacity['Va']/self.anchorCapacity['Vmax']
            disp_lat = abs(self.anchorCapacity.get('Lateral displacement', 0.0))
            disp_rot = abs(self.anchorCapacity.get('Rotational displacement', 0.0))
            limit_disp = 0.10*D0   # 10% of the pile diameter
            limit_rot = 10.0        # 10 deg 
            direction = 'shrink' if (UC_v <= 1.0 and disp_lat <= limit_disp and disp_rot <= limit_rot) else 'grow'
    
            max_iter = 200
            iter_count = 0
    
            if direction == 'shrink':
                for L in np.arange(L0, Lmin - 1e-6, -0.25):
                    self.dd['design']['L'] = L
                    for D in np.arange(Dmax, Dmin - 1e-6, -0.05):
                        if L/D > lambdap_con[1] or L/D < lambdap_con[0]:
                            continue
                        self.dd['design']['L'] = L
                        update_zlug()
                        self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                                               line_type=line_type, d=d, w=w, mass_update=True, plot=False, display=display)
                        UC_v = self.anchorCapacity['Va']/self.anchorCapacity['Vmax']
                        disp_lat = abs(self.anchorCapacity.get('Lateral displacement', 0.0))
                        disp_rot = abs(self.anchorCapacity.get('Rotational displacement', 0.0))
                        if self.display > 0: print(f'[Iter {iter_count}] L={L:.2f}, D={D:.2f}, UC_v={UC_v:.3f}, lat={disp_lat:.3f} m, rot={disp_rot:.3f} deg')
                        iter_count += 1
                        if not all(is_valid(v) for v in [UC_v, disp_lat, disp_rot]):
                            continue
                        if termination_condition_drilled():
                            print(f'\nTermination criteria met.')
                            print('Design:', self.dd['design'])
                            print('Capacity Results:', self.anchorCapacity)
                            return
            elif direction == 'grow':
                for L in np.arange(L0, Lmax + 1e-6, 0.25):
                    self.dd['design']['L'] = L
                    for D in np.arange(Dmin, Dmax + 1e-6, 0.05):
                        if L/D > lambdap_con[1] or L/D < lambdap_con[0]:
                            continue
                        self.dd['design']['L'] = L
                        update_zlug()
                        self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'],
                                               line_type=line_type, d=d, w=w, mass_update=True, plot=False, display=display)
                        UC_v = self.anchorCapacity['Va']/self.anchorCapacity['Vmax']
                        disp_lat = abs(self.anchorCapacity.get('Lateral displacement', 0.0))
                        disp_rot = abs(self.anchorCapacity.get('Rotational displacement', 0.0))
                        if self.display > 0: print(f'[Iter {iter_count}] L={L:.2f}, D={D:.2f}, UC_v={UC_v:.3f}, lat={disp_lat:.3f} m, rot={disp_rot:.3f} deg')
                        iter_count += 1
                        status = termination_condition_drilled()
                        if status == 'terminate':
                            print(f'Termination criteria met.')
                            print('Design:', self.dd['design'])
                            print('Capacity Results:', self.anchorCapacity)
                            return
                        elif status == 'continue':
                            continue
                    status = termination_condition_drilled()
                    if status == 'terminate':
                        print(f'\nTermination criteria met.')
                        print('Design:', self.dd['design'])
                        print('Capacity Results:', self.anchorCapacity)
                        return
            else:
                raise ValueError(f"Unknown optimization direction: {direction}")
    
            if self.display > 0: print('[Warning] While-loop search reached bounds without meeting criteria.')
    
        else:
            raise ValueError(f"Anchor type '{anchType_clean}' not supported for safety factor input.")
                   
   
    def getSafetyFactor(self):
        '''
        Calculate the safety factor based on the unity checks stored in capacity results.

        Returns
        -------
        dict
            Dictionary containing safety factors.
        '''

        anchType_clean = self.anchType.lower().replace(' ', '')

        if anchType_clean in ['helical', 'driven', 'drilled']:
            UC_v = self.anchorCapacity.get('Unity check (vertical)', None)
            UC_h = self.anchorCapacity.get('Unity check (horizontal)', None)

            if UC_v is None or UC_h is None:
                print("Warning: Vertical or horizontal unity check (UC) not found in capacity results. Returning NaN.")
                return {'SF_vertical': np.nan, 'SF_horizontal': np.nan}

            SF_v = 1.0/UC_v if UC_v != 0 else np.inf
            SF_h = 1.0/UC_h if UC_h != 0 else np.inf

            return {'SF_vertical': SF_v, 'SF_horizontal': SF_h}

        else:
            UC = self.anchorCapacity.get('UC', None)

            if UC is None:
                print("Warning: Unity check (UC) not found in capacity results. Returning NaN.")
                return {'SF_combined': np.nan}

            SF = 1.0/UC if UC != 0 else np.inf

            return {'SF_combined': SF} 
                            
    def getCost(self, ms=None, mass_update=True):
        '''
        Assign material cost using a Point object and getCost_and_MBL().
        
        Parameters
        ----------
        ms : MoorPy System, optional
            The mooring system to which the anchor point belongs. If None, a new one is created.
        mass_update : bool, optional
            If True, update mpAnchor mass from self.mass.
            If False, preserve existing mpAnchor.m if already set.
        '''

        # Create or use existing MoorPy system
        if ms is None:
            ms = mp.System()

            # Create MoorPy Point using makeMoorPyAnchor
            self.makeMoorPyAnchor(ms)

        # Assign self.mass if missing
        if self.mass is None or mass_update:
            if 'Weight pile' in self.anchorCapacity:
                self.mass = self.anchorCapacity['Weight pile']/self.g
            elif 'Weight plate' in self.anchorCapacity:
                self.mass = self.anchorCapacity['Weight plate']/self.g
            else:
                raise KeyError("Missing 'Weight pile' or 'Weight plate' in anchorCapacity. \
                Run getCapacityAnchor() before getCostAnchor(), or define self.mass explicitly.")
        
        # Assign mass to MoorPy point
        self.mpAnchor.m = self.mass

        cost, MBL, info = self.mpAnchor.getCost_and_MBL()

        # Store results
        self.cost = {
            'Material cost': cost,
            #'MBL': MBL,
            #'unit_cost': cost/self.mpAnchor.m 
            }

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




           
        
    # def getSuctionSize(self,D,L,loads=None,minfs={'Ha':1.6,'Va':2},LD_con=[4,8]):
    #     '''
        

    #     Parameters
    #     ----------
    #     D : float
    #         Diameter of suction bucket
    #     L : float
    #         Length of suction bucket
    #     loads : dict, optional
    #         Dictionary of maximum anchor loads in horizontal and vertical directions. The default is None.
    #     minfs : dict,optoinal
    #         Minimum factors of safety in horizontal and vertical directions
    #     LD_con : float
    #         Constraint for L/D parameter

    #     Returns
    #     -------
    #     None.

    #     '''
    #     from scipy.optimize import minimize
    #     anchType = self.dd['type']
    #     if not loads:
    #         loads = self.loads
            
    #     if not 'Ha' in loads:
    #         loads = self.getLugForces(mudloads=loads)
           
    #     loads['Ha'] = minfs['Ha']*loads['Ha']
    #     loads['Va'] = minfs['Va']*loads['Va']
                   
    #     if not 'zlug' in self.dd['design']:
    #         self.dd['design']['zlug'] = (2/3)*L
            
    #     # Define the objective function: Minimize |UC - 1| (aim for UC to be 1)
    #     def objective(vars):
    #         D, L = vars
    #         self.dd['design']['D'] = D
    #         self.dd['design']['L'] = L
    #         self.dd['design']['zlug'] = (2/3)*L
    #         results = self.getAnchorCapacity(plot=False)
    #         return abs(results['UC'] - 1)  
        
    #     def conFun(vars,LD_con):
    #         D, L = vars
    #         if L/D >= LD_con[0] and L/D <= LD_con[1]:
    #             conval = 1 
    #         else:
    #             conval = -1 
            
    #         return(conval)
        
    #     # Initial guess for D and L
    #     initial_guess = [D, L]       # Input values for D and L
        
    #     # Bounds for D and L (adjust as needed)
    #     bounds = [(1, 5), (5, 50)]   # Bounds for D and L
        
    #     # constraints
    #     constraints = [{'type':'ineq','fun':conFun,'args':(LD_con,)}]
        
    #     # Run the optimization to find D and L that satisfy UC close to 1
    #     solution = minimize(objective, initial_guess, bounds=bounds,method="COBYLA",
    #                         constraints=constraints,options={'rhobeg':0.1, 'catol':0.001})
        
    #     # Extract the optimized values of D and L
    #     self.dd['design']['D'], self.dd['design']['L'] = solution.x
    #     self.dd['design']['zlug'] = (2/3)*self.dd['design']['L'] 
    #     results = self.getAnchorCapacity(plot=False)    
   

    def getSize(self, geom, geomKeys, geomBounds=None, loads=None, minfs={'Ha':1.6,'Va':2}, 
                LD_con=[4,8], fix_zlug=False, FSdiff_max=None, plot=False):
        '''
        
    
        Parameters
        ----------
        geom: list
            starting guess geometry values
        geomKeys : list
            List of keys that match the geom list values i.e. 'L','D','zlug'
        geomBounds : list,optional
            List of upper and lower bounds for each geometry value. 
            Each entry should be a tuple of upper and lower bounds for each geometry i.e. [(5,10),(10,20)]
        loads : dict, optional
            Dictionary of maximum anchor loads in horizontal and vertical directions (not including factor of safety). The default is None.
        minfs : dict,optional
            Minimum factors of safety in horizontal and vertical directions
        LD_con : float
            Constraint for L/D parameter
        fix_zlug : bool
            Boolean to decide if zlug should be altered as geometric values are altered. 
            True = fixed zlug, False = zlug may be changed
        plot : bool
            Boolean controls if capacity plots are generated or not for the final configuration
    
        Returns
        -------
        None.
    
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
        
        def conFun_Drilled(vars, geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs):

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
        
        elif 'drilled' in anchType:
            constraints = [{'type':'ineq','fun':conFun_LD,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)},
                           {'type':'ineq','fun':conFun_Drilled,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)},
                           {'type':'ineq','fun':conFunH,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)},
                           {'type':'ineq','fun':conFunV,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)},
                           {'type':'ineq','fun':conBounds,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)}]
        
        else:
            constraints = [{'type':'ineq','fun':conFunH,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)},
                           {'type':'ineq','fun':conFunV,'args':(geomKeys, input_loads, fix_zlug, LD_con, geomBounds, minfs)}]
        
        # Run the optimization to find sizing that satisfy UC close to 1
        print('optimizing anchor size')
        
        if 'suction' in anchType or 'drilled' in anchType:
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
                if 'suction' in anchType or 'drilled' in anchType:
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
