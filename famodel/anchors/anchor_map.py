"""Anchor class for FAModel, containing information and key methods for anchors of mooring lines
    Work in progress
"""
import moorpy as mp
import numpy as np
from scipy.optimize import minimize
from famodel.famodel_base import Node
from famodel.mooring.mooring import Mooring
import famodel.platform.platform 
import matplotlib.pyplot as plt

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

        self.anchType = dd.get('type') if dd else None
        self.soil_type = None
        self.soil_profile = None
        self.profile_name = None
        self.soil_type_list = []

        self.mpAnchor = None
        self.capacity_format = None
        self.mass = dd.get('design', {}).get('mass', None) if dd else None

        self.anchorCapacity = {}
        self.cost = {}
        self.loads = {}
        self.soilProps = {}
        self.failure_probability = {}
        self.env_impact = {}

        # --- Assign soil profile if map is provided ---
        if profile_map is not None:
            self.setSoilProfile(profile_map)

    def setSoilProfile(self, profile_map):
        '''
        Assign a bilinearly interpolated soil profile from the 4 nearest CPTs.
    
        Parameters
        ----------
        profile_map : list of dict
            Each CPT must have keys: 'x', 'y', and 'layers'
    
        Returns
        -------
        None
        '''
        import numpy as np
    
        x_anchor, y_anchor = self.r[0], self.r[1]
    
        # Sort all CPTs by distance
        distances = [np.hypot(p['x'] - x_anchor, p['y'] - y_anchor) for p in profile_map]
        idx_sorted = np.argsort(distances)
        CPTs = [profile_map[i] for i in idx_sorted[:4]]
       
        # Extract positions and weights (inverse distance squared)
        x = np.array([cpt['x'] for cpt in CPTs])
        y = np.array([cpt['y'] for cpt in CPTs])
        dx = x - x_anchor
        dy = y - y_anchor
        d = np.hypot(dx, dy)
        w = 1/np.maximum(d, 1e-3)**2
        w /= np.sum(w)
    
        # Interpolate layer-by-layer
        layers_list = [cpt['layers'] for cpt in CPTs]
        n_layers = len(layers_list[0])
        interpolated_layers = []
    
        for i in range(n_layers):
            layer = {'soil_type': layers_list[0][i]['soil_type']}
            keys = layers_list[0][i].keys()
    
            for key in keys:
                if key == 'soil_type':
                    continue
                if all(key in l[i] for l in layers_list):
                    vals = [l[i][key] for l in layers_list]
                    layer[key] = np.dot(w, vals)
    
            interpolated_layers.append(layer)
    
        self.soil_profile = interpolated_layers
        self.profile_name = f'Interpolated_2D'
    
        # Assign soil type
        soil_types = [layer['soil_type'] for layer in self.soil_profile]
        self.soil_type_list = list(set(soil_types))
        self.soil_type = soil_types[0] if len(self.soil_type_list) == 1 else 'mixed'
    
        print(f"[Anchor] Assigned interpolated soil profile: {self.profile_name} weighting with soil types {self.soil_type_list}")

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
        import moorpy as mp

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

        # Set the point as an anchor entity
        self.mpAnchor.entity = {'type': 'anchor'}
        if 'type' in self.dd:
            self.mpAnchor.entity['anchor_type'] = self.dd['type']

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

    def getMudlineForces(self, lines_only=False, seabed=True, xyz=False):
        '''
        Find forces on anchor at mudline using the platform.getWatchCircle method or MoorPy Point.getForces method.

        Parameters
        ----------
        lines_only : boolean, optional
            Calculate forces from just mooring lines (True) or not (False). Default is False.
        seabed : boolean, optional
            Include effect of seabed pushing up the anchor (True) or not (False). Default is True.
        xyz : boolean, optional
            Return forces in x,y,z DOFs (True) or only the enabled DOFs (False). Default is False.

        Returns
        -------
        dict
            Dictionary containing mudline forces.
        '''
        loads = self.mpAnchor.getForces(lines_only=lines_only, seabed=seabed, xyz=xyz)

        self.loads = {
            'Hm': np.sqrt(loads[0]**2 + loads[1]**2),
            'Vm': loads[2],
            'thetam': np.degrees(np.arctan2(loads[2], np.sqrt(loads[0]**2 + loads[1]**2))),
            'method': 'static',
            'mudline_load_type': 'current_state'
        }

        return self.loads
       
    def getLugForces(self, Hm, Vm, zlug, line_type=None, d=None, w=None, plot=False):
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
        from famodel.anchors.anchors_famodel_map.capacity_load_map import getTransferLoad

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
                        w = 5000.0

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

        return layers, Ha, Va

    def getCapacityAnchor(self, Hm, Vm, zlug, line_type=None, d=None, w=None, plot=False):
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
        plot : bool, optional
            Whether to plot the load transfer and pile geometry
    
        Returns
        -------
        results : dict
            Capacity results dictionary from the selected capacity function.
        '''
        from famodel.anchors.anchors_famodel_map.capacity_plate_map import getCapacityPlate
        from famodel.anchors.anchors_famodel_map.capacity_suction_map import getCapacitySuction
        from famodel.anchors.anchors_famodel_map.capacity_torpedo_map import getCapacityTorpedo
        from famodel.anchors.anchors_famodel_map.capacity_helical_map import getCapacityHelical
        from famodel.anchors.anchors_famodel_map.capacity_driven_map import getCapacityDriven
        from famodel.anchors.anchors_famodel_map.capacity_dandg_map import getCapacityDandG
        from famodel.anchors.anchors_famodel_map.capacity_load_map import getTransferLoad
        from famodel.anchors.anchors_famodel_map.capacity_plots_map import plot_load
        # import numpy as np
    
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
                        w = 5000.0
    
            else:
                layers, resultsLoad = getTransferLoad(
                    profile_map=[{'layers': soil_profile}],
                    Tm=np.sqrt(Hm**2 + Vm**2),
                    thetam=np.degrees(np.arctan2(Vm, Hm)),
                    zlug=zlug,
                    line_type=line_type,
                    d=d,
                    w=w,
                    plot=False
                )
                if plot:
                    plot_load(
                        layers,
                        resultsLoad['drag_values'],
                        resultsLoad['depth_values'],
                        resultsLoad['Tm'],
                        resultsLoad['thetam'],
                        resultsLoad['Ta'],
                        resultsLoad['thetaa'],
                        zlug=zlug
                    )
    
                Ta = resultsLoad['Ta']
                thetaa = resultsLoad['thetaa']
                Ha = Ta*np.cos(np.deg2rad(thetaa))
                Va = Ta*np.sin(np.deg2rad(thetaa))
                
                print(f'Input Hm = {Hm}, Vm = {Vm}, zlug = {zlug}')
                print(f'Input Ha = {Ha}, Va = {Va}, zlug = {zlug}')
                print(f'Input Ta = {Ta}, thetaa = {(thetaa)}')
        else:
            Ha = Hm
            Va = Vm
    
        # --- Call the appropriate capacity function ---
        if anchType_clean in ['sepla', 'dea', 'depla', 'vla', 'plate']:
            self.capacity_format = 'plate'
            B = self.dd['design']['B']
            L = self.dd['design']['L']
            beta = self.dd['design'].get('beta', 0.0)
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
        self.capacity_results = {
            'Hmax': results.get('Horizontal max.', np.nan),
            'Vmax': results.get('Vertical max.', np.nan),
            'UC': results.get('Unity check', np.nan),
            'Ha': Ha,
            'Va': Va,
            'zlug': zlug,
            'z0': z0
        }
    
        if 'Weight pile' in results:
            self.capacity_results['Weight pile'] = results['Weight pile']
        if 'Weight plate' in results:
            self.capacity_results['Weight plate'] = results['Weight plate']
    
        if anchType_clean in ['driven', 'dandg']:
            self.capacity_results['Lateral displacement'] = results.get('Lateral displacement', np.nan)
            self.capacity_results['Rotational displacement'] = results.get('Rotational displacement', np.nan)
    
        return results
                   
    def getSafetyFactor(self):
        '''
        Calculate the safety factor based on the unity checks stored in capacity results.

        Returns
        -------
        dict
            Dictionary containing safety factors.
        '''
        anchType_clean = self.anchType.lower().replace(' ', '')

        if anchType_clean in ['helical', 'driven', 'dandg']:
            UC_v = self.capacity_results.get('Unity check (vertical)', None)
            UC_h = self.capacity_results.get('Unity check (horizontal)', None)

            if UC_v is None or UC_h is None:
                print("Warning: Vertical or horizontal unity check (UC) not found in capacity results. Returning NaN.")
                return {'SF_vertical': np.nan, 'SF_horizontal': np.nan}

            SF_v = 1.0/UC_v if UC_v != 0 else np.inf
            SF_h = 1.0/UC_h if UC_h != 0 else np.inf

            return {'SF_vertical': SF_v, 'SF_horizontal': SF_h}

        else:
            UC = self.capacity_results.get('UC', None)

            if UC is None:
                print("Warning: Unity check (UC) not found in capacity results. Returning NaN.")
                return {'SF_combined': np.nan}

            SF = 1.0/UC if UC != 0 else np.inf

            return {'SF_combined': SF} 

    def getCostAnchor(self, costDict='default'):
        '''
        Calculate the cost of the anchor based on material, installation, and decommissioning costs.

        Parameters
        ----------
        costDict : str or dict, optional
            If 'default', uses mean values from Task 49 Design Basis ranges.
            If dict or yaml path, loads user-defined cost dictionaries.

        Returns
        -------
        float
            Total cost of the anchor.
        '''
        if isinstance(costDict, str) and costDict != 'default':
            import yaml
            costDict = yaml.load(costDict, Loader=yaml.FullLoader)

        anchType = self.dd['type']

        if costDict == 'default':
            matCostDict = {
                'suction_pile': 4.435,
                'DEA': 5.705,
                'SEPLA': 5.705,
                'DEPLA': 5.705,
                'VLA': 5.705,
                'torpedo_pile': 5.0,
                'helical_pile': 6.0,
                'driven_pile': 4.0,
                'dandg_pile': 5.5
            }
            instCostDict = {
                'suction_pile': 2.0,
                'DEA': 1.5,
                'SEPLA': 1.5,
                'DEPLA': 1.5,
                'VLA': 1.5,
                'torpedo_pile': 2.5,
                'helical_pile': 3.0,
                'driven_pile': 2.0,
                'dandg_pile': 2.2
            }
            decomCostDict = {
                'suction_pile': 1.0,
                'DEA': 0.8,
                'SEPLA': 0.8,
                'DEPLA': 0.8,
                'VLA': 0.8,
                'torpedo_pile': 1.2,
                'helical_pile': 1.5,
                'driven_pile': 1.0,
                'dandg_pile': 1.1
            }
        else:
            matCostDict = costDict.get('material', {})
            instCostDict = costDict.get('install', {})
            decomCostDict = costDict.get('decom', {})

        keyFail = True

        # Ensure mass is available
        if self.mass is None or self.mass == 0:
            # Try to extract from capacity_results if already available
            if 'Weight pile' in self.capacity_results:
                self.mass = self.capacity_results['Weight pile']/self.g
            elif 'Weight plate' in self.capacity_results:
                self.mass = self.capacity_results['Weight plate']/self.g
            else:
                # If capacity_results missing, attempt to calculate capacity to retrieve weight
                if 'soil_properties' in self.dd:
                    self.getAnchorCapacity(plot=False)
                    if 'Weight pile' in self.capacity_results:
                        self.mass = self.capacity_results['Weight pile']/self.g
                    elif 'Weight plate' in self.capacity_results:
                        self.mass = self.capacity_results['Weight plate']/self.g
                    else:
                        print('Warning: Weight not found after capacity calculation, setting mass to 0.')
                        self.mass = 0
                else:
                    print('Soil properties needed to calculate anchor mass for cost. Setting mass to 0.')
                    self.mass = 0

        # Calculate material cost based on mass
        if anchType in matCostDict:
            self.cost['Material Cost'] = matCostDict[anchType]*self.mass
            keyFail = False
        else:
            raise KeyError(f'Anchor type {anchType} not found in material cost dictionary.')

        # Install and decom costs if available
        self.cost['Installation Cost'] = instCostDict.get(anchType, 0.0)
        self.cost['Decommissioning Cost'] = decomCostDict.get(anchType, 0.0)

        # Total cost
        self.cost['Total Cost'] = (self.cost['Material Cost'] +
                                   self.cost['Installation Cost'] +
                                   self.cost['Decommissioning Cost'])

        return sum(self.cost.values())
        
    def getSizeAnchor(self, geom, geomKeys, geomBounds=None, loads=None,
                      minfs={'Ha': 1.6, 'Va': 2.0}, lambdap_con=[4, 8], zlug_fix=False, plot=False):
        '''
        Generalized optimization method for all anchor types.
        '''
        
        anchType_clean = self.dd['type'].lower().replace(' ', '')
        
        if loads is None:
            loads = self.loads
    
        Hm = loads['Hm']
        Vm = loads['Vm']
    
        line_type = getattr(self, 'line_type', 'chain')
        d = getattr(self, 'd', 0.16)
        w = getattr(self, 'w', 5000.0)
        
        def update_zlug_if_suction():
            if anchType_clean == 'suction' and not zlug_fix and 'zlug' not in geomKeys:
                self.dd['design']['zlug'] = (2/3)*self.dd['design']['L']   

        # --- Stage 1: Safety Optimization to reach UC <= 1 ---
        def safety_objective(vars):
            for i, key in enumerate(geomKeys):
                self.dd['design'][key] = vars[i]
            update_zlug_if_suction()
    
            _, Ha, Va = self.getLugForces(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'], line_type=line_type, d=d, w=w, plot=False)
            
            self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'], line_type=line_type, d=d, w=w, plot=False)
    
            if self.capacity_format == 'envelope':
                UC = self.capacity_results.get('UC', 2.0)
            elif self.capacity_format == 'component':
                UC = max(
                    self.capacity_results.get('Unity check (horizontal)', 2.0),
                    self.capacity_results.get('Unity check (vertical)', 2.0)
                )
            elif self.capacity_format == 'plate':
                UC = self.capacity_results.get('UC', 2.0)
            else:
                UC = 2.0
            return (UC - 1.0)**2
    
        minimize(
            safety_objective,
            geom,
            method='COBYLA',
            bounds=geomBounds if geomBounds else None,
            options={'rhobeg': 0.02, 'catol': 0.001, 'maxiter': 1500}
        )
    
        # --- Stage 2: Weight Minimization with constraints ---
        if anchType_clean != 'torpedo':    
            def weight_objective(vars):
                for i, key in enumerate(geomKeys):
                    self.dd['design'][key] = vars[i]
                update_zlug_if_suction()
    
                _, Ha, Va = self.getLugForces(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'], line_type=line_type, d=d, w=w, plot=False)
        
                self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'], line_type=line_type, d=d, w=w, plot=False)
        
                return self.capacity_results.get('Weight pile',
                       self.capacity_results.get('Weight plate',
                       self.capacity_results.get('Weight', 1e9)))
        
            def constraint_uc_envelope(vars):
                return 1.0 - self.capacity_results.get('UC', 2.0)
        
            def constraint_uc_component(vars):
                return 1.0 - max(
                    self.capacity_results.get('Unity check (horizontal)', 2.0),
                    self.capacity_results.get('Unity check (vertical)', 2.0)
                )
        
            def constraint_fs_horizontal(vars):
                return (self.capacity_results.get('Horizontal max.', 0)/self.capacity_results.get('Ha', 1)) - minfs['Ha']
        
            def constraint_fs_vertical(vars):
                return (self.capacity_results.get('Vertical max.', 0)/self.capacity_results.get('Va', 1)) - minfs['Va']
        
            def constraint_lambda_min(vars):
                anchType_clean = self.dd['type'].lower().replace(' ', '')
                
                if anchType_clean == 'torpedo':
                    L = self.dd['design']['L1'] + self.dd['design']['L2']
                    A_wing = (self.dd['design']['D1'] - self.dd['design']['D2'])*self.dd['design']['L1']
                    A_shaft = self.dd['design']['D2']*L
                    D = (A_wing + A_shaft)/L
                elif anchType_clean in ['driven', 'dandg', 'helical', 'suction']:
                    L = self.dd['design']['L']
                    D = self.dd['design']['D']
                elif anchType_clean in ['plate', 'sepla', 'dea', 'depla', 'vla']:
                    L = self.dd['design']['L']
                    D = self.dd['design']['B']  
                else:
                    raise ValueError(f'lambda constraints not defined for anchor type: {anchType_clean}')
                return (L/D) - lambdap_con[0]
        
            def constraint_lambda_max(vars):
                anchType_clean = self.dd['type'].lower().replace(' ', '')
                
                if anchType_clean == 'torpedo':
                    L = self.dd['design']['L1'] + self.dd['design']['L2']
                    A_wing = (self.dd['design']['D1'] - self.dd['design']['D2'])*self.dd['design']['L1']
                    A_shaft = self.dd['design']['D2']*L
                    D = (A_wing + A_shaft)/L
                elif anchType_clean in ['driven', 'dandg', 'helical', 'suction']:
                    L = self.dd['design']['L']
                    D = self.dd['design']['D']
                elif anchType_clean in ['plate', 'sepla', 'dea', 'depla', 'vla']:
                    L = self.dd['design']['L']
                    D = self.dd['design']['B']  # use plate width
                else:
                    raise ValueError(f'lambda constraints not defined for anchor type: {anchType_clean}')
                return lambdap_con[1] - (L/D)
        
            if self.capacity_format == 'envelope':
                constraints = [
                    {'type': 'ineq', 'fun': constraint_uc_envelope},
                    {'type': 'ineq', 'fun': constraint_fs_horizontal},
                    {'type': 'ineq', 'fun': constraint_fs_vertical},
                    {'type': 'ineq', 'fun': constraint_lambda_min},
                    {'type': 'ineq', 'fun': constraint_lambda_max},
                ]
            elif self.capacity_format == 'component':
                constraints = [
                    {'type': 'ineq', 'fun': constraint_uc_component},
                    {'type': 'ineq', 'fun': constraint_lambda_min},
                    {'type': 'ineq', 'fun': constraint_lambda_max},
                ]
            elif self.capacity_format == 'plate':
                constraints = [
                    {'type': 'ineq', 'fun': constraint_uc_envelope}
                ]
            else:
                raise ValueError(f"Unknown capacity_format: {self.capacity_format}")
        
            result = minimize(
                weight_objective,
                [self.dd['design'][key] for key in geomKeys],
                method='COBYLA',
                constraints=constraints,
                bounds=geomBounds if geomBounds else None,
                options={'rhobeg': 0.5, 'catol': 0.01, 'maxiter': 100}
            )
        
            endGeom = dict(zip(geomKeys, result.x))
            print('Optimized geometry:', endGeom)
            self.dd['design'].update(endGeom)
        
            update_zlug_if_suction()
        
            self.getCapacityAnchor(Hm=Hm, Vm=Vm, zlug=self.dd['design']['zlug'], line_type=line_type, d=d, w=w, plot=plot)
        
            print('\nFinal Optimized Anchor:')
            print('Design:', self.dd['design'])
            print('Capacity Results:', self.capacity_results)


    def getCombinedPlot(self):
        '''
        Create a plot showing the suction pile and the inverse catenary overlay in the same coordinate system.
        '''
        from famodel.anchors.anchors_famodel_map.capacity_load_map import getTransferLoad
        from famodel.anchors.anchors_famodel_map.capacity_plots_map import plot_suction

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
