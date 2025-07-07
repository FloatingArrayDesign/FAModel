"""Anchor class for FAModel, containing information and key methods for anchors of mooring lines
    Work in progress
"""
import moorpy as mp
import numpy as np
from scipy.optimize import minimize
from famodel.famodel_base import Node
from famodel.mooring.mooring import Mooring
import famodel.platform.platform 

class Anchor(Node):
    
    def __init__(self, dd=None, ms=None, r=[0,0,0], aNum=None, id=None, g=9.81, rho=1025):
        '''
        Initialize an Anchor object.

        Parameters
        ----------
        dd : dict
            Design dictionary containing all information on the anchor:
            {
                type : str
                    Anchor type ('plate', 'suction_pile', 'torpedo_pile', 'helical_pile', 'driven_pile', 'dandg_pile')
                design : dict
                    Geometric properties (e.g., A, D, D1, D2, d, L, L1, L2, zlug, beta)
                cost : dict
                    Cost breakdown (matCost, instCost, decomCost)
            }
        ms : MoorPy system object
            The MoorPy system instance the anchor is added to.
        r : list of float
            Location of the anchor in (x, y, z) coordinates (m).
        aNum : int, optional
            Entry number for anchor within the project anchorList dictionary.
        id : str or int, optional
            Unique identifier for the anchor object.
        g : float, optional
            Gravitational acceleration (m/s²). Default is 9.81.
        rho : float, optional
            Water density (kg/m³). Default is 1025.
        '''

        from famodel.famodel_base import Node
        Node.__init__(self, id)

        self.dd = dd
        self.ms = ms
        self.r = r
        self.aNum = aNum
        self.g = g
        self.rho = rho

        # Initialize anchor type and soil type
        self.anchType = dd.get('type') if dd else None
        self.soil_type = None

        # Initialize MoorPy anchor object
        self.mpAnchor = None

        # Extract mass if available
        self.mass = dd.get('design', {}).get('mass', None) if dd else None

        # Initialize other dictionaries
        self.anchorCapacity = {}
        self.cost = {}
        self.loads = {}
        self.soilProps = {}
        self.failure_probability = {}
        self.env_impact = {}
        
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
            Type of mooring line ('chain' or 'wire').
        d : float
            Nominal diameter (m).
        w : float
            Unit weight (N/m).
        nolugload : bool
            True if no lug load transfer should be applied.
        '''
        for att in self.attachments.values():
            if isinstance(att['obj'], Mooring):
                mtype = att['obj'].dd['sections'][0]['type']['material'].lower()
                if 'chain' not in mtype:
                    print('No chain on seafloor, setting Ta=Tm (no load transfer).')
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
       
    def getLugForces(self, ground_conds, Hm, Vm, thetam, zlug, line_type=None, d=None, w=None, plot=False):
        '''
        Calculate the lug forces Ha and Va based on mudline loads.

        Parameters
        ----------
        ground_conds : dict
            Dictionary of ground conditions where keys are soil types.
        Hm : float
            Horizontal mudline load (N).
        Vm : float
            Vertical mudline load (N).
        thetam : float
            Mudline load angle (deg).
        zlug : float
            Padeye embedment depth (m).
        line_type : str, optional
            Type of mooring line ('chain' or 'wire').
        d : float, optional
            Mooring line diameter (m).
        w : float, optional
            Mooring line unit weight (N/m).
        plot : bool, optional
            Whether to plot the load transfer profile.

        Returns
        -------
        Ha : float
            Horizontal load at lug (N).
        Va : float
            Vertical load at lug (N).
        '''
        from famodel.anchors.anchors_famodel_profile.capacity_load import getTransferLoad

        if self.soil_type is None:
            self.soil_type = self.dd.get('design', {}).get('soil_type')
        
        soil_profile = self.dd.get('soil_properties', {}).get(self.soil_type)

        # Determine mudline depth
        z0 = soil_profile[0][0]

        # Load transfer if padeye is embedded
        if zlug > z0:
            # Fallback mechanism for line properties
            if line_type is None or d is None or w is None:
                try:
                    line_type, d, w, _ = self.getLineProperties()
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
        
            loads = getTransferLoad(
                soil_profile, self.soil_type, np.sqrt(Hm**2 + Vm**2),
                thetam, zlug, line_type, d, w, plot=plot
            )
            Ta = loads['Ta']
            thetaa = loads['thetaa']
            Ha = Ta*np.cos(thetaa)
            Va = Ta*np.sin(thetaa)
        
        else:
            Ha = Hm
            Va = Vm
        
        return Ha, Va


    def getCapacityAnchor(self, ground_conds, Hm, Vm, thetam, zlug, line_type=None, d=None, w=None, plot=False):
        '''
        Calculate anchor capacity based on anchor type and ground conditions.

        Parameters
        ----------
        ground_conds : dict
            Dictionary of ground conditions where keys are soil types.
        Hm : float
            Horizontal mudline load (N).
        Vm : float
            Vertical mudline load (N).
        thetam : float
            Mudline load angle (deg).
        zlug : float
            Padeye embedment depth (m).
        line_type : str, optional
            Type of mooring line ('chain' or 'wire').
        d : float, optional
            Mooring line diameter (m).
        w : float, optional
            Mooring line unit weight (N/m).
        plot : bool, optional
            Whether to plot the results.

        Returns
        -------
        None
            Updates anchor object with capacity results.
        '''
        from famodel.anchors.anchors_famodel_profile.capacity_plate import getCapacityPlate
        from famodel.anchors.anchors_famodel_profile.capacity_suction import getCapacitySuction
        from famodel.anchors.anchors_famodel_profile.capacity_torpedo import getCapacityTorpedo
        from famodel.anchors.anchors_famodel_profile.capacity_helical import getCapacityHelical
        from famodel.anchors.anchors_famodel_profile.capacity_driven import getCapacityDriven
        from famodel.anchors.anchors_famodel_profile.capacity_dandg import getCapacityDandG
        from famodel.anchors.anchors_famodel_profile.capacity_load import getTransferLoad
        import numpy as np

        # --- Dispatch dictionary ---
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

        # Normalize anchor type
        anchType_clean = self.anchType.lower().replace(' ', '')

        # Find function
        capacity_func = capacity_dispatch.get(anchType_clean)
        if capacity_func is None:
            raise ValueError(f"Unknown anchor type '{self.anchType}' for anchor capacity calculation.")

        # Get ground conditions
        soil_profile = ground_conds.get(self.soil_type)
        if soil_profile is None:
            raise ValueError(f"Ground condition '{self.soil_type}' not found in provided ground_conds.")

        # Determine if load transfer is needed
        z0 = soil_profile[0][0]
        
        # Load transfer if padeye is embedded
        if zlug > z0:
            if line_type is None or d is None or w is None:
                try:
                    line_type, d, w, nolugload = self.getLineProperties()
                except ValueError:
                    print('[Warning] No mooring attachment found. Trying anchor-level line properties...')
                    line_type = getattr(self, 'line_type', None)
                    d = getattr(self, 'd', None)
                    w = getattr(self, 'w', None)
                    nolugload = False
                
                    if line_type is None or d is None or w is None:
                        print('[Fallback] Using default chain properties.')
                        line_type = 'chain'
                        d = 0.16
                        w = 5000.0

                if nolugload:
                    Ha, Va = Hm, Vm
                else:
                    loads = getTransferLoad(soil_profile, self.soil_type, Tm=np.sqrt(Hm**2 + Vm**2), thetam=thetam, zlug=zlug, line_type=line_type, d=d, w=w, plot=plot)
                    Ta = loads['Ta']
                    thetaa = loads['thetaa']
                    Ha = Ta*np.cos(thetaa)
                    Va = Ta*np.sin(thetaa)
            else:
                loads = getTransferLoad(soil_profile, self.soil_type, Tm=np.sqrt(Hm**2 + Vm**2), thetam=thetam, zlug=zlug, line_type=line_type, d=d, w=w, plot=plot)
                Ta = loads['Ta']
                thetaa = loads['thetaa']
                Ha = Ta*np.cos(thetaa)
                Va = Ta*np.sin(thetaa)
        else:
            Ha = Hm
            Va = Vm
            
        # Call the capacity function based on anchor type
        if anchType_clean == 'suction':
            D = self.dd['design']['D']
            L = self.dd['design']['L']
            zlug = self.dd['design']['zlug']
            
            results = capacity_func(soil_profile, self.soil_type, D, L, zlug, Ha, Va, plot=plot)

        elif anchType_clean in ['sepla', 'dea', 'depla', 'vla', 'plate']:
            results = capacity_func(soil_profile, self.soil_type, self.B, self.L, zlug, self.beta, Ha, Va, plot=plot)

        elif anchType_clean == 'torpedo':
            results = capacity_func(soil_profile, self.soil_type, self.D1, self.D2, self.L1, self.L2, zlug, self.ballast, Ha, Va, plot=plot)

        elif anchType_clean == 'helical':
            results = capacity_func(soil_profile, self.soil_type, self.D, self.L, self.d, zlug, Ha, Va, plot=plot)

        elif anchType_clean == 'driven':
            y, z, results = capacity_func(soil_profile, self.soil_type, self.L, self.D, zlug, Ha, Va, plot=plot)

        elif anchType_clean == 'dandg':
            y, z, results = capacity_func(soil_profile, self.soil_type, self.L, self.D, zlug, Ha, Va, plot=plot)

        else:
            raise ValueError(f"Anchor type '{self.anchType}' not supported.")

        # --- Standardize and store capacity results ---
        self.capacity_results = {
            'Hmax': results.get('Horizontal max.', np.nan),
            'Vmax': results.get('Vertical max.', np.nan),
            'UC': results.get('Unity check', np.nan),
            'Ha': Ha,
            'Va': Va,
            'zlug': zlug,
            'z0': z0
        }

        # Add mass if weight is available
        if 'Weight pile' in results:
            self.capacity_results['Weight pile'] = results['Weight pile']
        if 'Weight plate' in results:
            self.capacity_results['Weight plate'] = results['Weight plate']

        # Special cases for displacement-based anchors
        if anchType_clean in ['driven_pile', 'dandg_pile']:
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
            UCv = self.capacity_results.get('Unity check (vertical)', None)
            UCh = self.capacity_results.get('Unity check (horizontal)', None)

            if UCv is None or UCh is None:
                print("Warning: Vertical or Horizontal Unity Check not found in capacity results. Returning NaN.")
                return {'SF_vertical': np.nan, 'SF_horizontal': np.nan}

            SFv = 1.0/UCv if UCv != 0 else np.inf
            SFh = 1.0/UCh if UCh != 0 else np.inf

            return {'SF_vertical': SFv, 'SF_horizontal': SFh}

        else:
            UC = self.capacity_results.get('UC', None)

            if UC is None:
                print("Warning: Unity Check (UC) not found in capacity results. Returning NaN.")
                return {'SF_combined': np.nan}

            SF = 1.0/UC if UC != 0 else np.inf

            return {'SF_combined': SF} 

    def getCost(self, costDict='default'):
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


    def getSize(self, geom, geomKeys, geomBounds=None, loads=None, minfs={'Ha':1.0,'Va':1.0}, 
                 lambdap_con=[4,6], zlug_fix=False, FSdiff_max=None, plot=False):
        '''
        Resize the anchor dimensions to meet the target safety factor and geometric constraints.
    
        Parameters
        ----------
        geom : list
            Starting guess geometry values.
        geomKeys : list
            List of keys matching the geom list values (e.g., 'L', 'D', 'zlug').
        geomBounds : list, optional
            List of tuples of upper and lower bounds for each geometry value.
        loads : dict, optional
            Dictionary of maximum anchor loads.
        minfs : dict, optional
            Minimum factors of safety in horizontal and vertical directions.
        lambdap_con : list, optional
            Constraint for L/D parameter as [min, max].
        zlug_fix : bool, optional
            Whether zlug should be fixed (True) or updated (False).
        FSdiff_max : dict, optional
            Maximum allowable difference between achieved FS and target FS.
        plot : bool, optional
            Whether to plot results.
    
        Returns
        -------
        None
        '''
        from scipy.optimize import minimize
        import numpy as np
    
        anchType = self.dd['type']
    
        if loads is None:
            loads = self.loads
    
        # Compute thetam internally from Hm and Vm
        Hm = loads['Hm']
        Vm = loads['Vm']
        thetam = np.degrees(np.arctan2(Vm, Hm))
        zlug = loads['zlug']
        
        # Read mooring properties from anchor attributes
        line_type = self.line_type
        d = self.d
        w = self.w
    
        Ha, Va = self.getLugForces(
            ground_conds=self.dd.get('soil_properties'),
            Hm=Hm,
            Vm=Vm,
            thetam=thetam,
            zlug=zlug,
            line_type=line_type,
            d=d,
            w=w,
            plot=plot
        )
    
        ground_conds = self.dd.get('soil_properties')
    
        input_loads = {'Ha': Ha*minfs['Ha'], 'Va': Va*minfs['Va']}
    
        # Objective: minimize weight
        def objective_(vars, geomKeys, Ha, Va, ground_conds, thetam, zlug, plot):
            newGeom = dict(zip(geomKeys, vars))
            self.dd['design'].update(newGeom)
            if 'suction' in self.dd['type'] and not zlug_fix:
                self.dd['design']['zlug'] = (2/3)*newGeom['L']
            results = self.getCapacityAnchor(ground_conds, Ha, Va, thetam, zlug, plot=plot)
            return results.get('Weight pile', results.get('Weight plate', 1e6))
    
        # Constraints
        def conFun_lambdap_(vars, lambdap_con, geomKeys):
            newGeom = dict(zip(geomKeys, vars))
            lambdap = newGeom['L']/newGeom['D']
            return min(lambdap - lambdap_con[0], lambdap_con[1] - lambdap)
    
        def conFun_UC_(vars, Ha, Va, ground_conds, thetam, zlug, plot):
            results = self.getCapacityAnchor(ground_conds, Ha, Va, thetam, zlug, plot=plot)
            return results.get('Unity check', 0) - 1
    
        def conFun_H_(vars, Ha, Va, ground_conds, thetam, zlug, plot):
            results = self.getCapacityAnchor(ground_conds, Ha, Va, thetam, zlug, plot=plot)
            FS = self.getFS()
            return FS.get('SF_horizontal', FS.get('SF_combined', 0)) - 1
    
        def conFun_V_(vars, minfs, Ha, Va, ground_conds, thetam, zlug, plot):
            results = self.getCapacityAnchor(ground_conds, Ha, Va, thetam, zlug, plot=plot)
            FS = self.getFS()
            if minfs['Va'] == 0:
                return 1
            return FS.get('SF_vertical', FS.get('SF_combined', 0)) - 1
    
        # Initial geometry setup
        startGeom = dict(zip(geomKeys, geom))
        self.dd['design'].update(startGeom)
    
        if not 'zlug' in self.dd['design']:
            if 'suction' in anchType and not zlug_fix:
                self.dd['design']['zlug'] = (2/3)*startGeom['L']
            else:
                self.dd['design']['zlug'] = 0
    
        if zlug_fix and 'zlug' in geomKeys:
            zlug_loc = geomKeys.index('zlug')
            geomKeys.pop(zlug_loc)
            geom.pop(zlug_loc)
            if geomBounds:
                geomBounds.pop(zlug_loc)
    
        initial_guess = geom
    
        # Setup constraints
        constraints = []
    
        if 'suction' in anchType:
            constraints.append({'type': 'ineq', 'fun': conFun_lambdap_, 'args': (lambdap_con, geomKeys)})
    
        if 'torpedo' in anchType or 'suction' in anchType:
            constraints.append({'type': 'ineq', 'fun': conFun_UC_, 'args': (Ha, Va, ground_conds, thetam, zlug, plot)})
        else:
            constraints.append({'type': 'ineq', 'fun': conFun_H_, 'args': (Ha, Va, ground_conds, thetam, zlug, plot)})
            constraints.append({'type': 'ineq', 'fun': conFun_V_, 'args': (minfs, Ha, Va, ground_conds, thetam, zlug, plot)})
    
        print('Starting optimization of anchor size')
    
        if geomBounds is None:
            solution = minimize(objective_, initial_guess, args=(geomKeys, Ha, Va, ground_conds, thetam, zlug, plot),
                                 method='COBYLA', constraints=constraints,
                                 options={'rhobeg': 0.1, 'catol': 0.001})
        else:
            solution = minimize(objective_, initial_guess, args=(geomKeys, Ha, Va, ground_conds, thetam, zlug, plot),
                                 method='COBYLA', constraints=constraints, bounds=geomBounds,
                                 options={'rhobeg': 0.1, 'catol': 0.001})
    
        # Update final geometry
        endGeom = dict(zip(geomKeys, solution.x))
        print('Optimized geometry: ', endGeom)
        self.dd['design'].update(endGeom)
    
        if 'suction' in anchType and not zlug_fix:
            self.dd['design']['zlug'] = (2/3)*self.dd['design']['L']
            zlug = self.dd['design']['zlug']  # update local zlug
    
        self.getCapacityAnchor(ground_conds, Ha, Va, thetam, zlug, plot=plot)
        
    def getSizeSuction(self, geom, geomKeys, geomBounds=None, loads=None, minfs={'Ha':1.6,'Va':2}, 
                 lambdap_con=[4,8], zlug_fix=False, plot=False):
        '''
        Two-stage optimization:
        Stage 1 - Grow anchor to satisfy UC <= 1.
        Stage 2 - Minimize weight while keeping UC <= 1 and satisfying L/D constraints.
        '''
    
        anchType = self.dd['type']
    
        if loads is None:
            loads = self.loads
    
        Hm = loads['Hm']
        Vm = loads['Vm']
        zlug = self.dd['design']['zlug']
        thetam = np.degrees(np.arctan2(Vm, Hm))
    
        line_type = self.line_type
        d = self.d
        w = self.w
    
        initial_guess = [self.dd['design']['L'], self.dd['design']['D']]
        bounds = geomBounds if geomBounds else [(5.0, 30.0), (2.0, 5.0)]

        ground_conds = self.dd.get('soil_properties')
    
        # --- Stage 1: Safety First ---
        def safety_objective(vars):
            L, D = vars
            self.dd['design']['L'] = L
            self.dd['design']['D'] = D
            self.dd['design']['zlug'] = (2/3) * L
    
            Ha, Va = self.getLugForces(
                ground_conds=ground_conds,
                Hm=Hm,
                Vm=Vm,
                thetam=thetam,
                zlug=self.dd['design']['zlug'],
                line_type=line_type,
                d=d,
                w=w,
                plot=False
            )
    
            self.getCapacityAnchor(
                ground_conds=ground_conds,
                Hm=Hm,
                Vm=Vm,
                thetam=thetam,
                zlug=self.dd['design']['zlug'],
                line_type=line_type,
                d=d,
                w=w,
                plot=False
            )
    
            UC = self.capacity_results.get('UC', 2.0)
            return max(0.0, UC - 1.0)**2
    
        minimize(
            safety_objective,
            initial_guess,
            method='COBYLA',
            bounds=bounds,
            options={'rhobeg': 0.1, 'catol': 0.001, 'maxiter': 300}
        )
    
        # --- Stage 2: Weight Minimization ---
        def weight_objective(vars):
            L, D = vars
            self.dd['design']['L'] = L
            self.dd['design']['D'] = D
            self.dd['design']['zlug'] = (2/3) * L
    
            Ha, Va = self.getLugForces(
                ground_conds=ground_conds,
                Hm=Hm,
                Vm=Vm,
                thetam=thetam,
                zlug=self.dd['design']['zlug'],
                line_type=line_type,
                d=d,
                w=w,
                plot=False
            )
    
            self.getCapacityAnchor(
                ground_conds=ground_conds,
                Hm=Hm,
                Vm=Vm,
                thetam=thetam,
                zlug=self.dd['design']['zlug'],
                line_type=line_type,
                d=d,
                w=w,
                plot=False
            )
    
            return self.capacity_results.get('Weight pile', 1e9)
    
        def constraint_uc(vars):
            L, D = vars
            return 1.0 - self.capacity_results.get('UC', 2.0)
        
        def constraint_fs_horizontal(vars):
            L, D = vars
            return (self.capacity_results.get('Hmax', 0) / self.capacity_results.get('Ha', 1)) - minfs['Ha']
        
        def constraint_fs_vertical(vars):
            L, D = vars
            return (self.capacity_results.get('Vmax', 0) / self.capacity_results.get('Va', 1)) - minfs['Va']

        def constraint_lambda_min(vars):
            L, D = vars
            return (L/D) - lambdap_con[0]
    
        def constraint_lambda_max(vars):
            L, D = vars
            return lambdap_con[1] - (L/D)
    
        result = minimize(
            weight_objective,
            [self.dd['design']['L'], self.dd['design']['D']],
            method='COBYLA',
            constraints=[
                {'type': 'ineq', 'fun': constraint_fs_horizontal},
                {'type': 'ineq', 'fun': constraint_fs_vertical},
                {'type': 'ineq', 'fun': constraint_lambda_min},
                {'type': 'ineq', 'fun': constraint_lambda_max}
            ],
            bounds=bounds,
            options={'rhobeg': 0.5, 'catol': 0.01, 'maxiter': 100}
        )
    
        # Update final geometry
        endGeom = dict(zip(geomKeys, result.x))
        print('Optimized geometry:', endGeom)
        self.dd['design'].update(endGeom)
    
        if 'suction' in anchType and not zlug_fix:
            self.dd['design']['zlug'] = (2/3) * self.dd['design']['L']
    
        self.getCapacityAnchor(
            ground_conds=ground_conds,
            Hm=Hm,
            Vm=Vm,
            thetam=thetam,
            zlug=self.dd['design']['zlug'],
            line_type=line_type,
            d=d,
            w=w,
            plot=plot
        )
    
        print('\nFinal Optimized Anchor:')
        print('Design:', self.dd['design'])
        print('Capacity Results:', self.capacity_results)

    def getCombinedPlot(self):
        '''
        Create a single plot showing the suction pile and the inverse catenary overlay in the same coordinate system.
        '''
        from famodel.anchors.anchors_famodel_profile.capacity_load import getTransferLoad
        from famodel.anchors.anchors_famodel_profile.capacity_plots import plot_suction
        import numpy as np
        import matplotlib.pyplot as plt
    
        if self.anchType.lower() != 'suction':
            raise NotImplementedError("getCombinedPlot only supports suction piles for now.")
    
        # Extract key parameters
        design = self.dd['design']
        D = design['D']
        L = design['L']
        zlug = design['zlug']
        soil_type = self.soil_type
        soil_profile = self.dd['soil_properties'][soil_type]
        z0 = soil_profile[0][0]
    
        Hm = self.loads['Hm']
        Vm = self.loads['Vm']
        thetam = self.loads.get('thetam', np.degrees(np.arctan2(Vm, Hm)))
    
        line_type = getattr(self, 'line_type', 'chain')
        d = getattr(self, 'd', 0.16)
        w = getattr(self, 'w', 5000.0)
    
        # Get inverse catenary path
        result = getTransferLoad(
            soil_profile, soil_type, np.sqrt(Hm**2 + Vm**2), thetam, zlug, line_type, d, w, plot=False
        )
        drag_values  = np.array(result['drag_values'])
        depth_values = np.array(result['depth_values'])
        depth_values = -depth_values[::-1]
        Tm = result['Tm']; thetam = result['thetam']
        Ta = result['Ta']; thetaa = result['thetaa']        
    
        # Transform to suction pile coordinate system
        x_start = D/2 + drag_values[0]
        z_start = zlug
        drag_transformed  = x_start - drag_values
        depth_transformed = z_start + (depth_values - depth_values[0])
    
        # Set up plot
        fig, ax = plt.subplots(figsize=(5, 5))
    
        # Plot suction pile
        plot_suction(soil_profile, soil_type, L, D, zlug=zlug, title='', ax=ax)
    
        # Overlay inverse catenary
        ax.plot(drag_transformed, depth_transformed, color='b', lw=2.0, label='Inverse Catenary')
    
        # Annotate ends
        ax.plot(drag_transformed[-1], depth_transformed[-1], 'ro', label='Mudline End')
        ax.plot(drag_transformed[0], depth_transformed[0], 'go', label='Embedded End')
        
        n = 2e6
        # Add load vectors
        ax.arrow(drag_transformed[-1], depth_transformed[-1], Tm*np.cos(np.deg2rad(thetam))/n, -Tm*np.sin(np.deg2rad(thetam))/n,
                 head_width=0.25, head_length=0.5, color='r', label='Mudline Load')
        ax.arrow(drag_transformed[0], depth_transformed[0], Ta*np.cos(thetaa)/n, -Ta*np.sin(thetaa)/n,
                 head_width=0.25, head_length=0.5, color='g', label='Padeye Load')
    
        # Finalize plot
        xmax = max(drag_transformed[-1] + D, 2*D)
        ax.set_xlim(-D, xmax)
        ax.set_title('Suction Pile with Inverse Catenary')
        ax.legend()
        ax.grid(True)
        plt.tight_layout()
        plt.show()
