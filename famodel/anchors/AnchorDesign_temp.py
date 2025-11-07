# -*- coding: utf-8 -*-
"""
temp storage of different anchor sizing functions that use different 
optimization methods. These were built-in methods to Anchor class.

Eventually can be converted into a full AnchorDesign class...
"""

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
               
def getSizeAnchor2(self, geom, geomBounds=None, loads=None, lambdap_con=[3, 6], 
                   zlug_fix=True, safety_factor={}, plot=False):
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

    sf_uc = safety_factor.get('SF_combined', 1.0)
    sf_Hm = safety_factor.get('Hm', 1.0)
    sf_Vm = safety_factor.get('Vm', 1.0)

    Hm = loads['Hm']*sf_Hm
    Vm = loads['Vm']*sf_Vm

    line_type = getattr(self, 'line_type', 'chain')
    d = getattr(self, 'd', 0.16)
    w = getattr(self, 'w', 5000.0)

    if anchType_clean not in ['suction', 'torpedo', 'plate']:
        raise ValueError(f"Grid-based getSizeAnchor only supports envelope anchors, not '{anchType_clean}'")

    UC_target = 1.0/sf_uc

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
                plot=False)
            
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