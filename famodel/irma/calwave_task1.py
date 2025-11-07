# calwave_task1.py
# Build CalWave Task 1 (Anchor installation) following the theory flow:
# 1) addAction → structure only (type, name, objects, deps)
# 2) evaluateAssets → assign vessels/roles (+ durations/costs)
# 3) (schedule/plot handled by your existing tooling)
import matplotlib.pyplot as plt
from famodel.project import Project
from calwave_irma import Scenario
import calwave_chart as chart
# from calwave_task import Task  # calwave_task module (Felipe)
from task import Task            # generic Task module ( Rudy )

import matplotlib.pyplot as plt

sc = Scenario()  # now sc exists in *this* session

# ---------- Core builder ----------
def build_task1_calwave(sc: Scenario, project: Project):
    """
    Creates Task 1 actions + dependencies (no scheduling/plotting here).
    """

    # --- Pre-ops ---
    mob_sd = sc.addAction('mobilize', 'mobilize_SanDiego')
    linehaul_convoy = sc.addAction(
        'transit_linehaul_tug', 'linehaul_to_site_convoy',
        dependencies=[mob_sd])
    
    mob_by = sc.addAction(
        'mobilize', 'mobilize_Beyster', 
        #dependencies=[mob_sd]
        )
    linehaul_by  = sc.addAction(
        'transit_linehaul_self', 'linehaul_to_site_Beyster',
        dependencies=[mob_sd]) 
    
    # --- Compute anchor centroid (x,y) for first onsite leg start ---
    anchors_all = list(project.anchorList.values())
    rs = [getattr(a, 'r', None) for a in anchors_all if getattr(a, 'r', None) is not None]
    xs = [float(r[0]) for r in rs]
    ys = [float(r[1]) for r in rs]
    anchor_centroid = (sum(xs)/len(xs), sum(ys)/len(ys)) if xs and ys else None
    try:
        print('[task1] anchor_centroid =', anchor_centroid)
    except Exception:
        pass

    # --- On-site (domain objects REQUIRED) ---
    installs, onsite_tug, onsite_by, monitors = [], [], [], []
    
    # first convoy leg starts after the linehaul convoy reaches site
    prev_for_next_tug = linehaul_convoy
    # Beyster’s first in-field leg starts after her own linehaul
    prev_for_by       = linehaul_by
    
    for i, (key, anchor) in enumerate(project.anchorList.items(), start=1):
        # 1) Onsite convoy (tug + barge) to this anchor
        a_tug = sc.addAction(
            'transit_onsite_tug', f'transit_convoy-{key}',
            objects=[anchor],
            dependencies=[prev_for_next_tug]      # first = linehaul_convoy; then = previous install
        )
        
        # 2) Beyster to this anchor (after previous monitor), independent of tug
        a_by = sc.addAction(
            'transit_onsite_self', f'transit_Beyster-{key}',
            objects=[anchor],
            dependencies=[prev_for_by, prev_for_next_tug]
        )
        
        # Inject centroid for the FIRST onsite legs only (centroid → first anchor)
        if i == 1 and anchor_centroid is not None:
            a_by.meta = getattr(a_by, 'meta', {}) or {}
            a_by.meta['anchor_centroid'] = anchor_centroid
            a_tug.meta = getattr(a_tug, 'meta', {}) or {}
            a_tug.meta['anchor_centroid'] = anchor_centroid        
    
        # 3) Install at this anchor (wait for both tug+barge and Beyster on station)
        a_inst = sc.addAction(
            'install_anchor', f'install_anchor-{key}',
            objects=[anchor],
            dependencies=[a_tug, a_by]
        )
        
        # 4) Monitor at this anchor (while anchor is installed)
        a_mon = sc.addAction(
            'monitor_installation', f'monitor_installation-{key}',
            objects=[anchor],
            dependencies=[a_tug]
        )
    
        # collect handles
        onsite_tug.append(a_tug)
        installs.append(a_inst)
        onsite_by.append(a_by)
        monitors.append(a_mon)
    
        # chain next legs:
        prev_for_next_tug = a_inst   # next convoy starts from this installed anchor
        prev_for_by       = a_mon    # or set to a_inst if you want Beyster to move immediately post-install


    # --- Post-ops (objectless) ---
    linehome_convoy  = sc.addAction(
        'transit_linehaul_tug', 'linehaul_to_home_convoy',
        dependencies=monitors)
        
    linehome_by  = sc.addAction(
        'transit_linehaul_self', 'transit_to_home_Beyster',
        dependencies=monitors)
    
    # --- Post-ops ---
    demob_sd  = sc.addAction(
        'demobilize', 'demobilize_SanDiego',
        dependencies=[linehome_convoy])
    
    demob_by  = sc.addAction(
        'demobilize', 'demobilize_Beyster',
        dependencies=[linehome_by])

    # Return a simple list for downstream evaluate/schedule/plot steps
    return {
        'mobilize': [mob_sd, mob_by],
        'linehaul_to_site': [linehaul_convoy, linehaul_by],
        'install': installs,
        'onsite_tug': onsite_tug,
        'onsite_by': onsite_by,
        'monitor': monitors,
        'linehaul_to_home': [linehome_convoy, linehome_by],
        'demobilize': [demob_sd, demob_by]}

# ---------- Assignment step (assign vessels & durations) ----------
def assign_actions(sc: Scenario, actions: dict):
    """
    Assign vessels/roles and set durations where the evaluator doesn't.
    Keeps creation and evaluation clearly separated.
    """
    V = sc.vessels  # shorthand

    # Mobilize
    actions['mobilize'][0].assignAssets({'operator': V['San_Diego']})
    actions['mobilize'][1].assignAssets({'operator': V['Beyster']})
    
    # Transit to site
    convoy_to_site, beyster_to_site = actions['linehaul_to_site']
    convoy_to_site.assignAssets({'carrier': V['Jag'], 'operator': V['San_Diego']})
    beyster_to_site.assignAssets({'vessel': V['Beyster']})

    # Onsite convoy (tug+barge)
    for a_tug in actions['onsite_tug']:
        a_tug.assignAssets({'carrier': V['Jag'], 'operator': V['San_Diego']})
        
    # Install (Jag carries, San_Diego operates the install)    
    for a_inst in actions['install']:
        a_inst.assignAssets({'carrier': V['Jag'], 'operator': V['San_Diego']})

    # Onsite self-propelled (Beyster)    
    for a_by in actions['onsite_by']:
        a_by.assignAssets({'vessel': V['Beyster']})

    # Monitor (Beyster as support)
    for a_mon in actions['monitor']:
        a_mon.assignAssets({'support': V['Beyster']})

    # Transit to home
    convoy_to_home, beyster_to_home = actions['linehaul_to_home']
    convoy_to_home.assignAssets({'carrier': V['Jag'], 'operator': V['San_Diego']})
    beyster_to_home.assignAssets({'vessel': V['Beyster']})
    
    # Demobilize
    actions['demobilize'][0].assignAssets({'operator': V['San_Diego']})
    actions['demobilize'][1].assignAssets({'operator': V['Beyster']})
        

if __name__ == '__main__':
    # 1) Load ontology that mirrors the sample schema (mooring_systems + mooring_line_configs)
    project = Project(file='calwave_ontology.yaml', raft=False)
    project.getMoorPyArray(cables=1)
    
    # 2) Scenario with CalWave catalogs
    sc = Scenario()
    
    # 3) Build (structure only)
    actions = build_task1_calwave(sc, project)
    
    # 4) Assign (assign vessels/roles)
    assign_actions(sc, actions)
    
    # # 5) schedule once, in the Task
    # calwave_task1 = Task.from_scenario(
    #     sc,
    #     name='calwave_task1',
    #     strategy='earliest',               # 'earliest' or 'levels'
    #     enforce_resources=False,         # keep single-resource blocking if you want it
    #     resource_roles=('vessel', 'carrier', 'operator'))


    # 5) Build Task
    task1 = Task(name='calwave_task1', actions=sc.actions, action_sequence='dependencies')

    # Check assets
    # task1.checkAssets(sc.vessels)

    # task1.updateStartTime(newStart=0)

    # 6) Build the Gantt chart
    task1.GanttChart(color_by='asset')
    plt.show()

    # Old chart building code:
    # 7) build the chart input directly from the Task and plot  #TODO: Rudy / Improve this later (maybe include it in Task.py/Scenario and let it plot the absolute time instead of relative time)
    chart_view = chart.view_from_task(task1, sc, title='CalWave Task 1 - Anchor installation plan')
    chart.plot_task(chart_view, outpath='calwave_task1_chart.png')



