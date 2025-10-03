# calwave_task1.py
# Build CalWave Task 1 (Anchor installation) following the theory flow:
# 1) addAction → structure only (type, name, objects, deps)
# 2) evaluateAssets → assign vessels/roles (+ durations/costs)
# 3) (schedule/plot handled by your existing tooling)

from famodel.project import Project
from calwave_irma import Scenario
from calwave_chart2 import Bubble, VesselTimeline, Task, plot_task

sc = Scenario()  # now sc exists in *this* session

# list vessel keys
print(list(sc.vessels.keys()))

# get a vessel name (dict-of-dicts access)
print(sc.vessels['San_Diego'])

# ------- Transit durations from plan [h] (one-way); tune as needed -------
SAN_DIEGO_TO_SITE = 4.5
JAG_TO_SITE       = 4.5
BEYSTER_TO_SITE   = 1.8

AT_SITE_SUPPORT_BLOCK = 7.5   # Jag on-station support window (fixed block)

SAN_DIEGO_TO_HOME = SAN_DIEGO_TO_SITE
JAG_TO_HOME       = JAG_TO_SITE
BEYSTER_TO_HOME   = BEYSTER_TO_SITE

# ---------- Helpers: map & normalize Project objects ----------
def map_project_objects(project):
    """
    Return A, M dicts with stable keys A1..A4 / M1..M4 mapped from project lists.
    Ensures each object has .type ('anchor'/'mooring') and .name for nice labels.
    """
    # Choose a stable order (sorted by key)
    a_keys = sorted(project.anchorList.keys())    # e.g., ['weca','wecb','wecc','wecd']
    m_keys = sorted(project.mooringList.keys())

    A = {f'A{i+1}': project.anchorList[k]  for i, k in enumerate(a_keys)}
    M = {f'M{i+1}': project.mooringList[k] for i, k in enumerate(m_keys)}

    # Normalize .type and .name (some libs don't set these)
    for k, obj in A.items():
        if getattr(obj, 'type', None) != 'anchor':
            setattr(obj, 'type', 'anchor')
        if not hasattr(obj, 'name'):
            setattr(obj, 'name', k)
    for k, obj in M.items():
        if getattr(obj, 'type', None) != 'mooring':
            setattr(obj, 'type', 'mooring')
        if not hasattr(obj, 'name'):
            setattr(obj, 'name', k)
    return A, M

def eval_set(a, roles, duration=None, **params):
    """
    Convenience: call evaluateAssets with roles/params and optionally set .duration.
    Always stores assigned_assets for plotting/scheduling attribution.
    """
    # Your Action.evaluateAssets may return (duration, cost); we still set explicit duration if passed.
    res = a.evaluateAssets(roles | params)
    if duration is not None:
        a.duration = float(duration)
    elif isinstance(res, tuple) and len(res) > 0 and res[0] is not None:
        try:
            a.duration = float(res[0])
        except Exception:
            pass
    # keep roles visible on the action
    a.assigned_assets = roles
    return a

# ---------- Core builder ----------
def build_task1_calwave(sc: Scenario, project: Project):
    """
    Creates Task 1 actions + dependencies (no scheduling/plotting here).
    Generic vessel actions are objectless; domain actions carry domain objects.
    """
    # Real domain instances
    A, M = map_project_objects(project)

    # --- Pre-ops (objectless) ---
    mob_sd  = sc.addAction('mobilize', 'mobilize_SanDiego')

    tr_sd   = sc.addAction('transit_tug', 'transit_site_SanDiego',
                           dependencies=[mob_sd])
    tr_jag  = sc.addAction('transit', 'transit_site_Jag',
                           dependencies=[mob_sd]) 
    
    mob_by  = sc.addAction('mobilize', 'mobilize_Beyster',
                           dependencies=[mob_sd])
    tr_by   = sc.addAction('transit', 'transit_site_Beyster',
                           dependencies=[mob_by])

    # Jag support window (objectless)
    sup_jag = sc.addAction('at_site_support', 'at_site_support_Jag',
                           dependencies=[tr_jag])

    # --- On-site (domain objects REQUIRED) ---
    inst, trans_tug, trans, mon = [], [], [], []
    for i in range(1, 5):
        ak, sk = f'A{i}', f'S{i}'
        a_inst = sc.addAction('install_anchor', f'install_anchor-{ak}',
                              objects=[A[ak]],
                              dependencies=[tr_sd])  # SD on site + Jag supporting
        a_trans_tug  = sc.addAction('transit_tug', f'transit_tug-{ak}',
                              # objects=[],
                              dependencies=[a_inst])
        a_trans  = sc.addAction('transit', f'transit-{ak}',
                              # objects=[],
                              dependencies=[a_inst])
        a_mon = sc.addAction('monitor_installation', f'monitor_installation-{sk}',
                             objects=[A[ak]],
                             dependencies=[a_trans_tug])
        inst.append(a_inst)
        trans_tug.append(a_trans_tug)
        trans.append(a_trans)
        mon.append(a_mon)
        
    mon_last = mon[-1]

    # --- Post-ops (objectless) ---
    home_sd  = sc.addAction('transit_tug', 'transit_homeport_SanDiego',
                            dependencies=mon)
    home_jag = sc.addAction('transit', 'transit_homeport_Jag',
                            dependencies=mon)
    home_by  = sc.addAction('transit', 'transit_homeport_Beyster',
                            dependencies=mon)
    
    # --- Pre-ops (objectless) ---
    demob_sd  = sc.addAction('demobilize', 'demobilize_SanDiego',
                             dependencies=[home_sd])
    demob_by  = sc.addAction('demobilize', 'demobilize_Beyster',
                             dependencies=[home_by])

    # Return a simple list for downstream evaluate/schedule/plot steps
    return {
        'mobilize': [mob_sd, mob_by],
        'transit_site': [tr_sd, tr_jag, tr_by],
        'support': [sup_jag],
        'install': inst,
        'transit_tug': trans_tug,
        'transit': trans,
        'monitor': mon,
        'transit_homeport': [home_sd, home_jag, home_by],
        'demobilize': [demob_sd, demob_by]
    }

# ---------- Evaluation step (assign vessels & durations) ----------
def evaluate_task1(sc: Scenario, actions: dict):
    """
    Assign vessels/roles and set durations where the evaluator doesn't.
    Keeps creation and evaluation clearly separated.
    """
    V = sc.vessels  # shorthand

    # Mobilize
    eval_set(actions['mobilize'][0], {'operator': V['San_Diego']}, duration=2.0)
    eval_set(actions['mobilize'][1], {'operator': V['Beyster']},   duration=1.0)

    # Transit to site
    tr_sd, tr_jag, tr_by = actions['transit_site']
    eval_set(tr_sd,  {'carrier': V['Jag'], 'operator': V['San_Diego']}, duration=SAN_DIEGO_TO_SITE)
    eval_set(tr_jag, {'carrier': V['Jag']},                             duration=JAG_TO_SITE)
    eval_set(tr_by,  {'carrier': V['Beyster']},                         duration=BEYSTER_TO_SITE)

    # Jag support block (fixed)
    eval_set(actions['support'][0], {'operator': V['Jag']}, duration=AT_SITE_SUPPORT_BLOCK)

    # Install / Lay (San Diego operates; durations can come from evaluator or set defaults)
    for a_inst in actions['install']:
        eval_set(a_inst, {'carrier': V['Jag'], 'operator': V['San_Diego']})
        if not getattr(a_inst, 'duration', 0):
            a_inst.duration = 1.2  # h, placeholder if evaluator didn’t set it
    for a_trans_tug in actions['transit_tug']:
        eval_set(a_trans_tug, {'carrier': V['Jag'], 'operator': V['San_Diego']})
        if not getattr(a_trans_tug, 'duration', 0):
            a_trans_tug.duration = 2.0  # h, placeholder
    for a_trans in actions['transit']:
        eval_set(a_trans, {'carrier': V['Beyster']})
        if not getattr(a_trans, 'duration', 0):
            a_trans.duration = 1.0  # h, placeholder

    # Monitor (Beyster)
    for a_mon in actions['monitor']:
        eval_set(a_mon, {'support': V['Beyster']})
        if not getattr(a_mon, 'duration', 0):
            a_mon.duration = 2.0  # h, placeholder

    # Transit home
    home_sd, home_jag, home_by = actions['transit_homeport']
    eval_set(home_sd,  {'carrier': V['Jag'], 'operator': V['San_Diego']}, duration=SAN_DIEGO_TO_HOME)
    eval_set(home_jag, {'carrier': V['Jag']},                             duration=JAG_TO_HOME)
    eval_set(home_by,  {'carrier': V['Beyster']},                         duration=BEYSTER_TO_HOME)
    
    # Demobilize
    eval_set(actions['demobilize'][0], {'operator': V['San_Diego']}, duration=2.0)
    eval_set(actions['demobilize'][1], {'operator': V['Beyster']},   duration=1.0)
    
def _action_resources(a) -> set[str]:
    """
    Return the set of resource keys (vessel names) this action occupies.
    We look in assigned_assets for any vessel-like roles.
    If nothing is assigned, we put the action into an 'unknown' pool (no blocking effect).
    """
    aa = getattr(a, 'assigned_assets', {}) or {}
    keys = []
    for role in ('vessel', 'carrier', 'operator'):
        v = aa.get(role)
        if v is not None:
            keys.append(getattr(v, 'name', str(v)))
    return set(keys) if keys else {'unknown'}

def schedule_actions(actions_by_name: dict[str, object]) -> dict[str, float]:
    """
    Compute earliest-start times for all actions given durations and dependencies,
    with single-resource constraints per vessel (i.e., a vessel can't overlap itself).
    Returns: {action_name: start_time_hours}
    """
    # Build dependency maps
    deps = {name: [d if isinstance(d, str) else getattr(d, 'name', str(d))
                    for d in getattr(a, 'dependencies', [])]
            for name, a in actions_by_name.items()}
    indeg = {name: len(dlist) for name, dlist in deps.items()}
    children = {name: [] for name in actions_by_name}
    for child, dlist in deps.items():
        for parent in dlist:
            if parent not in children:
                children[parent] = []
            children[parent].append(child)

    # Ready queue
    ready = sorted([n for n, k in indeg.items() if k == 0])

    start, finish = {}, {}
    avail = {}  # vessel_name -> time available

    scheduled = []

    while ready:
        name = ready.pop(0)
        a = actions_by_name[name]
        scheduled.append(name)

        # Dependency readiness
        dep_ready = 0.0
        for d in deps[name]:
            if d not in finish:
                raise RuntimeError(f"Dependency '{d}' of '{name}' has no finish time; check graph.")
            dep_ready = max(dep_ready, finish[d])

        # Resource readiness (all vessels the action occupies)
        res_keys = _action_resources(a)
        res_ready = max(avail.get(r, 0.0) for r in res_keys) if res_keys else 0.0

        # Start at the latest of dependency- and resource-readiness
        s = max(dep_ready, res_ready)
        d = float(getattr(a, 'duration', 0.0) or 0.0)
        f = s + d

        start[name] = s
        finish[name] = f

        # Block all involved resources until 'f'
        for r in res_keys:
            avail[r] = f

        # Release children
        for c in children.get(name, []):
            indeg[c] -= 1
            if indeg[c] == 0:
                ready.append(c)
        ready.sort()

    if len(scheduled) != len(actions_by_name):
        missing = [n for n in actions_by_name if n not in scheduled]
        raise RuntimeError(f"Cycle or missing predecessors detected; unscheduled: {missing}")

    return start

# --------------------------- Scenario indexes ------------------------------

def build_indexes(sc):
    '''Create quick lookups from Scenario content.
    Returns a dict with:
      - vessel_keys: set of vessel names as stored in sc.vessels
      - type_to_vessel: map of unique vessel type -> vessel key (only if unique)
      - action_to_vessel: map of action base-name -> vessel key declared in YAML
    '''
    vessel_keys = set((getattr(sc, 'vessels', {}) or {}).keys())

    # type -> vessel (only if unique across vessels)
    type_count = {}
    type_first = {}
    for vkey, vdesc in (getattr(sc, 'vessels', {}) or {}).items():
        vtype = vdesc.get('type')
        if isinstance(vtype, str) and vtype:
            type_count[vtype] = type_count.get(vtype, 0) + 1
            type_first.setdefault(vtype, vkey)
    type_to_vessel = {t: v for t, v in type_first.items() if type_count.get(t, 0) == 1}

    # actions listed under each vessel in YAML (if present)
    action_to_vessel = {}
    for vkey, vdesc in (getattr(sc, 'vessels', {}) or {}).items():
        v_actions = vdesc.get('actions') or {}
        for aname in v_actions.keys():
            action_to_vessel.setdefault(aname, vkey)

    return {
        'vessel_keys': vessel_keys,
        'type_to_vessel': type_to_vessel,
        'action_to_vessel': action_to_vessel,
    }

# --------------------------- Lane resolution --------------------------------

def lane_from_action_name(aname: str, vessel_keys: set[str]) -> str | None:
    # direct match of any vessel key as a token or suffix
    # examples: 'mobilize_San_Diego', 'transit_site_Beyster', 'foo-Jag'
    tokens = [aname]
    for sep in ('_', '-', ':'):
        tokens.extend(aname.split(sep))
    tokens = [t for t in tokens if t]
    for vk in vessel_keys:
        if vk in tokens or aname.endswith(vk) or aname.replace('_', '').endswith(vk.replace('_', '')):
            return vk
    return None

def lane_from_assigned_assets(action, vessel_keys: set[str], type_to_vessel: dict[str, str]) -> str | None:
    aa = getattr(action, 'assigned_assets', {}) or {}

    def pick(asset):
        # direct vessel key string
        if isinstance(asset, str) and asset in vessel_keys:
            return asset
        # object with .name equal to a vessel key
        nm = getattr(asset, 'name', None)
        if isinstance(nm, str) and nm in vessel_keys:
            return nm
        # dict with an explicit key equal to a vessel key
        if isinstance(asset, dict):
            for k in ('key', 'name', 'vessel', 'vessel_name', 'display_name', 'id', 'ID'):
                v = asset.get(k)
                if isinstance(v, str) and v in vessel_keys:
                    return v
            # dict with a type that uniquely maps to a vessel
            vtype = asset.get('type')
            if isinstance(vtype, str) and vtype in type_to_vessel:
                return type_to_vessel[vtype]
        # object with .type that uniquely maps
        vtype = getattr(asset, 'type', None)
        if isinstance(vtype, str) and vtype in type_to_vessel:
            return type_to_vessel[vtype]
        return None

    # try roles in a reasonable priority
    for role in ('vessel', 'operator', 'carrier', 'carrier1', 'carrier2', 'support'):
        lane = pick(aa.get(role))
        if lane:
            return lane
    return None

# --------------------------- Task builder -----------------------------------

def task_from_scenario(sc, start_times: dict, title: str, show_unknown: bool = False) -> Task:
    idx = build_indexes(sc)
    vessel_keys = idx['vessel_keys']
    type_to_vessel = idx['type_to_vessel']
    action_to_vessel = idx['action_to_vessel']

    buckets = {}
    warnings = []

    for a in (getattr(sc, 'actions', {}) or {}).values():
        aname = a.name
        t0 = float(start_times.get(aname, 0.0))
        dur = float(getattr(a, 'duration', 0.0) or 0.0)

        lane = (
            lane_from_action_name(aname, vessel_keys)
            or action_to_vessel.get(aname)
            or lane_from_assigned_assets(a, vessel_keys, type_to_vessel)
        )
        if not lane:
            if show_unknown:
                lane = 'unknown'
            else:
                warnings.append(f'No lane for action {aname}; dropping from plot')
                continue

        b = Bubble(
            action=aname,
            duration_hr=dur,
            label_time=f'{dur:.1f}',
            capabilities=[],
            period=(t0, t0 + dur)
        )
        buckets.setdefault(lane, []).append(b)

    # assemble rows and keep only vessels that exist in Scenario
    lanes = []
    for vname, blist in buckets.items():
        if vname != 'unknown' and vname not in vessel_keys:
            warnings.append(f'Lane {vname} not in Scenario.vessels; skipping')
            continue
        blist.sort(key=lambda b: b.period[0])
        lanes.append(VesselTimeline(vessel=vname, bubbles=blist))

    # order rows according to the order in Scenario (if you prefer a fixed order, hardcode it)
    order_map = {vk: i for i, vk in enumerate(sc.vessels.keys())}
    lanes.sort(key=lambda vt: order_map.get(vt.vessel, 999))

    return Task(name=title, vessels=lanes)

if __name__ == '__main__':
    # 1) Load ontology that mirrors the sample schema (mooring_systems + mooring_line_configs)
    project = Project(file='calwave_ontology.yaml', raft=False)
    project.getMoorPyArray(cables=1)
    
    # 2) Scenario with CalWave catalogs
    sc = Scenario()
    
    # 3) Build (structure only)
    actions = build_task1_calwave(sc, project)
    
    # 4) Evaluate (assign vessels/roles + durations)
    evaluate_task1(sc, actions)
       
    # 5) Schedule (replace with proper scheduler call)
    start_times = schedule_actions(sc.actions)  # <- dict {action_name: t_star}
    
    # 6) plot with the calwave_chart2 visual
    task = task_from_scenario(sc, start_times, title='CalWave Task 1', show_unknown=False)
    plot_task(task)
