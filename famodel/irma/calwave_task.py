"""Enhanced Task class for CalWave scheduling

- Adds earliest-start (critical-path) scheduler with single-resource constraints
- Keeps legacy level-based checkpoint scheduler (via getSequenceGraph)

Style: single quotes, spaces around + and -, no spaces around * or /
"""

from collections import defaultdict

class Task:
    def __init__(self, actions, action_sequence, **kwargs):
        '''
        Create a Task from a list of actions and a dependency map.

        Parameters
        ----------
        actions : list
            All Action objects that are part of this task.
        action_sequence : dict or None
            {action_name: [predecessor_name, ...]}.
            If None, dependencies are inferred from Action.dependencies.
        kwargs :
            Optional tuning:
              • strategy='earliest' | 'levels'
              • enforce_resources=True | False
              • resource_roles=('vessel','carrier','operator')
        '''
        # ---- options with sensible defaults (all via kwargs) ----
        strategy = kwargs.get('strategy', 'earliest')
        enforce_resources = kwargs.get('enforce_resources', True)
        resource_roles = kwargs.get('resource_roles', ('vessel', 'carrier', 'operator'))

        # ---- core storage ----
        self.actions = {a.name: a for a in actions}
        # allow None → infer solely from Action.dependencies
        self.action_sequence = {k: list(v) for k, v in (action_sequence or {}).items()}
        self.actions_ti = {}
        self.duration = 0.0
        self.cost = 0.0
        self.ti = 0.0
        self.tf = 0.0
        self.resource_roles = tuple(resource_roles)

        # ---- scheduling ----
        if strategy == 'levels':
            self._schedule_by_levels()
        else:
            self._schedule_by_earliest(enforce_resources=enforce_resources)

        # ---- roll-ups ----
        self.cost = sum(float(getattr(a, 'cost', 0.0) or 0.0) for a in self.actions.values())

    # -------- Convenience constructors / helpers (build deps inside the class) --------

    @staticmethod
    def _names_from_dependencies(a):
        deps = []
        for d in list(getattr(a, 'dependencies', []) or []):
            deps.append(d if isinstance(d, str) else getattr(d, 'name', str(d)))
        seen = set()
        clean = []
        for dn in deps:
            if dn != a.name and dn not in seen:
                clean.append(dn); seen.add(dn)
        return clean

    @classmethod
    def from_scenario(cls, sc, **kwargs):
        actions = list(sc.actions.values())
        base = {a.name: cls._names_from_dependencies(a) for a in actions}
        extra = kwargs.pop('extra_dependencies', None) or {}
        for k, v in extra.items():
            base.setdefault(k, [])
            for d in v:
                if d != k and d not in base[k]:
                    base[k].append(d)
        return cls(actions=actions, action_sequence=base, **kwargs)

    # --------------------------- Resource & Scheduling ---------------------------

    def _action_resources(self, a):
        '''Return set of resource keys (e.g., vessel names) this action occupies.
        Looks into a.assigned_assets for roles in self.resource_roles.
        If none assigned, returns {'unknown'} to avoid blocking anything real.
        '''
        aa = getattr(a, 'assigned_assets', {}) or {}
        keys = []
        for role in self.resource_roles:
            v = aa.get(role)
            if v is not None:
                keys.append(getattr(v, 'name', str(v)))
        return set(keys) if keys else {'unknown'}

    def _schedule_by_earliest(self, enforce_resources=True):
        '''Earliest-start (critical-path) schedule with optional single-resource blocking.'''
        # Merge dependencies from action attributes and explicit action_sequence
        deps = {}
        for name, a in self.actions.items():
            dlist = []
            # from Action.dependencies (may be objects or names)
            for d in list(getattr(a, 'dependencies', []) or []):
                dlist.append(d if isinstance(d, str) else getattr(d, 'name', str(d)))
            # from explicit dict
            dlist.extend(self.action_sequence.get(name, []))
            # hygiene
            seen = set()
            clean = []
            for d in dlist:
                if d != name and d not in seen:
                    clean.append(d)
                    seen.add(d)
            deps[name] = clean

        # ensure all nodes present
        for name in self.actions.keys():
            deps.setdefault(name, [])

        # Build children and indegrees
        children = {n: [] for n in self.actions}
        indeg = {n: len(dl) for n, dl in deps.items()}
        for child, dlist in deps.items():
            for parent in dlist:
                children.setdefault(parent, []).append(child)

        # Ready queue (roots)
        ready = sorted([n for n, k in indeg.items() if k == 0])

        start, finish = {}, {}
        avail = {}  # resource -> available time
        scheduled = []

        while ready:
            name = ready.pop(0)
            a = self.actions[name]
            scheduled.append(name)

            dep_ready = 0.0
            for d in deps[name]:
                if d not in finish:
                    raise RuntimeError(f"Dependency '{d}' of '{name}' missing finish time.")
                dep_ready = max(dep_ready, finish[d])

            if enforce_resources:
                res_keys = self._action_resources(a)
                res_ready = max(avail.get(r, 0.0) for r in res_keys) if res_keys else 0.0
            else:
                res_ready = 0.0

            s = max(dep_ready, res_ready)
            dur = float(getattr(a, 'duration', 0.0) or 0.0)
            f = s + dur

            start[name] = s
            finish[name] = f

            if enforce_resources:
                for r in self._action_resources(a):
                    avail[r] = f

            # release children
            for c in children.get(name, []):
                indeg[c] -= 1
                if indeg[c] == 0:
                    ready.append(c)
            ready.sort()

        if len(scheduled) != len(self.actions):
            missing = [n for n in self.actions if n not in scheduled]
            raise RuntimeError(f'Cycle or missing predecessors; unscheduled: {missing}')

        # Stamp fields on actions and compute task duration
        self.actions_ti = start
        for a in self.actions.values():
            a.start_hr = start[a.name]
            dur = float(getattr(a, 'duration', 0.0) or 0.0)
            a.end_hr = a.start_hr + dur
            a.period = (a.start_hr, a.end_hr)
            a.label_time = f'{dur:.1f}'
        self.duration = max((finish[n] for n in finish), default=0.0)
        self.tf = self.ti + self.duration

    def _schedule_by_levels(self):
        '''Wrapper that reuses the legacy level-based sequence graph to set starts.'''
        # Build levels using provided sequence dict
        levels = {}

        def level_of(a, path):
            if a in levels:
                return levels[a]
            if a in path:
                raise ValueError(f"Cycle detected at '{a}'.")
            path.add(a)
            pres = self.action_sequence.get(a, [])
            if not pres:
                lv = 1
            else:
                lv = 1 + max(level_of(p, path) if p in self.action_sequence else 1 for p in pres)
            levels[a] = lv
            return lv

        for name in self.action_sequence:
            level_of(name, set())

        max_level = max(levels.values(), default=1)
        groups = defaultdict(list)
        for n, lv in levels.items():
            groups[lv].append(n)
        level_dur = {lv: max(self.actions[a].duration for a in acts) for lv, acts in groups.items()}

        t = 0.0
        starts = {}
        for lv in range(1, max_level + 1):
            for n in groups.get(lv, []):
                starts[n] = t
            t += level_dur.get(lv, 0.0)

        self.actions_ti = starts
        self.duration = sum(level_dur.values())
        for a in self.actions.values():
            a.start_hr = starts.get(a.name, 0.0)
            dur = float(getattr(a, 'duration', 0.0) or 0.0)
            a.end_hr = a.start_hr + dur
            a.period = (a.start_hr, a.end_hr)
            a.label_time = f'{dur:.1f}'
        self.tf = self.ti + self.duration
