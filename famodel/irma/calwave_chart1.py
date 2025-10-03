
import math
from dataclasses import dataclass
from typing import List, Optional, Dict, Any, Union, Tuple
import matplotlib.pyplot as plt

# ===============================
# Data structures
# ===============================
@dataclass
class Bubble:
    action: str                      # action name (matches YAML if desired)
    duration_hr: float               # used to space bubbles proportionally along the row
    label_time: Union[str, float]    # text inside the bubble (e.g., 0.0, 0.1, 'A1')
    capabilities: Union[List[str], Dict[str, List[str]]]  # shown below bubble; list or dict-by-role
    period: Optional[Tuple[float, float]] = None  # (start_time, end_time)

@dataclass
class VesselTimeline:
    vessel: str
    bubbles: List[Bubble]

@dataclass
class Task:
    name: str
    vessels: List[VesselTimeline]

# ===============================
# Helper: format capabilities nicely (supports roles or flat list)
# ===============================

def _capabilities_to_text(capabilities: Union[List[str], Dict[str, List[str]]]) -> str:
    if isinstance(capabilities, dict):
        parts = []
        for role, caps in capabilities.items():
            if not caps:
                parts.append(f'{role}: (none)')
            else:
                parts.append(f"{role}: " + ', '.join(caps))
        return '\n'.join(parts)
    if isinstance(capabilities, list):
        return ', '.join(capabilities) if capabilities else '(none)'
    return str(capabilities)

# ===============================
# Core plotters
# ===============================

def _accumulate_starts(durations: List[float]) -> List[float]:
    starts = [0.0]
    for d in durations[:-1]:
        starts.append(starts[-1] + d)
    return starts

def plot_task(task: Task, outpath: Optional[str] = None, dpi: int = 200,
              show_title: bool = True) -> None:
    """
    Render a Gantt-like chart for a single Task with one timeline per vessel.
    • One axes with a horizontal lane per vessel (vessel names as y-tick labels)
    • Bubble per action: title above, time inside, capabilities below
    • Horizontal placement uses b.period when available, otherwise cumulative within vessel
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    # --- figure geometry (stable/sane) ---
    nrows = max(1, len(task.vessels))
    fig_h = max(3.0, 1.2 + 1.6 * nrows)
    fig_w = 16.0

    plt.rcdefaults()       # avoid stray global styles from other modules
    plt.close('all')
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)

    # --- build y lanes (top -> bottom) ---
    vessels_top_to_bottom = task.vessels  # keep given order
    y_positions = list(range(nrows))[::-1]
    name_to_y = {vt.vessel: y_positions[i] for i, vt in enumerate(vessels_top_to_bottom[::-1])}

    ax.set_yticks(y_positions)
    ax.set_yticklabels([vt.vessel for vt in vessels_top_to_bottom[::-1]])
    ax.tick_params(axis='y', labelrotation=0)

    # baselines that span the axes width (independent of x limits)
    for yi in y_positions:
        ax.plot([0, 1], [yi, yi], transform=ax.get_xaxis_transform(), color='k', lw=1)

    if show_title:
        ax.set_title(task.name, loc='left', fontsize=12, pad=8)

    # --- gather periods, compute x-range ---
    x_min, x_max = 0.0, 0.0
    per_row = {vt.vessel: [] for vt in task.vessels}

    for vt in task.vessels:
        t_cursor = 0.0
        for b in vt.bubbles:
            if b.period:
                s, e = float(b.period[0]), float(b.period[1])
            else:
                # fall back to cumulative placement within the row
                s = t_cursor
                e = s + float(b.duration_hr or 0.0)
            per_row[vt.vessel].append((s, e, b))
            x_min = min(x_min, s)
            x_max = max(x_max, e)
            t_cursor = e

    # --- drawing helpers (always in DATA coords) ---
    def _cap_text(caps: dict) -> str:
        return _capabilities_to_text(caps) if isinstance(caps, dict) else ""

    def _draw_bubble(x0, x1, y, b):
        xc = 0.5 * (x0 + x1)
        w  = max(0.001, (x1 - x0))
        radius = 0.22 + 0.02 * w  # visual tweak similar to the newer plotter

        circ = Circle((xc, y), radius=radius, facecolor='#1976d2', edgecolor='k', lw=1.0, zorder=3)
        ax.add_patch(circ)

        # label INSIDE bubble (duration/time)
        ax.text(xc, y, f"{b.label_time}", ha='center', va='center',
                fontsize=9, color='white', weight='bold', transform=ax.transData, clip_on=False)

        # title ABOVE bubble
        ax.text(xc, y + (radius + 0.18), b.action, ha='center', va='bottom',
                fontsize=9, transform=ax.transData, clip_on=False)

        # capabilities BELOW bubble
        caps_txt = _cap_text(b.capabilities)
        if caps_txt:
            ax.text(xc, y - (radius + 0.24), caps_txt, ha='center', va='top',
                    fontsize=7, transform=ax.transData, clip_on=False)

    # --- draw bubbles, stagger titles a bit when very close ---
    for v in task.vessels:
        items = sorted(per_row[v.vessel], key=lambda t: t[0])
        y = name_to_y[v.vessel]
        last_x = -1e9
        alt = 0
        for s, e, b in items:
            _draw_bubble(s, e, y, b)

            # if adjacent centers are too close, nudge the title slightly up/down
            xc = 0.5 * (s + e)
            if abs(xc - last_x) < 0.6:   # hours; tweak if your labels still collide
                alt ^= 1
                bump = 0.28 if alt else -0.28
                ax.text(xc, y + bump, b.action, ha='center', va='bottom',
                        fontsize=9, transform=ax.transData, clip_on=False)
            last_x = xc

    # --- axes cosmetics & limits ---
    if x_max <= x_min:
        x_max = x_min + 1.0
    pad = 0.3 * (x_max - x_min)
    ax.set_xlim(x_min - pad, x_max + pad)
    ax.set_ylim(y_positions[-1] - 1, y_positions[0] + 1)

    ax.set_xlabel('Hours (proportional)')
    ax.grid(False)
    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)

    fig = ax.figure
    fig.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.14)

    if outpath:
        fig.savefig(outpath, dpi=dpi)
    else:
        plt.show()
    
# ========= Adapters from Scenario → chart =========

def _vessel_name_from_assigned(action) -> str:
    """Pick a vessel label from the roles you assigned in evaluateAssets."""
    aa = getattr(action, 'assigned_assets', {}) or {}
    # order of preference (tweak if your roles differ)
    for key in ('vessel', 'carrier', 'operator'):
        v = aa.get(key)
        if v is not None:
            return getattr(v, 'name', str(v))
    # fallback: unknown bucket
    return 'unknown'

def _caps_from_type_spec(action) -> dict:
    """Turn the action type's roles spec into a {role: [caps]} dict for display."""
    spec = getattr(action, 'type_spec', {}) or {}
    roles = spec.get('roles') or {}
    # be defensive: ensure it's a dict[str, list[str]]
    if not isinstance(roles, dict):
        return {'roles': []}
    clean = {}
    for r, caps in roles.items():
        clean[r] = list(caps or [])
    return clean

def _label_for(action) -> str:
    """Label inside the bubble. You can change to action.type or object name etc."""
    # show duration with one decimal if present; else use the action's short name
    dur = getattr(action, 'duration', None)
    if isinstance(dur, (int, float)) and dur >= 0:
        return f"{dur:.1f}"
    return action.name

def scenario_to_chart_task(sc, start_times: dict[str, float], title: str):
    """
    Convert Scenario actions + a start-time map into a calwave_chart.Task for plotting.
      - start_times: {action_name: t0} from your scheduler
      - action.duration must be set (via evaluateAssets or by you)
    """
    # 1) bucket actions by vessel label
    buckets: dict[str, list[Bubble]] = {}

    for a in sc.actions.values():
        # period
        t0 = start_times.get(a.name, None)
        dur = float(getattr(a, 'duration', 0.0) or 0.0)
        period = None
        if t0 is not None:
            period = (float(t0), float(t0) + dur)

        bubble = Bubble(
            action=a.name,
            duration_hr=dur,
            label_time=_label_for(a),
            capabilities=_caps_from_type_spec(a),  # roles/caps from the type spec
            period=period
        )
        vessel_label = _vessel_name_from_assigned(a)
        buckets.setdefault(vessel_label, []).append(bubble)

    # 2) sort bubbles per vessel by start time (or keep input order if no schedule)
    vessels = []
    for vname, bubbles in buckets.items():
        bubbles_sorted = sorted(
            bubbles,
            key=lambda b: (9999.0 if b.period is None else b.period[0])
        )
        vessels.append(VesselTimeline(vessel=vname, bubbles=bubbles_sorted))

    # 3) stable vessel ordering: San Diego, Jag, Beyster, then others
    order_hint = {'San_Diego': 0, 'San Diego': 0, 'Jag': 1, 'Beyster': 2}
    vessels.sort(key=lambda vt: order_hint.get(vt.vessel, 10))

    return Task(name=title, vessels=vessels)
       
# Optional convenience: do everything after scheduling
def stage_and_plot(sc, start_times: dict[str, float], title: str, outpath: str | None = None, dpi: int = 200):
    t = scenario_to_chart_task(sc, start_times, title)
    plot_task(t, outpath=outpath, dpi=dpi, show_title=True)


if __name__ == '__main__':
    # Support vessel monitors 4 anchors
    Support = VesselTimeline(
    vessel='Beyster',
    bubbles=[
    Bubble(action='mobilize', duration_hr=1.0, label_time='1.0', 
           capabilities={'vessel': []},
           period=(5.0, 6.0)),
    Bubble(action='transit_site A2', duration_hr=0.5, label_time='0.5',
           capabilities={'carrier': []},
           period=(6.5, 7.0)),
    # Monitor each anchor install (x4)
    Bubble(action='site_survey A2', duration_hr=1.0, label_time='1.0',
           capabilities=[], 
           period=(7.0, 8.0)),
    Bubble(action='transit_site A1', duration_hr=0.5, label_time='0.5',
           capabilities=[], 
           period=(8.0, 8.5)),
    Bubble(action='site_survey A1', duration_hr=1.0, label_time='1.0',
           capabilities=[], 
           period=(8.5, 9.5)),
    Bubble(action='transit_site A4', duration_hr=0.5, label_time='0.5',
           capabilities=[], 
           period=(9.5, 10.0)),
    Bubble(action='site_survey A4', duration_hr=1.0, label_time='1.0',
           capabilities=[], 
           period=(10.0, 11.0)),
    Bubble(action='transit_site A3', duration_hr=0.5, label_time='0.5',
           capabilities=[], 
           period=(11.0, 11.5)),
    Bubble(action='site_survey A3', duration_hr=1.0, label_time='1.0',
           capabilities=[], 
           period=(11.5, 12.5)),
    Bubble(action='transit_homeport', duration_hr=0.75, label_time='0.75', 
           capabilities=[],
           period=(12.5, 13.25)),
    Bubble(action='mobilize', duration_hr=1.0, label_time='1.0', 
           capabilities=[],
           period=(13.25, 14.25)),
    ]
    )
      
    # Tug (Jar) stages/load moorings; lay/install handled by San_Diego per your note
    Tug = VesselTimeline(
    vessel='Jar',
    bubbles=[
    Bubble(action='transit_site', duration_hr=0.5, label_time='4.5', 
           capabilities={'carrier1': [], 'carrier2': [], 'operator': []}, 
           period=(2.0, 6.5)),
    Bubble(action='at_site_support', duration_hr=5.5, label_time='5.5', 
           capabilities=[], period=(6.5, 12.0)),
    Bubble(action='transit_homeport', duration_hr=4.5, label_time='4.5', 
           capabilities={'carrier1': [], 'carrier2': [], 'operator': []}, 
           period=(12.0, 16.5)),
    ]
    )
    
    # Barge performs lay_mooring and install_anchor for 4 anchors
    Barge = VesselTimeline(
    vessel='San_Diego',
    bubbles=[
    Bubble(action='mobilize', duration_hr=2.0, label_time='2.0', 
           capabilities={'carrier1': [], 'carrier2': ['deck_space', 'winch', 'positioning_system'], 'operator': ['crane']},
           period =(0.0, 2.0)),
    # Anchor 2
    Bubble(action='transit_homeport', duration_hr=4.5, label_time='4.5', 
           capabilities={'carrier1': [], 'carrier2': ['deck_space', 'winch', 'positioning_system'], 'operator': ['crane']},
           period =(2.0, 6.5)),
    Bubble(action='install_anchor A2', duration_hr=1.0, label_time='1.0',
           capabilities={'carrier': ['deck_space'], 'operator': []},
           period=(6.5, 7.5)),
    # Anchor 1
    Bubble(action='transit_site A1', duration_hr=0.5, label_time='0.5',
           capabilities=[], 
           period=(7.5, 8.0)),
    Bubble(action='install_anchor A1', duration_hr=1.0, label_time='1.5',
           capabilities={'carrier': ['deck_space'], 'operator': []},
           period=(8.0, 9.0)),
    # Anchor 4
    Bubble(action='transit_site A4', duration_hr=0.5, label_time='0.5',
           capabilities=[], 
           period=(9.0, 9.5)),
    Bubble(action='install_anchor A4', duration_hr=1.0, label_time='1.0',
           capabilities={'carrier': ['deck_space'], 'operator': []},
           period=(9.5, 10.5)),
    # Anchor 3
    Bubble(action='transit_site A3', duration_hr=0.5, label_time='0.5',
           capabilities=[], 
           period=(10.5, 11.0)),
    Bubble(action='install_anchor A3', duration_hr=1.0, label_time='1.0',
           capabilities={'carrier': ['deck_space'], 'operator': []},
           period=(11.0, 12.0)),
    Bubble(action='transit_homeport', duration_hr=4.5, label_time='4.5', 
           capabilities={'carrier1': [], 'carrier2': ['deck_space', 'winch', 'positioning_system'], 'operator': ['crane']},
           period =(12.0, 16.5)),
    Bubble(action='mobilize', duration_hr=1.0, label_time='1.0', 
           capabilities={'carrier1': [], 'carrier2': ['deck_space', 'winch', 'positioning_system'], 'operator': ['crane']},
           period =(16.5, 17.0)),
    ]
    )
     
    t = Task(name='Task 1 — Anchor installation plan (x4 anchors)', vessels=[Support, Tug, Barge])
    plot_task(t, outpath=None)