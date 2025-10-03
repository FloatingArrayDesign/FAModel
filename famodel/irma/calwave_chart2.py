
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
# Core plotter (single-axes, multiple lanes)
# ===============================

def plot_task(task: Task, outpath: Optional[str] = None, dpi: int = 200,
              show_title: bool = True) -> None:
    """
    Render a Gantt-like chart for a single Task with one axes and one horizontal lane per vessel.
    • Vessel names as y-tick labels (structure like calwave_chart1)
    • Visual styling aligned with calwave_chart: baseline arrows, light span bars, circle bubbles
      with time inside, title above, capabilities below, and consistent font sizes.
    • Horizontal placement uses Bubble.period when available; otherwise cumulative within vessel.
    """
    from matplotlib.patches import FancyArrow

    # --- figure geometry ---
    nrows = max(1, len(task.vessels))
    fig_h = max(3.0, 1.2 + 1.6*nrows)
    fig_w = 16.0

    plt.rcdefaults()
    plt.close('all')
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)

    # --- y lanes (top -> bottom keeps given order) ---
    vessels_top_to_bottom = task.vessels
    y_positions = list(range(nrows))[::-1]
    name_to_y = {vt.vessel: y_positions[i] for i, vt in enumerate(vessels_top_to_bottom[::-1])}

    ax.set_yticks(y_positions)
    ax.set_yticklabels([vt.vessel for vt in vessels_top_to_bottom[::-1]])
    ax.tick_params(axis='y', labelrotation=0)

    if show_title:
        ax.set_title(task.name, loc='left', fontsize=12, pad=8)

    # --- gather periods, compute x-range ---
    x_min, x_max = 0.0, 0.0
    per_row: Dict[str, List[Tuple[float, float, Bubble]]] = {vt.vessel: [] for vt in task.vessels}

    for vt in task.vessels:
        t_cursor = 0.0
        for b in vt.bubbles:
            if b.period:
                s, e = float(b.period[0]), float(b.period[1])
            else:
                s = t_cursor
                e = s + float(b.duration_hr or 0.0)
            per_row[vt.vessel].append((s, e, b))
            x_min = min(x_min, s)
            x_max = max(x_max, e)
            t_cursor = e

    # --- drawing helpers ---
    def _draw_lane_baseline(y_val: float):
        # Baseline with arrow (like calwave_chart) spanning current x-lims later
        ax.annotate('', xy=(x_max, y_val), xytext=(x_min, y_val),
                    arrowprops=dict(arrowstyle='-|>', lw=2))

    def _draw_span_hint(s: float, e: float, y_val: float):
        ax.plot([s, e], [y_val, y_val], lw=6, alpha=0.15, color='k')

    def _draw_bubble(s: float, e: float, y_val: float, b: Bubble, i_in_row: int):
        xc = 0.5*(s + e)
        # Bubble marker (match calwave_chart size)
        ax.plot(xc, y_val, 'o', ms=45, color=plt.rcParams['axes.prop_cycle'].by_key()['color'][0], zorder=3)
        # Time/label inside bubble
        ax.text(xc, y_val, f'{b.label_time}', ha='center', va='center', fontsize=20,
                color='white', weight='bold')
        # Title above (alternate small offset pattern as in calwave_chart)
        title_offset = 0.28 if (i_in_row % 2) else 0.20
        ax.text(xc, y_val + title_offset, b.action, ha='center', va='bottom', fontsize=10)
        # Capabilities below
        caps_txt = _capabilities_to_text(b.capabilities)
        if caps_txt:
            ax.text(xc, y_val - 0.26, caps_txt, ha='center', va='top', fontsize=8, wrap=True)

    # --- draw per lane ---
    for vt in task.vessels:
        y = name_to_y[vt.vessel]
        items = sorted(per_row[vt.vessel], key=lambda t: t[0])

        # Lane baseline w/ arrow
        _draw_lane_baseline(y)

        # Span hints and bubbles
        for j, (s, e, b) in enumerate(items):
            _draw_span_hint(s, e, y)
            _draw_bubble(s, e, y, b, j)

    # --- axes cosmetics & limits ---
    if x_max <= x_min:
        x_max = x_min + 1.0
    pad = 0.02*(x_max - x_min) if (x_max - x_min) > 0 else 0.02
    ax.set_xlim(x_min - pad, x_max + pad)

    ax.set_xlabel('Hours (proportional)')
    ax.grid(False)
    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)

    # y limits and a little margin
    ax.set_ylim(min(y_positions) - 0.5, max(y_positions) + 0.5)

    fig = ax.figure
    fig.subplots_adjust(left=0.10, right=0.98, top=0.90, bottom=0.14)

    if outpath:
        fig.savefig(outpath, dpi=dpi, bbox_inches='tight')
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
    if not isinstance(roles, dict):
        return {'roles': []}
    clean: Dict[str, List[str]] = {}
    for r, caps in roles.items():
        clean[r] = list(caps or [])
    return clean


def _label_for(action) -> str:
    """Label inside the bubble. You can change to action.type or object name etc."""
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
    buckets: dict[str, list[Bubble]] = {}

    for a in sc.actions.values():
        t0 = start_times.get(a.name, None)
        dur = float(getattr(a, 'duration', 0.0) or 0.0)
        period = None
        if t0 is not None:
            period = (float(t0), float(t0) + dur)

        bubble = Bubble(
            action=a.name,
            duration_hr=dur,
            label_time=_label_for(a),
            capabilities=_caps_from_type_spec(a),
            period=period
        )
        vessel_label = _vessel_name_from_assigned(a)
        buckets.setdefault(vessel_label, []).append(bubble)

    vessels: List[VesselTimeline] = []
    for vname, bubbles in buckets.items():
        bubbles_sorted = sorted(
            bubbles,
            key=lambda b: (9999.0 if b.period is None else b.period[0])
        )
        vessels.append(VesselTimeline(vessel=vname, bubbles=bubbles_sorted))

    order_hint = {'San_Diego': 0, 'San Diego': 0, 'Jag': 1, 'Beyster': 2}
    vessels.sort(key=lambda vt: order_hint.get(vt.vessel, 10))

    return Task(name=title, vessels=vessels)

def stage_and_plot(sc, start_times: dict[str, float], title: str,
                   outpath: str | None = None, dpi: int = 200):
    t = scenario_to_chart_task(sc, start_times, title)
    plot_task(t, outpath=outpath, dpi=dpi, show_title=True)


if __name__ == '__main__':
    # Demo scene (unchanged structure)
    Support = VesselTimeline(
        vessel='Beyster',
        bubbles=[
            Bubble(action='mobilize', duration_hr=1.0, label_time='1.0', capabilities={'vessel': []}, period=(5.0, 6.0)),
            Bubble(action='transit_A2', duration_hr=0.5, label_time='0.5', capabilities={'carrier': []}, period=(6.5, 7.0)),
            Bubble(action='site_survey_A2', duration_hr=1.0, label_time='1.0', capabilities=[], period=(7.0, 8.0)),
            Bubble(action='transit_A1', duration_hr=0.5, label_time='0.5', capabilities=[], period=(8.0, 8.5)),
            Bubble(action='site_survey_A1', duration_hr=1.0, label_time='1.0', capabilities=[], period=(8.5, 9.5)),
            Bubble(action='transit_A4', duration_hr=0.5, label_time='0.5', capabilities=[], period=(9.5, 10.0)),
            Bubble(action='site_survey_A4', duration_hr=1.0, label_time='1.0', capabilities=[], period=(10.0, 11.0)),
            Bubble(action='transit_A3', duration_hr=0.5, label_time='0.5', capabilities=[], period=(11.0, 11.5)),
            Bubble(action='site_survey_A3', duration_hr=1.0, label_time='1.0', capabilities=[], period=(11.5, 12.5)),
            Bubble(action='transit', duration_hr=0.75, label_time='0.75', capabilities=[], period=(12.5, 13.25)),
            Bubble(action='mobilize', duration_hr=1.0, label_time='1.0', capabilities=[], period=(13.25, 14.25)),
        ]
    )

    Tug = VesselTimeline(
        vessel='Jar',
        bubbles=[
            Bubble(action='transit', duration_hr=0.5, label_time='4.5', capabilities={'carrier1': [], 'carrier2': [], 'operator': []}, period=(2.0, 6.5)),
            Bubble(action='at_site_support', duration_hr=5.5, label_time='5.5', capabilities=[], period=(6.5, 12.0)),
            Bubble(action='transit', duration_hr=4.5, label_time='4.5', capabilities={'carrier1': [], 'carrier2': [], 'operator': []}, period=(12.0, 16.5)),
        ]
    )

    Barge = VesselTimeline(
        vessel='San_Diego',
        bubbles=[
            Bubble(action='mobilize', duration_hr=2.0, label_time='2.0', capabilities={'carrier1': [], 'carrier2': ['deck_space', 'winch', 'positioning_system'], 'operator': ['crane']}, period=(0.0, 2.0)),
            Bubble(action='transit_tug', duration_hr=4.5, label_time='4.5', capabilities={'carrier1': [], 'carrier2': ['deck_space', 'winch', 'positioning_system'], 'operator': ['crane']}, period=(2.0, 6.5)),
            Bubble(action='install_anchor A2', duration_hr=1.0, label_time='1.0', capabilities={'carrier': ['deck_space'], 'operator': []}, period=(6.5, 7.5)),
            Bubble(action='transit_site A1', duration_hr=0.5, label_time='0.5', capabilities=[], period=(7.5, 8.0)),
            Bubble(action='install_anchor A1', duration_hr=1.0, label_time='1.5', capabilities={'carrier': ['deck_space'], 'operator': []}, period=(8.0, 9.0)),
            Bubble(action='transit_site A4', duration_hr=0.5, label_time='0.5', capabilities=[], period=(9.0, 9.5)),
            Bubble(action='install_anchor A4', duration_hr=1.0, label_time='1.0', capabilities={'carrier': ['deck_space'], 'operator': []}, period=(9.5, 10.5)),
            Bubble(action='transit_site A3', duration_hr=0.5, label_time='0.5', capabilities=[], period=(10.5, 11.0)),
            Bubble(action='install_anchor A3', duration_hr=1.0, label_time='1.0', capabilities={'carrier': ['deck_space'], 'operator': []}, period=(11.0, 12.0)),
            Bubble(action='transit_tug', duration_hr=4.5, label_time='4.5', capabilities={'carrier1': [], 'carrier2': ['deck_space', 'winch', 'positioning_system'], 'operator': ['crane']}, period=(12.0, 16.5)),
            Bubble(action='mobilize', duration_hr=1.0, label_time='1.0', capabilities={'carrier1': [], 'carrier2': ['deck_space', 'winch', 'positioning_system'], 'operator': ['crane']}, period=(16.5, 17.0)),
        ]
    )

    t = Task(name='Task 1 — Anchor installation plan (x4 anchors)', vessels=[Support, Tug, Barge])
    plot_task(t, outpath=None)
