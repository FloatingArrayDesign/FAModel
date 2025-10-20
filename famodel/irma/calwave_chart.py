
from dataclasses import dataclass
from typing import List, Optional, Dict, Tuple
import matplotlib.pyplot as plt

# ===============================
# Data structures
# ===============================

@dataclass
class Bubble:
    action: str
    duration_hr: float
    label_time: str
    period: Optional[Tuple[float, float]] = None
    category: Optional[str] = None  # new: action category for coloring

@dataclass
class VesselTimeline:
    vessel: str
    bubbles: List[Bubble]

@dataclass
class Task:
    name: str
    vessels: List[VesselTimeline]

# ===============================
# Color palette + categorization
# ===============================

# User-requested color scheme
ACTION_TYPE_COLORS: Dict[str, str] = {
    'Mobilization': '#d62728',                # red
    'Towing & Transport': '#2ca02c',          # green
    'Mooring & Anchors': '#0056d6',           # blue
    'Heavy Lift & Installation': '#ffdd00',   # yellow
    'Cable Operations': '#9467bd',            # purple
    'Survey & Monitoring': '#ff7f0e',         # orange
    'Other': '#1f77b4'}                       # fallback color (matplotlib default)


# Keyword buckets → chart categories
CAT_KEYS = [
    ('Mobilization', ('mobilize', 'demobilize')),
    ('Towing & Transport', ('transit', 'towing', 'tow', 'convoy', 'linehaul')),
    ('Mooring & Anchors', ('anchor', 'mooring', 'pretension', 'pre-tension')),
    ('Survey & Monitoring', ('monitor', 'survey', 'inspection', 'rov', 'divers')),
    ('Heavy Lift & Installation', ('install_wec', 'install device', 'install', 'heavy-lift', 'lift', 'lower', 'recover_wec', 'recover device')),
    ('Cable Operations', ('cable', 'umbilical', 'splice', 'connect', 'wet-mate', 'dry-mate'))]


def view_from_task(sched_task, sc, title: str | None = None):
    """
    Minimal map: scheduler Task -> chart view Task
    Show an action on multiple lanes if it uses multiple assets.

    Rules per role value:
      • If str and in sc.vessels → use as key.
      • Else if object → resolve by identity to sc.vessels.
      • Else if dict → try ['name'] as key; else if ['type'] is unique in sc.vessels, use that key.
      • Add the bubble to every resolved lane (deduped).
      • Skip actions with dur<=0 or with no resolvable lanes.
    """
    # reverse lookup for identity → key
    id2key = {id(obj): key for key, obj in sc.vessels.items()}

    # unique type → key (used only if type is unique in catalog)
    type_counts = {}
    for k, obj in sc.vessels.items():
        t = obj.get('type') if isinstance(obj, dict) else getattr(obj, 'type', None)
        if t:
            type_counts[t] = type_counts.get(t, 0) + 1
    unique_type2key = {}
    for k, obj in sc.vessels.items():
        t = obj.get('type') if isinstance(obj, dict) else getattr(obj, 'type', None)
        if t and type_counts.get(t) == 1:
            unique_type2key[t] = k

    buckets = {}

    for a in sched_task.actions.values():
        dur = float(getattr(a, 'duration', 0.0) or 0.0)
        if dur <= 0.0:
            continue

        aa = getattr(a, 'assigned_assets', {}) or {}

        # collect ALL candidate roles → multiple lanes allowed
        lane_keys = set()
        for role in ('vessel', 'carrier', 'operator', 'support'):
            if role not in aa:
                continue
            v = aa[role]

            # resolve lane key
            lane = None
            if isinstance(v, str):
                lane = v if v in sc.vessels else None
            elif v is not None:
                lane = id2key.get(id(v))
                if lane is None and isinstance(v, dict):
                    nm = v.get('name')
                    if isinstance(nm, str) and nm in sc.vessels:
                        lane = nm
                    else:
                        t = v.get('type')
                        if t in unique_type2key:
                            lane = unique_type2key[t]
            if lane:
                lane_keys.add(lane)

        if not lane_keys:
            continue

        t0 = float(getattr(a, 'start_hr', 0.0) or 0.0)
        t1 = float(getattr(a, 'end_hr', t0) or 0.0)
        
        # Color code for action categories based on CAT_KEYS
        def cat_for(act):
            s = f"{getattr(act, 'type', '')} {getattr(act, 'name', '')}".lower().replace('_', ' ')
            for cat, keys in CAT_KEYS:
                if any(k in s for k in keys):
                    return cat
            return 'Other'
        
        # one bubble per lane (same fields)
        for lane in lane_keys:
            b = Bubble(
                action=a.name,
                duration_hr=dur,
                label_time=getattr(a, 'label_time', f'{dur:.1f}'),
                period=(t0, t1),
                category=cat_for(a))
            
            buckets.setdefault(lane, []).append(b)

    # preserve sc.vessels order; only include lanes with content
    lanes = []
    for vname in sc.vessels.keys():
        blist = sorted(buckets.get(vname, []), key=lambda b: b.period[0])
        if blist:
            lanes.append(VesselTimeline(vessel=vname, bubbles=blist))

    return Task(name=title or getattr(sched_task, 'name', 'Task'), vessels=lanes)

# ===============================
# Core plotter (single-axes, multiple lanes)
# ===============================

def plot_task(task: Task, outpath: Optional[str] = None, dpi: int = 200,
              show_title: bool = True) -> None:
    """
    Render a Gantt-like chart for a single Task with one axes and one horizontal lane per vessel.
    • Vessel names as y-tick labels
    • Baseline arrows, light span bars, circle bubbles with time inside, title above,
      and consistent font sizes.
    • Horizontal placement uses Bubble.period when available; otherwise cumulative within vessel.
    • Bubbles are colored by Bubble.category (legend added).
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Circle

    # --- figure geometry ---
    nrows = max(1, len(task.vessels))
    fig_h = max(3.0, 1.2 + 1.6*nrows)
    fig_w = 16.0

    plt.rcdefaults()
    plt.close('all')
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)

    # --- y lanes (top -> bottom keeps given order) ---
    vessels_top_to_bottom = task.vessels
    nrows = max(1, len(task.vessels))
    y_positions = list(range(nrows))[::-1]
    name_to_y = {vt.vessel: y_positions[i] for i, vt in enumerate(vessels_top_to_bottom[::-1])}
    
    ax.set_yticks(y_positions)
    ax.set_yticklabels([])          
    ax.tick_params(axis='y', labelrotation=0)   
    
    if show_title:
        ax.set_title(task.name, loc='left', fontsize=16, pad=12)

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
        ax.annotate('', xy=(x_max, y_val), xytext=(x_min, y_val),
                    arrowprops=dict(arrowstyle='-|>', lw=2))

    def _draw_span_hint(s: float, e: float, y_val: float):
        ax.plot([s, e], [y_val, y_val], lw=6, alpha=0.15, color='k')

    def _bubble_face_color(b: Bubble) -> str:
        cat = b.category or 'Other'
        return ACTION_TYPE_COLORS.get(cat, ACTION_TYPE_COLORS['Other'])

    def _text_color_for_face(face: str) -> str:
        return 'black' if face.lower() in ('#ffdd00', 'yellow') else 'white'

    def _draw_bubble(s: float, e: float, y_val: float, b: Bubble, i_in_row: int):
        xc = 0.5*(s + e)
        face = _bubble_face_color(b)
        txtc = _text_color_for_face(face)
        ax.plot(xc, y_val, 'o', ms=45, color=face, zorder=3)
        ax.text(xc, y_val, f'{b.label_time}', ha='center', va='center', fontsize=20,
                color=txtc, weight='bold')
        title_offset = 0.30 if (i_in_row % 2) else 0.20
        ax.text(xc, y_val + title_offset, b.action, ha='center', va='bottom', fontsize=10)
        # caps_txt = _capabilities_to_text(b.capabilities)
        # if caps_txt:
        #    ax.text(xc, y_val - 0.26, caps_txt, ha='center', va='top', fontsize=8, wrap=True)

    # --- draw per lane ---
    seen_cats: set[str] = set()
    for vt in task.vessels:
        y = name_to_y[vt.vessel]
        items = sorted(per_row[vt.vessel], key=lambda t: t[0])
        _draw_lane_baseline(y)
        for j, (s, e, b) in enumerate(items):
            _draw_span_hint(s, e, y)
            _draw_bubble(s, e, y, b, j)
            seen_cats.add(b.category or 'Other')

    # --- legend ---
    handles = []
    legend_cats = [c for c in ACTION_TYPE_COLORS.keys() if c in seen_cats]
    # if you prefer to always show all categories, replace the line above with: legend_cats = list(ACTION_TYPE_COLORS.keys())
    for cat in legend_cats:
        handles.append(Line2D([0], [0], marker='o', linestyle='none', markersize=12,
                              markerfacecolor=ACTION_TYPE_COLORS[cat], markeredgecolor='none', label=cat))
    if handles:
        # Place the legend below the x-axis label (bottom center)
        fig_ = ax.figure
        fig_.legend(handles=handles,
                    loc='lower center',
                    bbox_to_anchor=(0.5, -0.12), # move below the axis label
                    ncol=3,
                    title='Action Types',
                    frameon=False)

    # --- axes cosmetics & limits ---
    if x_max <= x_min:
        x_max = x_min + 1.0
    pad = 0.02*(x_max - x_min) if (x_max - x_min) > 0 else 0.5
    ax.set_xlim(x_min - pad, x_max + pad)
    
    # Draw circled vessel names at the same y positions
    x_name = x_min - 3*pad       # small left offset inside the axes
    
    # After you have vessels_top_to_bottom, name_to_y, x_min/x_max, pad, left_extra, x_name...
    max_len = max(len(vt.vessel) for vt in vessels_top_to_bottom)  # longest label
    
    # make the circle tighter/looser:
    circle_pad = 0.18   
    
    for vt in vessels_top_to_bottom[::-1]:
        y = name_to_y[vt.vessel]
        fixed_text = vt.vessel.center(max_len)  # pad with spaces to max length
        ax.text(
            x_name, y, fixed_text,
            ha='center', va='center', zorder=6, clip_on=False,
            fontsize=12, color='black', fontfamily='monospace',  # <- key: monospace
            bbox=dict(boxstyle='circle,pad={:.2f}'.format(circle_pad),
                      facecolor='lightgrey', edgecolor='tomato', linewidth=6))

        ax.set_xlabel('Timeline (h)')
        ax.grid(False)
        for spine in ['top', 'right', 'left']:
            ax.spines[spine].set_visible(False)
    
        ax.set_ylim(min(y_positions) - 0.5, max(y_positions) + 0.5)

    fig = ax.figure
    # Add extra bottom margin to make space for the legend below the x-axis label
    fig.subplots_adjust(left=0.10, right=0.98, top=0.90, bottom=0.15)

    if outpath:
        fig.savefig(outpath, dpi=dpi, bbox_inches='tight')
    else:
        plt.show()
