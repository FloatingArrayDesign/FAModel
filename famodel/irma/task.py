"""Action base class"""

import numpy as np
import matplotlib.pyplot as plt

import moorpy as mp
from moorpy.helpers import set_axes_equal
from moorpy import helpers
import yaml
from copy import deepcopy
import networkx as nx

# Import select required helper functions
from famodel.helpers import (check_headings, head_adjust, getCableDD, getDynamicCables, 
                            getMoorings, getAnchors, getFromDict, cleanDataTypes, 
                            getStaticCables, getCableDesign, m2nm, loadYAML, 
                            configureAdjuster, route_around_anchors)


class Task():
    '''
    A Task is a general representation of a set of marine operations
    that follow a predefined sequence/strategy. There can be multiple
    tasks that achieve the same end, each providing an alternative strategy.
    Each Task consists of a set of Actions with internal dependencies.
    
    
    For now, we'll assume each Task must be port-to-port,
    i.e. its vessel(s) must start and end at port over the course of the task.
    
    
    
    '''
    
    def __init__(self, name, actions, action_sequence=None, **kwargs):
        '''Create an action object...
        It must be given a name and a list of actions.
        The action list should be by default coherent with actionTypes dictionary.
        
        Parameters
        ----------
        name : string
            A name for the action. It may be appended with numbers if there
            are duplicate names.
        actions : list
            A list of all actions that are part of this task.
        action_sequence : dict, optional
            A dictionary where each key is the name of each action, and the values are
            each a list of which actions (by name) must be completed before the current
            one. If None, the action_sequence will be built by calling self.stageActions(from_deps=True)
            [building from the dependencies of each action].
        kwargs 
            Additional arguments may depend on the task type.
        
        '''
    
        
        

        self.name = name

        if isinstance(actions, dict):
            self.actions = actions
        elif isinstance(actions, list):
            self.actions = {a.name: a for a in actions}


        if action_sequence is None:
            self.stageActions(from_deps=True)
        else:
            self.action_sequence = {k: list(v) for k, v in action_sequence.items()}


        self.status = 0  # 0, waiting;  1=running;  2=finished
        self.actions_ti = {}  # relative start time of each action [h]
        
        self.duration = 0.0  # duration must be calculated based on lengths of actions
        self.cost     = 0.0  # cost must be calculated based on the cost of individual actions.
        self.ti       = 0.0  # task start time [h?]
        self.tf       = 0.0  # task end time [h?]
        
        # Calculate duration and cost
        self.calcDuration()  # organizes actions and calculates duration
        self.calcCost()

        print(f"---------------------- Initializing Task '{self.name} ----------------------")        
        print(f"Task '{self.name}' initialized with duration = {self.duration:.2f} h.")
        print(f"Task '{self.name}' initialized with cost     = ${self.cost:.2f} ")

    def getSequenceGraph(self, action_sequence=None, plot=True):
        '''Generate a multi-directed graph that visalizes action sequencing within the task.
                Build a MultiDiGraph with nodes:
                Start -> CP1 -> CP2 -> ... -> End
        
        Checkpoints are computed from action "levels":
          level(a) = 1 if no prerequisites.
          level(a) = 1 + max(level(p) for p in prerequisites)      1 + the largest level among a’s prerequisites.
          Number of checkpoints = max(level) - 1.           
        '''
        if action_sequence is None:
            action_sequence = self.action_sequence
        # Compute levels
        levels: dict[str, int] = {}
        def level_of(a: str, b: set[str]) -> int:
            '''Return the level of action a. b is the set of actions currently being explored'''
                    
            # If we have already computed the level, return it
            if a in levels:
                return levels[a]

            if a in b:
                raise ValueError(f"Cycle detected in action sequence at '{a}' in task '{self.name}'. The action cannot be its own prerequisite.")

            b.add(a)

            # Look up prerequisites for action a.
            pres = action_sequence.get(a, [])
            if not pres:
                lv = 1  # No prerequisites, level 1
            else:
                # If a prerequisites name is not in the dict, treat it as a root (level 1)
                lv = 1 + max(level_of(p, b) if p in action_sequence else 1 for p in pres)
            # b.remove(a)  # if you want to unmark a from the explored dictionary, b, uncomment this line. 
            levels[a] = lv
            return lv

        for a in action_sequence:
            level_of(a, set())
        
        max_level = max(levels.values(), default=1)
        num_cps = max(0, max_level - 1)

        H = nx.MultiDiGraph()

        # Add the Start -> [checkpoints] -> End nodes
        H.add_node("Start")
        for i in range(1, num_cps + 1):
            H.add_node(f"CP{i}")
        H.add_node("End")

        shells = [["Start"]]
        if num_cps > 0:
            # Middle shells
            cps = [f"CP{i}" for i in range(1, num_cps + 1)]
            shells.append(cps)
        shells.append(["End"])

        pos = nx.shell_layout(H, nlist=shells)

        xmin, xmax = -2.0, 2.0  # maybe would need to change those later on.
        pos["Start"] = (xmin, 0)
        pos["End"]   = (xmax, 0)

        # Add action edges
        # Convention:
        #   level 1 actions: Start -> CP1 (or Start -> End if no CPs)
        #   level L actions (2 <= L < max_level): CP{L-1} -> CP{L}
        #   level == max_level actions: CP{num_cps} -> End
        for action, lv in levels.items():
            action = self.actions[action]
            if num_cps == 0:
                # No checkpoints: all actions from Start to End
                H.add_edge("Start", "End", key=action, duration=action.duration, cost=action.cost)
            else:
                if lv == 1:
                    H.add_edge("Start", "CP1", key=action, duration=action.duration, cost=action.cost)
                elif lv < max_level:
                    H.add_edge(f"CP{lv-1}", f"CP{lv}", key=action, duration=action.duration, cost=action.cost)
                else:  # lv == max_level
                    H.add_edge(f"CP{num_cps}", "End", key=action, duration=action.duration, cost=action.cost)                                                    
        # 3. Compute cumulative start time for each level
        level_groups = {}
        for action, lv in levels.items():
            level_groups.setdefault(lv, []).append(action)
        
        level_durations = {lv: max(self.actions[a].duration for a in acts)
                        for lv, acts in level_groups.items()}   

        
        task_duration = sum(level_durations.values())
        level_start_time = {}
        elapsed = 0.0
        cp_string = []
        for lv in range(1, max_level + 1):
            level_start_time[lv] = elapsed
            elapsed += level_durations.get(lv, 0.0)
            # also collect all actions at this level for title
            acts = [a for a, l in levels.items() if l == lv]
            if acts and lv <= num_cps:
                cp_string.append(f"CP{lv}: {', '.join(acts)}")
            elif acts and lv > num_cps:
                cp_string.append(f"End: {', '.join(acts)}")        
        # Assign to self:
        self.duration = task_duration
        self.actions_ti = {a: level_start_time[lv] for a, lv in levels.items()}        
        self.sequence_graph = H
        title_str = f"Task {self.name}. Duration {self.duration:.2f} : " + " | ".join(cp_string)
        if plot:
            fig, ax = plt.subplots()
            # pos = nx.shell_layout(G)
            nx.draw(H, pos, with_labels=True, node_size=500, node_color="lightblue", edge_color='white')

            label_positions = {}  # to store label positions for each edge
            # Group edges by unique (u, v) pairs
            for (u, v) in set((u, v) for u, v, _ in H.edges(keys=True)):
                # get all edges between u and v (dict keyed by edge key)
                edge_dict = H.get_edge_data(u, v)  # {key: {attrdict}, ...}
                n = len(edge_dict)

                # curvature values spread between -0.3 and +0.3 [helpful to visualize multiple edges]
                if n==1:
                    rads = [0]
                    offsets = [0.5]
                else:
                    rads = np.linspace(-0.3, 0.3, n)
                    offsets = np.linspace(0.2, 0.8, n)
                
                # draw each edge
                durations = [d.get("duration", 0.0) for d in edge_dict.values()]
                scale = max(max(durations), 0.0001)  # avoid div by zero
                width_scale = 4.0 / scale  # normalize largest to ~4px

                for rad, offset, (k, d) in zip(rads, offsets, edge_dict.items()):
                    nx.draw_networkx_edges(
                        H, pos, edgelist=[(u, v)], ax=ax,
                        connectionstyle=f"arc3,rad={rad}",
                        arrows=True, arrowstyle="-|>",
                        edge_color="gray",
                        width=max(0.5, d.get("duration", []) * width_scale),
                    )
                    label_positions[(u, v, k)] = offset  # store position for edge label

            ax.set_title(title_str, fontsize=12, fontweight="bold")
            ax.axis("off")
            plt.tight_layout()

        return H
                
                    

    def stageActions(self, from_deps=True):
        '''
        This method stages the action_sequence for a proper execution order.

        Parameters
        ----------
        from_deps : bool
            If True, builds the action_sequence from the dependencies of each action.
            More options will be added in the future.
        '''
        if from_deps:
            # build from dependencies
            def getDeps(action):
                deps = []
                for dep in action.dependencies:
                    deps.append(dep)
                return deps
            
            self.action_sequence = {self.actions[name].name: getDeps(self.actions[name]) for name in self.actions}
            

    def calcDuration(self):
        '''Organizes the actions to be done by this task into the proper order
        based on the action_sequence. This is used to fill out 
        self.actions_ti, self.ti, and self.tf. This method assumes that action.duration
        have already been evaluated for each action in self.actions.
        '''
        # Initialize dictionaries to hold start and finish times
        starts = {}
        finishes = {}

        # Iterate through actions in the sequence
        for action, dep_actions in self.action_sequence.items():
            # Calculate start time as the max finish time of dependencies
            starts[action] = max((finishes[dep] for dep in dep_actions), default=0)

            # get duration from actions
            duration = self.actions[action].duration  # in hours

            # Calculate finish time
            finishes[action] = starts[action] + duration

        # Update self.actions_ti with relative start times
        self.actions_ti = starts

        # Set task start time and finish time
        self.ti = min(starts.values(), default=0)
        self.tf = max(finishes.values(), default=0)

        # Task duration
        self.duration = self.tf - self.ti

    def calcCost(self):
        '''Calculates the total cost of the task based on the costs of individual actions.
        Updates self.cost accordingly. This method assumes that action.cost has 
        already been evaluated for each action in self.actions.
        '''
        total_cost = 0.0
        for action in self.actions.values():
            total_cost += action.cost
        self.cost = total_cost
        return self.cost

    def updateTaskTime(self, newStart=0.0):
        '''Update the start time of all actions based on a new task start time.

        Parameters
        ----------
        newStart : float
            The new start time for the task. All action start times will be adjusted accordingly.
        '''
        # Calculate the time shift
        time_shift = newStart - self.ti

        # Update task start and finish times
        self.ti = newStart
        self.tf += time_shift

        # Update action start times
        for action in self.actions_ti:
            self.actions_ti[action] += self.ti
    



    
    def get_row(self, assets):
        '''Get a matrix of (cost, duration) tuples for each asset to perform this task. Will be a row in the task_asset matrix.
        
        Parameters
        ----------
        assets : list
            A list of all assets available to perform the task.
        
        Returns
        -------
        matrix : array-like
            A 2D array of (cost, duration) tuples indicating the cost and duration for each asset to perform this task.
            Must be 2x len(assets).
        
        '''

        matrix = np.zeros((len(assets), 2))
        # TODO: build row of matrix that holds (cost, duration) tuple of asset / assets (?) to complete task

        # Could look something like...
        ''' 
        for i, asset in enumerate(assets):
            for action in self.actions: # can we do this without the double loop?
                if asset in action.roles:
                    action = self.actions[asset.name]
                    matrix[i, 0] = action.cost
                    matrix[i, 1] = action.duration
                else:
                    matrix[i, 0] = -1 # negative cost/duration means asset cannot perform task
                    matrix[i, 1] = -1
        '''

        return np.zeros((len(assets), 2)) # placeholder, replace with actual matrix


    def GanttChart(self, start_at_zero=True, color_by=None):
        '''Generate a Gantt chart for the task showing the schedule of actions.
        
        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure object containing the Gantt chart.
        ax : matplotlib.axes.Axes
            The axes object containing the Gantt chart.
        '''

        # --- color palette ---
        colors = [
            "lime", "orange", "magenta", "blue",
            "red", "yellow", "cyan", "purple"
        ]        

        fig, ax = plt.subplots(figsize=(10, 10))

        # Prepare data for Gantt chart
        action_names = list(self.actions.keys())
        start_times = [self.actions_ti[name] for name in action_names]
        durations = [self.actions[name].duration for name in action_names]
        
        # Get asset information from action.assets
        all_assets = set()
        all_roles  = set()
        for action in self.actions.values():
            for role, asset in action.assets.items():
                all_assets.add(asset['name'])
                all_roles.add(role)

        # Assign colors
        if color_by == 'asset':
            asset_list = list(all_assets)
            color_dict = {asset: colors[i] for i, asset in enumerate(asset_list)}
        elif color_by == 'role':
            # Flip the colors
            colors = colors[::-1]
            role_list = list(all_roles)
            color_dict = {role: colors[i] for i, role in enumerate(role_list)}
        
        # Generate vertical lines to indicate the start and finish of the whole task
        ax.axvline(x=self.ti, ymin=0, ymax=len(action_names), color='black', linestyle='-', linewidth=2.0)
        ax.axvline(x=self.tf, ymin=0, ymax=len(action_names), color='black', linestyle='-', linewidth=2.0)

        # Create bars for each action
        ht = 0.4
        for i, (name, start, duration) in enumerate(zip(action_names, start_times, durations)):
            opp_i = len(action_names) - i - 1  # to have first action on top
            action = self.actions[name]
            assets = list({asset['name'] for asset in action.assets.values()})
            roles  = list({role for role in action.assets.keys()})
            
            assets = list(set(assets))  # Remove duplicates from assets

            n_assets = len(assets)
            n_roles  = len(roles)
   
            if color_by is None:
                ax.barh(opp_i, duration, color='cyan', left=start, height=ht, align='center')
            elif color_by == 'asset':
                # Compute vertical offsets if multiple assets
                if n_assets == 0:
                    # No assets info
                    ax.barh(i, duration, left=start, height=ht, color='cyan', align='center')
                else:
                    sub_ht = ht / n_assets
                    for j, asset in enumerate(assets):
                        bottom = opp_i - ht/2 + j * sub_ht
                        color = color_dict.get(asset, 'gray')
                        ax.barh(bottom + sub_ht/2, duration, left=start, height=sub_ht * 0.9, 
                                color=color, edgecolor='k', linewidth=0.3, align='center')                             
            elif color_by == 'role':
                # Compute vertical offsets if multiple roles
                if n_roles == 0:
                    # No roles info
                    ax.barh(opp_i, duration, left=start, height=ht, color='cyan', align='center')
                else:
                    sub_ht = ht / n_roles
                    for j, role in enumerate(roles):
                        bottom = opp_i - ht/2 + j * sub_ht
                        color = color_dict.get(role, 'gray')
                        ax.barh(bottom + sub_ht/2, duration, left=start, height=sub_ht * 0.9, 
                                color=color, edgecolor='k', linewidth=0.3, align='center')       
            else:
                color_by = None
                raise Warning(f"color_by option '{color_by}' not recognized. Use 'asset', 'role'. None will be used")

            ax.text(self.ti, opp_i, f'  {name}', va='center', ha='left', color='black')
            ax.axhline(y=opp_i - ht/2, xmin=0, xmax=self.tf, color='gray', linestyle='--', linewidth=0.5)
            ax.axhline(y=opp_i + ht/2, xmin=0, xmax=self.tf, color='gray', linestyle='--', linewidth=0.5)
            ax.axvline(x=start, ymin=0, ymax=len(action_names), color='gray', linestyle='--', linewidth=0.5)
        
        # Set y-ticks and labels
        ax.set_yticks(range(len(action_names)))
        ax.set_yticklabels([])

        # Set labels and title
        ax.set_xlabel('time (hrs.)')
        ax.set_title(f'Gantt Chart for Task: {self.name}')

        if color_by == 'asset':
            handles = [plt.Rectangle((0, 0), 1, 1, color=color_dict[a]) for a in all_assets]
            ax.legend(handles, all_assets, title='Assets', bbox_to_anchor=(1.02, 1), loc='upper right')
        elif color_by == 'role':
            handles = [plt.Rectangle((0, 0), 1, 1, color=color_dict[a]) for a in all_roles]
            ax.legend(handles, all_roles, title='Roles', bbox_to_anchor=(1.02, 1), loc='upper right')            

        if start_at_zero:
            ax.set_xlim(0, self.tf + 1)
        # Create a grid and adjust layout
        # ax.grid(True)
        plt.tight_layout()
        return fig, ax

    def chart(self, start_at_zero=True):
        '''Generate a chart grouped by asset showing when each asset is active across all actions.

        Parameters
        ----------
        start_at_zero : bool, optional
            If True, the x-axis starts at zero. Defaults to True.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure object containing the Gantt chart.
        ax : matplotlib.axes.Axes
            The axes object containing the Gantt chart.
        '''
        pass
        