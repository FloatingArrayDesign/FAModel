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
    
    def __init__(self, actions, action_sequence, name, **kwargs):
        '''Create an action object...
        It must be given a name and a list of actions.
        The action list should be by default coherent with actionTypes dictionary.
        
        Parameters
        ----------
        actions : list
            A list of all actions that are part of this task.
        action_sequence : dict
            A dictionary where each key is the name of each action, and the values are
            each a list of which actions (by name) must be completed before the current
            one.
        name : string
            A name for the action. It may be appended with numbers if there
            are duplicate names.
        kwargs 
            Additional arguments may depend on the task type.
        
        '''
        
        # Make a dict by name of all actions that are carried out in this task
        self.actions = {}
        for act in actions:
            self.actions[act.name] = act
        
        
        # Create a graph of the sequence of actions in this task based on action_sequence
        self.getSequenceGraph(action_sequence)

        self.name = name
        
        self.status = 0  # 0, waiting;  1=running;  2=finished
        
        self.actions_ti = {}  # relative start time of each action [h]
        self.t_actions = {}  # timing of task's actions, relative to t1 [h]
        # t_actions is a dict with keys same as action names, and entries of [t1, t2]
        
        self.duration = 0  # duration must be calculated based on lengths of actions
        self.cost     = 0  # cost must be calculated based on the cost of individual actions.
        self.ti =0  # task start time [h?]
        self.tf =0  # task end time [h?]
        
        # what else do we need to initialize the task?
        
        # internal graph of the actions within this task.
        self.G = self.getTaskGraph()
    
    
    def organizeActions(self):
        '''Organizes the actions to be done by this task into the proper order
        based on the strategy of this type of task...
        '''
        
        if self.type == 'parallel_anchor_install':
            
            pass
            # make a graph that reflects this strategy?
    
    
    def calcDuration(self):
        '''Calculates the duration of the task based on the durations of the
        individual actions and their order of operation.'''
        
        # Does Rudy have graph-based code that can do this?
        
        
        
        
        
        
    def getSequenceGraph(self, action_sequence):
        '''Generate a multi-directed graph that visalizes action sequencing within the task.
                Build a MultiDiGraph with nodes:
                Start -> CP1 -> CP2 -> ... -> End
        
        Checkpoints are computed from action "levels":
          level(a) = 1 if no prerequisites.
          level(a) = 1 + max(level(p) for p in prerequisites)      1 + the largest level among a’s prerequisites.
          Number of checkpoints = max(level) - 1.           
        '''

        # Compute levels
        levels: dict[str, int] = {}
        def level_of(a: str, b: set[str]) -> int:
            '''Return the level of action a. b is the set of actions currently being explored'''

            # If we have already computed the level, return it
            if a in levels:
                return levels[a]
            
            # The action cannot be its own prerequisite 
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

        xmin, xmax = -2.0, 2.0
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

        fig, ax = plt.subplots()
        # pos = nx.shell_layout(G)
        nx.draw(H, pos, with_labels=True, node_size=500, node_color="lightblue", edge_color='white')

        # Group edges by unique (u, v) pairs
        for (u, v) in set((u, v) for u, v, _ in H.edges(keys=True)):
            # get all edges between u and v (dict keyed by edge key)
            edge_dict = H.get_edge_data(u, v)  # {key: {attrdict}, ...}
            n = len(edge_dict)

            # curvature values spread between -0.3 and +0.3
            if n==1:
                rads = [0]
            else:
                rads = np.linspace(-0.3, 0.3, n)
            
            # draw each edge
            durations = [d.get("duration", 0.0) for d in edge_dict.values()]
            scale = max(max(durations), 0.0001)
            width_scale = 4.0 / scale  # normalize largest to ~4px

            for rad, (k, d) in zip(rads, edge_dict.items()):
                nx.draw_networkx_edges(
                    H, pos, edgelist=[(u, v)], ax=ax,
                    connectionstyle=f"arc3,rad={rad}",
                    arrows=True, arrowstyle="-|>",
                    edge_color="gray",
                    width=max(0.5, d.get("duration", []) * width_scale),
                )

        # --- after drawing edges ---
        edge_labels = {}
        for u, v, k, d in H.edges(keys=True, data=True):
            # each edge may have a unique key; include it in the label if desired
            label = k.name
            edge_labels[(u, v, k)] = label

        nx.draw_networkx_edge_labels(
            H,
            pos,
            edge_labels=edge_labels,
            font_size=8,
            label_pos=0.5,   # position along edge (0=start, 0.5=middle, 1=end)
            rotate=False     # keep labels horizontal
        )

        ax.axis("off")
        plt.tight_layout()
        plt.show()
        
        self.sequence_graph = H
        return H
                
        

    def getTaskGraph(self, plot=True):
        '''Generate a graph of the action dependencies.
        '''
        
        # Create the graph
        G = nx.DiGraph()
        for item, data in self.actions.items():
            for dep in data.dependencies:
                G.add_edge(dep, item, duration=data.duration)  # Store duration as edge attribute

        # Compute longest path & total duration
        longest_path = nx.dag_longest_path(G, weight='duration')
        longest_path_edges = list(zip(longest_path, longest_path[1:]))  # Convert path into edge pairs

        total_duration = sum(self.actions[node].duration for node in longest_path) 
        return G  

    
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

