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
    that follow a predefined sequency/strategy. There can be multiple
    tasks that achieve the same end, each providing an alternative strategy.
    Each Task consists of a set of Actions with internal dependencies.
    
    
    For now, we'll assume each Task must be port-to-port,
    i.e. its vessel(s) must start and end at port over the course of the task.
    
    
    
    '''
    
    def __init__(self, actionList, name, **kwargs):
        '''Create an action object...
        It must be given a name and a list of actions.
        The action list should be by default coherent with actionTypes dictionary.
        
        Parameters
        ----------
        actionList : list
            A list of all actions that are part of this task.
        name : string
            A name for the action. It may be appended with numbers if there
            are duplicate names.
        kwargs 
            Additional arguments may depend on the task type.
        
        '''
        
        self.actionList   = actionList  # all actions that are carried out in this task
        self.name = name
        self.status = 0  # 0, waiting;  1=running;  2=finished
        
        self.duration = 0  # duration must be calculated based on lengths of actions
        self.cost     = 0  # cost must be calculated based on the cost of individual actions.

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

    def getTaskGraph(self):
        '''Generate a graph of the action dependencies.
        '''
        
        # Create the graph
        G = nx.DiGraph()
        for item, data in self.actionList.items():
            for dep in data.dependencies:
                G.add_edge(dep, item, duration=data.duration)  # Store duration as edge attribute

        # Compute longest path & total duration
        longest_path = nx.dag_longest_path(G, weight='duration')
        longest_path_edges = list(zip(longest_path, longest_path[1:]))  # Convert path into edge pairs
        total_duration = sum(self.actionList[node].duration for node in longest_path)
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
            for action in self.actionList: # can we do this without the double loop?
                if asset in action.roles:
                    action = self.actionList[asset.name]
                    matrix[i, 0] = action.cost
                    matrix[i, 1] = action.duration
                else:
                    matrix[i, 0] = -1 # negative cost/duration means asset cannot perform task
                    matrix[i, 1] = -1
        '''

        return np.zeros((len(assets), 2)) # placeholder, replace with actual matrix