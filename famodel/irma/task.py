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
        self.actions = {a.name: a for a in actions.values()}

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
            self.actions_ti[action] += time_shift
    



    
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

