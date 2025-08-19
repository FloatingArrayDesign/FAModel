"""Action base class"""

import numpy as np
import matplotlib.pyplot as plt

import moorpy as mp
from moorpy.helpers import set_axes_equal
from moorpy import helpers
import yaml
from copy import deepcopy


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
    
    def __init__(self, taskType, name, **kwargs):
        '''Create an action object...
        It must be given a name.
        The remaining parameters should correspond to items in the actionType dict...
        
        Parameters
        ----------
        taskType : dict
            Dictionary defining the action type (typically taken from a yaml).
        name : string
            A name for the action. It may be appended with numbers if there
            are duplicate names.
        kwargs 
            Additional arguments may depend on the action type and typically
            include a list of FAModel objects that are acted upon, or
            a list of dependencies (other action names/objects).
        
        '''
        
        # list of things that will be controlled during this action
        self.vesselList   = []  # all vessels required for the action
        self.objectList   = []  # all objects that are acted on
        self.actionList   = []  # all actions that are carried out in this task
        self.dependencies = {}  # list of other tasks this one depends on
        
        self.type = getFromDict(taskType, 'type', dtype=str)
        self.name = name
        self.status = 0  # 0, waiting;  1=running;  2=finished
        
        self.duration = 0  # duration must be calculated based on lengths of actions
        
        # what else do we need to initialize the task?
        
        # internal graph of the actions within this task?
    
    
    def organizeActions(self, actions):
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
        