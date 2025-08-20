"""Action base class"""

import numpy as np
import matplotlib.pyplot as plt

import moorpy as mp
from moorpy.helpers import set_axes_equal
from moorpy import helpers
import yaml
from copy import deepcopy

#from shapely.geometry import Point, Polygon, LineString
from famodel.seabed import seabed_tools as sbt
from famodel.mooring.mooring import Mooring
from famodel.platform.platform import Platform
from famodel.anchors.anchor import Anchor
from famodel.mooring.connector import Connector
from famodel.substation.substation import Substation
from famodel.cables.cable import Cable
from famodel.cables.dynamic_cable import DynamicCable
from famodel.cables.static_cable import StaticCable
from famodel.cables.cable_properties import getCableProps, getBuoyProps, loadCableProps,loadBuoyProps
from famodel.cables.components import Joint
from famodel.turbine.turbine import Turbine
from famodel.famodel_base import Node

# Import select required helper functions
from famodel.helpers import (check_headings, head_adjust, getCableDD, getDynamicCables, 
                            getMoorings, getAnchors, getFromDict, cleanDataTypes, 
                            getStaticCables, getCableDesign, m2nm, loadYAML, 
                            configureAdjuster, route_around_anchors)


def incrementer(text):
    split_text = text.split()[::-1]
    for ind, spl in enumerate(split_text):
        try:
            split_text[ind] = str(int(spl) + 1)
            break
        except ValueError:
            continue
    return " ".join(split_text[::-1])


def increment_name(name):
    '''Increments an end integer after a dash in a name'''
    name_parts = name.split(sep='-')
    
    # if no numeric suffix yet, add one
    if len(name_parts) == 1 or not name_parts[-1].isdigit():  
        name = name+'-0'
    # otherwise there must be a suffix, so increment it
    else:
        name_parts[-1] = str( 1 + int(name_parts[-1]))
        
        name = '-'.join(name_parts) # reassemble name string
        
    return name


class Action():
    '''
    An Action is a general representation of a marine operations action
    that involves manipulating a system/design/structure using assets/
    equipment. The Action base class contains generic routines and parameters.
    Specialized routines for performing each action should be set up in
    subclasses.
    '''
    
    def __init__(self, actionType, name, **kwargs):
        '''Create an action object...
        It must be given a name.
        The remaining parameters should correspond to items in the actionType dict...
        
        Parameters
        ----------
        actionType : dict
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
        self.assets = {}  # dict of named roles for the vessel(s) or port required to perform the action
        self.requirements = {}  # the capabilities required of each asset role assets (same keys as self.assets)
        self.objectList   = []  # all objects that could be acted on
        self.dependencies = {}  # list of other actions this one depends on
        
        self.type = getFromDict(actionType, 'type', dtype=str)
        self.name = name
        self.status = 0  # 0, waiting;  1=running;  2=finished
        
        self.duration = getFromDict(actionType, 'duration', default=3)
        
        '''
        # Create a dictionary of supported object types (with empty entries)
        if 'objects' in actionType:  #objs = getFromDict(actionType, 'objects', shape=-1, default={})
            for obj in actionType['objects']:  # go through keys in objects dictionary
                self.objectList[obj] = None    # make blank entries with the same names
        
        
        # Process objects according to the action type
        if 'objects' in kwargs:  #objects = getFromDict(kwargs, objects, default=[])
            for obj in kwargs['objects']:
                objType = obj.__class__.__name__.lower()
                if objType in self.objectList:
                    self.objectList[objType] = obj
                else:
                    raise Exception(f"Object type '{objType}' is not in the action's supported list.")
        '''
        
        # Process objects to be acted upon
        # make list of supported object type names
        if 'objects' in actionType:  
            if isinstance(actionType['objects'], list):
                supported_objects = actionType['objects']
            elif isinstance(actionType['objects'], dict):
                supported_objects = list(actionType['objects'].keys())
        else:
            supported_objects = []
        
        # Add objects to the action's object list as long as they're supported
        if 'objects' in kwargs:
            for obj in kwargs['objects']:
                objType = obj.__class__.__name__.lower() # object class name
                if objType in supported_objects:
                    self.objectList.append(obj)
                else:
                    raise Exception(f"Object type '{objType}' is not in the action's supported list.")
        
        # Create placeholders for asset roles based on the "requirements"
        if 'roles' in actionType:
            for asset, caplist in actionType['roles'].items():
                self.assets[asset] = None  # each asset role holds a None value until assigned
                self.requirements[asset] = {}
                for cap in caplist:
                    self.requirements[asset][cap] = 0  # fill in each required metric with zero to start with?
                    
        
        # Process dependencies
        if 'dependencies' in kwargs:
            for dep in kwargs['dependencies']:
                self.dependencies[dep.name] = dep
        
        # Process some optional kwargs depending on the action type
        
        
    
        
    def addDependency(self, dep):
        '''Registers other action as a dependency of this one.'''
        self.dependencies[dep.name] = dep
        # could see if already a dependency and raise a warning if so...
    
    
    def checkAsset(self, role_name, asset):
        '''Checks if a specified asset has sufficient capabilities to fulfil
        a specified role in this action.
        '''
        
        # Make sure role_name is valid for this action
        if not role_name in self.assets.keys():
            raise Exception(f"The specified role name '{role_name}' is not a named asset role in this action.")
        
        for cap, req in self.requirements['role_name']:
            if asset.capabilities[cap] >= req:  # <<< this is pseudocode. Needs to look at capability numbers! <<<
                pass
                # so far so good
            else:
                return False  # at least on capability is not met, so return False
                
        return True
    
    
    def assignAsset(self, role_name, asset):
        '''Assigns a vessel or port to a certain role in the action.
        
        Parameters
        ----------
        role_name : string
            Name of the asset role being filled (must be in the action's list)
        asset : Vessel or Port object
            The asset to be registered with the class.
        '''
        
        # Make sure role_name is valid for this action
        if not role_name in self.assets.keys():
            raise Exception(f"The specified role name '{role_name}' is not a named asset role in this action.")
        
        self.assets[role_name] = asset
        
    
    def calcDurationAndCost(self):
        pass
    
    
    def evaluateAssets(self, assets):
        '''Check whether an asset can perform the task, and if so calculate
        the time and cost associated with using those assets.
        
        assets : dict
            Dictionary of role_name, asset object pairs for assignment to the action.
            Each asset is a vessel or port object.
        
        '''
        
        # Assign each specified asset to its respective role
        for akey, aval in assets.items():
            self.assignAsset(akey, aval)
            
        self.calcDurationAndCost()
        
        
        # can store the cost and duration in self as well
        
    
    
    # ----- Below are drafts of methods for use by the engine -----
    """
    def begin(self):
        '''Take control of all objects'''
        for vessel in self.vesselList:
            vessel._attach_to(self)
        for object in self.objectList:
            object._attach_to(self)
    
    
    def end(self):
        '''Release all objects'''
        for vessel in self.vesselList:
            vessel._detach_from()
        for object in self.objectList:
            object._detach_from()
    """
    
    def timestep(self):
        '''Advance the simulation of this action forward one step in time.'''
        
        # (this is just documenting an idea for possible future implementation)
        # Perform the hourly action of the task
        
        if self.type == 'tow':
            # controller - make sure things are going in right direction...
            # (switch mode if need be)
            if self.mode == 0 :  # gathering vessels
                ves = self.assets['vessel']
                dr = self.r_start - ves.r
                ves.setCourse(dr)  # sets vessel velocity
                
                # if vessel is stopped (at the object), time to move
                if np.linalg.norm(ves.v) == 0:
                    self.mode = 1
                    
            if self.mode == 1:  # towing
                ves = self.assets['vessel']
                dr = self.r_finish - ves.r
                ves.setCourse(dr)  # sets vessel velocity
                
                # if all vessels are stopped (at the final location), time to end
                if np.linalg.norm(ves.v) == 0:
                    self.mode = 2
            
            if self.mode == 2:  # finished
                self.end()
        

