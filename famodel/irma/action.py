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
    '''
    Increments the last integer found in a string.

    Inputs
    ------
    `text` : `str`
        The input string to increment.

    Returns
    -------
    `str`
        The incremented string.
    '''
    split_text = text.split()[::-1]
    for ind, spl in enumerate(split_text):
        try:
            split_text[ind] = str(int(spl) + 1)
            break
        except ValueError:
            continue
    return " ".join(split_text[::-1])


def increment_name(name):
    '''
    Increments an end integer after a dash in a name.

    Inputs
    ------
    `name` : `str`
        The input name string.

    Returns
    -------
    `str`
        The incremented name string.
    '''
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
        
        Inputs
        ----------
        `actionType` : `dict`
            Dictionary defining the action type (typically taken from a yaml).
        `name` : `string`
            A name for the action. It may be appended with numbers if there
            are duplicate names.
        `kwargs`
            Additional arguments may depend on the action type and typically
            include a list of FAModel objects that are acted upon, or
            a list of dependencies (other action names/objects).

        Returns
        -------
        `None`
        '''
        
        # list of things that will be controlled during this action
        self.assets = {}  # dict of named roles for the vessel(s) or port required to perform the action
        self.requirements = {}  # the capabilities required of each role (same keys as self.assets)
        self.objectList   = []  # all objects that could be acted on
        self.dependencies = {}  # list of other actions this one depends on
        
        self.type = getFromDict(actionType, 'type', dtype=str)
        self.name = name
        self.status = 0  # 0, waiting;  1=running;  2=finished
        
        self.duration = getFromDict(actionType, 'duration', default=0) # this will be overwritten by calcDurationAndCost. TODO: or should it overwrite any duration calculation?
        self.cost = 0 # this will be overwritten by calcDurationAndCost

        self.supported_objects = [] # list of FAModel object types supported by the action
        
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
                self.supported_objects = actionType['objects']
            elif isinstance(actionType['objects'], dict):
                self.supported_objects = list(actionType['objects'].keys())
        
        # Add objects to the action's object list as long as they're supported
        if 'objects' in kwargs:
            for obj in kwargs['objects']:
                objType = obj.__class__.__name__.lower() # object class name
                if objType in self.supported_objects:
                    self.objectList.append(obj)
                else:
                    raise Exception(f"Object type '{objType}' is not in the action's supported list.")
        
        # Create placeholders for asset roles based on the "requirements"
        if 'roles' in actionType:
            for role, caplist in actionType['roles'].items():
                self.requirements[role] = {key: None for key in caplist}  # each role requirment holds a dict of capabilities with values set to None for now
                for cap in caplist:
                    # self.requirements[role][cap] = {}  # fill in each required capacity with empty dict
                    self.requirements[role][cap] = {'area_m2': 1000, 'max_load_t': 1600}  # dummy values for now, just larger than MPSV_01 values to trigger failure

                self.assets[role] = None  # placeholder for the asset assigned to this role

        # Process dependencies
        if 'dependencies' in kwargs:
            for dep in kwargs['dependencies']:
                self.dependencies[dep.name] = dep
        
        # Process some optional kwargs depending on the action type
    
        
    def addDependency(self, dep):
        '''
        Registers other action as a dependency of this one.

        Inputs
        ------
        `dep` : `Action`
            The action to be added as a dependency.

        Returns
        -------
        `None`
        '''
        self.dependencies[dep.name] = dep
        # could see if already a dependency and raise a warning if so...
    

    def assignObjects(self, objects):
        '''
        Adds a list of objects to the actions objects list, 
        checking they are valid for the actions supported objects.

        Inputs
        ------
        `objects` : `list`
            A list of FAModel objects to be added to the action.

        Returns
        -------
        `None`
        ''' 

        for obj in objects:
            objType = obj.__class__.__name__.lower() # object class name
            if objType in self.supported_objects:
                if obj in self.objectList:
                    print(f"Warning: Object '{obj}' is already in the action's object list.")
                self.objectList.append(obj)
            else:
                raise Exception(f"Object type '{objType}' is not in the action's supported list.")



    # def setUpCapability(self):
    #     '''
    #     Example of what needs to happen to create a metric.
    #
    #     Inputs
    #     ------
    #     `None`
    #
    #     Returns
    #     -------
    #     `None`
    #     '''
    #     # WIP: example of what needs to happen to create a metric

    #     # figure out how to assign required metrics to capabilies in the roles based on the objects
    #     for role, caps in self.requirements.items():
    #         for cap, metrics in caps.items():
    #             for obj in self.objectList:
    #                 # this is for the deck_space capability
    #                 metrics = {'area_m2': obj.area, 'max_load_t': obj.mass / 1000} # / 1000 to convert kg to T
    #                 metrics.update(obj.get_capability_metrics(cap))
    #     pass

    
    def checkAsset(self, role_name, asset):
        '''
        Checks if a specified asset has sufficient capabilities to fulfil
        a specified role in this action.

        Inputs
        ------
        `role_name` : `string`
            The name of the role to check.
        `asset` : `dict`
            The asset to check against the role's requirements. 

        Returns
        -------
        `bool`
            True if the asset meets the role's requirements, False otherwise.
        `str`
            A message providing additional information about the check.
        '''        

        # Make sure role_name is valid for this action
        if not role_name in self.assets.keys():
            raise Exception(f"The specified role '{role_name}' is not a named in this action.")
        
        if self.assets[role_name] is not None: 
            return False, f"Role '{role_name}' is already filled in action '{self.name}'."

        for capability in self.requirements[role_name].keys():

            if capability in asset['capabilities'].keys(): # check capability is in asset
                
                # TODO: does this work if there are no metrics in a capability? This should be possible, as not all capabilities will require a constraint. 
                for metric in self.requirements[role_name][capability].keys(): # loop over the capacity requirements for the capability (if more than one)
                    
                    if metric not in asset['capabilities'][capability].keys(): # value error because capabilities are defined in capabilities.yaml. This should only be triggered if something has gone wrong (i.e. overwriting values somewhere)
                        raise ValueError(f"The '{capability}' capability does not have metric: '{metric}'.")

                    if self.requirements[role_name][capability][metric] > asset['capabilities'][capability][metric]: # check requirement is met
                        return False, f"The asset does not have sufficient '{metric}' for '{capability}' capability in '{role_name}' role of '{self.name}' action."

                return True, 'All capabilities in role met'
            
            else:
                return False, f"The asset does not have the '{capability}' capability for '{role_name}' role of '{self.name}' action."  # a capability is not met
        
    
    def calcDurationAndCost(self):
        '''
        Calculates duration and cost for the action. The structure here is dependent on `actions.yaml`.
        TODO: finish description

        Inputs
        ------
        `None`

        Returns
        -------
        `None`
        '''
        
        # Check that all roles in the action are filled
        for role_name in self.requirements.keys():
            if self.assets[role_name] is None:
                raise Exception(f"Role '{role_name}' is not filled in action '{self.name}'. Cannot calculate duration and cost.")

        # --- Towing & Transport ---
        if self.type == 'tow':
            pass
        elif self.type == 'transport_components':
            pass

        # --- Mooring & Anchors ---
        elif self.type == 'install_anchor':
            pass
        elif self.type == 'retrieve_anchor':
            pass
        elif self.type == 'load_mooring':
            pass
        elif self.type == 'lay_mooring':
            pass
        elif self.type == 'mooring_hookup':
            pass

        # --- Heavy Lift & Installation ---
        elif self.type == 'install_wec':
            pass
        elif self.type == 'install_semisub':
            pass
        elif self.type == 'install_spar':
            pass
        elif self.type == 'install_tlp':
            pass
        elif self.type == 'install_wtg':
            pass

        # --- Cable Operations ---
        elif self.type == 'lay_cable':
            pass
        elif self.type == 'retrieve_cable':
            pass
        elif self.type == 'lay_and_bury_cable':
            pass
        elif self.type == 'backfill_rockdump':
            pass

        # --- Survey & Monitoring ---
        elif self.type == 'site_survey':
            pass
        elif self.type == 'monitor_installation':
            pass
        else:
            raise ValueError(f"Action type '{self.type}' not recognized.")
        
        return self.duration, self.cost


    def evaluateAssets(self, assets):
        '''
        Checks assets for all the roles in the action. This calls `checkAsset()`
        for each role/asset pair and then calculates the duration and
        cost for the action as if the assets were assigned. Does not assign 
        the asset(s) to the action. WARNING: this function will clear the values 
        (but not keys) in `self.assets`.

        Inputs
        ------
        `assets` : `dict`
            Dictionary of {role_name: asset} pairs for assignment of the 
            assets to the roles in the action.

        Returns
        -------
        `cost` : `float`
            Estimated cost of using the asset.
        `duration` : `float`
            Estimated duration of the action when performed by asset.
        '''
        
        # Check each specified asset for its respective role
        for role_name, asset in assets.items():
            assignable, message = self.checkAsset(role_name, asset)
            if assignable:
                self.assets[role_name] = asset # Assignment required for calcDurationAndCost(), will be cleared later
            else:
                print('INFO: '+message+' Action cannot be completed by provided asset list.')
                return -1, -1 # return negative values to indicate incompatibility. Loop is terminated becasue assets not compatible for roles. 
        
        # Check that all roles in the action are filled
        for role_name in self.requirements.keys():
            if self.assets[role_name] is None:
                
                for role_name in assets.keys(): # Clear the assets dictionary
                    assets[role_name] = None
                raise Exception(f"Role '{role_name}' is not filled in action '{self.name}'. Cannot calculate duration and cost.") # possibly just a warning and not an exception?


        duration, cost = self.calcDurationAndCost()

        for role_name in assets.keys(): # Clear the assets dictionary
            assets[role_name] = None

        return duration, cost # values returned here rather than set because will be used to check compatibility and not set properties of action
    

    def assignAsset(self, role_name, asset):
        '''
        Checks if asset can be assigned to an action. 
        If yes, assigns asset to role in the action.

        Inputs
        ------
        `role_name` : `str`
            The name of the role to which the asset will be assigned.
        `asset` : `dict`
            The asset to be assigned to the role.

        Returns
        -------
        `None`
        '''
        # Make sure role_name is valid for this action
        if not role_name in self.assets.keys():
            raise Exception(f"The specified role name '{role_name}' is not in this action.")

        if self.assets[role_name] is not None:
            raise Exception(f"Role '{role_name}' is already filled in action '{self.name}'.")

        assignable, message = self.checkAsset(role_name, asset)
        if assignable:
            self.assets[role_name] = asset
        else:
            raise Exception(message) # throw error message

    def assignAssets(self, assets):
        '''
        Assigns assets to all the roles in the action. This calls
        `assignAsset()` for each role/asset pair and then calculates the
        duration and cost for the action. Similar to `evaluateAssets()`
        however here assets are assigned and duration and cost are 
        set after evaluation.

        Inputs
        ------
        `assets` : `dict`
            Dictionary of {role_name: asset} pairs for assignment of the 
            assets to the roles in the action.

        Returns
        -------
        `None`
        '''
        
        # Assign each specified asset to its respective role
        for role_name, asset in assets.items():
            self.assignAsset(role_name, asset)
        
        # Check that all roles in the action are filled
        for role_name in self.requirements.keys():
            if self.assets[role_name] is None:
                raise Exception(f"Role '{role_name}' is not filled in action '{self.name}'. Cannot calculate duration and cost.") # possibly just a warning and not an exception?

        self.calcDurationAndCost()
    
    
    # ----- Below are drafts of methods for use by the engine -----
    """
    def begin(self):
        '''
        Take control of all objects.

        Inputs
        ------
        `None`

        Returns
        -------
        `None`
        '''
        for vessel in self.vesselList:
            vessel._attach_to(self)
        for object in self.objectList:
            object._attach_to(self)
    
    
    def end(self):
        '''
        Release all objects.

        Inputs
        ------
        `None`

        Returns
        -------
        `None`
        '''
        for vessel in self.vesselList:
            vessel._detach_from()
        for object in self.objectList:
            object._detach_from()
    """
    
    def timestep(self):
        '''
        Advance the simulation of this action forward one step in time.

        Inputs
        ------
        `None`

        Returns
        -------
        `None`
        '''
        
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
        

