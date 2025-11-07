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
        self.assets = {}        # dict of named roles for the vessel(s) or port required to perform the action
        self.requirements = {}  # capabilities required of each role (same keys as self.assets)
        self.objectList   = []  # all objects that could be acted on
        self.dependencies = {}  # list of other actions this one depends on
        
        self.actionType = actionType  # <— keep the YAML dict on the instance
        
        self.type = getFromDict(actionType, 'type', dtype=str)
        self.name = name
        self.status = 0  # 0, waiting;  1=running;  2=finished
        
        self.duration = getFromDict(actionType, 'duration', default=0) # this will be overwritten by calcDurationAndCost. TODO: or should it overwrite any duration calculation?
        self.cost = 0 # this will be overwritten by calcDurationAndCost
        self.ti = 0  # action start time [h?]
        self.tf = 0  # action end time [h?]

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
        
        # Create placeholders for asset roles based on the "requirements"
        if 'roles' in actionType:
            for role, caplist in actionType['roles'].items():
                self.requirements[role] = {key: {} for key in caplist}  # each role requirment holds a dict of capabilities with each capability containing a dict of metrics and values, metrics dict set to empty for now. 
                self.assets[role] = None  # placeholder for the asset assigned to this role

        # Process objects to be acted upon. NOTE: must occur after requirements and assets placeholders have been assigned. 
        # make list of supported object type names
        if 'objects' in actionType:  
            if isinstance(actionType['objects'], list):
                self.supported_objects = actionType['objects']
            elif isinstance(actionType['objects'], dict):
                self.supported_objects = list(actionType['objects'].keys())
        
        # Add objects to the action's object list as long as they're supported
        if 'objects' in kwargs:
            self.assignObjects(kwargs['objects'])

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
    

    def getMetrics(self, cap, met, obj):
        '''
        Retrieves the minimum metric(s) for a given capability required to act on target object.
        A metric is the number(s) associated with a capability. A capability is what an action 
        role requires and an asset has. 

        These minimum metrics are assigned to capabilities in the action's role in `assignObjects`.

        Inputs
        ------
        `cap` : `str`
            The capability for which the metric is to be retrieved.
        `met` : `dict`
            The metrics dictionary containing any existing metrics for the capability.
        `obj` : FAModel object
            The target object on which the capability is to be acted upon.

        Returns
        -------
        `metrics` : `dict`
            The metrics and values for the specified capability and object.

        '''

        metrics = met  # metrics dict with following form: {metric_1 : required_value_1, ...}. met is assigned here in case values have already been assigned
        objType = obj.__class__.__name__.lower()

        """
        Note to devs:
        This function contains hard-coded evaluations of all the possible combinations of capabilities and objects. 
        The intent is we generate the minimum required <metric> of a given <capability> to work with the object. An
        example would be minimum bollard pull required to tow out a platform. The capabilities (and their metrics) 
        are from capabilities.yaml and the objects are from objects.yaml. There is a decent ammount of assumptions 
        made here so it is important to document sources where possible. 

        Some good preliminary work on this is in https://github.com/FloatingArrayDesign/FAModel/blob/IOandM_development/famodel/installation/03_step1_materialItems.py 

        ### Code Explanation ###
        This function has the following structure 

        ```
        if cap == <target_cap>:
            # some comments

            if objType == 'mooring':
                metric_value = calc based on obj
            elif objType == 'platform':
                metric_value = calc based on obj
            elif objType == 'anchor':
                metric_value = calc based on obj
            elif objType == 'component':
                metric_value = calc based on obj
            elif objType == 'turbine':
                metric_value = calc based on obj
            elif objType == 'cable':
                metric_value = calc based on obj
            else:
                metric_value = -1
                
            # Assign the capabilties metrics (keep existing metrics already in dict if larger than calc'ed value)
            metrics[<target_cap_metric>] = metric_value if metric_value > metrics.get(<target_cap_metric>) else metrics.get(<target_cap_metric>)
        ```
        
        Some of the logic for checking object types can be omitted if it doesnt make sense. For example, the chain_locker capability 
        only needs to be checked against the Mooring object. The comment `# object logic checked` shows that the logic in that capability
        has been thought through. 

        A metric_value of -1 indicates the object is not compatible with the capability. This is indicated by a warning printed at the end. 

        A completed example of what this can look like is the line_reel capability.
        """


        if cap == 'deck_space':
            # logic for deck_space capability (platforms and sites not compatible)
            # TODO: how do we account for an action like load_mooring (which has two roles, 
            # representing vessels to be loaded). The combined deck space of the carriers
            # should be the required deck space for the action. Right now I believe it is
            # set up that only one asset can fulfill the capability minimum. 

            # object logic checked
            if objType == 'mooring':
                pass
            elif objType == 'platform':
                pass
            elif objType == 'anchor':
                pass
            elif objType == 'component':
                pass
            elif objType == 'turbine':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['area_m2'] = None if None > metrics.get('area_m2') else metrics.get('area_m2')
            # metrics['max_load_t'] = None if None > metrics.get('max_load_t') else metrics.get('max_load_t')

        elif cap == 'chain_locker':
            # logic for chain_locker capability (only mooring objects compatible)
            # object logic checked

            if objType == 'mooring':

                # set baseline values for summation
                vol = 0
                length = 0

                for i, sec in enumerate(obj.dd['sections']): # add up the volume and length of all chain in the object
                    if sec['type']['chain']:
                        diam = sec['type']['d_nom']           # diameter [m]
                        vol += 0.0 # TODO: calculate chain_locker volume from sec['L'] and diam. Use Delmar data from Rudy. Can we make function of chain diam?
                        length += sec['L']                     # length [m]

            else:
                vol = -1

            # Assign the capabilties metrics
            metrics['volume_m3'] = vol if vol > metrics.get('volume_m3') else metrics.get('volume_m3')

        elif cap == 'line_reel':
            # logic for line_reel capability (only mooring objects compatible)
            # object logic checked, complete

            if objType == 'mooring':

                # set baseline values for summation
                vol = 0
                length = 0

                for i, sec in enumerate(obj.dd['sections']): # add up the volume and length of all non_chain line in the object
                    if not sec['type']['chain']: # any line type thats not chain
                        vol += sec['L'] * np.pi * (sec['type']['d_nom'] / 2) ** 2 # volume [m^3]
                        length += sec['L']                     # length [m]

            else:
                vol = -1
                length = -1

            # Assign the capabilties metrics
            metrics['volume_m3'] = vol if vol > metrics.get('volume_m3') else metrics.get('volume_m3')
            metrics['rope_capacity_m'] = length if length > metrics.get('rope_capacity_m') else metrics.get('rope_capacity_m')

        elif cap == 'cable_reel':
            # logic for cable_reel capability (only cable objects compatible)
            # object logic checked
            vol = 0
            length = 0
            '''
            if objType == 'cable':
                for cable in cables: # TODO: figure out this iteration
                    if cable is cable and not other thing in cables object: # TODO figure out how to only check cables, not j-tubes or any other parts
                        vol += cable['L'] * np.pi * (cable['type']['d_nom'] / 2) ** 2
                        length += cable['L']                     # length [m]
            else:
                vol = -1
                length = -1
            '''
            # Assign the capabilties metrics
            metrics['volume_m3'] = vol if vol > metrics.get('volume_m3') else metrics.get('volume_m3')
            metrics['cable_capacity_m'] = length if length > metrics.get('cable_capacity_m') else metrics.get('cable_capacity_m')

        elif cap == 'winch':
            # logic for winch capability
            if objType == 'mooring':
                pass
            elif objType == 'platform':
                pass
            elif objType == 'anchor':
                pass
            elif objType == 'component':
                pass
            elif objType == 'turbine':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # # Assign the capabilties metrics
            # metrics['max_line_pull_t'] = None if None > metrics.get('max_line_pull_t') else metrics.get('max_line_pull_t')
            # metrics['brake_load_t'] = None if None > metrics.get('brake_load_t') else metrics.get('brake_load_t')
            # metrics['speed_mpm'] = None if None > metrics.get('speed_mpm') else metrics.get('speed_mpm')

        elif cap == 'bollard_pull':
            # per calwave install report (section 7.2): bollard pull can be described as function of vessel speed and load
            
            # logic for bollard_pull capability
            if objType == 'mooring':
                pass
            elif objType == 'platform':
                pass
            elif objType == 'anchor':
                pass
            elif objType == 'component':
                pass
            elif objType == 'turbine':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['max_force_t'] = None if None > metrics.get('max_force_t') else metrics.get('max_force_t')

        elif cap == 'crane':
            # logic for deck_space capability (all compatible)
            # object logic checked
            if objType == 'mooring':
                pass
            elif objType == 'platform':
                pass
            elif objType == 'anchor':
                pass
            elif objType == 'component':
                pass
            elif objType == 'turbine':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['capacity_t'] = None if None > metrics.get('capacity_t') else metrics.get('capacity_t')
            # metrics['hook_height_m'] = None if None > metrics.get('hook_height_m') else metrics.get('hook_height_m')

        elif cap == 'station_keeping':
            # logic for station_keeping capability
            if objType == 'mooring':
                pass
            elif objType == 'platform':
                pass
            elif objType == 'anchor':
                pass
            elif objType == 'component':
                pass
            elif objType == 'turbine':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['type'] = None if None > metrics.get('type') else metrics.get('type')

        elif cap == 'mooring_work':
            # logic for mooring_work capability
            if objType == 'mooring':
                pass
            elif objType == 'platform':
                pass
            elif objType == 'anchor':
                pass
            elif objType == 'component':
                pass
            elif objType == 'turbine':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['line_types'] = None if None > metrics.get('line_types') else metrics.get('line_types')
            # metrics['stern_roller'] = None if None > metrics.get('stern_roller') else metrics.get('stern_roller')
            # metrics['shark_jaws'] = None if None > metrics.get('shark_jaws') else metrics.get('shark_jaws')
            # metrics['towing_pin_rating_t'] = None if None > metrics.get('towing_pin_rating_t') else metrics.get('towing_pin_rating_t')

        elif cap == 'pump_surface':
            # logic for pump_surface capability
            if objType == 'mooring':
                pass
            elif objType == 'platform':
                pass
            elif objType == 'anchor':
                pass
            elif objType == 'component':
                pass
            elif objType == 'turbine':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['power_kW'] = None if None > metrics.get('power_kW') else metrics.get('power_kW')
            # metrics['pressure_bar'] = None if None > metrics.get('pressure_bar') else metrics.get('pressure_bar')
            # metrics['weight_t'] = None if None > metrics.get('weight_t') else metrics.get('weight_t')
            # metrics['dimensions_m'] = None if None > metrics.get('dimensions_m') else metrics.get('dimensions_m')

        elif cap == 'pump_subsea':
            # logic for pump_subsea capability
            if objType == 'mooring':
                pass
            elif objType == 'platform':
                pass
            elif objType == 'anchor':
                pass
            elif objType == 'component':
                pass
            elif objType == 'turbine':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['power_kW'] = None if None > metrics.get('power_kW') else metrics.get('power_kW')
            # metrics['pressure_bar'] = None if None > metrics.get('pressure_bar') else metrics.get('pressure_bar')
            # metrics['weight_t'] = None if None > metrics.get('weight_t') else metrics.get('weight_t')
            # metrics['dimensions_m'] = None if None > metrics.get('dimensions_m') else metrics.get('dimensions_m')

        elif cap == 'pump_grout':
            # logic for pump_grout capability
            if objType == 'mooring':
                pass
            elif objType == 'platform':
                pass
            elif objType == 'anchor':
                pass
            elif objType == 'component':
                pass
            elif objType == 'turbine':
                pass
            elif objType == 'cable':
                pass
            else:
                pass
            
            # Assign the capabilties metrics
            # metrics['power_kW'] = None if None > metrics.get('power_kW') else metrics.get('power_kW')
            # metrics['flow_rate_m3hr'] = None if None > metrics.get('flow_rate_m3hr') else metrics.get('flow_rate_m3hr')
            # metrics['pressure_bar'] = None if None > metrics.get('pressure_bar') else metrics.get('pressure_bar')
            # metrics['weight_t'] = None if None > metrics.get('weight_t') else metrics.get('weight_t')
            # metrics['dimensions_m'] = None if None > metrics.get('dimensions_m') else metrics.get('dimensions_m')

        elif cap == 'hydraulic_hammer':
            # logic for hydraulic_hammer capability (only platform and anchor objects compatible)
            # object logic checked
            if objType == 'platform':
                pass
            elif objType == 'anchor': # for fixed bottom installations
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['power_kW'] = None if None > metrics.get('power_kW') else metrics.get('power_kW')
            # metrics['energy_per_blow_kJ'] = None if None > metrics.get('energy_per_blow_kJ') else metrics.get('energy_per_blow_kJ')
            # metrics['weight_t'] = None if None > metrics.get('weight_t') else metrics.get('weight_t')
            # metrics['dimensions_m'] = None if None > metrics.get('dimensions_m') else metrics.get('dimensions_m')

        elif cap == 'vibro_hammer':
            # logic for vibro_hammer capability (only platform and anchor objects compatible)
            # object logic checked
            if objType == 'platform':
                pass
            elif objType == 'anchor': # for fixed bottom installations
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['power_kW'] = None if None > metrics.get('power_kW') else metrics.get('power_kW')
            # metrics['centrifugal_force_kN'] = None if None > metrics.get('centrifugal_force_kN') else metrics.get('centrifugal_force_kN')
            # metrics['weight_t'] = None if None > metrics.get('weight_t') else metrics.get('weight_t')
            # metrics['dimensions_m'] = None if None > metrics.get('dimensions_m') else metrics.get('dimensions_m')

        elif cap == 'drilling_machine':
            # logic for drilling_machine capability (only platform, anchor, and cable objects compatible)
            # Considering drilling both for export cables, interarray, and anchor/fixed platform install
            # object logic checked
            if objType == 'platform':
                pass
            elif objType == 'anchor':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['power_kW'] = None if None > metrics.get('power_kW') else metrics.get('power_kW')
            # metrics['weight_t'] = None if None > metrics.get('weight_t') else metrics.get('weight_t')
            # metrics['dimensions_m'] = None if None > metrics.get('dimensions_m') else metrics.get('dimensions_m')

        elif cap == 'torque_machine':
            # logic for torque_machine capability
            if objType == 'mooring':
                pass
            elif objType == 'platform':
                pass
            elif objType == 'anchor':
                pass
            elif objType == 'component':
                pass
            elif objType == 'turbine':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['power_kW'] = None if None > metrics.get('power_kW') else metrics.get('power_kW')
            # metrics['torque_kNm'] = None if None > metrics.get('torque_kNm') else metrics.get('torque_kNm')
            # metrics['weight_t'] = None if None > metrics.get('weight_t') else metrics.get('weight_t')
            # metrics['dimensions_m'] = None if None > metrics.get('dimensions_m') else metrics.get('dimensions_m')

        elif cap == 'cable_plough':
            # logic for cable_plough capability (only cable objects compatible)
            # object logic checked
            if  objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['power_kW'] = None if None > metrics.get('power_kW') else metrics.get('power_kW')
            # metrics['weight_t'] = None if None > metrics.get('weight_t') else metrics.get('weight_t')
            # metrics['dimensions_m'] = None if None > metrics.get('dimensions_m') else metrics.get('dimensions_m')

        elif cap == 'rock_placement':
            # logic for rock_placement capability (only platform, anchor, and cable objects compatible)
            # object logic checked
            if objType == 'platform':
                pass
            elif objType == 'anchor':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['placement_method'] = None if None > metrics.get('placement_method') else metrics.get('placement_method')
            # metrics['max_depth_m'] = None if None > metrics.get('max_depth_m') else metrics.get('max_depth_m')
            # metrics['accuracy_m'] = None if None > metrics.get('accuracy_m') else metrics.get('accuracy_m')
            # metrics['rock_size_range_mm'] = None if None > metrics.get('rock_size_range_mm') else metrics.get('rock_size_range_mm')

        elif cap == 'container':
            # logic for container capability (only platform, turbine, and cable objects compatible)
            # object logic checked
            if objType == 'wec':
                pass
            elif objType == 'turbine':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['weight_t'] = None if None > metrics.get('weight_t') else metrics.get('weight_t')
            # metrics['dimensions_m'] = None if None > metrics.get('dimensions_m') else metrics.get('dimensions_m')

        elif cap == 'rov':
            # logic for rov capability (all compatible)
            # object logic checked
            if objType == 'mooring':
                pass
            elif objType == 'platform':
                pass
            elif objType == 'anchor':
                pass
            elif objType == 'component':
                pass
            elif objType == 'turbine':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['class'] = None if None > metrics.get('class') else metrics.get('class')
            # metrics['depth_rating_m'] = None if None > metrics.get('depth_rating_m') else metrics.get('depth_rating_m')
            # metrics['weight_t'] = None if None > metrics.get('weight_t') else metrics.get('weight_t')
            # metrics['dimensions_m'] = None if None > metrics.get('dimensions_m') else metrics.get('dimensions_m')

        elif cap == 'positioning_system':
            # logic for positioning_system capability (only platform, anchor, and cable objects compatible)
            # object logic checked
            if objType == 'platform':
                pass
            elif objType == 'anchor':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['accuracy_m'] = None if None > metrics.get('accuracy_m') else metrics.get('accuracy_m')
            # metrics['methods'] = None if None > metrics.get('methods') else metrics.get('methods')

        elif cap == 'monitoring_system':
            # logic for monitoring_system capability
            if objType == 'mooring':
                pass
            elif objType == 'platform':
                pass
            elif objType == 'anchor':
                pass
            elif objType == 'component':
                pass
            elif objType == 'turbine':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['metrics'] = None if None > metrics.get('metrics') else metrics.get('metrics')
            # metrics['sampling_rate_hz'] = None if None > metrics.get('sampling_rate_hz') else metrics.get('sampling_rate_hz')

        elif cap == 'sonar_survey':
            # logic for sonar_survey capability (only anchor and cable objects compatible)
            # object logic checked
            if objType == 'anchor':
                pass
            elif objType == 'cable':
                pass
            else:
                pass

            # Assign the capabilties metrics
            # metrics['types'] = None if None > metrics.get('types') else metrics.get('types')
            # metrics['resolution_m'] = None if None > metrics.get('resolution_m') else metrics.get('resolution_m')

        else:
            raise Exception(f"Unsupported capability '{cap}'.")
        
        for met in metrics.keys():
            if metrics[met] == -1:
                print(f"WARNING: No metrics assigned for '{met}' metric in '{cap}' capability based on object type '{objType}'.")
        

        return metrics # return the dict of metrics and required values for the capability


    def assignObjects(self, objects):
        '''
        Adds a list of objects to the actions objects list and 
        calculates the required capability metrics, checking objects 
        are valid for the actions supported objects.

        The minimum capability metrics are used by when checking for 
        compatibility and assinging assets to the action in `assignAsset`.
        Thus this function should only be called in the intialization 
        process of an action.

        Inputs
        ------
        `objects` : `list`
            A list of FAModel objects to be added to the action.

        Returns
        -------
        `None`
        ''' 

        for obj in objects:

            # Check compatibility, set capability metrics based on object, and assign object to action

            objType = obj.__class__.__name__.lower() # object class name
            if objType not in self.supported_objects:
                raise Exception(f"Object type '{objType}' is not in the action's supported list.")
            else:
                if obj in self.objectList:
                    print(f"Warning: Object '{obj}' is already in the action's object list. Capabilities will be overwritten.")
                '''
                # Set capability requirements based on object
                for role, caplist in self.requirements.items():
                    for cap in caplist:
                        metrics = self.getMetrics(cap, caplist[cap], obj) # pass in the metrics dict for the cap and the obj

                        self.requirements[role][cap] = metrics  # assign metric of capability cap based on value required by obj
                # MH: commenting our for now just so the code will run, but it may be better to make the above a separate step anyway
                '''
                self.objectList.append(obj)

    
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
            
        # Initialize cost and duration
        self.cost = 0.0 # [$]
        self.duration = 0.0 # [h]
        
        """
        Note to devs:
        The code here calculates the cost and duration of an action. Each action in the actions.yaml has a hardcoded 'model' 
        here that is used to evaluate the action based on the assets assigned to it. 

        This is where a majority of assumptions about the action's behavior are made, so it is key to cite references behind
        any abnormal approaches. 

        Some good preliminary work on this is in https://github.com/FloatingArrayDesign/FAModel/blob/IOandM_development/famodel/installation/
        and in assets.py
        """
        
        # --- Mobilization ---
        if self.type == 'mobilize':
            # Hard-coded example of mobilization times based on vessel type
            durations = {
                'crane_barge': 3.0,
                'research_vessel': 1.0
            }
            for role_name, vessel in self.assets.items():
                vessel_type = vessel['type'].lower()
                for key, duration in durations.items():
                    if key in vessel_type:
                        self.duration += duration
                        break

        elif self.type == 'demobilize':
            # Hard-coded example of demobilization times based on vessel type
            durations = {
                'crane_barge': 3.0,
                'research_vessel': 1.0
            }
            for role_name, vessel in self.assets.items():
                vessel_type = vessel['type'].lower()
                for key, duration in durations.items():
                    if key in vessel_type:
                        self.duration += duration
        elif self.type == 'load_cargo':
            pass

        # --- Towing & Transport ---
        elif self.type == 'tow':
            pass
        
        elif self.type == 'transit_linehaul_self':
            # YAML override
            try:
                v = getFromDict(self.actionType, 'duration_h', dtype=float); self.duration += v
            except ValueError:
                try:
                    v = getFromDict(self.actionType, 'default_duration_h', dtype=float); self.duration += v
                except ValueError:
                    vessel = self.assets.get('vessel') or self.assets.get('operator') or self.assets.get('carrier')
                    if vessel is None:
                        raise ValueError('transit_linehaul_self: no vessel assigned.')
        
                    tr = vessel['transport']
                    
                    # distance
                    dist_m = float(tr['route_length_m'])
                    
                    # speed: linehaul uses transport.cruise_speed_mps
                    speed_mps = float(tr['cruise_speed_mps'])

                    dur_h = dist_m/speed_mps/3600.0
                    self.duration += dur_h
            # cost
            rate_per_hour = 0.0
            for _, asset in self.assets.items():
                rate_per_hour += float(asset['day_rate'])/24.0
            self.cost += self.duration*rate_per_hour
            return self.duration, self.cost


        elif self.type == 'transit_linehaul_tug':
            # YAML override
            try:
                v = getFromDict(self.actionType, 'duration_h', dtype=float); self.duration += v
            except ValueError:
                try:
                    v = getFromDict(self.actionType, 'default_duration_h', dtype=float); self.duration += v
                except ValueError:
                    tug   = self.assets.get('operator') or self.assets.get('vessel')
                    barge = self.assets.get('carrier')
                    if tug is None or barge is None:
                        raise ValueError('transit_linehaul_tug: need tug (operator) and barge (carrier).')
        
                    tr_b = barge.get('transport', {})
                    tr_t = tug.get('transport', {})
                    
                    # distance: prefer barge’s transport
                    dist_m = float(tr_b.get('route_length_m', tr_t['route_length_m']))
                    
                    # speed for convoy linehaul: barge (operator) cruise speed
                    operator = self.assets.get('operator') or self.assets.get('vessel')
                    if operator is None:
                        raise ValueError('transit_linehaul_tug: operator (barge) missing.')
                    
                    speed_mps = float(operator['transport']['cruise_speed_mps'])

                    dur_h = dist_m/speed_mps/3600.0
                    self.duration += dur_h

            # cost
            rate_per_hour = 0.0
            for _, asset in self.assets.items():
                rate_per_hour += float(asset['day_rate'])/24.0
            self.cost += self.duration*rate_per_hour
            return self.duration, self.cost

        elif self.type == 'transit_onsite_self':
            # YAML override
            try:
                v = getFromDict(self.actionType, 'duration_h', dtype=float); self.duration += v
            except ValueError:
                try:
                    v = getFromDict(self.actionType, 'default_duration_h', dtype=float); self.duration += v
                except ValueError:
                    # vessel (Beyster) required
                    vessel = self.assets.get('vessel') or self.assets.get('operator') or self.assets.get('carrier')
                    if vessel is None:
                        raise ValueError('transit_onsite_self: no vessel assigned.')
        
                    # NEW: quick vessel print
                    try:
                        print(f"[onsite_self] {self.name}: vessel={vessel.get('type')}")
                    except Exception:
                        pass
        
                    # destination anchor from objects (required)
                    if not self.objectList:
                        raise ValueError('transit_onsite_self: destination anchor missing in objects.')
                    dest = self.objectList[0]
                    r_dest = getattr(dest, 'r', None)
        
                    # NEW: print dest
                    try:
                        print(f"[onsite_self] {self.name}: r_dest={r_dest}")
                    except Exception:
                        pass
        
                    # infer start from dependency chain (BFS up to depth 3)
                    r_start = None
                    from collections import deque
                    q, seen = deque(), set()
                    for dep in self.dependencies.values():
                        q.append((dep, 0)); seen.add(id(dep))
                    while q:
                        node, depth = q.popleft()
                        if node.objectList and hasattr(node.objectList[0], 'r'):
                            r_start = node.objectList[0].r
                            break
                        # if depth < 3:
                        #     for nxt in node.dependencies.values():
                        #         if id(nxt) in seen: continue
                        #         seen.add(id(nxt)); q.append((nxt, depth+1))
        
                    # NEW: print BFS result
                    try:
                        print(f"[onsite_self] {self.name}: r_start(BFS)={r_start}")
                    except Exception:
                        pass
        
                    # CHANGED: fallback for first onsite leg → try centroid, else keep old zero-distance fallback
                    if r_start is None and r_dest is not None:
                        # NEW: centroid read (linehaul_to_site should set it on this action)
                        cent = (getattr(self, 'meta', {}) or {}).get('anchor_centroid')
                        if cent is None:
                            cent = (getattr(self, 'params', {}) or {}).get('anchor_centroid')
                        if cent is not None and len(cent) >= 2:
                            r_start = (float(cent[0]), float(cent[1]))
                            try:
                                print(f"[onsite_self] {self.name}: using centroid as r_start={r_start}")
                            except Exception:
                                pass
                        else:
                            # ORIGINAL behavior: assume zero in-field distance
                            r_start = r_dest
                            try:
                                print(f"[warn] {self.name}: could not infer start from deps; assuming zero in-field distance.")
                            except Exception:
                                pass
        
                    # 2D distance [m]
                    from math import hypot
                    dx = float(r_dest[0]) - float(r_start[0])
                    dy = float(r_dest[1]) - float(r_start[1])
                    dist_m = hypot(dx, dy)
        
                    # NEW: print distance
                    try:
                        print(f"[onsite_self] {self.name}: dist_m={dist_m:.1f} (start={r_start} → dest={r_dest})")
                    except Exception:
                        pass
        
                    # onsite speed from capabilities.engine (SI)
                    cap_eng = vessel.get('capabilities', {}).get('engine', {})
                    speed_mps = float(cap_eng['site_speed_mps'])

                    self.duration += dist_m/speed_mps/3600.0
        
                    # NEW: print duration increment
                    try:
                        print(f"[onsite_self] {self.name}: speed_mps={speed_mps:.3f}, dT_h={dist_m/speed_mps/3600.0:.3f}, total={self.duration:.3f}")
                    except Exception:
                        pass
        
            # cost
            rate_per_hour = 0.0
            for _, asset in self.assets.items():
                rate_per_hour += float(asset['day_rate'])/24.0
            self.cost += self.duration*rate_per_hour
            return self.duration, self.cost

        elif self.type == 'transit_onsite_tug':
            # YAML override
            try:
                v = getFromDict(self.actionType, 'duration_h', dtype=float); self.duration += v
            except ValueError:
                try:
                    v = getFromDict(self.actionType, 'default_duration_h', dtype=float); self.duration += v
                except ValueError:
                    # assets required (operator = San_Diego tug; carrier = Jag barge)
                    operator = self.assets.get('operator') or self.assets.get('vessel')
                    carrier  = self.assets.get('carrier')
                    if operator is None and carrier is None:
                        raise ValueError('transit_onsite_tug: no operator/carrier assigned.')
        
                    # quick prints
                    try:
                        op_name = operator.get('type') if operator else None
                        ca_name = carrier.get('type')  if carrier  else None
                        print(f"[onsite_tug] {self.name}: operator={op_name} carrier={ca_name}")
                    except Exception:
                        pass
        
                    # destination anchor from objects (required)
                    if not self.objectList:
                        raise ValueError('transit_onsite_tug: destination anchor missing in objects.')
                    dest = self.objectList[0]
                    r_dest = getattr(dest, 'r', None)
        
                    try:
                        print(f"[onsite_tug] {self.name}: r_dest={r_dest}")
                    except Exception:
                        pass
        
                    # infer start from dependency chain (BFS up to depth 3)
                    r_start = None
                    from collections import deque
                    q, seen = deque(), set()
                    for dep in self.dependencies.values():
                        q.append((dep, 0)); seen.add(id(dep))
                    while q:
                        node, depth = q.popleft()
                        if node.objectList and hasattr(node.objectList[0], 'r'):
                            r_start = node.objectList[0].r
                            break
                        # if depth < 3:
                        #     for nxt in node.dependencies.values():
                        #         if id(nxt) in seen: continue
                        #         seen.add(id(nxt)); q.append((nxt, depth+1))
        
                    try:
                        print(f"[onsite_tug] {self.name}: r_start(BFS)={r_start}")
                    except Exception:
                        pass
        
                    # fallback for first onsite leg: use centroid if present, else zero-distance fallback
                    if r_start is None and r_dest is not None:
                        cent = (getattr(self, 'meta', {}) or {}).get('anchor_centroid')
                        if cent is None:
                            cent = (getattr(self, 'params', {}) or {}).get('anchor_centroid')
                        if cent is not None and len(cent) >= 2:
                            r_start = (float(cent[0]), float(cent[1]))
                            try:
                                print(f"[onsite_tug] {self.name}: using centroid as r_start={r_start}")
                            except Exception:
                                pass
                        else:
                            r_start = r_dest
                            try:
                                print(f"[warn] {self.name}: could not infer start from deps; assuming zero in-field distance.")
                            except Exception:
                                pass
        
                    # 2D distance [m]
                    from math import hypot
                    dx = float(r_dest[0]) - float(r_start[0])
                    dy = float(r_dest[1]) - float(r_start[1])
                    dist_m = hypot(dx, dy)
        
                    try:
                        print(f"[onsite_tug] {self.name}: dist_m={dist_m:.1f} (start={r_start} → dest={r_dest})")
                    except Exception:
                        pass
        
                    # speed for convoy onsite: barge (operator) site speed
                    operator = self.assets.get('operator') or self.assets.get('vessel')
                    if operator is None:
                        raise ValueError('transit_onsite_tug: operator (barge) missing.')
                    
                    cap_eng = operator.get('capabilities', {}).get('bollard_pull', {})
                    speed_mps = float(cap_eng['site_speed_mps'])

                    self.duration += dist_m/speed_mps/3600.0
        
                    try:
                        print(f"[onsite_tug] {self.name}: speed_mps={speed_mps:.3f}, dT_h={dist_m/speed_mps/3600.0:.3f}, total={self.duration:.3f}")
                    except Exception:
                        pass
        
            # cost (unchanged)
            rate_per_hour = 0.0
            for _, asset in self.assets.items():
                rate_per_hour += float(asset['day_rate'])/24.0
            self.cost += self.duration*rate_per_hour
            return self.duration, self.cost

        elif self.type == 'at_site_support':
            pass
        elif self.type == 'transport_components':
            pass

        # --- Mooring & Anchors ---
        elif self.type == 'install_anchor':
            # YAML override (no model if present)
            default_duration = None
            try:
                default_duration = getFromDict(self.actionType, 'duration_h', dtype=float)
            except ValueError:
                default_duration = None
        
            if default_duration is not None:
                computed_duration_h = default_duration
       
            else:
                # Expect an anchor object in self.objectList
                if not self.objectList:
                    raise ValueError("install_anchor: no anchor object provided in 'objects'.")
                                
                # 1) Relevant metrics for cost and duration
                anchor = self.objectList[0]
                L = anchor.dd['design']['L']
                depth_m = abs(float(anchor.r[2]))
                
                # 2) Winch vertical speed [mps]
                v_mpm = float(self.assets['carrier']['capabilities']['winch']['speed_mpm'])
                t_lower_min = depth_m/v_mpm
                
                # 3) Penetration time ~ proportional to L 
                rate_pen = 15. # [min] per [m]
                t_pen_min = L*rate_pen
                
                # 4) Connection / release (fixed)
                t_ops_min = 15
                
                duration_min = t_lower_min + t_pen_min + t_ops_min
                computed_duration_h = duration_min/60.0 # [h]
                
            # print(f'[install_anchor] yaml_duration={yaml_duration} -> used={computed_duration_h} h')
            
            # Duration addition
            self.duration += computed_duration_h 
            
            # Cost assessment
            rate_per_hour = 0.0
            for _, asset in self.assets.items():
                rate_per_hour += float(asset['day_rate'])/24.0
            
            self.cost += self.duration*rate_per_hour

        elif self.type == 'retrieve_anchor':
            pass       
        elif self.type == 'install_mooring':
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
        elif self.type == 'install_turbine':
            pass

        # --- Cable Operations ---
        elif self.type == 'lay_cable':
            pass
        elif self.type == 'cable_hookup':
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
            # 1) YAML override first
            try:
                v = getFromDict(self.actionType, 'duration_h', dtype=float); self.duration += v
            except ValueError:
                try:
                    v = getFromDict(self.actionType, 'default_duration_h', dtype=float); self.duration += v
                except ValueError:
                    # --- find the paired install ---
                    ref_install = getattr(self, 'paired_install', None)
        
                    # fallback: BFS through deps to find an install on the same anchor
                    if ref_install is None:
                        anchor_obj = self.objectList[0] if self.objectList else None
                        from collections import deque
                        q, seen = deque(), set()
                        for dep in self.dependencies.values():
                            q.append((dep, 0)); seen.add(id(dep))
                        while q:
                            node, depth = q.popleft()
                            if getattr(node, 'type', None) == 'install_anchor':
                                if anchor_obj and node.objectList and node.objectList[0] is anchor_obj:
                                    ref_install = node
                                    break
                                if ref_install is None:
                                    ref_install = node
                            if depth < 3:
                                for nxt in node.dependencies.values():
                                    if id(nxt) in seen: continue
                                    seen.add(id(nxt)); q.append((nxt, depth+1))
        
                    # --- get install duration, compute-on-demand if needed (no side effects) ---
                    inst_dur = 0.0
                    if ref_install is not None:
                        inst_dur = float(getattr(ref_install, 'duration', 0.0) or 0.0)
        
                        # if not computed yet, safely compute and restore
                        if inst_dur <= 0.0 and not getattr(ref_install, '_in_monitor_pull', False):
                            try:
                                ref_install._in_monitor_pull = True  # guard re-entrancy
                                prev_cost = ref_install.cost
                                prev_dur  = ref_install.duration
                                d, _ = ref_install.calcDurationAndCost()
                                inst_dur = float(d) if d is not None else 0.0
                                # restore to avoid double counting later
                                ref_install.cost = prev_cost
                                ref_install.duration = prev_dur
                            finally:
                                ref_install._in_monitor_pull = False
        
                    self.duration += inst_dur
        
            # cost (same pattern you use elsewhere)
            rate_per_hour = 0.0
            for _, asset in self.assets.items():
                rate_per_hour += float(asset['day_rate'])/24.0
            self.cost += self.duration * rate_per_hour
            return self.duration, self.cost

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
        

