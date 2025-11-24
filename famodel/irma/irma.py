"""Core code for setting up a IO&M scenario"""

import os
import numpy as np
import matplotlib.pyplot as plt

import moorpy as mp
from moorpy.helpers import set_axes_equal
from moorpy import helpers
import yaml
from copy import deepcopy
import string
try: 
    import raft as RAFT
except:
    pass

#from shapely.geometry import Point, Polygon, LineString
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

import networkx as nx
from action import Action, increment_name
from task import Task

from assets import Vessel, Port



def loadYAMLtoDict(info, already_dict=False):
    '''Reads a list or YAML file and prepares a dictionary'''
    
    if isinstance(info, str):
    
        with open(info) as file:
            data = yaml.load(file, Loader=yaml.FullLoader)
            if not data:
                raise Exception(f'File {info} does not exist or cannot be read. Please check filename.')
    elif isinstance(info, list):
        data = info
    else:
        raise Exception('loadYAMLtoDict must be passed a filename or list')
    
    # Go through contents and product the dictionary
    info_dict = {}
    
    if already_dict:
        # assuming it's already a dict 
        info_dict.update(data)
        
    else: # a list of dicts with name parameters
        # So we will convert into a dict based on those names
        for entry in data:
            if not 'name' in entry:
                print(entry)
                raise Exception('This entry does not have a required name field.')
                
            if entry['name'] in info_dict:
                print(entry)
                raise Exception('This entry has the same name as an existing entry.')
                
            info_dict[entry['name']] = entry  # could make this a copy operation if worried
    
    return info_dict
    

#def storeState(project,...):


#def applyState():


def unifyUnits(d):
    '''Converts any capability specification/metric in supported non-SI units
    to be in SI units. Converts the key names as well.'''
    
    # >>> not working yet <<<
    
    
    # load conversion data from YAML (eventually may want to store this in a class)
    with open('spec_conversions.yaml') as file:
        data = yaml.load(file, Loader=yaml.FullLoader)
    
    keys1 = []
    facts = []  # conversion factors
    keys2 = []
    
    for line in data:
        keys1.append(line[0])
        facts.append(line[1])
        keys2.append(line[2])
        
    # >>> dcopy = deepcopy(d)
    
    for asset in d.values():  # loop through each asset's dict
        
        capabilities = {}  # new dict of capabilities to built up
        
        for cap_key, cap_val in asset['capabilities'].items():
            
            # make the capability type sub-dictionary
            capabilities[cap_key] = {}
        
            for key, val in cap_val.items():  # look at each capability metric
                try:
                    i = keys1.index(key)  # find if key is on the list to convert
                    
                    
                    if keys2[i] in cap_val.keys():
                        raise Exception(f"Specification '{keys2[i]}' already exists")
                    
                    capabilities[cap_key][keys2[i]] = val * facts[i]  # make converted entry
                    #capability[keys2[i]] = val * facts[i]  # create a new SI entry
                    #del capability[keys1[i]]  # remove the original?
                    
                except:
                    
                    capabilities[cap_key][key] = val  # copy over original form
                 

class Scenario():

    def __init__(self):
        '''Initialize a scenario object that can be used for IO&M modeling of
        of an offshore energy system. Eventually it will accept user-specified 
        settings files.
        '''
        
        # ----- Load database of supported things -----
        
        actionTypes = loadYAMLtoDict('actions.yaml', already_dict=True)  # Descriptions of actions that can be done
        requirements = loadYAMLtoDict('requirements.yaml', already_dict=True)  # Descriptions of requirements that can be done
        capabilities = loadYAMLtoDict('capabilities.yaml', already_dict=True)
        vessels = loadYAMLtoDict('vessels.yaml', already_dict=True)
        objects = loadYAMLtoDict('objects.yaml', already_dict=True)
        
        unifyUnits(vessels)  # (function doesn't work yet!) <<<
        
        # ----- Validate internal cross references -----
        
        # Make sure vessels don't use nonexistent capabilities or actions
        for key, ves in vessels.items():
            
            #if key != ves['name']:
            #    raise Exception(f"Vessel key ({key}) contradicts its name ({ves['name']})")
            
            # Check capabilities
            if not 'capabilities' in ves:
                raise Exception(f"Vessel '{key}' is missing a capabilities list.")
                
            for capname, cap in ves['capabilities'].items():
                if not capname in capabilities:
                    raise Exception(f"Vessel '{key}' capability '{capname}' is not in the global capability list.")
                
                # Could also check the sub-parameters of the capability
                for cap_param in cap:
                    if not cap_param in capabilities[capname]:
                        raise Exception(f"Vessel '{key}' capability '{capname}' parameter '{cap_param}' is not in the global capability's parameter list.")
            
            # Check actions
            if not 'actions' in ves:
                raise Exception(f"Vessel '{key}' is missing an actions list.")
            
            for act in ves['actions']:
                if not act in actionTypes:
                    raise Exception(f"Vessel '{key}' action '{act}' is not in the global action list.")
        
        
        # Make sure actions refer to supported object types/properties and capabilities
        for key, act in actionTypes.items():
            
            act['type'] = key
            
            #if key != act['name']:
            #    raise Exception(f"Action key ({key}) contradicts its name ({act['name']})")
            
            # Check capabilities
            #if 'capabilities' in act:
            #    raise Exception(f"Action '{key}' is missing a capabilities list.")
            
            if 'capabilities' in act:
                
              for cap in act['capabilities']:
                if not cap in capabilities:
                    raise Exception(f"Action '{key}' capability '{cap}' is not in the global capability list.")
                
                # Could also check the sub-parameters of the capability
                #for cap_param in cap:
                #    if not cap_param in capabilities[cap['name']]:
                #        raise Exception(f"Action '{key}' capability '{cap['name']}' parameter '{cap_param}' is not in the global capability's parameter list.")
            
            if 'roles' in act:   # look through capabilities listed under each role
                for caps in act['roles'].values():
                    for cap in caps:
                        if not cap in capabilities:
                            raise Exception(f"Action '{key}' capability '{cap}' is not in the global capability list.")
                        
            
            # Check objects
            if not 'objects' in act:
                raise Exception(f"Action '{key}' is missing an objects list.")
            
            for obj in act['objects']:
                if not obj in objects:
                    raise Exception(f"Action '{key}' object '{obj}' is not in the global objects list.")
        
                # Could also check the sub-parameters of the object
                if isinstance(act['objects'], dict):  # if the object
                    for obj_param in act['objects'][obj]:
                        if not obj_param in objects[obj]:
                            raise Exception(f"Action '{key}' object '{obj}' parameter '{obj_param}' is not in the global object's parameter list.")
            
        
        # Store some things
        self.actionTypes  = actionTypes
        
        self.requirements  = requirements
        self.capabilities = capabilities  
        self.vessels      = vessels
        self.objects      = objects
        
        
        # Initialize some things
        self.actions = {}
        self.tasks = {}
        
        
    def registerAction(self, action):
        '''Registers an already created action'''
        
        # this also handles creation of unique dictionary keys
        
        if action.name in self.actions:  # check if there is already a key with the same name
            raise Warning(f"Action '{action.name}' is already registered.")
            print(f"Action name '{action.name}' is in the actions list so incrementing it...")
            action.name = increment_name(action.name)

        # What about handling of dependencies?? <<< done in the action object, 
        # but could check that each is in the list already...
        for dep in action.dependencies.values():
            if not dep in self.actions.values():
                raise Exception(f"New action '{action.name}' has a dependency '{dep.name}' this is not in the action list.")        
        
        # Check that all the requirements of all actions conform to the
        # options in requirements.yaml.
        for reqname, req in action.requirements.items():
            if req['base'] in self.requirements:  # ensure this requirement is listed
                for cap in req['capabilities']:
                    if not cap in self.capabilities:
                        raise Exception(f"Requirement '{reqname}' capability '{cap}' is not in the global capability list.")
            else:
                raise Exception(f"Action {action.name} requirement {req['base']} is not in requirements.yaml")

        # Add it to the actions dictionary
        self.actions[action.name] = action
        
        
    def addAction(self, action_type_name, action_name, **kwargs):
        '''Creates an action and adds it to the register'''
        
        if not action_type_name in self.actionTypes:
            raise Exception(f"Specified action type name {'action_type_name'} is not in the list of loaded action types.")
        
        # Get dictionary of action type information
        action_type = self.actionTypes[action_type_name]
        
        # Initialize full zero-valued dictionary of possible required capability specs
        reqs = {}  # Start a dictionary to hold the requirements -> capabilities -> specs
        for req in action_type['requirements']:
            
            # make sure it's a valid requirement
            if '-in' in req:  
                req_dir = 1  # this means the req is for storage and storage is being filled
                req_base = req[:-3]  # this is the name of the req as in requirements.yaml, no suffix
            elif '-out' in req:  
                req_dir = -1
                req_base = req[:-4]
            else:
                req_dir = 0
                req_base = req
            
            # Make sure the requirement and its direction are supported in the requirements yaml
            if not req_base in self.requirements:
                raise Exception(f"Requirement '{req_base}' is not in the requirements yaml.")
            if abs(req_dir) > 0 and ('directions' not in self.requirements[req_base] 
                                     or req_dir not in self.requirements[req_base]['directions']):
                raise Exception(f"Requirement '{req_base}' direction '{req_dir}' is not supported in the requirements yaml.")
            
            # Make the new requirements entry
            reqs[req] = {'base':req_base, 'direction':req_dir, 'capabilities':{}}
            
            # add the caps of the req
            for cap in self.requirements[req_base]['capabilities']: 
                reqs[req]['capabilities'][cap] = {}
                #print(f'   {cap}')
                # add the specs of the capability
                for spec in self.capabilities[cap]:
                    reqs[req]['capabilities'][cap][spec] = 0
                    #print(f'     {spec} = 0')
        # swap in the filled-out dict
        action_type['requirements'] = reqs
        
        # Create the action
        act = Action(action_type, action_name, **kwargs)        
        
        # Register the action
        self.registerAction(act)
        
        return act   # return the newly created action object, or its name?

    
    def addActionDependencies(self, action, dependencies):  
        '''Adds dependencies to an action, provided those dependencies have
        already been registered in the action list.
        '''
        
        if not isinstance(dependencies, list):
            dependencies = [dependencies]  # get into list form if singular
                
        for dep in dependencies:
            # Make sure the dependency is already registered
            if dep in self.actions.values():  
                action.addDependency(dep)
            else:
                raise Exception(f"New action '{action.name}' has a dependency '{dep.name}' this is not in the action list.")                  


    def visualizeActions(self):
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
        if len(longest_path)>=1:
            last_node = longest_path[-1]  # Identify last node of the longest path
            # Define layout
            pos = nx.shell_layout(G)
            # Draw all nodes and edges (default gray)
            nx.draw(G, pos, with_labels=True, node_size=500, 
                    node_color='skyblue', font_size=10, font_weight='bold', 
                    font_color='black', edge_color='gray')

            # Highlight longest path in red
            nx.draw_networkx_edges(G, pos, edgelist=longest_path_edges, edge_color='red', width=2)

            # Annotate last node with total duration in red
            plt.text(pos[last_node][0], pos[last_node][1] - 0.1, f"{total_duration:.2f} hr", fontsize=12, color='red', fontweight='bold', ha='center')          
        else:
            pass
        plt.axis('equal')

        # Color first node (without dependencies) green
        i = 0
        for node in G.nodes():
            if G.in_degree(node) == 0:  # Check if the node has no incoming edges
                nx.draw_networkx_nodes(G, pos, nodelist=[node], node_color='green', node_size=500, label='Action starters' if i==0 else None)
                i += 1
        plt.legend()
        return G    
    
    
    def registerTask(self, task):
        '''Registers an already created task'''
        
        # this also handles creation of unique dictionary keys
        
        if task.name in self.tasks:  # check if there is already a key with the same name
            raise Warning(f"Action '{task.name}' is already registered.")
            print(f"Task name '{task.name}' is in the tasks list so incrementing it...")
            task.name = increment_name(task.name)

        # Add it to the actions dictionary
        self.tasks[task.name] = task
        
        
    def addTask(self, task_name, actions, action_sequence, **kwargs):
        '''Creates a task and adds it to the register'''
        
        # Create the action
        task = Task(task_name, actions, action_sequence, **kwargs)        
        
        # Register the action
        self.registerTask(task)
        
        return task
    
    
    
    def findCompatibleVessels(self):
        '''Go through actions and identify which vessels have the required
        capabilities (could be based on capability presence, or quantitative.
        '''
        
        pass


    def figureOutTaskRelationships(self):
        '''Calculate time constraints
        between tasks.
        '''
        
        # Figure out task durations (for a given set of asset assignments?)
        #for task in self.tasks.values():
            #task.calcTiming()
        
        # Figure out timing constraints between tasks based on action dependencies
        n = len(self.tasks)
        dt_min = np.zeros((n,n))  # matrix of required time offsets between tasks
        
        for i1, task1 in enumerate(self.tasks.values()):
            for i2, task2 in enumerate(self.tasks.values()):
                # look at all action dependencies from tasks 1 to 2 and
                # identify the limiting case (the largest time offset)...
                dt_min_1_2, dt_min_2_1 = findTaskDependencies(task1, task2)
                
                # for now, just look in one direction
                dt_min[i1, i2] = dt_min_1_2

        return dt_min
    

def findTaskDependencies(task1, task2):
    '''Finds any time dependency between the actions of two tasks.
    Returns the minimum time separation required from task 1 to task 2,
    and from task 2 to task 1. I
    '''
    
    time_1_to_2 = []
    time_2_to_1 = []
    
    # Look for any dependencies where act2 depends on act1:
    #for i1, act1 in enumerate(task1.actions.values()):
    #    for i2, act2 in enumerate(task2.actions.values()):
    for a1, act1 in task1.actions.items():
        for a2, act2 in task2.actions.items():
        
            if a1 in act2.dependencies:  # if act2 depends on act1
                time_1_to_2.append(task1.actions_ti[a1] + act1.duration
                                   - task2.actions_ti[a2])
        
            if a2 in act1.dependencies:  # if act2 depends on act1
                time_2_to_1.append(task2.actions_ti[a2] + act2.duration
                                   - task1.actions_ti[a1])
    
    print(time_1_to_2)
    print(time_2_to_1)
    
    # TODO: provide cleaner handling of whether or not there is a time constraint in either direction <<<
    
    dt_min_1_2 = min(time_1_to_2) if time_1_to_2 else -np.inf  # minimum time required from t1 start to t2 start
    dt_min_2_1 = min(time_2_to_1) if time_2_to_1 else -np.inf  # minimum time required from t2 start to t1 start
    
    if dt_min_1_2 + dt_min_2_1 > 0:
        print(f"The timing between these two tasks seems to be impossible...")
    
    #breakpoint()
    return dt_min_1_2, dt_min_2_1


def implementStrategy_staged(sc):
    '''This sets up Tasks in a way that implements a staged installation
    strategy where all of one thing is done before all of a next thing.
    '''
    
    # ----- Create a Task for all the anchor installs -----
    
    # gather the relevant actions
    acts = []
    for action in sc.actions.values():
        if action.type == 'install_anchor':
            acts.append(action)

    # create the task, passing in the sequence of actions
    sc.addTask('install_all_anchors', acts, action_sequence='series')
    
    
    # ----- Create a Task for all the mooring installs -----
    
    # gather the relevant actions
    acts = []
    # first load each mooring
    for action in sc.actions.values():
        if action.type == 'load_mooring':
            acts.append(action)
    # next lay each mooring (eventually route logic could be added)
    for action in sc.actions.values():
        if action.type == 'lay_mooring':
            acts.append(action)

    # create the task, passing in the sequence of actions
    sc.addTask('install_all_moorings', acts, action_sequence='series')
    
    
    # ----- Create a Task for the platform tow-out and hookup -----
    
    # gather the relevant actions
    acts = []
    # first tow out the platform
    acts.append(sc.actions['tow'])
    # next hook up each mooring
    for action in sc.actions.values():
        if action.type == 'mooring_hookup':
            acts.append(action)

    # create the task, passing in the sequence of actions
    sc.addTask('tow_and_hookup', acts, action_sequence='series')
    


if __name__ == '__main__':
    '''This is currently a script to explore how some of the workflow could
    work. Can move things into functions/methods as they solidify.
    '''
    
    # ----- Load up a Project -----    
    
    from famodel.project import Project


    #%% Section 1: Project without RAFT
    print('Creating project without RAFT\n')
    print(os.getcwd())
    # create project object
    # project = Project(file='C:/Code/FAModel/examples/OntologySample200m_1turb.yaml', raft=False) # for Windows
    project = Project(file='../../examples/OntologySample200m_1turb.yaml', raft=False) # for Mac
    # create moorpy system of the array, include cables in the system
    project.getMoorPyArray(cables=1)
    # plot in 3d, using moorpy system for the mooring and cable plots
    # project.plot2d()
    # project.plot3d()

    '''
    # project.arrayWatchCircle(ang_spacing=20)
    # save envelopes from watch circle information for each mooring line
    for moor in project.mooringList.values():
        moor.getEnvelope()

    # plot motion envelopes with 2d plot
    project.plot2d(save=True, plot_bathymetry=False)
    '''
    
    # Tally up some object properties (eventually make this built-in Project stuff)
    for mooring in project.mooringList.values():
        # sum up mooring quantities of interest
        L = 0  # length
        m = 0  # mass
        V = 0  # volume

        for sec in mooring.sections(): # add up the length of all sections in the mooring
            L += sec['L']
            m += sec['L'] * sec['type']['m']
            V += sec['L'] * np.pi/4 * sec['type']['d_vol']**2
        
        mooring.props = {}
        mooring.props['length'] = L
        mooring.props['pretension'] = 0   # <<< get this from MoorPy once this is moved into Mooring class?
        mooring.props['weight'] = 9.8*(m - 1025*V)
        mooring.props['mass'] = m
        mooring.props['volume'] = V
        
    print("should do the same for platforms and anchors...") # <<<
    
    sc = Scenario()  # class instance holding most of the info
                
    
    # ----- Create the interrelated actions (including their individual requirements) -----
    print("===== Create Actions =====")
    # When an action is created, its requirements will be calculated based on
    # the nature of the action and the objects involved.
    
    for akey, anchor in project.anchorList.items():
        
        ## Test action.py for anchor install   

        # add and register anchor install action(s)
        a1 = sc.addAction('install_anchor', f'install_anchor-{akey}', objects=[anchor])
        #duration, cost = a1.evaluateAssets({'carrier' : sc.vessels["MPSV_01"], 'operator':sc.vessels["AHTS_alpha"]})
        #print(f'Anchor install action {a1.name} duration: {duration:.2f} days, cost: ${cost:,.0f}')
        
        # register the actions as necessary for the anchor <<< do this for all objects??
        anchor.install_dependencies = [a1]
    
    
    hookups = []  # list of hookup actions
    
    for mkey, mooring in project.mooringList.items():
        
        # note origin and destination
        
        # --- lay out all the mooring's actions (and their links)

        ## Test action.py for mooring load

        # create load vessel action
        a2 = sc.addAction('load_mooring', f'load_mooring-{mkey}', objects=[mooring])
        #duration, cost = a2.evaluateAssets({'carrier2' : sc.vessels["HL_Giant"], 'carrier1' : sc.vessels["Barge_squid"], 'operator' : sc.vessels["HL_Giant"]})
        #print(f'Mooring load action {a2.name} duration: {duration:.2f} days, cost: ${cost:,.0f}')

        # create ship out mooring action
        
        # create lay mooring action
        a3 = sc.addAction('lay_mooring', f'lay_mooring-{mkey}', objects=[mooring], dependencies=[a2])
        sc.addActionDependencies(a3, mooring.attached_to[0].install_dependencies) # in case of shared anchor
        #print(f'Lay mooring action {a3.name} duration: {duration:.2f} days, cost: ${cost:,.0f}')
        
        # mooring could be attached to anchor here - or could be lowered with anchor!!
        #(r=r_anch, mooring=mooring, anchor=mooring.anchor...)
        # the action creator can record any dependencies related to actions of the anchor
        
        # create hookup action
        a4 = sc.addAction('mooring_hookup', f'mooring_hookup-{mkey}', 
                          objects=[mooring, mooring.attached_to[1]], dependencies=[a2, a3])
        #(r=r, mooring=mooring, platform=platform, depends_on=[a4])
        # the action creator can record any dependencies related to actions of the platform
        
        hookups.append(a4)
        
    
    # add the FOWT install action
    a5 = sc.addAction('tow', 'tow', objects=[list(project.platformList.values())[0]])
    for a in hookups:
        sc.addActionDependencies(a, [a5])  # make each hookup action dependent on the FOWT being towed out
    

    # ----- Do some graph analysis -----
    
    #G = sc.visualizeActions()
    
    # ----- Generate tasks (sequences of Actions following specific strategies) -----
    print('Generating tasks')
    # Call one of the task strategy implementers, which will create the tasks
    # (The created tasks also contain information about their summed requirements)
    implementStrategy_staged(sc)
    
    
    # ----- Try assigning assets to the tasks -----
    print('Trying to assign assets to tasks')
    for task in sc.tasks.values():
        print(f"--- Looking at task {task.name} ---")
        if task.checkAssets([sc.vessels['AHTS_alpha']], display=1)[0]:
            print('Assigned AHTS')
            task.assignAssets([sc.vessels['AHTS_alpha']])
        elif task.checkAssets([sc.vessels['CSV_A']], display=1)[0]:
            print('Assigned CSV_A')
            task.assignAssets([sc.vessels['CSV_A']])
        else:
            task.checkAssets([sc.vessels['AHTS_alpha'], sc.vessels['HL_Giant'], sc.vessels['CSV_A']], display=1)
            print('assigning the kitchen sink')
            task.assignAssets([sc.vessels['AHTS_alpha'], sc.vessels['HL_Giant'], sc.vessels['CSV_A']])
        
        # Calculation durations of the actions, and then of the task
        for a in task.actions.values():
            a.calcDurationAndCost()
        task.calcDuration()
    
   
    # Example task time adjustment and plot
    sc.tasks['tow_and_hookup'].setStartTime(5)
    sc.tasks['tow_and_hookup'].chart()
    
    #dt_min = sc.figureOutTaskRelationships()
    
    '''
    # inputs for scheduler
    offsets_min = {}        # min: 'taskA->taskB': offset     max: 'taskA->taskB': (offset, 'exact')
    for taskA_index, taskA_name in enumerate(sc.tasks.keys()):
        for taskB_index, taskB_name in enumerate(sc.tasks.keys()):
            if dt_min[taskA_index, taskB_index] != 0:
                offsets_min[f'{taskA_name}->{taskB_name}'] = dt_min[taskA_index, taskB_index]
    '''
    
    # ----- Check tasks for suitable vessels and the associated costs/times -----
    '''
    # preliminary/temporary test of anchor install asset suitability
    for akey, anchor in project.anchorList.items():
        for a in anchor.install_dependencies:  # go through required actions (should just be the anchor install)
            a.evaluateAssets([sc.vessels["MPSV_01"]])  # see if this example vessel can do it


    # ----- Generate the task_asset_matrix for scheduler -----
    # UNUSED FOR NOW
    task_asset_matrix = np.zeros((len(sc.tasks), len(sc.vessels), 2))
    for i, task in enumerate(sc.tasks.values()):
        row = task.get_row(sc.vessels)
        if row.shape != (len(sc.vessels), 2):
            raise Exception(f"Task '{task.name}' get_row output has wrong shape {row.shape}, should be {(2, len(sc.vessels))}")
        task_asset_matrix[i, :] = row
    '''
    # ----- Call the scheduler -----
    # for timing with weather windows and vessel assignments 
    
    records = []
    for task in sc.tasks.values():
        print('')
        print(task.name)
        for act in task.actions.values():
            print(f" {act.name}:  duration: {act.duration:8.2f}   start time: {task.actions_ti[act.name]:8.2f}")
            # start = float(task.actions_ti[name])      # start time [hr]
            # dur   = float(act.duration)               # duration [hr]
            # end   = start + dur
    
            # records.append({
            #     'task'       : task.name,
            #     'action'     : name,
            #     'duration_hr': dur,
            #     'time_label' : f'{start:.1f}–{end:.1f} hr',
            #     'periods'    : [(start, end)],        # ready for future split periods
            #     'start_hr'   : start,                 # optional but handy
            #     'end_hr'     : end
            # })
    
    # Example: 
    # for r in records:
        # print(f"{r['task']} :: {r['action']}  duration_hr={r['duration_hr']:.1f}  "
        #       f"start={r['start_hr']:.1f}  label='{r['time_label']}'  periods={r['periods']}")

    
    # ----- Run the simulation -----
    '''
    for t in np.arange(8760):
        
        # run the actions - these will set the modes and velocities of things...
        for a in actionList:
            if a.status == 0:
                pass
                #check if the event should be initiated
            elif a.status == 1:
                a.timestep()  # advance the action
            # if status == 2: finished, then nothing to do
            
        # run the time integrator to update the states of things...
        for v in self.vesselList:
            v.timestep()
       
        # log the state of everything...
    '''
        
    plt.show()
    