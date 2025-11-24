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
    
    def __init__(self, name, actions, action_sequence='series', **kwargs):
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
        action_sequence : string or dict, optional
            If a dictionary, each key is the name of each action, and the values are
            each a list of which actions (by name) must be completed before the current
            one. 
            If a string, indicates which approach is used for automatically
            setting the sequence of actions:
            'series': one after the other based on the order in actions (default),
            'dependencies': based on the dependencies of each action.
        kwargs 
            Additional arguments may depend on the task type.
        
        '''
    

        self.name = name
        print(f" Initializing Task '{self.name}")     
        
        # Save the task's dictionary of actions
        if isinstance(actions, dict):
            self.actions = actions
        elif isinstance(actions, list):  # turn list into a dict based on name
            self.actions = {a.name: a for a in actions}

        # --- Set up the sequence of actions ---
        # key is action name, value is a list of what action names are to be completed before it
        
        
        if isinstance(action_sequence, dict):  # Full dict provided (use directly)
            self.action_sequence = {k: list(v) for k, v in action_sequence.items()}
        
        elif isinstance(action_sequence, str):
            self.action_sequence = {}
            
            if action_sequence == 'series':  # Puts the actions in linear sequence
                actions = list(self.actions.values())
                for i in range(len(actions)):
                    if i==0:  # first action has no dependencies
                        self.action_sequence[actions[i].name] = []
                    else:  # previous action must be done first
                        self.action_sequence[actions[i].name] = [ actions[i-1].name ] 
                
            elif action_sequence == 'dependencies':  # Sequences based on the dependencies of each action
            
                def getDeps(action):
                    deps = []
                    for dep in action.dependencies:
                        deps.append(dep)
                    return deps
            
                self.action_sequence = {self.actions[name].name: getDeps(self.actions[name]) for name in self.actions}
            else:
                raise Exception("Action_sequence must be either 'series' or 'dependencies', or a dict.")
        else:
            raise Exception("Action_sequence must be either a string or dict.")
        
        
        # Initialize some task variables
        self.status = 0  # 0, waiting;  1=running;  2=finished
        self.actions_ti = {}  # relative start time of each action [h]
        
        self.duration = 0.0  # duration must be calculated based on lengths of actions
        self.cost     = 0.0  # cost must be calculated based on the cost of individual actions.
        self.ti       = 0.0  # task start time [h?]
        self.tf       = 0.0  # task end time [h?]
        '''
        # Calculate duration and cost
        self.calcDuration()  # organizes actions and calculates duration
        self.calcCost()

        print(f"---------------------- Initializing Task '{self.name} ----------------------")        
        print(f"Task '{self.name}' initialized with duration = {self.duration:.2f} h.")
        print(f"Task '{self.name}' initialized with cost     = ${self.cost:.2f} ")
        '''
        
        # --- Make a list that conveys the action sequence (similar format as Mooring subcomponents)
        act_sequence = dependenciesToSequence(self.action_sequence)
        # (contents of sequence or action names/keys)
        
        # >>> temporarily hard coded here >>>
        # Here's a list of specs we might want to take the max of instead of sum: Add more as needed
        specs_to_max = ['hook_height_m', 'depth_rating_m',
                        'max_depth_m', 'accuracy_m',
                        'speed_mpm', 'capacity_t']
        
        
        # ----- Get Task requirements -----
        
        # Go through each series step in the task action sequence and figure out its requirements
        # (storage capacities will add, for example)
        
        # A Task-level dependency dict that describes the overally requirements
        # when all actions are combined in sequence
        self.requirements = {}
        
        #req_bases = {}  # accumulation of all the requirements over the task's action sequence (maxes, etc.)
        #
        req_sequences = [[]]  # sequence of totalled up requirements at each step
        # capacity specs will add/subtract, while others will be instantaneous
        # req_sequences is as nested list-list-dict-dict-dict of breakdown -> i_step -> req -> capacity* -> spec
        # Whenever there are multiply capacities in a req, the number of 
        # breakdowns multiplies (branches) by the number of capacities.
        # For storage-y reqs, there will only be one capacity listed in the final data.
        
        for i, step in enumerate(act_sequence):  # go through each step in the action sequence
            
            #print(f" === Step {i} ==== ({len(req_sequences)} breakdowns)")
            
            for j in range(len(req_sequences)):
                #if i == 0:  # first step, start with a blank dict
                req_sequences[j].append({}) # for the reqs to log at this step
                #else:  # copy over the last step's requirements as a starting point
                #    req_sequences[j].append(deepcopy(req_sequences[j][i-1]))
            
            # ----- Parallel actions case -----
            if isinstance(step, list):  # parallel actions
                # Currently, this approach just sums up the requirements/capabilities across
                # parallel actions. This has the implication that one requirement must be 
                # fulfilled by only one type of capability. I.e. chain storage can't be 
                # divided between both deck space and chain locker.
                
                # A constraint to be considered is for parallel actions to be performed
                # by separate vessels. That would require more thought, then implementation.
                
                for step_act in step:
                
                    act = self.actions[step_act]
                
                    # Go through requirements of the single action at this step
                    for req in act.requirements.values():
                    
                        #print(f"   Requirement {req['base']}")
        
                        # is this requirement related to storage?
                        storey = req['base'] in ['storage','chain_storage','rope_storage','cable_storage']
                        
                        nold = len(req_sequences)  # number of possible breakdowns SO FAR
                                                
                        # Iterate for each possible requirements breakdown
                        for j in range(nold):
                            #print(f"   j={j}")
                            
                            # Add an entry for this requirement if one doesn't already exist from last step
                            if not req['base'] in req_sequences[j][i]:
                                req_sequences[j][i][req['base']] = {}  # add entry for this requirement
                            
                            ncaps = len(req['capabilities'])
                            
                            for k, cap in enumerate(req['capabilities']):  # go through capabilities in the req
                                                                
                                #print(f"     k={k} - capability is {cap}")
                                
                                # force to use the same storage as established previously if unloading
                                if storey and req['direction'] == -1:  # if unloading
                                    # look through prevous actions...
                                    doNothing = True
                                    for iprev in range(i-1, -1, -1):
                                        if isinstance(act_sequence[iprev], list): # if parallel actions here
                                            for act2_name in act_sequence[iprev]:
                                                act2 = self.actions[act_sequence[iprev]]
                                                if act.dependsOn(act2):  # if dependency, look for related action
                                                    for req2 in act2.requirements.values():
                                                        # if the same storage requirement gets added to or loaded
                                                        if req['base'] == req2['base'] and req2['direction']==1:
                                                            # Check if the current capability is what was loaded to
                                                            if cap in req_sequences[j][iprev][req['base']]:
                                                                doNothing = False   # flag that this is the case to keep
                                                                break
                                        else:
                                            act2 = self.actions[act_sequence[iprev]] 
                                            if act.dependsOn(act2):  # if dependency, look for related action
                                                for req2 in act2.requirements.values():
                                                    # if the same storage requirement gets added to or loaded
                                                    if req['base'] == req2['base'] and req2['direction']==1:
                                                        # Check if the current capability is what was loaded to
                                                        if cap in req_sequences[j][iprev][req['base']]:
                                                            doNothing = False   # flag that this is the case to keep
                                                            break
                                    
                                    # this must mean we aren't unloading from a prevoiusly loaded capacity in this
                                    # particular loop, so skip it
                                    if doNothing:
                                        continue  # skip the rest of this
                                
                                # make a copy of things if it's a storage-y requirement and being added to
                                # (there are only options when adding to storage, not when removing)
                                if k < ncaps-1 and storey and req['direction'] == 1: 
                                # I guess we need to make a copy of everything that happened before this, 
                                    this_req_sequence = deepcopy(req_sequences[j])
                                
                                else:  # otherwise (i.e. on the last one) work with the current sequence
                                    this_req_sequence = req_sequences[j] 
                                
                                # If this capacity isn't already stored at this req in this step (i.e. this is parallel action 0)
                                #if not cap in this_req_sequence[i][req['base']]:
                                new_cap = {}  # add an entry for the capacity's specs
                                #else:  # otherwise use the existing one 
                                
                                
                                if i==0:  # first step (starts off the values)
                                    
                                    for spec, val in req['capabilities'][cap].items():
                                        new_cap[spec] = val  # add the specs
                                    
                                else:  # subsequent steps (accumulates some things)
                                    
                                    # -- add/subtract capability specs depending on direction --
                                    
                                    last_specs = {}
                                    
                                    if storey:  # if it's a storage related spec, make sure to work with prevous values
                                        for iprev in range(i-1, -1, -1): # look for last time this requirement's capacity came up
                                            if req['base'] in req_sequences[j][iprev] and cap in req_sequences[j][iprev][req['base']]:  
                                                last_specs = req_sequences[j][iprev][req['base']][cap] # cap value in previous step
                                                break
                                    
                                    for spec, val in req['capabilities'][cap].items():  # go through specs of the capability
                                    
                                        if spec in this_req_sequence[i][req['base']][cap]:  # check if it's already here (from a parallel action)
                                            last_val = this_req_sequence[i][req['base']][cap][spec]
                                        elif spec in last_specs:  # otherwise use the previous value if available, so that we add to it
                                            last_val = last_specs[spec]
                                        else:
                                            last_val = 0
                                    
                                        # note: this logic should be good for storagey reqs, but unsure for others, e.g. cranes
                                    
                                        if req['direction'] == 0 or spec in specs_to_max:
                                            new_cap[spec] = max(last_val, val)
                                            
                                        elif req['direction'] == 1:  # add to the previous value
                                            new_cap[spec] = last_val + val  # add to previous value
                                            
                                        elif req['direction'] == -1:  # subtract from the previous value
                                            new_cap[spec] = last_val - val  # subtract from previous value
                                        
                                        else:
                                            raise Exception("Invalid direction value (must be 0, 1, or -1).")
                                        
                                                           
                                this_req_sequence[i][req['base']][cap] = new_cap  # add this req's info (may overwrite in parallel case)
                                
                                # Append this as a new possible sequence
                                if k < ncaps-1 and storey and req['direction'] == 1:
                                    req_sequences.append(this_req_sequence)
                                    # Note: if k==0 then the existing req sequence has already been adjusted
            
            
            # ----- normal case, single action -----
            else:
                act = self.actions[step]
            
                # Go through requirements of the single action at this step
                for req in act.requirements.values():
                
                    #print(f"   Requirement {req['base']}")
    
                    # is this requirement related to storage?
                    storey = req['base'] in ['storage','chain_storage','rope_storage','cable_storage']
                    
                    nold = len(req_sequences)  # number of possible breakdowns SO FAR
                    
                    # >>> bifurcate the current branch of the req_sequences, 
                    # adding n-1 new branches where n is the number of capabilities
                    # (each which represents one possibility for satisfying the req)
                    #n = len(req['capabilities'])
                    
                    # Iterate for each possible requirements breakdown
                    for j in range(nold):
                        #print(f"   j={j}")
                        
                        # Add an entry for this requirement if one doesn't already exist from last step
                        #if not req['base'] in req_sequences[j][i]:
                        req_sequences[j][i][req['base']] = {}  # add entry for this requirement
                        
                        ncaps = len(req['capabilities'])
                        
                        for k, cap in enumerate(req['capabilities']):  # go through capabilities in the req
                            #   for k in range(len(req['capabilities'])-1, -1, -1):  # go through capabilities in the req
                            #cap req['capabilities'].keys()
                            
                            #print(f"     k={k} - capability is {cap}")
                            
                            # force to use the same storage as established previously if unloading:
                            if storey and req['direction'] == -1:  # if unloading
                                keepThisCapability = False
                                for iprev in range(i-1, -1, -1):  # look through prevous actions...
                                    act2 = self.actions[act_sequence[iprev]]
                                    if act.dependsOn(act2):  # do something special here?
                                        
                                        for req2 in act2.requirements.values():
                                            # if the same storage requirement gets added to or loaded
                                            if req['base'] == req2['base'] and req2['direction']==1:
                                                # Check if the current capability is what was loaded to
                                                if cap in req_sequences[j][iprev][req['base']]:
                                                    #if i==4 and j==1:
                                                    #    breakpoint()
                                                    keepThisCapability = True   # flag that this is the case to keep
                                                    break
                                
                                # this must mean we aren't unloading from a prevoiusly loaded capacity in this
                                # particular loop, so kip it
                                if not keepThisCapability:
                                    #if act.name=='lay_mooring-fowt0a':
                                    #    breakpoint()  
                                    #print(f"WARNING - action {act.name} involves unloading storage but the prior load action was not found")
                                    continue  # skip adding this capability
                                    
                                    # >>> still need to add support for parallel actions <<<
                            
                            # make a copy of things if it's a storage-y requirement and being added to
                            # (there are only options when adding to storage, not when removing)
                            #if k < ncaps-1 and storey and req['direction'] == 1:  <<< old one
                            if k < ncaps-1 and not (storey and req['direction'] == -1): 
                            # I guess we need to make a copy of everything that happened before this, 
                                if j>20000:
                                    breakpoint()
                                this_req_sequence = deepcopy(req_sequences[j])
                            
                            else:  # otherwise (i.e. on the last one) work with the current sequence
                                this_req_sequence = req_sequences[j] 
                            
                            new_cap = {}  # add an entry for the capacity's specs
                                
                            if i==0:  # first step (starts off the values)
                                
                                for spec, val in req['capabilities'][cap].items():
                                    new_cap[spec] = val  # add the specs
                                
                            else:  # subsequent steps (accumulates some things)
                                
                                # -- add/subtract capability specs depending on direction --
                                # if this req and cap already exist
                                
                                last_specs = {}
                                
                                if storey:  # if it's a storage related spec, make sure to work with prevous values
                                    for iprev in range(i-1, -1, -1): # look for last time this requirement's capacity came up
                                        if req['base'] in req_sequences[j][iprev] and cap in req_sequences[j][iprev][req['base']]:  
                                            last_specs = req_sequences[j][iprev][req['base']][cap] # cap value in previous step
                                            break
                                
                                #if not cap in req_sequences[j][i][req['base']]:  # if capacity doesn't exist in past
                                #    req_sequences[j][i][req['base']][cap] = {}  # add a blank for it
                                
                                for spec, val in req['capabilities'][cap].items():  # go through specs of the capability
                                
                                    if spec in last_specs:
                                        last_val = last_specs[spec]
                                    #if spec in req_sequences[j][i][req['base']][cap]:
                                    #    last_val = req_sequences[j][i][req['base']][cap][spec]
                                    else:
                                        last_val = 0
                                
                                    if req['direction'] == 0 or spec in specs_to_max:
                                        new_cap[spec] = max(last_val, val)
                                        
                                    elif req['direction'] == 1:  # add to the previous value
                                        new_cap[spec] = last_val + val  # add to previous value
                                        
                                    elif req['direction'] == -1:  # subtract from the previous value
                                        new_cap[spec] = last_val - val  # subtract from previous value
                                    
                                    else:
                                        raise Exception("Invalid direction value (must be 0, 1, or -1).")
                                    
                                    #print(f" {act.name}   {req['base']}  {cap}  {spec} ")
                                    
                                    #... also check if a spec value is going to go below zero, leave at zero ^^^
                                    # also distinguish between stock and flow specs, e.g. some to max vs add/subtract ^^^
                            
                            #if act.name == 'mooring_hookup-fowt0a':
                            #    breakpoint()
                            
                            this_req_sequence[i][req['base']][cap] = new_cap  # add this req's info
                            
                            #if j > 40:
                            #    breakpoint()
                            
                            # Append this as a new possible sequence
                            #if k < ncaps-1 and storey and req['direction'] == 1:
                            if k < ncaps-1 and not (storey and req['direction'] == -1): 
                                req_sequences.append(this_req_sequence)
                                # Note: if k==0 then the existing req sequence has already been adjusted
                    
        print(f"Task requirements processed. There are {len(req_sequences)} possible combinations.")                    
        
        
        # Go through the requirements sequence and find the maximum values
        # These become the overall requirements of the task.
        task_reqs = []
        for j in range(len(req_sequences)):
            task_reqs.append({})  # An empty dictionary of requirements for this breakdown
            
            for i, rs in enumerate(req_sequences[j]):
                
                for req, caps in rs.items():
                    
                    # if req not already in the list, add it
                    if not req in task_reqs[j]:  
                        task_reqs[j][req] = {}
                
                    # go through req capabilities
                    for cap, specs in caps.items():
                        if not cap in task_reqs[j][req]:  # if cap not in the list,
                            task_reqs[j][req][cap] = {}  # add it
                        
                        # go through capability specs and take the maxes
                        for spec, val in specs.items():
                            if spec in task_reqs[j][req][cap]:
                                last_val = task_reqs[j][req][cap][spec]
                            else:
                                last_val = 0
                            
                            # Retain the max value of the spec
                            task_reqs[j][req][cap][spec] = max(last_val, val)
        
        if len(req_sequences) > 20000:
            breakpoint()
            print("there's a lot of options")
        
        # Save things
        self.act_sequence = act_sequence
        self.req_sequences = req_sequences
        self.task_reqs = task_reqs
    
    
    
    def checkAssets(self, assets, display=0):
        '''
        Checks if a specified set of assets has sufficient capabilities and 
        specs to fulfill all requirements of this task.
        
        Parameters
        ----------
        assets : list of assets
        '''
        
        # this should evaluate the assets w.r.t. self.task_reqs
        
        
        
        # Sum up the asset capabilities and their specs (not sure this is useful/valid)
        
        # Here's a list of specs we might want to take the max of instead of sum: Add more as needed
        '''
        specs_to_max = ['hook_height_m', 'depth_rating_m',
                        'max_depth_m', 'accuracy_m',
                        'speed_mpm', 'capacity_t']  # capacity_t is here because it doesn't make sense to have two cranes to lift a single anchor. 
        '''
        asset_caps = combineCapabilities(assets, display=display-1)
        
        # <<< maybe instead of all this we should do an approach that looks by asset
        #     because that could then also be used to decide asset assignment
        #     to each requirement >>>
        
        
        # See if summed asset capabilities satisfy any of the n task_req breakdowns
        # .>>> an output of this could also be assigning assets to action requirements!! >>>
        
        requirements_met = []
        assignable = []
        
        for i in range(len(self.task_reqs)):
        
            if display > 2:  print(f"Task {self.name} requirements breakdown #{i}:")
            
            requirements_met.append({})
            
            requirements_met[i] = doCapsMeetRequirements(asset_caps, self.task_reqs[i])
            '''
            
            for req, caps in self.task_reqs[i].items():  # go through each requirement
                                
                requirements_met[i][req] = False  # start assume it is not met

                # Let's check if each capability is sufficiently provided for
                capable = True  # starting with optimism...
                
                for cap, specs in caps.items():  # go throuch each capability of the requirement
                    
                    if cap not in asset_caps: # assets don't have this capability, fail
                        capable = False
                        if display > 2: print(f"Warning: capability '{cap}' is missing from the assets.")
                        break
                    
                    for key, val in specs.items():  # go through each spec for this capability
                        
                        if val == 0:  # if zero value, no spec required, move on
                            continue
                        if key not in asset_caps[cap]:  # if the spec is missing, fail
                            capable = False
                            if display > 2: print(f"Warning: capability '{cap}' does not have spec '{key}'.")
                            break
                        if asset_caps[cap][key] < val: # if spec is too small, fail
                            capable = False
                            if display > 2: print(f"Warning: capability '{cap}' does not meet spec '{key}' requirement of {val:.2f} (has {asset_caps[cap][key]:.2f}).")
                            break
                
                # Final call on whether requirement can be met 
                if capable:
                    requirements_met[i][req] = True 
                else:
                    requirements_met[i][req] = False
                    if display > 1: print(f"Requirement '{req}' is not met by asset(s):")
                    if display > 2: print(f"{assets}.")
                    break  
            '''
            # Check if all requirements are met by the assets for this breakdown
            assignable.append(all(requirements_met[i].values()))
            if display > 1: print(f" Suitability is {assignable[i]}.")
                
        if display > 0: 
            print(f"Finished checking assets. {sum(assignable)} of {len(assignable)} requirement breakdowns are feasible.")
        '''
        if self.name =='install_all_anchors':
            for i in range(len(self.task_reqs)):
            
            
                print(doCapsMeetRequirements(asset_caps, self.task_reqs[i]))
                
                if not 'divers' in self.task_reqs[i]['anchor_orienting']:
                    print(i)
                    printStruct(self.task_reqs[i])
        
            breakpoint()
        '''
        '''  (Older method that looks at any capability being satisfied)
        requirements_met = {}
        for req, vals in self.requirements.items():  # go through each requirement
            
            caps = vals['capabilities']
            dir = vals['direction']
            
            # The following logic should mark a requirement as met if any one of 
            # the requirement's needed capabilities has all of its specs by the
            # combined spec values of the assets
            
            requirements_met[req] = False  # start assume it is not met

            for cap, specs in caps.items():  # go throuch capability of the requirement
                if cap not in asset_caps: # assets don't have this capability, move on
                    continue
                
                # Let's check if this capability is sufficient
                capable = True
                for key, val in specs.items():  # go through each spec for this capability
                    
                    if val == 0:  # if zero value, no spec required, move on
                        continue
                    if key not in asset_caps[cap]:  # if the spec is missing, fail
                        capable = False
                        print(f"Warning: capability '{cap}' does not have metric '{key}'.")
                        break
                    if asset_caps[cap][key] < val: # if spec is too small, fail
                        # note: may need to add handling for lists/strings, or standardize specs more
                        capable = False
                        print(f"Warning: capability '{cap}' does not meet metric '{key}' requirement of {val:.2f} (has {asset_caps[cap][key]:.2f}).")
                        break
                    
                if capable:
                    requirements_met[req] = True  # one capability fully satisfies the requirement
                    break                         # no need to check other capabilities for this requirement
            
            if not requirements_met[req]:
                print(f"Requirement '{req}' is not met by asset(s): {assets}.")
        
        assignable = all(requirements_met.values())
        
        # message:
        if assignable:
            message = "Asset meets all required capabilities."
        else:
            unmet = [req for req, met in requirements_met.items() if not met]
            detailed = []
            for req in unmet:
                expected = [cap for cap in self.requirements[req].keys()]
                detailed.append(f"- {req}: {expected}.")
                detailed_msg = "\n".join(detailed)
            
            message = "Asset does not meet the following required capabilities:\n" + detailed_msg        
        '''
        
        
        # return bool of if any req breakdowns can be satisfied, and a list of which ones
        return any(assignable), assignable
    
    
    def assignAssets(self, assets, display=0):
        '''Figures out an assignment of the asset capabilities to the task's
        steps' requirements, including each action's requirements.'''
        
        doable, indices = self.checkAssets(assets)
        
        if doable:
            
            # sum up combined asset capabilities
            asset_caps = combineCapabilities(assets)
                        
            # take the first requirement breakdown that works
            ind = indices.index(True)  # get the index of the first true value
            
            # Select that breakdown (update Task's/actions' requirements) <<< can be turned into method
            # Go through and delete any requirements in the actions that don't correspond to this breakdown
            # traverse action sequence
            for i, step in enumerate(self.act_sequence):  # go through each step in the action sequence
            
                if isinstance(step, list): # Parallel actions case 
                    for j in len(step):  # each parallel action
                        pass  # (we don't actually know how to handle this yet) <<<
                        
                else:  # normal case (one action at a time)
                    
                    action = self.actions[self.act_sequence[i]]  # this step's action
                    reqs = self.req_sequences[ind][i] # altered/active requirements at this step
                                        
                    for areq in action.requirements.values():
                        if not areq['base'] in reqs: 
                            raise Exception(f"Action {action.name} somehow has a req that isn't in the task's req_sequence")
                        
                        # Create selected_capabilities (or clear it if it already exists)
                        areq['selected_capability'] = {}
                        
                        for acap, aspecs in areq['capabilities'].items():  # cycle through action's reqs capability keys
                            if acap in reqs[areq['base']]:  # if this capability is listed, it means we plan to use it
                                areq['selected_capability'][acap] = aspecs  # so copy it over
                                # (there should only be one capability selected per requirement)
                                
                                # Note which asset(s) are planned to fill this req
                                for ass in assets:
                                    met = checkCapability(areq['selected_capability'], [ass], acap)
                                    if met:
                                        areq['assigned_assets'] = [ass]
                                        break
                                
                                if not met:
                                    met = checkCapability(areq['selected_capability'], assets, acap)
                                    if met:
                                        areq['assigned_assets'] = assets
                                    else:
                                        raise Exception(f"Task {self.name} could not satisfy action {action.name} capability {acap} with the available assets.")
                        
                        
                    action.assignAssets(assets)
            
            self.assetList = assets
            
            if display > 0:
                print(f"For the task {self.name}, assigned the assets:")
                print([asset['name'] for asset in assets])
        else:
            print("This asset assignment is not feasible for the task.")
    
    
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
            # (set as zero if the action does not depend on other actions in the task)
            starts[action] = max((finishes[dep] for dep in dep_actions), default=0)

            # get duration from actions
            duration = self.actions[action].duration  # in hours

            # Calculate finish time
            finishes[action] = starts[action] + duration

        # Update self.actions_ti with relative start times
        self.actions_ti = starts

        # Task duration
        self.duration = max(finishes.values())

        # Update task finish time
        self.tf = self.ti + self.duration

    
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


    def setStartTime(self, start_time):
        '''Update the start time of all actions based on a new task start time.
        This requires that the task's duration and relative action start times are
        already calculated.

        Parameters
        ----------
        newStart : float
            The new start time for the task. All action start times will be adjusted accordingly.
        '''

        # Update task start and finish times
        self.ti = start_time
        self.tf = start_time + self.duration

        # Update action times
        for name, action in self.actions.items():
            action.setStartTime(start_time + self.actions_ti[name])
    

    def clearAssets(self):
        '''
        Clear all assigned assets from all actions in the task.
        This resets the asset assignments and re-evaluates the actions.
        '''
        for action in self.actions.values():
            action.clearAssets()
        
        # Reinitialize duration and cost after clearing assets.
        self.duration = 0
        self.cost     = 0

    
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

        fig, ax = plt.subplots(figsize=(6, 6))

        # Prepare data for Gantt chart
        action_names = list(self.actions.keys())
        start_times = [self.actions_ti[name] for name in action_names]
        durations = [self.actions[name].duration for name in action_names]
        
        # Get asset information from action.assets
        all_assets = [asset['name'] for asset in self.assetList]  # list of asset names
        '''
        all_assets = set()
        #all_roles  = set()
        for action in self.actions.values():
            for asset in action.assetList:
                all_assets.add(asset)
                #all_roles.add(role)
        '''
        # Assign colors
        if color_by == 'asset':
            color_dict = {asset: colors[i] for i, asset in enumerate(all_assets)}
        '''
        elif color_by == 'role':
            # Flip the colors
            colors = colors[::-1]
            role_list = list(all_roles)
            color_dict = {role: colors[i] for i, role in enumerate(role_list)}
        '''
        # Generate vertical lines to indicate the start and finish of the whole task
        ax.axvline(x=self.ti, ymin=0, ymax=len(action_names), color='black', linestyle='-', linewidth=2.0)
        ax.axvline(x=self.tf, ymin=0, ymax=len(action_names), color='black', linestyle='-', linewidth=2.0)

        # Create bars for each action
        ht = 0.4
        for i, (name, start, duration) in enumerate(zip(action_names, start_times, durations)):
            opp_i = len(action_names) - i - 1  # to have first action on top
            action = self.actions[name]
            assets = list({asset['name'] for asset in action.assetList})
            #roles  = list({role for role in action.assets.keys()})
            
            assets = list(set(assets))  # Remove duplicates from assets

            n_assets = len(assets)
            #n_roles  = len(roles)
   
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
                '''
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
                '''
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
        '''
        elif color_by == 'role':
            handles = [plt.Rectangle((0, 0), 1, 1, color=color_dict[a]) for a in all_roles]
            ax.legend(handles, all_roles, title='Roles', bbox_to_anchor=(1.02, 1), loc='upper right')            
        '''
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
    
    
    def chart(self, title='', outpath='', dpi=200):
        '''
        Render a Gantt-like chart for a single ChartTask with one axes and one horizontal lane per vessel.
        • Vessel names as y-tick labels
        • Baseline arrows, light span bars, circle bubbles with time inside, title above,
          and consistent font sizes.
        • Horizontal placement uses Bubble.period when available; otherwise cumulative within vessel.
        • Bubbles are colored by Bubble.category (legend added).
        
        Show an action on multiple lanes if it uses multiple assets.
        Skip actions with dur<=0 or with no resolvable lanes.
        '''
        
        # MH: unsure how much of this up-front stuff is needed

        from dataclasses import dataclass
        from typing import List, Optional, Dict, Tuple
        import matplotlib.pyplot as plt

        # Data structures

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

        # Color palette + categorization

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
                
        
        
        # MH: making a vessels dict to fit with what this was looking for (quick temporary solution)
        vessels = { ves['name'] : ves for ves in self.assetList }
        
        # reverse lookup for identity → key
        id2key = {id(obj): key for key, obj in vessels.items()}

        # unique type → key (used only if type is unique in catalog)
        type_counts = {}
        for k, obj in vessels.items():
            t = obj.get('type') if isinstance(obj, dict) else getattr(obj, 'type', None)
            if t:
                type_counts[t] = type_counts.get(t, 0) + 1
        unique_type2key = {}
        for k, obj in vessels.items():
            t = obj.get('type') if isinstance(obj, dict) else getattr(obj, 'type', None)
            if t and type_counts.get(t) == 1:
                unique_type2key[t] = k

        buckets = {}

        for a in self.actions.values():
            if a.duration <= 0.0:
                continue  # skip if no duration

            aa = getattr(a, 'assets', {}) or {}

            # collect ALL candidate roles → multiple lanes allowed
            lane_keys = set()
            for v in a.assetList:

                # resolve lane key
                lane = id2key.get(id(v))
                if lane is None and isinstance(v, dict):
                    nm = v.get('name')
                    if isinstance(nm, str) and nm in vessels:
                        lane = nm
                    else:
                        t = v.get('type')
                        if t in unique_type2key:
                            lane = unique_type2key[t]
                if lane:
                    lane_keys.add(lane)

            if not lane_keys:
                continue
            
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
                    duration_hr=a.duration,
                    label_time=getattr(a, 'label_time', f'{a.duration:.1f}'),
                    period=( a.ti, a.tf ),
                    category=cat_for(a))
                
                buckets.setdefault(lane, []).append(b)
                #breakpoint()
                #print('hi')

        # preserve sc.vessels order; only include lanes with content
        lanes = []
        for vname in vessels.keys():
            blist = sorted(buckets.get(vname, []), key=lambda b: b.period[0])
            if blist:
                lanes.append(VesselTimeline(vessel=vname, bubbles=blist))

        # Core plotter (single-axes, multiple lanes)

        from matplotlib.lines import Line2D
        from matplotlib.patches import Circle

        # --- figure geometry ---
        nrows = max(1, len(lanes))
        fig_h = max(3.0, 0.8 + 1*nrows)
        fig_w = 5

        plt.rcdefaults()
        plt.close('all')
        fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)

        # --- y lanes (top -> bottom keeps given order) ---
        vessels_top_to_bottom = lanes   # <<< this seems to just make a duplicate reference for the same data
        nrows = max(1, len(lanes))
        y_positions = list(range(nrows))[::-1]
        name_to_y = {vt.vessel: y_positions[i] for i, vt in enumerate(vessels_top_to_bottom[::-1])}
        
        ax.set_yticks(y_positions)
        ax.set_yticklabels([])          
        ax.tick_params(axis='y', labelrotation=0)   
        
        if len(title) > 0:
            ax.set_title(title, loc='left', pad=6)

        # --- gather periods, compute x-range ---
        x_min, x_max = 0.0, 0.0
        per_row: Dict[str, List[Tuple[float, float, Bubble]]] = {vt.vessel: [] for vt in lanes}

        for vt in lanes:
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
            ax.plot([s, e], [y_val, y_val], lw=4, alpha=0.15, color='k')

        def _bubble_face_color(b: Bubble) -> str:
            cat = b.category or 'Other'
            return ACTION_TYPE_COLORS.get(cat, ACTION_TYPE_COLORS['Other'])

        def _text_color_for_face(face: str) -> str:
            return 'black' if face.lower() in ('#ffdd00', 'yellow') else 'white'

        def _draw_bubble(s: float, e: float, y_val: float, b: Bubble, i_in_row: int):
            xc = 0.5*(s + e)
            face = _bubble_face_color(b)
            txtc = _text_color_for_face(face)
            ax.plot(xc, y_val, 'o', ms=22, color=face, zorder=3)
            ax.text(xc, y_val, f'{b.label_time}', ha='center', va='center', fontsize=10,
                    color=txtc, weight='bold')
            title_offset = 0.30 if (i_in_row % 2) else 0.20
            ax.text(xc, y_val + title_offset, b.action, ha='center', va='bottom', fontsize=6)
            # caps_txt = _capabilities_to_text(b.capabilities)
            # if caps_txt:
            #    ax.text(xc, y_val - 0.26, caps_txt, ha='center', va='top', fontsize=8, wrap=True)

        # --- draw per lane ---
        seen_cats: set[str] = set()
        for vt in lanes:
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
            handles.append(Line2D([0], [0], marker='o', linestyle='none', markersize=6,
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
                fontsize=6, color='black', fontfamily='monospace',  # <- key: monospace
                bbox=dict(boxstyle='circle,pad={:.2f}'.format(circle_pad),
                          facecolor='lightgrey', edgecolor='tomato', linewidth=3))

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






def dependenciesToSequence(dependencies):
    '''
    Receive a dictinoary of item dependencies that define a sequence,
    and generate a nested list that follows that sequence.
    
    Example:
        B    A    G  D  F   H
        C               E 
    
    dependencies = {'a': ['c'], 'b': [],
                    'c': [], 'd': ['g'],
                    'e': ['d'], 'f': ['d'],
                    'g': ['a'], 'h': ['e','f']}
    '''
    
    n = len(dependencies)  # number of actions in this task
    acts = list(dependencies.keys())  # list of action names
    deps = list(dependencies.values())  # dependencies of each action
    si = np.zeros(n)-1  # step index of action (-1 = TBD)
    
    sequence = [[]] # create first step slot in the sequence
    for i in range(n):  # go through action: dependencies
        if len(deps[i])==0:  # no dependency, it's in first (0) step
            si[i] = 0  # mark as being in the first step
            sequence[0].append(acts[i])

    for j in range(1,n):  # look for step j actions
        #print(f"Step {j} ----")
        sequence.append([]) # create next step slot in the sequence
        for i in range(n):  # go through action: dependencies
            #print(f"  Action {i}")
            if si[i] < 0:  # only look at actions that aren't yet sequenced
                if any([prev_act in deps[i] for prev_act in sequence[j-1]]): 
                    si[i] = j
                    sequence[j].append(acts[i])
    
    # Clean up the sequence
    clean_sequence = []
    for step in sequence:
        if len(step) == 1:
            clean_sequence.append(step[0])  # add single entry by itself (not a list)
        elif len(step) == 0:
            break  # if we've hit an empty step, we're at the end
        else:
            clean_sequence.append(step)
    
    return clean_sequence


def combineCapabilities(assets, display=0):
    '''Combines the capabilies across multiple assets.'''
    
    specs_to_max = ['hook_height_m', 'depth_rating_m',
                    'max_depth_m', 'accuracy_m',
                    'speed_mpm', 'capacity_t']
    
    asset_caps = {}
    for asset in assets:
        for cap, specs in asset['capabilities'].items():
            if not cap in asset_caps:  # add the capability entry if absent
                asset_caps[cap] = {}
            for key, val in specs.items():
                if key in asset_caps[cap]:
                    if key in specs_to_max:
                        asset_caps[cap][key] = max(asset_caps[cap][key], val)
                    else:
                        asset_caps[cap][key] += val  # add to the spec
                else:
                    asset_caps[cap][key] = val  # create the spec
    
    if display > 0: 
        print('Combined asset specs are as follows:')
        for cap, specs in asset_caps.items():
            print(f'  Capability {cap}')
            for key, val in specs.items():
                 print(f'    Total spec {key} = {val}')

    return asset_caps


def checkCapability(required_capability, assets, capability_name, display=0):
    '''Check if the required capability can be met by the combination
    of the assets specified.'''
    
    # required_capability is assumed tobe a dict of cap_name : specs, meaning capability_name is probably redundant <<<

    asset_caps = combineCapabilities(assets)
    
    
    # See if summed asset capabilities satisfy any of the n task_req breakdowns
    # .>>> an output of this could also be assigning assets to action requirements!! >>>
    
    requirements_met = []
    assignable = []

    # Let's check if each capability is sufficiently provided for
    capable = True  # starting with optimism...
    
    for cap, specs in required_capability.items():  # go throuch each capability of the requirement
    
        if not cap == capability_name:
            breakpoint()
            print('there is a contradiction...')
            
        
        if cap not in asset_caps: # assets don't have this capability, fail
            capable = False
            if display > 2: print(f"Warning: capability '{cap}' is missing from the assets.")
            break
        
        for key, val in specs.items():  # go through each spec for this capability
            
            if val == 0:  # if zero value, no spec required, move on
                continue
            if key not in asset_caps[cap]:  # if the spec is missing, fail
                capable = False
                if display > 2: print(f"Warning: capability '{cap}' does not have spec '{key}'.")
                break
            if asset_caps[cap][key] < val: # if spec is too small, fail
                capable = False
                if display > 2: print(f"Warning: capability '{cap}' does not meet spec '{key}' requirement of {val:.2f} (has {asset_caps[cap][key]:.2f}).")
                break
    
    # Final call on whether requirement can be met 
    if capable:
        return True 
    else:
        return False
    

def doCapsMeetRequirements(asset_caps, requirements, display=0):
    '''Checks if asset capabilities collectively can satisfy the listed
    requirements.'''
    
    requirements_met = {}  # dictionary of requirements being met True/False
    
    for req, caps in requirements.items():  # go through each requirement
                        
        requirements_met[req] = False  # start assume it is not met

        # Let's check if each capability is sufficiently provided for
        capable = True  # starting with optimism...
        
        for cap, specs in caps.items():  # go throuch each capability of the requirement
            
            if cap not in asset_caps: # assets don't have this capability, fail
                capable = False
                if display > 2: print(f"Warning: capability '{cap}' is missing from the assets.")
                break
            
            for key, val in specs.items():  # go through each spec for this capability
                
                if val == 0:  # if zero value, no spec required, move on
                    continue
                if key not in asset_caps[cap]:  # if the spec is missing, fail
                    capable = False
                    if display > 2: print(f"Warning: capability '{cap}' does not have spec '{key}'.")
                    break
                if asset_caps[cap][key] < val: # if spec is too small, fail
                    capable = False
                    if display > 2: print(f"Warning: capability '{cap}' does not meet spec '{key}' requirement of {val:.2f} (has {asset_caps[cap][key]:.2f}).")
                    break
        # Final call on whether requirement can be met 
        if capable:
            requirements_met[req] = True 
        else:
            requirements_met[req] = False
            if display > 1: print(f"Requirement '{req}' is not met by asset(s):")
            if display > 2: print(f"{assets}.")
    
    return requirements_met


def printStruct(t, s=0):

    if not isinstance(t,dict) and not isinstance(t,list):
        print(" "*s+str(t))
    else:
        for key in t:
            if isinstance(t,dict) and not isinstance(t[key],dict) and not isinstance(t[key],list):
                print(" "*s+str(key)+"  :  "+str(t[key]))
            else:
                print(" "*s+str(key))
                if not isinstance(t,list):
                    printStruct(t[key], s=s+2)


if __name__ == '__main__':
    pass
    
    
    
        