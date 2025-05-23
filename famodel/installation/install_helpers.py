import matplotlib.pyplot as plt
import networkx as nx

def tranportTo_actionItem(vessel, distance2port):
    """
    Creates an action item for transporting a vessel to a site.

    Parameters
    ----------
    vessel : dict
        The vessel to be transported.
    distance2port : float
        The distance to the port.

    Returns
    -------
    action : dict
        The action item for transporting the vessel.
    """
    action = {
        "transport_to_site": {
            "duration": distance2port/vessel['transport_specs']['transit_speed'],
            "dependencies": []
        }
    }    
    return action

def tranportFrom_actionItem(vessel, distance2port, empty_factor=1.0):
    """
    Creates an action item for transporting a vessel from a site.

    Parameters
    ----------
    vessel : dict
        The vessel to be transported.
    distance2port : float
        The distance to the port.
    empty_factor : float, optional
        The factor to account for empty return trip.

    Returns
    -------
    action : dict
        The action item for transporting the vessel.
    """
    action = {
        "transport_from_site": {
            "duration": empty_factor * distance2port/vessel['transport_specs']['transit_speed'],
            "dependencies": []
        }
    }    
    return action

def mobilizeV_actionItem(vessel):
    """
    Creates an action item for mobilizing a vessel.

    Parameters
    ----------
    vessel : dict
        The vessel to be mobilized.

    Returns
    -------
    action : dict
        The action item for mobilizing the vessel.
    """
    action = {
        "mobilize_vessel": {
            "duration": vessel['vessel_specs'],
            "dependencies": []
        }
    }
    return action

def mobilizeM_actionItem(vessel, pkg): # TODO: this seems like a replication of vessel.mobilize. Which one do we use? 
    """
    Creates an action item for mobilizing materials on a vessel.

    Parameters
    ----------
    vessel : dict
        The vessel to be mobilized.
    pkg : dict
        The package of materials to be mobilized.

    Returns
    -------
    action : dict
        The action item for mobilizing the materials.
    vessel : dict
        The updated vessel state.
    """
    winch_speed = vessel['storage_specs']['winch_speed']*60  # m/hr
    anchor_loading_speed = vessel['storage_specs']['anchor_loading_speed']
    
    action = {
        "load spooling": {
            "duration": 1,
            "dependencies": []
        },
        "load line": {
            "duration": 0,
            "dependencies": ["load spooling"]
        },
        "load anchor": {
            "duration": 0,
            "dependencies": []
        },
        "load gear": {
            "duration": 2,
            "dependencies": []
        },
        "seafasten": {
            "duration": 3,
            "dependencies": ["load spooling", "load line", "load anchor", "load gear"]
        }
    }
    
    for key, item in pkg.items():
        item['obj'].inst['mobilized'] = True
        if key.startswith("sec"):  # agnostic to line type
            action["load line"]["duration"] += item['length'] / winch_speed
            vessel['state']['remaining_spool_capacity'] -= item['length']
            
        elif key.startswith("anchor"):  # anchor
            if item['load'] > vessel['storage_specs']['max_deck_load']:
                raise ValueError(f"item {key} has a load higher than what the vessel can withstand.")
            
            action["load anchor"]["duration"] += anchor_loading_speed  # Assuming 1 anchor load = 1 * speed
            vessel['state']['remaining_deck_space'] -= item['space']
    
        vessel['state']['remaining_cargo'] -= item['mass']
        vessel['state']['assigned_materials'].append(item['obj'])
        

    
    return action, vessel

def install_actionItem(vessel, pkg):
    """
    Creates an action item for installing a materials package.

    Parameters
    ----------
    vessel : dict
        The vessel to be used for installation.
    pkg : dict
        The package of materials to be installed.

    Returns
    -------
    action : dict
        The action item for installing the materials.
    vessel : dict
        The updated vessel state.
    """
    action = {
        "position onsite": {"duration": 0, "dependencies": []},
        "site survey": {"duration": 0, "dependencies": ["position onsite"]},
        "install anchor": {"duration": 0, "dependencies": ["position onsite", "site survey"]},
        "rerig deck": {"duration": 0, "dependencies": ["position onsite", "install anchor"]},
        "install line": {"duration": 0, "dependencies": ["install anchor", "rerig deck"]},
    }

    def installItem(key):
        '''
        NOT A PUBLIC FUNCTION
        This function installs an item and its dependencies.
        It checks if the item is already installed and if not, it installs its dependencies first.
        Then, it installs the item itself and updates the vessel state.

        Parameters
        ----------
        key : str
            The key of the item to be installed.

        Returns
        -------
        None
        '''
        item = pkg.get(key)
        for dep in item['dependencies']:
            if not pkg[dep]['obj'].inst['installed']:
                installItem(dep)
        
        if key.startswith("anchor"):
            action["position onsite"]["duration"] = 2  # from PPI (only once per anchor)           
            action["site survey"]["duration"] = 2       # from PPI 
            if item['obj'].dd['design']['type']=='suction':
                pile_fixed = vessel["vessel_specs"]["pile_fixed_install_time"]
                pile_depth = 0.005 * abs(item['obj'].r[-1])

                action["install anchor"]["duration"] = pile_fixed + pile_depth
            else:
                # support for other anchor types
                pass
            
            vessel['state']['remaining_deck_space'] += item.get('space', 0)

        elif key.startswith("sec"):
            if action["install line"]["duration"]==0:
                # first line to install
                action["rerig deck"]["duration"] = vessel['storage_specs'].get('rerig_deck', 0)
            winch_speed = vessel['storage_specs']['winch_speed']*60  # m/hr
            line_fixed = vessel["vessel_specs"]["line_fixed_install_time"]
            line_winch = item['length']/winch_speed
            action["install line"]["duration"] += line_fixed + line_winch
            action["install line"]["dependencies"] = ["install anchor", "rerig deck"]

            vessel['state']['remaining_spool_capacity'] += item.get('length', 0)

        item['obj'].inst['installed'] = True
        vessel['state']['remaining_cargo'] += item['mass']
        vessel['state']['assigned_materials'].remove(item['obj'])

    for key in pkg.keys():
        installItem(key)
    
    return action, vessel

def visualizeAction(action):
    """
    Visualizes the action items as a directed graph.

    Parameters
    ----------
    action : dict
        The action items to be visualized.

    Returns
    -------
    None
    """
    # Create the graph
    G = nx.DiGraph()
    for item, data in action.items():
        for dep in data['dependencies']:
            G.add_edge(dep, item, duration=data['duration'])  # Store duration as edge attribute

    # Compute longest path & total duration
    longest_path = nx.dag_longest_path(G, weight='duration')
    longest_path_edges = list(zip(longest_path, longest_path[1:]))  # Convert path into edge pairs
    total_duration = sum(action[node]['duration'] for node in longest_path)
    if len(longest_path)>=1:
        last_node = longest_path[-1]  # Identify last node of the longest path
        # Define layout
        pos = nx.shell_layout(G)        
        # Draw all nodes and edges (default gray)
        nx.draw(G, pos, with_labels=True, node_size=500, node_color='skyblue', font_size=10, font_weight='bold', font_color='black', edge_color='gray')

        # Highlight longest path in red
        nx.draw_networkx_edges(G, pos, edgelist=longest_path_edges, edge_color='red', width=2)

        # Annotate last node with total duration in red
        plt.text(pos[last_node][0], pos[last_node][1] - 0.1, f"{total_duration:.2f} hr", fontsize=12, color='red', fontweight='bold', ha='center')          
    else:
        pass