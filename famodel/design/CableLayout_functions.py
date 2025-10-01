# -*- coding: utf-8 -*-
import os
import numpy as np
from sklearn.cluster import SpectralClustering
# from sklearn.cluster import SpectralClustering
from scipy.spatial.distance import cdist, pdist, squareform
import networkx as nx
import math
import pandas as pd
import matplotlib.pyplot as plt
from famodel.cables.cable_properties import *
from shapely.geometry import Point, LineString, MultiPoint
import shapely as sh
from copy import deepcopy


# TODO: rename and reorder inputs

def getCableLayout(turb_coords, subs_coords, conductor_sizes, 
                   cableProps_type, turb_rating_MW, turb_cluster_id=[], turb_subs_id=[], 
                   n_cluster_sub=0, n_tcmax=8, plot=False, oss_rerouting=False,
                   substation_id=None):
    ''' Function creating the cable layout of whole wind farm, including 
    estimation of cable conductor sizes. It currently supports a single 
    substation.
    
    Parameters
    ----------        
    turb_coords : 2D array
        Coordinates of each turbine, provided as an N-by-2 array of [x,y] values [m].
    subs_coords : list or array
        Substation [x,y] coordinates [m].
    conductor_sizes ; list
        Conductor sizes to be allowed when sizing cables [mm^2].
    cableProps_type : string
        Name of cable type in cableProps property scaling coefficients yaml.
    turb_rating_MW : float
        Turbine rated power [MW]
    turb_cluster_id : list (optional)
        The index of the cluster (integers starting from zero) that each 
        turbine belongs to. This is specified to determines the clusters.
    turb_subs_id : list (optional)
        The index of the substation (integers starting from zero) that each 
        turbine should feed to.
    n_cluster_sub : int (optional)
        Then number of clusters per substation to create if clustering automatically
        (turb_cluster_id should not be specified in this case).
    n_tcmax : int (optional)
        Then number of clusters to create if clustering automatically
        (turb_cluster_id should not be specified in this case).
    plot : bool (optional, default False)
        Displays a plot of the array cable network if True.
    
    Returns
    -------
    iac_dic : list of dicts
        List of an array cable information dictionary for each cable.
    '''                       

    # Handle if coordinates are inputted as lists rather than arrays
    if type(turb_coords) == list:
        turb_coords = np.array(turb_coords)
    n_turb = turb_coords.shape[0]  # number of turbines
    
    if type(subs_coords) == list:
        subs_coords = np.array(subs_coords)
        
    if subs_coords.shape == (2,):  # if just an x,y pair, put in a 2D array
        subs_coords = np.array([subs_coords])
    n_subs = subs_coords.shape[0]  # number of substations
        

    # Get cable properties
    iac_props = []
    for A in conductor_sizes:
        cprops = getCableProps(A, cableProps_type, cableProps=None, source='default', name="", rho=1025.0, g=9.81)
        iac_props.append(cprops)
    
    
    # ----- Divide turbines among substations -----
    
    if len(turb_subs_id) > 0: # if substation assignment indices are provided
        
        subs_labels_unique, subs_labels_counts = np.unique(turb_subs_id,return_counts=True)
        if not n_subs == len(subs_labels_unique):
            raise Exception("There are more unique entries in turb_subs_id than number of subs_coords provided.")
        turb_subs_id = np.array(turb_subs_id)
        
        # Check that substation labels are integers counting up from 0
        for i in range(n_subs):
            if not i in subs_labels_unique:
                raise Exception(f"provided substation assignment labels must be integers counting up from 0. Integer {i} was not found.")
        
    else:  # If no substation assignments are provided, divide the turbines by distance
        # determine max # of turbines allowed per substation
        max_turbs_per_substation = n_cluster_sub*n_tcmax + 5 # num clusters x num turbs per cluster + a few extra
        turb_subs_id = assignSubstationTurbines(turb_coords, subs_coords, max_turbs_per_substation)
    
    
    # ----- Handle turbine clustering -----
    
    if len(turb_cluster_id) > 0: # if cluster indices are provided
        
        cluster_labels_unique, cluster_labels_counts = np.unique(turb_cluster_id,return_counts=True)
        n_cluster = len(cluster_labels_unique)
        turb_cluster_id = np.array(turb_cluster_id)

        
        # Check that cluster labels are integers counting up from 0
        cluster_subs_id = []
        for i in range(n_cluster):
            cluster_subs_id.append(int(np.unique(turb_subs_id[turb_cluster_id==i])))
            if not i in cluster_labels_unique:
                raise Exception(f"provided cluster labels must be integers counting up from 0. Integer {i} was not found.")
                
        # TODO: figure out how to deal with inconsistencies between turbine cluster vs substation assignments
        
    else:  # If no clusters are provided, create clusters
        
        n_cluster = 0
        if n_cluster_sub == 0:  # if number of clusters (per substation) not specified, use default
            n_cluster_sub = int(np.ceil(n_turb/n_tcmax/n_subs))
    
        # cluster turbines (for each substation if multiple)
        turb_cluster_id = [None]*n_turb  # cluster ID of each turbine
        cluster_labels_counts = []  # the number of turbines in each cluster
        cluster_subs_id = []  # substation ID of each cluster
        for i in range(n_subs):
            turbs = np.where(turb_subs_id==i)[0]
            cluster_id, labels_counts = clusteringSpectral(turb_coords[turbs], 
                                        subs_coords[i,:], n_cluster_sub, n_tcmax)
        
            # Store each turbine's cluster ID (adjusting IDs for multiple substations)
            for ii,cid in enumerate(cluster_id):
                turb_cluster_id[turbs[ii]] = int(cid + n_cluster) 
            cluster_subs_id += list([int(x) for x in np.zeros(len(labels_counts)) + i])
                
            cluster_labels_counts += list(labels_counts)
            
            n_cluster += len(labels_counts)  # tally up actual number of clusters
    
    
    
    # ----- Figure out cable connections for each cluster -----    
    if not substation_id:
        substation_id = []
        for i in range(n_subs):
            substation_id.append(n_turb + i)
    
    index_map = []  # maps local turbine id within each cluster to the global turbine list index
    
    # The main outputs of this part of the code (one entry per cable)
    global_edge_list = []  # global end connection ids of each cable (a, b)
    upstreamturb_count = []  # number of turbines upstream of each cable
    cable_cluster_id = []  # id number of the cluster each cable belongs to
    
    cable_types = []  # list of the cable type dict for each cable
    
    for ic in range(n_cluster):   # for each cluster
        # Select indices of points per cluster
        cluster_ind = np.where(np.array(turb_cluster_id) == ic)[0] 
        
        # Index of the substation for this cluster
        #isubs = turb_subs_id[cluster_ind[0]]
        isubs = cluster_subs_id[ic]
        
        # ----- Make coordinate lists for each cluster, and index map -----
        
        # Make array of just the coordinates in the cluster
        cluster_coords = turb_coords[cluster_ind,:]

        #cluster_arrays.append(cluster_coords)
        # Make list of global turbine indicies that are within this cluster
        index_map.append(np.arange(n_turb)[cluster_ind])
        
        
        # Distances from substation to turbine locations for cluster
        distances = np.linalg.norm(cluster_coords - subs_coords[isubs,:], axis=1)
        # Find the index of the closest turbine to substation
        gate_index0 = np.argmin(distances)
        
        # Calculate minimum spanning tree for the cluster
        cluster_edge_list = minimum_spanning_tree(cluster_coords, gate_index0)
        # This is a list of [a, b] pairs of turbine indices where, within each
        # pair, the power flow is from b to a, and a is closer to the substation.
        
        # Get number of upstream turbines per turbine, counting th 
        iac_upstreamturb_count_ic = getUpstreamTurbines(cluster_edge_list)
        # iac_upstreamturb_count_ic is now a list giving the number of 
        # upstream turbines for each cable, with the same indexing as 
        # cluster_edge_list.
        
        # Convert cluster edge list into global turbine IDs
        for ia, ib in cluster_edge_list:
            global_edge_list.append([index_map[ic][ia],
                                     index_map[ic][ib]])
            cable_cluster_id.append(ic)
            
            upstreamturb_count.append(iac_upstreamturb_count_ic[ib] + 1)
            
        # determine which substation this cable goes to based on cluster to substation index mapping
        subid = substation_id[isubs]
           
        # Add the cable that goes from the substation to the cluster gate
        global_edge_list.append([subid, index_map[ic][gate_index0]])
        cable_cluster_id.append(ic)
        upstreamturb_count.append(cluster_labels_counts[ic])  # (cable to substation)
        
        # Get cable id and assign cable to turbine
        #iac_cab2turb_ic2 = getCableID(cluster_coords, gate_coords[ic], 
        #                          cluster_edge_list, iac_upstreamturb_count_ic)
        # iac_cab2turb_ic = [[el[1], i, el[0]] for i, el in enumerate(cluster_edge_list)] 
        # Above is no longer used <<<
    
    
    # ----- Size cables and generate dictionary of cable information -----
    
    # results of the previous stage are stored in 
    # - global_edge_list
    # - upstreamturb_count
    # - cable_cluster_id
    
    # combine coordinates for easy plotting of everything
    coords = np.vstack([turb_coords, subs_coords])
        
    iac_dic = []  # list of dictionaries for each cable's information
    
    # loop through ALL cables
    for i in range(len(global_edge_list)):
        
        # Size cable to support cumulative power up to this point
        required_rated_power = turb_rating_MW * upstreamturb_count[i]
        selected_cable = selectCable(required_rated_power, iac_props)
        
        cable_types.append(selected_cable)
        
        # note: turb_id_A/B is currently opposite of cluster_edge_list [a,b] <<<
        turb_id_A = global_edge_list[i][1]
        turb_id_B = global_edge_list[i][0]
            
        coordinates = [[coords[turb_id_A][0], coords[turb_id_A][1]],
                       [coords[turb_id_B][0], coords[turb_id_B][1]]]
              
        iac_dic.append({'cluster_id': cable_cluster_id[i],
                     'turbineA_glob_id': turb_id_A,  # row_id_A,
                     'turbineB_glob_id': turb_id_B,  # row_id_B,
                     'cable_id': i,  # this is the global id
                     'upstream_turb_count': upstreamturb_count[i],
                     '2Dlength': np.linalg.norm(coords[turb_id_A] - coords[turb_id_B]),
                     'coordinates': coordinates,  # end/turbine coordinates: [[xA,yA],[xB,yB]]
                     'conductor_area': selected_cable['A'],
                     'cable_costpm': selected_cable['cost']})
    
    
    """
    
    # >>> This section has draft rerouting capability for cable to substation. <<<
    # oss_rerouting : cable rerouting to avoid intersections of cables between clusters and substation. True = on, False = off
    intersection_join = False  # ?
    
    # GATE ROUTING
    # Create a list to store the connections
    gate_connections = []
    gate_line = LineString(gate_coords)
    # Connect OSS coordinates to each point along the gate line
    for point in gate_line.coords:
        connection_line = LineString([subs_coords, point])
        gate_connections.append(connection_line)

    # Loop over each cluster
    for ic in range(n_cluster):

        # Check for intersection
        # Overwrite new path when there is an intersection
        # Define the first connection
        connection = gate_connections[ic]
        # Find the intersection between the first connection and the gate line
        intersection = connection.intersection(gate_line)

        # Check if there is an intersection
        # Multipoint means, there is another intersection, except the target gate  
        if intersection_join:
            if intersection.geom_type == "MultiPoint":
                # Create new path
                if ic == (len(gate_coords)) and oss_rerouting == 1:
                    # If last gate leads to an intersection
                    connection_new = [subs_coords, gate_coords[ic-1], gate_coords[ic]]   
    
                elif ic >= len(gate_coords) - 1:
                    connection_new = [subs_coords, gate_coords[ic-1], gate_coords[ic]] 
                
                else:    
                    connection_new = [subs_coords, gate_coords[ic+1], gate_coords[ic]]
                # Create a new LineString with the updated coordinates
                new_line = LineString(connection_new)
     
                '''
                # Second interation - check if new line is also intersecting
                intersection = connection.intersection(new_line)
                if intersection.geom_type == "MultiPoint":
                    # Create new path
                    if ic == range(len(gate_coords)):
                        # If last gate leads to an intersection
                        connection_new = [subs_coords, gate_coords[ic-2], gate_coords[ic-1], gate_coords[ic]]   
                    else:    
                        connection_new = [subs_coords, gate_coords[ic+2], gate_coords[ic+1], gate_coords[ic]]
                    # Create a new LineString with the updated coordinates
                    new_line = LineString(connection_new)
                '''
                
                
                # Replace gate connection with new line
                gate_connections[ic] = new_line
    """      
    

    # Make cable layout plot
    if plot == 1:
        plotCableLayout(iac_dic, turb_coords, subs_coords, save=False) 
    
    # cable_id = np.array([a['cable_id'] for a in iac_dic])
    # ia = np.array([a['turbineA_glob_id'] for a in iac_dic])
    # ib = np.array([a['turbineB_glob_id'] for a in iac_dic])
    # cid =np.array([a['cluster_id'] for a in iac_dic])
    
    return iac_dic, global_edge_list, cable_types

   
# ----- Cluster turbines -----
def clusteringSpectral(turb_coords, subs_coords, n_cluster, n_tcmax):
    ''' Clustering wind turbines based on their angles from a single 
    substation using Spectral Clustering.
    
    Input:
    self.turb_coords                : turbines coordinates
    self.subs_coords                 : offshore substation coordinates
    self.n_cluster                  : amount of clusters
    n_tcmax                         : max amount of turbines per cluster
    
    Output:
    self.cluster_arrays             : list with turbine coordinates per cluster
    self.turb_cluster_id             : array with cluster ID per turbine location
    
    https://scikit-learn.org/stable/modules/clustering.html#spectral-clustering      
    '''                       
    # ----- Clustering with Spectral clustering
    # Output: labels
    # Calculate vectors from root to each point
    vectors = turb_coords - subs_coords
    # Calculate angles (in radians) between vectors and x-axis
    angles = np.arctan2(vectors[:, 1], vectors[:, 0])
    # Rescale angles to [0, 2*pi]
    angles[angles < 0] += 2 * np.pi
    # Reshape angles to column vector for clustering
    angles = angles.reshape(-1, 1)
    
    # Calculate Euclidean distance from each point to the root
    # Clustering using spectral with angles as features
    spectral_clustering = SpectralClustering(n_clusters=n_cluster, random_state = 0, affinity='nearest_neighbors', n_neighbors=n_tcmax)
    spectral_clustering.fit(angles)
 
    # ----- Cluster labels
    turb_cluster_id = spectral_clustering.labels_
    # ----- Number of turbines per cluster
    cluster_labels_unique, cluster_labels_counts = np.unique(turb_cluster_id, return_counts=True)
    '''
    # -----  Cluster locations array
    cluster_arrays = []
    
    for name in cluster_labels_unique: 
        #name=0
        # Select indices of points per cluster
        cluster_ind = np.where(turb_cluster_id == name)[0] 
        cluster_points = turb_coords[cluster_ind]
        cluster_arrays.append(cluster_points) 
    '''
    return turb_cluster_id, cluster_labels_counts


def getclusterGates(turb_coords, subs_coords, turb_cluster_id):      
    ''' Get gates of turbines cluster, meaning the closest turbines to oss from each cluster.
    Input:
    turb_coords    : turbine coordinates - list of x,y pairs
    turb_cluster_id : cluster ID of each turbine
    subs_coords     : substation coordinates - list of x,y pairs
    
    Output:
    gate_coords    : list of gate coordinates per cluster 
    gate_index     : index of the turbine that is the gate per cluster
    '''         
    cluster_names = np.unique(turb_cluster_id)
    gate_coords = np.zeros((len(cluster_names),2))
    gate_index = np.zeros(len(cluster_names))
    
    for i in cluster_names :
        # Get locations in current cluster
        cluster_ind = np.where(turb_cluster_id == i)[0] 
        cluster_points = turb_coords[cluster_ind]
        # Calculate distances from OSS to Turb locations for cluster
        distances = np.linalg.norm(cluster_points - subs_coords, axis=1)
        # Find the index of the turbine with the minimum distance to OSS
        gate_index0 = np.argmin(distances)
        # Get the closest location to OSS for current cluster
        gate_coords[i,:] = cluster_points[gate_index0]
        gate_index[i] = gate_index0
    return gate_coords, gate_index


def minimum_spanning_tree(points, start_index):
    '''Find edges that form a minimum spanning tree of the provided node
    points, starting from a specified node. 
    X are edge weights of fully connected graph.
    This function is adapted from the 'Simplistic Minimum Spanning Tree in Numpy'
    from Andreas Mueller, 2012. 
    https://peekaboo-vision.blogspot.com/2012/02/simplistic-minimum-spanning-tree-in.html
    If only one point is provided, an empty list will be returned.
    
    Input:
    points           : List of turbine coordinate x,y pairs
    start_index      : index of which point to start at, which corresponds to
                       the turbine that will be attached to the substation.
    
    Output:
    spanning_edges : list of lists
        Collection of node pairs for each edge, where in each [a,b] pair, a 
        is the ID of the node closer to the substation.
    '''
    
    X = squareform(pdist(points))

    n_vertices = X.shape[0]
    spanning_edges = []
    
    # initialize with start_index:                                                                                         
    visited_vertices = [start_index]                                                                                            
    num_visited = 1
    # exclude self connections:
    diag_indices = np.arange(n_vertices)
    X[diag_indices, diag_indices] = np.inf  # set self-distances to infinite to exclude them
    
    while num_visited != n_vertices:
        # define new edge as shortest distance between visited vertices and others
        new_edge = np.argmin(X[visited_vertices], axis=None)
        # 2d encoding of new_edge from flat, get correct indices                                                      
        new_edge = divmod(new_edge, n_vertices)
        new_edge = [visited_vertices[new_edge[0]], new_edge[1]]                                                       
        # add edge to tree
        spanning_edges.append(new_edge)
        visited_vertices.append(new_edge[1])
        # remove all edges inside current tree so they aren't considered for the next new_edge
        X[tuple(visited_vertices), new_edge[1]] = np.inf
        X[new_edge[1], tuple(visited_vertices)] = np.inf
        num_visited += 1
        
    return spanning_edges


def selectCable(required_rated_power, cableTypes):
    '''Selected the cable type from a list that is the smallest option to
    exceed the required rated power.'''
    
    closest_rated_power = float('inf')  # Initialize with positive infinity to find the closest power
    selected_cable = None

    # Iterate through the list and find the closest power that is >= required_rated_power
    for cable_props_dict in cableTypes:
        if cable_props_dict['power'] >= required_rated_power and cable_props_dict['power'] < closest_rated_power:
            
            closest_rated_power = cable_props_dict['power']
            selected_cable = cable_props_dict
    
    if not selected_cable:
        raise Exception(f"No cable provided meets the required rated power of {required_rated_power}.")
        breakpoint()
        
    return selected_cable

def assignSubstationTurbines(turb_coords, sub_coords, max_turbines):
    '''
    Function to split turbines between substations based on which substation a turbine is closest to.
    
    Parameters
    ----------
    turb_coords : array
        Array of turbine x,y coordinates
    sub_coords : array
        Array of substation x,y coordinates
    max_turbines : int
        Maximum number of turbines allowed per substation
        
    Returns
    -------
    turb_subs_id : list
        The index of substation that each turbine should feed to
    '''
    turb_subs_id = np.zeros((len(turb_coords[:,0]))) # array of substations associated with each turbine
    turbs_for_oss = [] # list of turbine ids for each substation
    distlist = [] # list of distances for each turbine from each substation
    noss = len(sub_coords[:,0]) # number of substations

    # create list where each entry is an array of distances from turbine coords to a specific oss coord
    for oo in range(noss):
        turbs_for_oss.append([])
        distlist.append(np.linalg.norm(turb_coords - sub_coords[oo], axis=1))
    
    # find which oss is closest to each point
    for idx in range(len(distlist[0])):
        turb_subs_id[idx] = int(np.argmin([dist[idx] for dist in distlist]))
    # list of turbine ids broken out by substation
    turbs_for_oss = [list(np.where(turb_subs_id==subid)[0]) for subid in range(noss)]
    
    rturbs_for_oss = deepcopy(turbs_for_oss)
    # if an oss has too many turbines, need to switch some turbines to another oss
    overfilled_oss = [oo for oo in range(noss) if len(turbs_for_oss[oo])>max_turbines]
    if len(overfilled_oss)>0:
        # find oss with least number of turbines
        uoss = np.argmin([len(turbs_for_oss[oo]) for oo in range(noss)]) # underfilled oss

        # for each overfilled oss, switch some turbines to the underfilled oss
        for ooss in overfilled_oss:
            turbine_ids = np.array(turbs_for_oss[ooss]) # ids of turbines currently associated with overfilled oss
            # find difference in distance between each turbine and the over- and under-filled oss
            dist_disparity_margin = [distlist[uoss][tidx]-distlist[ooss][tidx] for tidx in turbine_ids]
            # sort list of indices by decreasing distance difference 
            sorted_dist_disp = np.flip(np.argsort(dist_disparity_margin)) 
            rturbs_for_oss[ooss] = list(turbine_ids[sorted_dist_disp[:max_turbines]]) # update overfilled oss turb list with turbines of largest distance disparity
            rturbs_for_oss[uoss].extend(list(turbine_ids[sorted_dist_disp[max_turbines:]])) # add remaining turbines to underfilled oss
            
        # update turb_subs_id
        for oo,ossid in enumerate(rturbs_for_oss):
            for tid in ossid: # tid is the turbine index/id number
                turb_subs_id[tid] = int(oo)
                
    # return vals
    return(turb_subs_id)


"""

# IN WORK => BACKLOG!
# ----- Advanced routing -----
def advancedCableRouting(iac_edges, cluster_arrays, exclusion_coords):
    '''Wrapping method to perform advanced cable routing, considering obstacles
    iac_edges                    : list of array with edge IDs
    cluster_arrays               : List of arrays with turbine coordinates per cluster
    exclusion_coords             : List of arrays with exclusion zone coordinates
    '''
    
    # Check cable intersection
    intersecting_lines, lines = checkCableIntersections(iac_edges, cluster_arrays, exclusion_coords)
    
    nearby_lines = getObstacles(intersecting_lines, lines, buffer_distance=2000)

    obstacles_list = [nearby_lines, lines, exclusion_polygons_sh]
 
    
 
def checkCableIntersections(iac_edges, cluster_arrays, exclusion_coords):
    '''Wrapping method to perform advanced cable routing, considering obstacles
    Input:
    iac_edges                    : list of array with edge IDs
    cluster_arrays               : List of arrays with turbine coordinates per cluster
    exclusion_coords             : List of arrays with exclusion zone coordinates
    
    Output:
    intersecting_indices         : list of array with edge IDs
    
    '''   
    # Exclusion zones
    exclusion = exclusion_coords 
    exclusion_polygons_sh = []  # List to store polygons   
    
    # Create exclusion zone polygons 
    for ie in range(len(exclusion)):
        exclusion_polygon = sh.Polygon(exclusion[ie])
        exclusion_polygons_sh.append(exclusion_polygon)    
               
    # Convert iac_edges and respective coordinates into Shapely LineString objects and identify intersecting lines
    intersecting_indices = []
    # Loop over clusters
    for ic in range(len(iac_edges)):
        edges = iac_edges[ic]
        coords = cluster_arrays[ic]
        
        # Loop over cables in cluster
        for ie in range(len(edges)):
            start, end = edges[ie]
            line = LineString([coords[start], coords[end]])
            
            # Check if the line intersects the exclusion polygon and get iac_edge index
            if line.intersects(exclusion_polygon):
                intersecting_indices.append((ic, ie))
    



    
    # Get insecting edges
    iac_edges[intersecting_indices[0][0]][intersecting_indices[0][1]]
    
    
    cluster_arrays[intersecting_indices[0][0],iac_edges[intersecting_indices[0][0]][intersecting_indices[0][1]]]
    
    # Convert iac_edges and respective coordinates into Shapely LineString objects
    lines = []
    for edges, coords in zip(iac_edges, cluster_arrays):
        for edge in edges:
            start, end = edge
            line = LineString([coords[start], coords[end]])
            lines.append(line)
    
    # Identify lines that intersect with the exclusion polygon
    intersecting_lines = [line for line in lines if line.intersects(exclusion_polygons_sh[0])]
   
    return intersecting_lines, lines 
    
    
def getObstacles(intersecting_lines, lines, buffer_distance):
    '''Wrapping method to perform advanced cable routing, considering obstacles
    intersecting_lines           : list of array shapely lines
    buffer_distance              : distance, integer
    '''        
    combined_buffer = intersecting_lines[0].buffer(buffer_distance)
    for intersecting_line in intersecting_lines[1:]:
        combined_buffer = combined_buffer.union(intersecting_line.buffer(buffer_distance))
    
    # Identify lines that intersect with the buffer
    # Currently cables only, later include mooring lines as well          
    nearby_lines = [line for line in lines if line.intersects(combined_buffer) and line not in intersecting_lines]
    
    return nearby_lines
 
    
    #x,y = nearby_lines[0].coords.xy
    
 
    # Plotting - with nearby lines
    plt.figure(figsize=(10, 10))
    
    # Plot all lines in blue
    for line in lines:
        x, y = line.xy
        plt.plot(x, y, marker='o', color='blue')
    
    # Plot intersecting lines in red
    for intersecting_line in intersecting_lines:
        x, y = intersecting_line.xy
        plt.plot(x, y, marker='o', color='red')
    
    # Plot nearby lines in orange
    for line in nearby_lines:
        x, y = line.xy
        plt.plot(x, y, marker='o', color='orange')
    
    # Plot the combined buffer
    #x, y = combined_buffer.exterior.xy
    #plt.plot(x, y, color='green', linestyle='--')
    
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Lines and Buffer Around Intersecting Line')
    plt.grid(True)
    plt.show()        
     
              
     
        
    # Plotting - different lines only
    plt.figure(figsize=(10, 10))
    for line in lines:
        x, y = line.xy
        plt.plot(x, y, marker='o')
    
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Shapely Lines from Edges and Coordinates')
    plt.grid(True)
    plt.show()

    # Display the results
    for line in lines:
        print(line)   
    
    # Plotting - insecting lines
    plt.figure(figsize=(10, 10))
    
    # Plot all lines
    for line in lines:
        x, y = line.xy
        plt.plot(x, y, marker='o', color='blue')
    
    # Plot intersecting lines in red
    for line in intersecting_lines:
        x, y = line.xy
        plt.plot(x, y, marker='o', color='red')
    
    # Plot the exclusion polygon
    for polygon in exclusion_polygons_sh:
        x, y = polygon.exterior.xy
        plt.plot(x, y, color='green')
    
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Shapely Lines and Exclusion Polygon')
    plt.grid(True)
    plt.show()    
    
"""    


def getUpstreamTurbines(edge_list):
    '''Calculate the number of turbines upstream of each turbine.
    Input:
    edge_list : list of list pairs
        List of the object ids at the ends of each cable [a ,b], where power
        flows from b to a.
    
    Output:
    self.iac_upstreamturb        : upstream turbines per cable
    self.iac_upstreamturb_count  : amount of upstream turbines per cable
    ''' 
    
    if len(edge_list) == 0:
        return []  # if there is only one turbine in the cluster
    
    # Create a directed graph from the iac edge list
    G = nx.DiGraph()
    G.add_edges_from(edge_list)

    # Initialize a list to store neighbors (turbines) of each point
    neighbors_list = []
    # Iterate over each point and find its neighbors until a point has no neighbors
    for point in range(np.max(edge_list) + 1):
        neighbors = bfs_neighbors(G, point)
        # Directly append neighbors which might be a set or None
        neighbors_list.append(neighbors if neighbors is not None else None)
    # Neighbor count
    neighbor_count = [len(neighbors) if neighbors is not None else 0 for neighbors in neighbors_list]  
    iac_upstreamturb_count = [nc for nc in neighbor_count]

    return iac_upstreamturb_count
     
     
# Function to perform BFS traversal
def bfs_neighbors(graph, start_point):
    '''Breadth-First Search. It's a algorithm for searching or traversing tree or graph data structures. 
    The algorithm starts at a chosen node of a graph and explores all of the neighbor nodes at the present 
    depth prior to moving on to the nodes at the next depth level.
    
    Input:
    graph           : network graph
    start_point     : start point
    
    Output:
    neighbors       : list of neighbors of each point       
    '''
    
    neighbors = set()  # Set to store neighbors
    visited = set()  # Set to store visited points
    queue = [start_point]  # Initialize queue with start point
    
    while queue:
        # Dequeue a point from the queue
        current_point = queue.pop(0)
        # Check if current point has neighbors
        if current_point not in visited:
            visited.add(current_point)
            current_neighbors = set(graph.neighbors(current_point))
            neighbors |= current_neighbors  # Union operation to add neighbors
            queue.extend(current_neighbors - visited)  # Add unvisited neighbors to the queue
  
    # Replace empty set with None
    if not neighbors:
        neighbors = None
    
    return neighbors


"""  
#Seems like the below function is equivalent to
#iac_cab2turb_ic = [[el[1], i, el[0]] for i, el in enumerate(cluster_edge_list)]
#but should double check if it has additional functionality:
def getCableID(coords, gate_coord, edge_list, iac_upstreamturb_count):     
    '''Identify cable (edge) number related to turbines.
    Input:
    self.iac_edges               : list of inter array cable edges per cluster 
    self.iac_upstreamturb_count  : amount of upstream turbines per cable
    gate_coord                 : list of gate coordinates per cluster 
    
    Output:
    self.iac_ID                  : Downflow cable ID at each wind turbine
    self.iac_cab2turb            : List with cables and turbines, without 999 (gates)
        
    '''

    # Identify cable (edge) number related to turbines
    # These cables are in flow direction of the respective wind turbine       
    #breakpoint()
    iac_ID = []
    TurbB_ID = []


    cab_id = np.zeros(len(iac_upstreamturb_count), dtype=int)
    turb_id_B = np.zeros(len(iac_upstreamturb_count), dtype=int) 
    
    gate_index = np.where((coords == gate_coord).all(axis=1))[0][0]
    
    # Iterate through the range of points
    for turb_id_A in range(np.min(edge_list), np.max(edge_list) + 1):
        # If gate point, then skip, because a gate point does not have a inner cluster cable
        if turb_id_A  == gate_index:
            cab_id[turb_id_A] = 999
            turb_id_B[turb_id_A] = 200 # Index for substation
        else:    
            connected_edges = []
            edge_neighbor_counts = []
            edge_points = []
        
            # Find edges connected to the point and their respective neighbor counts
            for edge_index, (start, end) in enumerate(edge_list):
                if start == turb_id_A  or end == turb_id_A:
                    connected_edges.append(edge_index)
     
                    # Add neighbor count for the opposite end of the edge
                    target_point = end if start == turb_id_A  else start
                    edge_neighbor_counts.append(iac_upstreamturb_count[target_point])
        
            # Select the edge (cable) that leads to the turbine with the most neighbors
            if edge_neighbor_counts:
                max_neighbors_index = np.argmax(edge_neighbor_counts)
                selected_edge_index = connected_edges[max_neighbors_index]
                cab_id[turb_id_A] = selected_edge_index
                
                # Get Index of turbine B
                cable = edge_list[selected_edge_index]
                turb_id_B[turb_id_A] = cable[cable != turb_id_A]

    #iac_cab2turb = relateCab2Turb(iac_ID, TurbB_ID)
    
    #iac_edges[iac_ID[ic]]
    array = np.column_stack((np.arange(len(cab_id)), cab_id, turb_id_B))
    iac_cab2turb = array[array[:, 1] != 999]
    
    return iac_cab2turb
"""


# ----- Plot wind farm layout -----
def plotCableLayout(iac_dic, turb_coords, subs_coords, gate_connections=[], exclusion_coords=[], save=False):
    '''Plot wind farm Cable layout.

    '''
    
    # combine coordinates for easy plotting of everything
    coords = np.vstack([turb_coords, subs_coords])
      
    # Exclusion zones
    if len(exclusion_coords) > 0:
        exclusion = exclusion_coords 
        exclusion_polygons_sh = []  # List to store polygons   
        
        # Create exclusion polygons 
        for ie in range(len(exclusion)):
            exclusion_polygon = sh.Polygon(exclusion[ie])
            exclusion_polygons_sh.append(exclusion_polygon)

 
    # Set font sizes
    #fsize_legend = 12    # Legend 
    #fsize_ax_label = 12  # Ax Label 
    #fsize_ax_ticks = 12  # Ax ticks
    #fsize_title = 16     # Title 
             
    # Create a colormap and a legend entry for each unique cable section
    # Find unique values
    # Convert dictionary into data frame
    iac_df=pd.DataFrame(iac_dic)

    unique_cables = np.unique([a['conductor_area'] for a in iac_dic])
    colors = plt.cm.viridis(np.linspace(0, 1, len(unique_cables)))  # Create a colormap based on the number of unique sections
    section_to_color = {sec: col for sec, col in zip(unique_cables, colors)}

    
    plt.figure(figsize=(10, 6))
    
    # ----- Lease area boundary
    #shape_polygon = sh.Polygon(self.boundary)
    #x, y =  self.boundary_sh.exterior.xy
    #plt.plot(x, y, label='Boundary', linestyle='dashed', color='black')
    
    # Plot Turbines
    plt.scatter(coords[:-1, 0], coords[:-1, 1], color='red', label='Turbines')
    
    # Annotate each point with its index
    for i in range(coords.shape[0]-1): #, point in enumerate(cluster_arrays[ic]):
        plt.annotate(str(i), coords[i,:], textcoords="offset points", xytext=(0, 10), ha='center')
      
    # Loop over edges / cable ids
    for i in range(len(iac_dic)):

        # Cable selection
        color = section_to_color[iac_dic[i]['conductor_area']]
        label = f"Section {int(iac_dic[i]['conductor_area'])} mm²" if int(iac_dic[i]['conductor_area']) not in plt.gca().get_legend_handles_labels()[1] else ""
        
        ia = iac_dic[i]['turbineA_glob_id']
        ib = iac_dic[i]['turbineB_glob_id']
        
        plt.plot( coords[[ia,ib], 0], coords[[ia,ib], 1], color=color, label=label)
        
        plt.text( np.mean(coords[[ia,ib], 0]), np.mean(coords[[ia,ib], 1]), str(i), fontsize=9, color='black')

        
        # Turbines
       # plt.scatter(cluster_arrays[ic][:, 0], cluster_arrays[ic][:, 1], color='red', label='Turbines')    
        # Plot gate as a diamond marker
        #plt.scatter(self.gate_coords[ic][0], self.gate_coords[ic][1], marker='D', color='green', label='Gate')
    
    """
    ## ----- Cables Gates to OSS
    # TODO: updated cable_id below >>>
    iac_oss = iac_df[iac_df['cable_id'] >= 100]
    iac_array_oss = iac_oss.values
 
    for i in range(n_cluster):
        cable_section_size = int(iac_array_oss[i, 9])  # Assuming cable section size is in the 7th column
        color = section_to_color.get(cable_section_size, 'black')  # Default to black if section size not found
        connection = gate_connections[i]
        x_connection, y_connection = connection.xy
        label = f'Section {cable_section_size} mm²' if cable_section_size not in plt.gca().get_legend_handles_labels()[1] else ""
        plt.plot(x_connection, y_connection, color=color, label=label)
        
        #plt.plot([gate_coords[i][0], x_oss], [gate_coords[i][1], y_oss], color=color, label=f'Section {cable_section_size} mm²' if cable_section_size not in plt.gca().get_legend_handles_labels()[1] else "")
    """
    
    if len(exclusion_coords) > 0:
        for ie in range(len(exclusion)):
            shape_polygon = exclusion_polygons_sh[ie]#sh.Polygon(self.exclusion[i])
            x, y = shape_polygon.exterior.xy
            plt.plot(x, y, linestyle='dashed', color='orange', label='Exclusion Zone')   
        #ax.plot([], [], linestyle='dashed', color='orange', label='Exclusion Zone')
        
    # turbine locations
    #ax.scatter(x0, y0, c='black', s=12, label='Turbines')
         

    # ----- OSS
    plt.scatter(subs_coords[:,0], subs_coords[:,1], label='substation', marker='*', color='black', s=100)


    # Set plot title and labels
    plt.title('Wind Turbine Cluster - Cable Conductor Sizes')
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
       
    # Create a custom legend for the unique cable sections
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # Removing duplicate labels
    plt.legend(by_label.values(), by_label.keys(),loc='upper left', fancybox=True, ncol=2)
    plt.gca().set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal

    # Create a custom legend for the unique cable sections
    #handles, labels = plt.gca().get_legend_handles_labels()
    #by_label = dict(zip(labels, handles))  # Removing duplicate labels
    #sorted_labels = sorted(by_label.keys())  # Sort the labels alphabetically
    #sorted_handles = [by_label[label] for label in sorted_labels]  # Get handles corresponding to sorted labels
    #plt.legend(sorted_handles, sorted_labels, loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, ncol=2)
    #plt.gca().set_aspect('equal', adjustable='box')  # Set aspect ratio to be equal
      
    
    plt.grid(True)
  
    # ----- Save plot with an incremented number if it already exists
    if save:
        counter = 1
        output_filename = f'wind farm layout_{counter}.png'
        while os.path.exists(output_filename):
            counter += 1
            output_filename = f'wind farm layout_{counter}.png'
        
        # Increase the resolution when saving the plot
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')  # Adjust the dpi as needed


# Test Script
if __name__ == '__main__':
    
    
    turb_coords = [[    0, 1000],
                   [    0, 2000],
                   [    0, 3000],
                   [    0, 4000],
                   [    0, 5000],
                   [ 1000,   0],
                   [ 1000, 1000],
                   [ 1000, 2000],
                   [ 1000, 3000],
                   [ 2000, 2000],
                   [ 2000, 3000],
                   [ 2000, 4000],
                   [ 2000, 5000]]
    
    cluster_id = [ 0, 0, 0, 0, 0, 
                   1, 1, 1, 1, 1, 
                   1, 1, 2]
    
    subs_coords = [ 1400, 200]
    
    conductor_sizes = np.array([300, 630, 1000])

    cableProps_type = 'dynamic_cable_66'
    turb_rating_MW = 15
    
    
    #iac_dic = getCableLayout(conductor_sizes, cableProps_type, turb_rating_MW, turb_coords, subs_coords, plot=1)
    iac_dic, connections, types = getCableLayout(turb_coords, subs_coords, conductor_sizes, cableProps_type, turb_rating_MW, turb_cluster_id=[], plot=1)
    
    cable_id = np.array([a['cable_id'] for a in iac_dic])
    ia = np.array([a['turbineA_glob_id'] for a in iac_dic])
    ib = np.array([a['turbineB_glob_id'] for a in iac_dic])
    cid =np.array([a['cluster_id'] for a in iac_dic])
    
    
    # set up a CableSystem!!
    from famodel.cables.cable_system import CableSystem
    cs = CableSystem(turb_coords)
    
    cs.update(connections, types, coords=turb_coords, 
              powers=[15]*len(turb_coords), 
              subcoords=subs_coords)
    
    cs.checkConnectivity()
    
    plt.show()
    