import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import scipy.stats
import copy
import math
import random

from failure.graphBuilder import *
from failure.searchMethods import *
from failure.twoTurbineCaseStudy import *

''' failureProbability -------------------------------------------------------------------------------------------------------------

            ****************************** Last Updated: 23 April 2024 ******************************

 Methods:
 1) transition_natrix: inputs adjacency matrix, probabilities --> output conditional probabilities, scaled probabilities

 2) conditional_probabilities_update: input start node, list of probabilities --> output conditional probabilities

 3) single_simulation: input start node, adjacency martrix, list of node nanes, update boolean --> output graph generated

 4) single_backward_simulation: input start node, adjacency martrix, list of node nanes, update boolean --> output graph generated

 5) monte_carlo_sim: number of iterations, plotting boolean, starting nodes, adjacency matrix, list of node names, random seed
 booleanm, midpoint boolean --> output average probabilities, similarity between average and conditional probabilities

 6) single_conditional_prob: inputs first probability, second probability, midpoint boolean --> output conditional probability

 7) multi_cond_prob: input parent nodes, current nodes, probabilities, midpoint boolean --> output conditional probability

 8) bayesian_table: input adjacency matrix, current node, midpoint boolean, node names, alternate probabilties boolean,
 alternate probabilities, multiple turbines boolean --> output parents list, probability distribution table

 9) write_bayesian_probabilities: input adjacency matrix, node names, probabilties, multiple turbines boolean, midpoint boolean calculation,
 number of iterations --> writes probability distribution tables to Excel file (nothing returned)

10) probability_over_time: input failure rate, time --> output probability of a failure at time t

11) bayesian_inference: input adjacency matrix, node names, indicator nodes, evidence nodes, hypothesis nodes, probabilties, hypothesis
boolean, printing boolean, multiple turbine boolean, parent or child indicator, twoTurbine boolean --> outputs inference probabilities

12) backward_bayesian_inference: input adjacency matrix, list of node names, array of nodes to start with, array of evidence nodes, 
 array of hypothesis nodes, probabilties, boolean for indexing into the right values of the array, boolean for multiple turbines,
 parent or child indicator, boolean for twoTurbine method --> output two arrays of inference probabilities (one normalized and one not)
 
------------------------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------------------------


    transition_natrix Documentation
 ----------------------
 This method inputs an adjacency matrix and vector of probabilities (corresponding to the probabilitiy that each
 failure will happen). Since the data we have does not tell us the probability that any two event will happen, we
 find the lower and upper bounds of what the conditional probability of every pair of related events (related indicated
 by the adjacency matrix). We then take the median of these values as the probability that the second failure will
 occur given that the first already occured. This method returns the array of these conditional probabilities and
 a scaled array of these probabilities (scaled so that the sum of the columns is 1).'''

def transition_matrix(arr, probabilities, mid_point = True):
    arr = make_binary(arr) # Array to show which failures lead to other failures
    nodes = diagonal_nodes(arr) # Matrix with numbers on the diagonal

    bounds = [] # Initialize list of upper and lower bounds
    transitions = np.zeros(arr.shape) # Initialize matrix of conditional probabilities

    for node in range(arr.shape[0]): # For each node...
        children_bool = arr[node] @ nodes # Find the children of the node
        children = children_bool[np.nonzero(children_bool)]

        parent_probability = probabilities[node][0] # Find the probabilitiy of the node occuring
        lower_bound = [] # Initialize lower bounds of all children of the node
        upper_bound = [] # Initialize upper bounds of all children of the node

        for child in children: # For each of the node's children...
            child_probability = probabilities[child - 1][0] # Find the probability of the child occuring

            # The lower bound is the probability of the child occuring (i.e. child and parent are completely indepedent evens)
            lb = child_probability 
             # The upper bound is the overlap of the child and parent probabilities (i.e. child and parent are completely dependent events)
            ub = min(parent_probability, child_probability)/parent_probability

            if mid_point: 
                tm_val = (lb + ub) / 2 # Calculate the midpoint of the two bounds
                # print("Midpoint = True") --> for debugging
            else: 
                rand_val = np.random.rand()
                tm_val = rand_val * (ub - lb) + lb
                # print(rand_val, tm_val, tm_val == ((lb + ub) / 2), lb, ub) --> for debugging
                # print("Midpoint = False") --> for debugging
            
            transitions[node][child - 1] =  tm_val
            lower_bound.append(lb) # Append the lower bound to the list
            upper_bound.append(ub) # Append the upper bound to the list
        
        lower_bound = np.array(lower_bound) # Convert lower and upper bounds to numpy arrays
        upper_bound = np.array(upper_bound)

        midpoint = (upper_bound + lower_bound) / 2 # Calculate the midpoints for all the children of the parent node

        bounds.append([lower_bound, upper_bound, midpoint]) # Append the lower, upper, and midpoints to the list

    transition_matrix = transitions / np.sum(transitions, axis = 0) # Scale the conditional probabilities
    transition_matrix[np.isnan(transition_matrix)] = 0 # Replace any NaN values with 0s
    return transitions, transition_matrix # Return conditional probabilities and scaled conditional probabilities

''' conditional_probabilities_update Documentation
 --------------------------------------------------
 This method inputs a start value (indicating the starting node) and a list of probabilities that each node will
 occur. We then compute conditional probabilities of all the nodes given the starting node. This method returns
 the list of calculated conditional probabilities.'''
    
def conditional_probabilities_update(start_arr, probabilities):
    # Create vector of the parent probabilities
    for start in start_arr:
        parent_probability = np.reshape(np.repeat(probabilities[start][0], probabilities.shape[0]), probabilities.shape) 

        # Calculate the conditional probabilities (using the midpoint method described in transition_matrix() method)
        conditional_probabilities = (probabilities + np.reshape(np.min(np.hstack((probabilities, parent_probability)), axis = 1), probabilities.shape)/parent_probability[0])/2
        probabilities = conditional_probabilities

    return conditional_probabilities # Return the calculated probabilities

''' bayesian_probabilities Documentation
 --------------------------------------------------
 This method inputs a start value (indicating the starting node) and a list of probabilities that each node will
 occur. We then compute conditional probabilities of all the nodes given the starting node. This method returns
 the list of calculated conditional probabilities.'''
    
def bayesian_probabilities(parents, child, probabilities):
    for i in range(len(parents)):
        rand_val = np.random.rand()



''' single_simulation Documentation
 ---------------------------------------------
 This method inputs a starting node, adjacency martrix, list of node nanes, and a boolean that either uses the method of
 updating conditional probabilities based on already visited nodes. We then compute the failure propagation and graph the
 propagation via Networkx. This method returns the graph generated.'''
    
def single_simulation(start_arr, arr, nodeNames, update = False, plot = False, rand_seed = True, mid_point = True):
    monte_carlo_array = np.zeros(arr.shape)
    probs = np.zeros((arr.shape[0], 1))

    # Use the probabilities from the COREWIND failure rates
    probabilities = np.array([0.0195, 0.0195, 0.013625, 0.0055, 0.0175, 0.2075, 0.001, 0.001, 0.001, 0.093185, 0.001, 0.001,
                        0.027310938, 0.033968125, 0.033968125, 0.01375, 0.01375, 0.01375, 0.01375, 0.01375, 0.01375,
                        0.0205, 0.0205, 0.02, 0.01, 0.01, 0.233, 0.288, 0.543374, 0.1285, 0.01, 0.01, 0.01, 0.015, 0.0155,
                        0.015, 0.0155, 0.015, 0.0155, 0.015, 0.33, 0.025, 0.025, 0.025, 0.025, 0.025, 0.105]) #0.01375, 
    probabilities = np.reshape(probabilities, (arr.shape[0], 1)) # Reshape these probabilities into a vector
    
    G = nx.DiGraph() # Initialize a graph and a list of modes and effects for plotting nodes in different colors
    effects = []
    modes = []

    if update: # If you want to update the probabilities, update the probabilities given the indicated starting node
        probabilities = conditional_probabilities_update(start, probabilities)
    transitions, tm = transition_matrix(arr, probabilities, mid_point) # Compute the probabilities for all linked nodes

    if rand_seed: random.seed(16) # Use a random seed to get the same results (yet random results) each time

    adj = make_binary(arr).astype(int) # Make a binary adjacency matrix to determine the children of each node
    nodes = diagonal_nodes(arr) # Diagonal array of node indices + 1

    nodeList = np.reshape(np.repeat(False, arr.shape[0]), (arr.shape[0], 1)) # Create a list ot determine which nodes have been visited
    gens = {}
    queue = []

    for start in start_arr:
        G.add_node(str(0) + ": " + str(nodeNames[start-1]))
        nodeList[start-1][0] = False
        queue.append([start, 0, 0])
        gens.update({str(0) + ": " + str(nodeNames[start-1]): {"layer": 0}})
        # Determine if the source node is an effect or a mode, and add it to the correct array
        '''if source < 49: modes.append(nodeNames[source - 1])
        else: effects.append(nodeNames[source - 1])'''
        if start < 27: effects.append(nodeNames[start - 1]) # Add the start node to either the effects or mode list
        else: modes.append(nodeNames[start - 1])

    while len(queue) > 0: # For each node in the queue...
        current = queue[0] # Get the node from the front of the queue
        children_bool = adj[current[0]-1] @ nodes # vector of zeros and child names (numerical names)
        children = children_bool[np.nonzero(children_bool)] #list of just the child names (numerical names)

        for child in children: # For each child of the current node...
            if nodeList[child - 1] == False or gens[nodeNames[child - 1]]['layer'] > current[1]: # If the node has not been visited...
                if np.random.rand() < transitions[current[0] - 1][child - 1]: # If the random value is less than the probability...
                    monte_carlo_array[current[0] - 1][child - 1] = 1
                    probs[child - 1] = 1
                    G.add_node(nodeNames[child-1]) # Add the child to the graph
                    if update: # If updateing is desired, update the probabilities with the addition of the child node
                        probabilities = conditional_probabilities_update(current[0]-1, probabilities)
                        transitions, tm = transition_matrix(arr, probabilities, mid_point)

                    # Determine if the child is an effect or mode and add it to the correct array
                    if child < 27: effects.append(nodeNames[child - 1])
                    else: modes.append(nodeNames[child - 1])

                    # Add an edge between the current node and child in the tree we are building
                    G.add_edge(nodeNames[current[0]-1], nodeNames[child-1])

                    queue.append([child, current[1]+1]) # Append the child to the queue
                    nodeList[child - 1] = True # Change the status of the child to say we have visited it
                    gens.update({nodeNames[child-1]: {"layer": current[1]+1}}) # Add the child to the layer dictionary
        queue = queue[1:] # Remove the current node from the queue

    nx.set_node_attributes(G, gens) # Set the layers of the graph

    # Plot the graph
    if plot:
        pos = nx.multipartite_layout(G, subset_key='layer')
        nx.draw_networkx_nodes(G, pos, nodelist=effects, node_color="#98c5ed", node_size=2700, edgecolors="#799dbd", node_shape="s")
        nx.draw_networkx_nodes(G, pos, nodelist=modes, node_color="#fabc98", node_size=2700, edgecolors="#c89679", node_shape="s")
        nx.draw_networkx_labels(G, pos, font_size=10, verticalalignment='center_baseline')
        nx.draw_networkx_edges(G, pos, arrowsize=60)
        plt.box(False)
        plt.show()

    # draw_bfs_multipartite(arr, nodeNames, start, "child") # Plot the graph when all probabilities are one for comparison
    return G, monte_carlo_array, probs # Return the graph




''' single_backward_simulation Documentation
 ---------------------------------------------
 This method inputs a starting node, adjacency martrix, list of node nanes, and a boolean that either uses the method of
 updating conditional probabilities based on already visited nodes. We then compute the backward failure propagation and graph the
 propagation via Networkx. This method returns the graph generated.'''
    
def single_backward_simulation(start_arr, arr, nodeNames, update = False):
    # Use the probabilities from the COREWIND failure rates
    probabilities = np.array([0.0195, 0.0195, 0.013625, 0.0055, 0.0175, 0.2075, 0.001, 0.001, 0.001, 0.093185, 0.001, 0.001,
                        0.027310938, 0.033968125, 0.033968125, 0.01375, 0.01375, 0.01375, 0.01375, 0.01375, 0.01375,
                        0.0205, 0.0205, 0.02, 0.01, 0.01, 0.233, 0.288, 0.543374, 0.1285, 0.01, 0.01, 0.01, 0.015, 0.0155,
                        0.015, 0.0155, 0.015, 0.0155, 0.015, 0.33, 0.025, 0.025, 0.025, 0.025, 0.025, 0.105]) #0.01375, 
    probabilities = np.reshape(probabilities, (arr.shape[0], 1)) # Reshape these probabilities into a vector
    
    G = nx.DiGraph() # Initialize a graph and a list of modes and effects for plotting nodes in different colors
    effects = []
    modes = []
    queue = []
    gens = {}

    if update: # If you want to update the probabilities, update the probabilities given the indicated starting node
        probs = conditional_probabilities_update(start, probabilities)
    else: probs = probabilities
    transitions, tm = transition_matrix(arr.T, probs) # Compute the probabilities for all linked nodes

    # random.seed(16) # Use a random seed to get the same results (yet random results) each time --> feel free to comment out

    adj = make_binary(arr).astype(int) # Make a binary adjacency matrix to determine the children of each node
    nodes = diagonal_nodes(adj) # Create a matrix with node names to determine the children of each node
    nodeList = np.reshape(np.repeat(False, arr.shape[0]), (arr.shape[0], 1)) # Create a list ot determine which nodes have been visited
    
    for start in start_arr:
        G.add_node(str(0) + ": " + str(nodeNames[start-1]))
        nodeList[start-1][0] = False
        queue.append([start, 0, 0])
        gens.update({str(0) + ": " + str(nodeNames[start-1]): {"layer": 0}})
        # Determine if the source node is an effect or a mode, and add it to the correct array
        '''if source < 49: modes.append(nodeNames[source - 1])
        else: effects.append(nodeNames[source - 1])'''
        if start < 27: effects.append(nodeNames[start - 1]) # Add the start node to either the effects or mode list
        else: modes.append(nodeNames[start - 1])

    while len(queue) > 0: # For each node in the queue...
        current = queue[0] # Get the node from the front of the queue
        children_bool = adj[current[0]-1] @ nodes # vector of zeros and child names (numerical names)
        children = children_bool[np.nonzero(children_bool)] #list of just the child names (numerical names)

        for child in children: # For each child of the current node...
            if nodeList[child - 1] == False or gens[nodeNames[child - 1]]['layer'] < current[1]: # If the node has not been visited...
                if np.random.rand() < transitions[current[0] - 1][child - 1]: # If the random value is less than the probability...
                    G.add_node(nodeNames[child-1]) # Add the child to the graph
                    if update: # If updateing is desired, update the probabilities with the addition of the child node
                        probs = conditional_probabilities_update(current[0]-1, probs)
                        transitions, tm = transition_matrix(arr.T, probs)

                    # Determine if the child is an effect or mode and add it to the correct array
                    if child < 27: effects.append(nodeNames[child - 1])
                    else: modes.append(nodeNames[child - 1])

                    # Add an edge between the current node and child in the tree we are building
                    G.add_edge(nodeNames[child-1], nodeNames[current[0]-1])

                    queue.append([child, current[1]-1]) # Append the child to the queue
                    nodeList[child - 1] = True # Change the status of the child to say we have visited it
                    gens.update({nodeNames[child-1]: {"layer": current[1]-1}}) # Add the child to the layer dictionary
        queue = queue[1:] # Remove the current node from the queue

    nx.set_node_attributes(G, gens) # Set the layers of the graph

    # Plot the graph
    # pos = nx.multipartite_layout(G, subset_key='layer')
    # nx.draw_networkx_nodes(G, pos, nodelist=effects, node_color="#98c5ed", node_size=2700, edgecolors="#799dbd", node_shape="s")
    # nx.draw_networkx_nodes(G, pos, nodelist=modes, node_color="#fabc98", node_size=2700, edgecolors="#c89679", node_shape="s")
    # nx.draw_networkx_labels(G, pos, font_size=10, verticalalignment='center_baseline')
    # nx.draw_networkx_edges(G, pos, arrowsize=60)
    # plt.box(False)
    # plt.show()
    return G, nx.to_numpy_array(G), probs # Return the graph




''' monte_carlo_sim Documentation
 ---------------------------------------------
 This method inputs the number of iterations, a boolean for plotting the average graph or not, an array of the nodes to start with,
 an adjacency matrix, list of nodeNames, boolean for a random seed, and a boolean for using a midpoint calculation. We then generate
 a graph with failure probabilities for the number of iterations and find the average graph and the probability that each node
 is in the graph (the number of times the node shows up divided by the number of iterations). Lastly, we calculate the similarity
 between the probability of the nodes in the graph (that we just calculated) compared to the estimated probability calculated via
 conditional probabilities. We return the first list of probabilities and cosine similarity between the two lists of probabilities.'''
    
def monte_carlo_sim(num_iterations, plot, start_arr, adjacency_matrix, nodeNames, rand_seed, mid_point):
    # Initialize the adjacency matrix to trrack which nodes are in each graph (which we will then calculate an average from)
    adj_matrices = np.zeros(adjacency_matrix.shape)

    # Initialize array for probabilities to track average probabilities for each node
    probs = np.zeros((adjacency_matrix.shape[0], 1))

    # List of probabilities for each failure happening
    probabilities = np.array([0.0195, 0.0195, 0.013625, 0.0055, 0.0175, 0.2075, 0.001, 0.001, 0.001, 0.093185, 0.001, 0.001,
                                0.027310938, 0.033968125, 0.033968125, 0.01375, 0.01375, 0.01375, 0.01375, 0.01375, 0.01375,
                                0.0205, 0.0205, 0.02, 0.01, 0.01, 0.233, 0.288, 0.543374, 0.1285, 0.01, 0.01, 0.01, 0.015, 0.0155,
                                0.015, 0.0155, 0.015, 0.0155, 0.015, 0.33, 0.025, 0.025, 0.025, 0.025, 0.025, 0.105]) #0.01375, 
    probabilities = np.reshape(probabilities, (adjacency_matrix.shape[0], 1))

    # Run each simulaiton
    for i in range(num_iterations):
        arr = copy.deepcopy(adjacency_matrix)
        G, adj_mat, prob = single_simulation(start_arr, arr, nodeNames, rand_seed = rand_seed, mid_point = mid_point)
        adj_matrices += adj_mat.astype(float)
        probs += prob
        # print(i+1) # Debugging --> feel free to uncomment

    # Calculate average graph and average probabilities
    adj_matrices = adj_matrices/num_iterations
    probs = probs/num_iterations
    K, abc = matrix2Graph(adj_matrices, nodeNames, True)

    # Calculate similarity between average and conditional probabilities
    v1 = conditional_probabilities_update(start_arr, probabilities)
    v2 = probs

    # Plot average graph
    if plot:
        draw_bfs_multipartite(adj_matrices, nodeNames, start_arr, "child")
    return v2, cosine_similarity(v1, v2) # Return average probabilities and similarity of average and conditional probabilities




''' single_conditional_prob Documentation
 ---------------------------------------------
 This method inputs the probability of first failure, probability of second failure, and boolean for using the midpoint calculation
 or not. We use these probabilities to estimate the conditional probability between the two failures. We return this conditional probability.'''
    
def single_conditional_prob(pa, pb, mid_point = True):
    lb = pb * pa # Lower bound for conditional probability
    ub = min(pa, pb) # Upper bound for conditional probability
    if mid_point: 
        overlap = (lb + ub) / 2 # Find midpoint between the lower and upper bounds
    else: 
        rand_val = np.random.rand()
        overlap = rand_val * (ub - lb) + lb # Find a random value between the lower and upper bounds
    return overlap/pa # Return the conditional probability




''' multi_cond_prob Documentation
 ---------------------------------------------
 This method inputs array of parent nodes, current nodes, array of probabilities, and boolean for using a midpoint calculation or not. We
 calculate the conditional probability when the node has multiple parents (i.e. probability of current event given the parent events). This
 method outputs the conditional probability.'''
    
def mult_cond_prob(parents, current, probabilities, midpoint = True):
    if len(parents) < 1:
        return probabilities[current - 1] # If no parents, return the probability of the current node
    
    # Create list of parent node probabilities, with the probability of the current node appended to the end
    parents_sub = [int(parent - 1) for parent in parents]
    parents_sub.append((current - 1))
    parents_sub = np.array(parents_sub)
    parent_probs = probabilities[parents_sub]
    # Initialize denomenator
    denomenator = 1

    # Calculate the single conditional probability for pairs of events until there are no events left
    while len(parent_probs) > 1:
        if (not isinstance(parent_probs[1][0], float)) and len(parent_probs[1][0]) > 1:
            parent_probs[1][0] = parent_probs[1][0][1]
        lb = parent_probs[0][0] * parent_probs[1][0] # Lower bound for conditional probability
        ub = min(parent_probs[0][0], parent_probs[1][0]) # Upper bound for conditional probability
        if midpoint: 
            val = (lb + ub) / 2 # Find midpoint between the lower and upper bounds
        else: 
            rand_val = np.random.rand()
            val = rand_val * (ub - lb) + lb
        parent_probs = np.vstack(([val], parent_probs[2:]))
        if len(parent_probs) == 2:
            denomenator = parent_probs[0][0]
    # print(parent_probs[0][0], denomenator, parent_probs[0][0]/denomenator)
    return parent_probs[0][0]/denomenator # Return the conditional probability




''' bayesian_table Documentation
 ---------------------------------------------
 This method inputs adjacency matrix, current node, boolean for midpoint calculation, list of node names, boolean for alternate probabilties,
 array of new probabilities, and boolean for calculations among multiple turbines. We calculate the probability distribution and format it
 into a table. This method outputs the list of parents and the probability distribution table.'''
    
def bayesian_table(adj, current, midpoint, nodeNames, pvs = False, prob_vals=None, mult_turbine = False):
    nodes = diagonal_nodes(adj) # Diagonal matrix of node's numerical names +1
    adj = make_binary(adj, 0.5) # Binarize adjacency matrix

    parents_bool = nodes @ adj[:, current-1] # vector of zeros and parent names (numerical names)
    parents = list(parents_bool[np.nonzero(parents_bool)]) # list of just the parent names (numerical names)
    prob_table = np.zeros((2 ** len(parents), len(parents) + 2)) # Initialize the array for the node's conditional probability distribution

    # *** Code in development for multiple turbine bayesian tables ***
    arr_nonbinary = adj.copy()
    num_parents = len(parents)

    '''if mult_turbine:
        for p in range(num_parents):
            if arr_nonbinary[parents[p] - 1][current - 1] > 1:
                print("Unlocked!", parents[p])
                parents.append(parents[p])'''
    # *** End of code in development *********************************

    for iteration in range(2 ** len(parents)): # For each combination of parents having either failed or not
        if pvs: # Determine the probabilities being used
            probabilities = prob_vals.copy()
        else:
            probabilities = np.array([0.0195, 0.0195, 0.013625, 0.0055, 0.0175, 0.2075, 0.001, 0.001, 0.001, 0.093185, 0.001, 0.001,
                                    0.027310938, 0.033968125, 0.033968125, 0.01375, 0.01375, 0.01375, 0.01375, 0.01375, 0.01375,
                                    0.0205, 0.0205, 0.02, 0.01, 0.01, 0.233, 0.288, 0.543374, 0.1285, 0.01, 0.01, 0.01, 0.015, 0.0155,
                                    0.015, 0.0155, 0.015, 0.0155, 0.015, 0.33, 0.025, 0.025, 0.025, 0.025, 0.025, 0.105]) #0.01375, 
        probabilities = np.reshape(probabilities, (len(probabilities), 1))
        # print(probabilities)

        true_false_array = [] # Iniitalize to record which parents have failed

        for p in range(len(parents)):
            probabilities[int(parents[p]-1)] = abs((int(iteration / (2 ** p)) % 2) - probabilities[int(parents[p]-1)]) # Determine parent probability (given if the parent has failed or not)
            true_false_array.append((int(iteration / (2 ** p)) % 2)) # Append whether or not the parent node has failed
            prob_table[iteration][p] = int(iteration / (2 ** p)) % 2 # Append whether or not the parent node has failed

            # if mult_turbine and nodeNames[int(current - 1)][:3] != nodeNames[int(parents[p] - 1)][:3]: #and p > 0) and parents[p] <= parents[p-1]:
                # probabilities[int(parents[p]-1)] *= 0.7 # If there is more than one turbine, then decrease the probability of failure by 30%
        prob = mult_cond_prob(parents, current, probabilities, midpoint) # Calculate the conditional probability of the node given if the parents have failed or not
        # print(true_false_array, prob)
        # print(probabilities)
        prob_table[iteration][-2] = prob # Add calculated probability to the array
        prob_table[iteration][-1] = 1 - prob # Add the probability of the node not failing to the array
    #print(prob_table)
    return parents, prob_table # Return a list of the parents and the probability distribution array




''' write_bayesian_probabilities Documentation
 ----------------------------------------------
 This method inputs adjacency matrix, list of node names, probabilties, boolean for calculations among multiple turbines, 
 boolean for midpoint calculation, and number of iterations. We calculate the probability distribution for each node in the graph and format it
 into a table. This method writes the probability distribution tables to an Excel file (nothing returned).'''
    
def write_bayesian_probabilities(adjacency_matrix, nodeNames, probabilities, mult_turbine, midpoint=True, num_iterations=1, twoTurbines = False, name=""):
    ps = probabilities.copy()
    with pd.ExcelWriter("bayesianProbabilities"+name+".xlsx") as writer: # Write to file
        for current in range(1, adjacency_matrix.shape[0]+1): # Repeat for each node in the graph
            if any(probabilities != ps):
                print("Unequal probabilities", current-1)
                break
            adjacency_matrix2 = adjacency_matrix.copy()
            parents, our_table = bayesian_table(adjacency_matrix2, current, midpoint, nodeNames, True, probabilities, mult_turbine) # Probability distrubution table
            if twoTurbines:
                parents, our_table = twoTurbine_bayesian_table(adjacency_matrix, adjacency_matrix, current, nodeNames, nodeNames) # Probability distrubution table
            parents = [int(parent) for parent in parents] # Format parent list

            # Find average probability distribution table for the number of iterations
            for i in range(num_iterations - 1):
                parents, our_table2 = bayesian_table(adjacency_matrix, current, midpoint, nodeNames = nodeNames, pvs=True, prob_vals=probabilities, mult_turbine=mult_turbine)
                our_table += our_table2

            # Update column titles for the probability distribution table
            parents = [nodeNames[parent - 1].replace("\n", " ") for parent in parents]
            parents = np.append(parents, "P(" + nodeNames[current - 1].replace("\n", " ") + " = True)")
            parents = np.append(parents, "P(" + nodeNames[current - 1].replace("\n", " ") + " = False)")
            df = pd.DataFrame(our_table, columns = parents)

            df.to_excel(writer, sheet_name=nodeNames[current - 1].replace("\n", " ").replace("/", " or ").replace(":", "")[:31]) # Write to file
            print(current, nodeNames[current - 1].replace("\n", " ")[:31]) # Print current node to keep track of which nodes have a probability distribution table




''' probability_over_time Documentation
 ----------------------------------------------
 This method inputs a failure rate, lamda, and a time, t. We return the probability of a failure at time t.'''
    
def probability_over_time(lamda, t):
    return 1 - math.exp((np.log(1-lamda)) * t) # Return the updated probability of failure given a certain time




''' bayesian_inference Documentation
 ----------------------------------------------
 This method inputs adjacency matrix, list of node names, array of nodes to start with for generating the graph of the Bayesian network, array of 
 evidence nodes, array of hypothesis nodes, probabilties, boolean for finding hte probability of failure of the hypothesis, boolean for printing 
 information about the calculations, boolean for multiple turbines, parent or child calculation (string), and boolean for using the twoTurbine case
 study method. We calculate Bayesian inference given the evidence and hypothesis variables. This method outputs an array of inference probabilities.'''
    
def bayesian_inference(arr, nodeNames, indicator, num_evidence, num_hypothesis, probabilities, tf = True, printing = False, multi= False, poc="parent", twoTurbine = False):
    # print("E", num_evidence)
    # print("H", num_hypothesis)
    a = arr.copy()
    non = nodeNames
    prblts = probabilities
    evidence = num_evidence
    hypothesis = num_hypothesis
    if not multi:
        K, a, g, e, m, non = breadth_first_multi(a, nodeNames, indicator, poc) # Generate tree for Bayesian network
        # draw_bfs_multipartite(arr, nodeNames, indicator, type = "multi-"+poc, multi_turbine = False)
        prblts = [] # Initialize array of node probabilities (in order of appearance in graph)
        for node in non:
            node_index = np.where(nodeNames == node)[0][0]
            prblts.append(probabilities[node_index]) # Add nodes to array of node probabilities
        prblts = np.array(prblts)

        evidence = []
        hypothesis = []
        for hypothesis_node in num_hypothesis: # Adjust hypothesis node so that index in tree matches up with index in original adjacency matrix
            if twoTurbine:
                non = np.array(non)
            if nodeNames[hypothesis_node - 1] not in non: # If node is not in tree, then there is a 0% probability of failure
                print("Probability table of", nodeNames[hypothesis_node - 1], "...", np.reshape(np.array([0,1]), (1,2)))
                return np.reshape(np.array([0,1]), (1,2)), np.reshape(np.array([0,1]), (1,2))
            hypothesis.append(np.where(non == nodeNames[hypothesis_node - 1])[0][0])

        for evidence_node in num_evidence: # Adjust evidence node so that index in tree matches up with index in original adjacency matrix
            if twoTurbine:
                non = np.array(non)
            if nodeNames[evidence_node - 1] not in non: # If node is not in tree, then this is an error!
                print("ERROR - evidence node, " + nodeNames[evidence_node - 1] + ", not in graph!")
            if evidence_node in num_hypothesis:# If node is in hypothesis, then probability of failure is 100%
                print("Probability table of", nodeNames[evidence_node - 1], "...", np.reshape(np.array([1,0]), (1,2)))
                return np.reshape(np.array([1,0]), (1,2)), np.reshape(np.array([1,0]), (1,2))
            evidence.append(np.where(non == nodeNames[evidence_node - 1])[0][0])

    probabilitiy_table = np.zeros((2, a.shape[0])) # Initialize table of inference probabilities
    nodes = diagonal_nodes(a) # Diagonal matrix of node names (numerical +1)
    a = make_binary(a, 0.5) # Binarize adjacency table
    hypothesis_nodes_array = np.zeros((len(hypothesis), 2))

    # Depending on tree (forward or backward propagation), either iterate through nodes forward or backwards
    if poc == "parent":
        this_range = reversed(range(a.shape[0]))
    elif poc == "child":
        this_range = range(a.shape[0])

    for node in this_range:
        pts_bool = nodes @ a[:, node] # vector of zeros and child names (numerical names)
        pts = pts_bool[np.nonzero(pts_bool)] #list of just the child names (numerical names)

        if len(pts) < 1: # If no parents, the probability is the initial probability
            probabilitiy_table[0][node] = prblts[node]
            probabilitiy_table[1][node] = 1 - prblts[node]
            continue

        parents, our_table = bayesian_table(a, node+1, True, nodeNames, True, prblts) # Calculate the probability distribution table

        # If only using two turbines (specific to case study), calculate using twoTurbine method
        if twoTurbine:
            parents, our_table = twoTurbine_bayesian_table(a, arr, node + 1, nodeNames, non) # Calculate the probability distribution table
        mlt_table = np.ones((our_table.shape[0],2)) # Initialize table for multiplying across rows of probability distribution table
        
        for i in range(our_table.shape[0]):
            for j in range(our_table.shape[1] - 2):
                parent = int(parents[j])

                # Replace boolean in probability distribution table with probability of this event
                if our_table[i,j] == 0:
                    # print(a.shape, probabilitiy_table.shape)
                    our_table[i,j] = probabilitiy_table[0][parent - 1]

                     # If the node is in the hypothesis or evidence array, do not include this 
                    if probabilitiy_table[0][parent - 1] == 0:
                        # print("indexing_error!!", parent, node)
                        # print(non[parent], non[node])
                        break
                    if (parent-1 in hypothesis and tf): # or (parent in evidence):
                        our_table[i,j] = 0
                        # print("in hypto and tf or in evidence", parent)
                else:
                    our_table[i,j] = probabilitiy_table[1][parent - 1]

                     # If the node is in the hypothesis or evidence array, do not include this 
                    if (parent-1 in hypothesis and not tf) or (parent-1 in evidence):
                        our_table[i,j] = 0
                        # print("in hypto and not tf or in evidence", parent)

                mlt_table[i,0] *= our_table[i,j] # Multiply the probabilities across the probability distribution table
            mlt_table[i,1] = mlt_table[i,0] * our_table[i, -1] # Multiple by the probability of event not failing given combination of parent failure
            mlt_table[i,0] *= our_table[i, -2] # Multiple by the probability of event failing given combination of parent failure

        sm_table = np.sum(mlt_table, axis = 0) #/np.sum(mlt_table) # Sum the products of probabilities across the columns

        if node in hypothesis:
            print("Probability table of", non[node].replace("\n", " "), "...", sm_table)
            hypothesis_nodes_array[np.where(np.array(hypothesis) == node)[0][0]][0] = sm_table[0]
            hypothesis_nodes_array[np.where(np.array(hypothesis) == node)[0][0]][1] = sm_table[1]
            if all(hypothesis_nodes_array != 0):
                return probabilitiy_table, hypothesis_nodes_array

        probabilitiy_table[0][node] = sm_table[0] # Update the inference probability table with the probabilites just calculated
        probabilitiy_table[1][node] = sm_table[1]

    if printing: # Print the inference probability, along with the indicators, evidence, and hypothesis conditionals
        print("Indicator:", nodeNames[np.array(indicator)])
        print("Evidence:", nodeNames[np.array(evidence)])
        print("Hypothesis:", nodeNames[np.array(hypothesis)])
        print("Probability of Failure:", probabilitiy_table[:, 0][0], probabilitiy_table[:, 0][0]/np.sum(probabilitiy_table[:, 0]))
        print("Probability of No Failure:", probabilitiy_table[:, 0][1], probabilitiy_table[:, 0][1]/np.sum(probabilitiy_table[:, 0]))
    return probabilitiy_table, hypothesis_nodes_array # Return array or inference probabilities




''' backward_bayesian_inference Documentation
 ----------------------------------------------
 This method inputs adjacency matrix, list of node names, array of nodes to start with for generating the graph of the Bayesian network, array of 
 evidence nodes, array of hypothesis nodes, probabilties, boolean for indexing into the right values of the array, boolean for multiple turbines, parent or child
 calculation type (string), and boolean for using twoTurbine case study method. We calculate Bayesian inference given the evidence and hypothesis variables. 
 This method outputs an array of inference probabilities (normalized and not).'''
    
def backward_bayesian_inference(adjacency_matrix, nodeNames, indicator, evidence, hypothesis, probabilities, start_bool = True, multi = False, poc="parent", twoTurbine = False):
    # Calculate Bayesian inference for hypothesis = True and hypothesis = False
    #print("--- First Run of Bayesian Inference ----------------")
    # print("Before bn", probabilities.shape)
    printing = False
    if len(evidence) > 0 and len(hypothesis)>0:
        printing = False

    pt1, hnd1 = bayesian_inference(adjacency_matrix, nodeNames, indicator, evidence, hypothesis, probabilities, True, printing=printing, multi = multi, poc=poc, twoTurbine=twoTurbine)
    #print()
    #print("--- Second Run of Bayesian Inference ---------------")
    if poc == "parent":
        pt2, hnd2 = bayesian_inference(adjacency_matrix, nodeNames, indicator, evidence, hypothesis, probabilities, False, printing=printing, multi = multi, poc=poc, twoTurbine=twoTurbine)
    # print("pt1", pt1[:, 0])
    # print("pt2", pt2[:, 0])
    # Choose correct indexing value
    index = 1
    if start_bool:
        index = 0

    if multi:
        return [hnd1[0][0]/(hnd1[0][0]+hnd1[0][1]), hnd1[0][1]/(hnd1[0][0]+hnd1[0][1])], [0,0]

    # Depending on what our evidence and hypothesis nodes are, index into the correct values
    if poc == "parent" and 0 in evidence:
        p_true = pt1[:, 0][index]
        p_false = pt2[:, 0][index]

    elif poc == "parent" and 0 in hypothesis:
        p_true = pt1[:, 0][0] + pt2[:, 0][0]
        p_false = pt1[:, 0][1] + pt2[:, 0][1]

    elif poc == "child":
        if hnd1.shape[0] == 1:
            p_true = hnd1[0][0]
            p_false = hnd1[0][1]

        else: # If multiple nodes in hypothesis, multiply all possibilities together to get normailized value
            print("Starting hypothesis for-loops...")
            p_false = 0
            for iteration in range(2 ** hnd1.shape[0]):
                p_false_mult = 1
                for node in range(hnd1.shape[0]):
                    p_false_mult *= hnd1[node][int(iteration / (2 ** node)) % 2]
                if iteration == 0:
                    p_true = p_false_mult
                else:
                    p_false += p_false_mult
                if (iteration) % 1000000000000 == 0:
                    print("Percent done:", iteration/(2 ** hnd1.shape[0]) * 100)
            print("Percent done:", 100)
    else:
        p_true = np.sum(pt1[:, 0])
        p_false = np.sum(pt2[:, 0])

    probability_distribution = [p_true/(p_true + p_false), p_false/(p_true + p_false)] # Normalize the probability distribution of the event happening

    return probability_distribution, [p_true, p_false] # Return the normalized probability distribution and unnormalized distribution



'''
# ----------- For running the code: feel free to uncomment -------------------
arr, nodeNames = excel2Matrix("failureData.xlsx", "bigMatrix")
#random.seed(18)
start = 30 #np.random.randint(0, arr.shape[0])
probabilities = np.array([0.0195, 0.0195, 0.013625, 0.0055, 0.0175, 0.2075, 0.001, 0.001, 0.001, 0.093185, 0.001, 0.001,
                            0.027310938, 0.033968125, 0.033968125, 0.01375, 0.01375, 0.01375, 0.01375, 0.01375, 0.01375,
                            0.0205, 0.0205, 0.02, 0.01, 0.01, 0.233, 0.288, 0.543374, 0.1285, 0.01, 0.01, 0.01, 0.015, 0.0155,
                            0.015, 0.0155, 0.015, 0.0155, 0.015, 0.33, 0.025, 0.025, 0.025, 0.025, 0.025, 0.105]) #0.01375, 
probabilities = np.reshape(probabilities, (arr.shape[0], 1))
start = 44 #np.random.randint(1, arr.shape[0]+1)
single_simulation(start, arr, nodeNames)
single_simulation(start, arr, nodeNames, update=True)

num_iterations = 100
plot = True
start = 11 #np.random.randint(1, arr.shape[0]+1)
monte_carlo_sim(num_iterations, plot, start, adjacency_matrix, nodeNames, rand_seed=False, mid_point=False)'''




'''# ----------- For running the code: feel free to uncomment -------------------
adjacency_matrix, nodeNames = excel2Matrix("ExcelFiles/failureData.xlsx", "bigMatrix")
probabilities = np.array([0.0195, 0.0195, 0.013625, 0.0055, 0.0175, 0.2075, 0.001, 0.001, 0.001, 0.093185, 0.001, 0.001,
                            0.027310938, 0.033968125, 0.033968125, 0.01375, 0.01375, 0.01375, 0.01375, 0.01375, 0.01375,
                            0.0205, 0.0205, 0.02, 0.01, 0.01, 0.233, 0.288, 0.543374, 0.1285, 0.01, 0.01, 0.01, 0.015, 0.0155,
                            0.015, 0.0155, 0.015, 0.0155, 0.015, 0.33, 0.025, 0.025, 0.025, 0.025, 0.025, 0.105]) #0.01375, 
probabilities = np.reshape(probabilities, (adjacency_matrix.shape[0], 1))
targets = [44]
starts = [16, 33, 35, 36]

# C, D = bayesian_table(adjacency_matrix, 16, True, pvs = True, prob_vals=probabilities, mult_turbine = False)
# A, B = backward_bayesian_inference(adjacency_matrix, nodeNames, [0], [0], [44], probabilities, start_bool = True)
# print(B)
write_bayesian_probabilities(adjacency_matrix, nodeNames, probabilities, mult_turbine = True, midpoint = True, num_iterations=1)'''