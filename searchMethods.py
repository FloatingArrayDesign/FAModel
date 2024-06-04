import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from graphBuilder import *

''' searchMethods--------------------------------------------------------------------------------------------------

            ****************************** Last Updated: 10 April 2024 ******************************

 Methods:
 1) breadth_first_child: inputs adjacency matrix, names of nodes, start node --> outputs tree, adjacency matrix, 
 dictionary of nodes and generations, list of 'effect' nodes, list of 'mode' nodes

 2) breadth_first_parent: inputs adjacency matrix, names of nodes, start node --> outputs tree, adjacency matrix, 
 dictionary of nodes and generations, list of 'effect' nodes, list of 'mode' nodes

 3) draw_bfs_multipartite: nputs adjacency matrix, list of node names, start node --> plots breadth first tree.
 Nothing is returned.

 4) breadth_first_multi: adjacency matrix, list of node names, starting node, type of breadth first search -->
  outputs tree, adjacency matrix, dictionary of nodes, arrays for effects and modes, and array of nodes
 
--------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------------


   breadth_first_child Documentation
 -----------------------------------
 This method takes in an adjacency matrix, list of node names, and an integer indicating te starting node. We 
 traverse through the graph. For each generation, starting from the source node, we find and add all the children
 of the nodes in this generation, unless the node has aldready been placed in the tree. As such, we create a tree. 
 This method returns this tree, the tree's adjacency matrix, a dictionary of nodes, and arrays for effects and modes.
 
 --> Note: a tree is a graph with no cycles or loops. It does not have to be directed, but it is acyclic.'''

def breadth_first_child(arr, nodeNames, source):
    # Initialize a new digraph for the tree we are going to make and add all the nodes (we will determine edges later)
    G = nx.DiGraph()

    # Initialize effects and modes arrays for plotting purposes
    effects = []
    modes = []

    # Add the source node to the graph
    G.add_node(nodeNames[source - 1]) 

    # Determine if the source node is an effect or a mode, and add it to the correct array
    if source < 27: effects.append(nodeNames[source - 1])
    elif source <= 47: modes.append(nodeNames[source - 1])

    '''if source < 49: modes.append(nodeNames[source - 1])
    else: effects.append(nodeNames[source - 1])'''

    # Note that we are changing the numerical names of the nodes so that the values of nodes are from 1 to 47 rather
    # than from 0 to 46. This is so we can find all the children of the nodes. Using the 0-46 range poses a problem
    # since we do not know if a zero indicates no relation or that the 0 node is a child of the node we are looking at.
    # We binarize the adjacency matrix so that relationships either exist or not (rather than having an intensity)
    adj = make_binary(arr, 0.5).astype(int)

    # Create an array of booleans such that the ith entry indicates if the node has already been visited. Nodes
    # that have been visited are set to TRUE.
    nodeList = np.reshape(np.repeat(False, arr.shape[0]), (arr.shape[0], 1))

    # Create a diagonal matrix such that the value of the i,i entry is 1+1, referencing the node with name i+1
    nodes = diagonal_nodes(adj)
    # print("nodes", nodes.shape)

    # Visit the source node
    nodeList[source-1] = True

    # Add the source node to the queue. We use a queue to iterate through the nodes. This allows us to search through
    # the graph generation by generation rather than following one specific path at a time.
    queue = [[source, 0]]

    # Initialize a dictionary that will tell us what generation each node is in. Label the generation of the source node
    # as zero.
    gens = {nodeNames[source-1]: {"layer": 0}}

    # Continue while there are still nodes in the queue (reference algorithms for a breadth first search for more info)
    while len(queue) > 0:
        # Set the node we are looking at equal to the node next in line in the queue
        current = queue[0]

        # Determine the children of the current node. 
        # Proof/Reason Why This Works --> Using the unweighted adjacency matrix, we can find the children
        # by multiplying the row of the adjacency matrix corresponding to the current node by the diagonal matrix of
        # node names. If the jth node is a child of the current node, then there is a 1 in the j-1 column of the
        # current row. When multiplied by the diagonal matrix, this 1 will be multiplied by the j-1 column of the.
        # since the only non zero entry in the j-1 column of the diagonal matrix is j (located on the diagonal), then
        # the entry in the j-1 column of the returned vector will be j. However, if the jth node is not a child of
        # the current node, then there is a zero in the j-1 column of the current row of the adjacency matrix. When
        # multiplied by the diagonal matrix, this zero is multiplied by every entry in the j-1 column of the diagonal
        # matrix. So, the j-1 column of the returned vector is zero.

        children_bool = adj[current[0]-1] @ nodes # vector of zeros and child names (numerical names)
        children = children_bool[np.nonzero(children_bool)] #list of just the child names (numerical names)

        for child in children: # For every child of the current node that was found above...
            # Check if the child has been visited or if it is in the same generation as child. If either not visited or
            # in the same generation, continue with the following:
            if nodeList[child - 1] == False or gens[nodeNames[child - 1]]['layer'] > current[1]:

                # Add the child to the graph
                G.add_node(nodeNames[child-1])

                # Determine if the child is an effect or mode and add it to the correct array
                if child < 27: effects.append(nodeNames[child - 1])
                elif child <= 47: modes.append(nodeNames[child - 1])
                '''if child < 49: modes.append(nodeNames[child - 1])
                else: effects.append(nodeNames[child - 1])'''

                # Add an edge between the current node and child in the tree we are building
                G.add_edge(nodeNames[current[0]-1], nodeNames[child-1])

                queue.append([child, current[1]+1]) # Append the child to the queue
                nodeList[child - 1] = True # Change the status of the child to say we have visited it
                gens.update({nodeNames[child-1]: {"layer": current[1]+1}}) # Add child to dictionary of nodes

        # Remove the current node from the queue
        queue = queue[1:]

    # Return the tree we made, its adjacency matrix, the dictionary of nodes, the effects array, and the modes arrat
    return G, nx.to_numpy_array(G), gens, effects, modes



''' breadth_first_parent Documentation
 ---------------------------------------
 This method takes in an adjacency matrix, list of node names, and an integer indicating te starting node. We 
 traverse through the graph. For each generation, starting from the source node, we find and add all the parents
 of the nodes in this generation, unless the node has aldready been placed in the tree. As such, we create a tree. 
 This method returns this tree, the tree's adjacency matrix, a dictionary of nodes, and arrays for effects and modes.'''

def breadth_first_parent(arr, nodeNames, source):
    # Initialize a new digraph for the tree we are going to make and add all the nodes (we will determine edges later)
        # Initialize a graph in Networkx
    G = nx.DiGraph()

    # Initialize effects and modes arrays for plotting purposes
    effects = []
    modes = []
    names_of_nodes = []

    # Add the source node to the graph
    G.add_node(nodeNames[source - 1]) 
    names_of_nodes.append(nodeNames[source - 1])

    # Determine if the source node is an effect or a mode, and add it to the correct array
    '''if source < 27: effects.append(nodeNames[source - 1])
    else: modes.append(nodeNames[source - 1])'''
    if source < 49: modes.append(nodeNames[source - 1])
    else: effects.append(nodeNames[source - 1])

    # Binarize the adjacency matrix
    arr2 = arr.copy()
    adj = make_binary(arr2).astype(int)

    # Create an array of booleans such that the ith entry indicates if the node has already been visited. Nodes
    # that have been visited are set to TRUE.
    nodeList = np.reshape(np.repeat(False, arr.shape[0]), (arr.shape[0], 1))

    # Create a diagonal matrix such that the value of the i,i entry is 1+1, referencing the node with name i+1
    nodes = diagonal_nodes(adj)

    # Visit the source node
    nodeList[source-1] = True

    # Add the source node to the queue. We use a queue to iterate through the nodes. This allows us to search through
    # the graph generation by generation rather than following one specific path at a time.
    queue = [[source, 100]]

    # Initialize a dictionary that will tell us what generation each node is in. Label the generation of the source node
    # as zero.
    gens = {nodeNames[source-1]: {"layer": 100}}

    # Continue while there are still nodes in the queue (reference algorithms for a breadth first search for more info)
    while len(queue) > 0:
    # Set the node we are looking at equal to the node next in line in the queue
        current = queue[0]

        # Determine the parents of the current node. 
        parents_bool = nodes @ adj[:, current[0]-1] # vector of zeros and parent names (numerical names)
        parents = parents_bool[np.nonzero(parents_bool)] #list of just the parent names (numerical names)

        for parent in parents: # For every child of the current node that was found above...
            # Check if the parent has been visited or if it is in the same generation as parent. If either not visited or
            # in the same generation, continue with the following:
            if nodeList[parent - 1] == False or gens[nodeNames[parent - 1]]['layer'] < current[1]:

                # Add the parent to the graph
                G.add_node(nodeNames[parent-1])

                # Determine if the parent is an effect or mode and add it to the correct array
                if parent < 27: effects.append(nodeNames[parent- 1])
                else: modes.append(nodeNames[parent - 1])
                '''if parent < 49: modes.append(nodeNames[parent- 1])
                else: effects.append(nodeNames[parent - 1])'''

                # Add an edge between the current node and parent in the tree we are building
                G.add_edge(nodeNames[parent-1], nodeNames[current[0]-1])

                queue.append([parent, current[1]-1]) # Append the parent to the queue
                nodeList[parent - 1] = True # Change the status of the parent to say we have visited it
                gens.update({nodeNames[parent-1]: {"layer": current[1]-1}}) # Add parent to dictionary of nodes

        # Remove the current node from the queue
        queue = queue[1:]

    # Return the tree we made, its adjacency matrix, the dictionary of nodes, the effects array, and the modes arrat
    return G, nx.to_numpy_array(G), gens, effects, modes



''' draw_bfs_multipartite Documentation
 ---------------------------------------
 This method inputs an adjacency matrix, list of node names, and a start node. We use the breadth first searh to 
 find generations of nodes starting from the start node. This method plots the resulting graph from the breadth
 first search. Nothing is returned.'''

def draw_bfs_multipartite(arr, nodeNames, start, type = "child", multi_turbine = False):
    # Obtain the graph, adjacency matrix, dictionary, effects, and modes from breadth_first_tree
    if type == "child":
        G, arr, gens, effects, modes = breadth_first_child(arr, nodeNames, start)
    elif type == "parent":
        G, arr, gens, effects, modes = breadth_first_parent(arr, nodeNames, start)
    elif type == "multi-child":
        G, arr, gens, effects, modes, non = breadth_first_multi(arr, nodeNames, start, "child")
    else:
        G, arr, gens, effects, modes, non = breadth_first_multi(arr, nodeNames, start, "parent")

    # Give each node the attribute of their generation
    nx.set_node_attributes(G, gens)
    # print(effects)

    if multi_turbine:
        effect_colors = ["#ffd6ed", "#ffb3ba", "#ffdfba", "#ffffba", "#baffc9", "#bae1ff", "#b1adff", "#e4adff", "#e5e5e5", "#e8d9c5"]
        mode_colors = ["#e5c0d5", "#e5a1a7", "#e5c8a7", "#e5e5a7", "#a7e5b4", "#a7cae5", "#9f9be5", "#cd9be5", "#cecece", "#d0c3b1"]
        pos = nx.multipartite_layout(G, subset_key='layer')
        for node in G.nodes:
            # print(int(str(node)[0]))
            if str(node) in effects:
                nx.draw_networkx_nodes(G, pos, nodelist=[node], node_color = effect_colors[int(str(node)[0])], node_size=750, node_shape="s")
            else:
                nx.draw_networkx_nodes(G, pos, nodelist=[node], node_color = effect_colors[int(str(node)[0])], node_size=750)
        nx.draw_networkx_labels(G, pos, font_size=5, verticalalignment='center_baseline')
        nx.draw_networkx_edges(G, pos, arrowsize=20)
        plt.box(False)
        plt.show()
    else:
        # Plot the graph
        pos = nx.multipartite_layout(G, subset_key='layer')
        nx.draw_networkx_nodes(G, pos, nodelist=effects, node_color="#98c5ed", node_size=2700, edgecolors="#799dbd", node_shape="s")
        nx.draw_networkx_nodes(G, pos, nodelist=modes, node_color="#fabc98", node_size=2700, edgecolors="#c89679", node_shape="s")
        nx.draw_networkx_labels(G, pos, font_size=10, verticalalignment='center_baseline')
        nx.draw_networkx_edges(G, pos, arrowsize=60)
        plt.box(False)
        plt.show()

    # Nothing is returned
    return



'''breadth_first_multi Documentation
 -----------------------------------
 This method takes in an adjacency matrix, list of node names, an integer indicating te starting node, and a string
 that indicates which type of breadth first search we are conducting (child or parent). We traverse through the graph. 
 For each generation, starting from the source nodes (handles multiple start nodes), we find and add all the children/parents
 of the nodes in this generation, unless the node has aldready been placed in the tree. As such, we create a tree. 
 This method returns this tree, the tree's adjacency matrix, a dictionary of nodes, arrays for effects and modes, and array
 of nodes (in the order they are added).'''

def breadth_first_multi(arr, nodeNames, sources, type_poc):
    # Function that determines how value of the layer changes, depending on the type of calculation
    def layer_fcn(layer, type_poc):
        if type_poc == "child":
            return layer + 1
        elif type_poc == "parent":
            return layer - 1
    # Determining if the current generation is after the former node's generation, depending on the type of calculation
    def same_layer(layer1, layer2, type_poc):
        if type_poc == "child":
            if layer1 > layer2:
                return True
            else:
                return False
        elif type_poc == "parent":
            if layer1 < layer2:
                return True
            else:
                return False

    # Initialize a new digraph for the tree we are going to make and add all the nodes (we will determine edges later)
    G = nx.DiGraph()

    # Initialize arrays for tracking nodes
    effects = []
    modes = []
    names_of_nodes = []
    adj = make_binary(arr).astype(int)
    nodeList = np.reshape(np.repeat(False, arr.shape[0]), (arr.shape[0], 1))
    nodes = diagonal_nodes(adj)
    queue = []
    gens = {}

    # Initialize value for the first layer created
    if type_poc == "child":
        layer_val = 0
    else:
        layer_val = 100

    # Add the source node to the graph
    for start in sources:
        G.add_node(nodeNames[start - 1]) 

        # Determine if the source node is an effect or a mode, and add it to the correct array
        if start < 27: effects.append(nodeNames[start - 1])
        else: modes.append(nodeNames[start - 1])

        '''if start < 49: modes.append(nodeNames[start - 1])
        else: effects.append(nodeNames[start - 1])'''

        # Visit the source node
        nodeList[start-1] = True
        names_of_nodes.append(nodeNames[start - 1])

        # Add the source node to the queue. We use a queue to iterate through the nodes. This allows us to search through
        # the graph generation by generation rather than following one specific path at a time.
        queue.append([start, 0])

        # Initialize a dictionary that will tell us what generation each node is in. Label the generation of the source node
        # as zero.
        gens.update({nodeNames[start-1]: {"layer": layer_val}})
        # print("non", names_of_nodes)
    
    # Continue while there are still nodes in the queue (reference algorithms for a breadth first search for more info)
    while len(queue) > 0:
        # Set the node we are looking at equal to the node next in line in the queue
        current = queue[0]

        if type_poc == "child":
        # Determine the children of the current node. 
            children_bool = adj[current[0]-1] @ nodes # vector of zeros and child names (numerical names)
            children = children_bool[np.nonzero(children_bool)] #list of just the child names (numerical names)
        else:
            children_bool = nodes @ adj[:, current[0]-1] # vector of zeros and parent names (numerical names)
            children = children_bool[np.nonzero(children_bool)] #list of just the parent names (numerical names)


        for child in children: # For every child of the current node that was found above...
            # Check if the child has been visited or if it is in the same generation as child. If either not visited or
            # in the same generation, continue with the following:
            # print("child", child)

            if nodeList[child - 1] == False or same_layer(gens[nodeNames[child - 1]]['layer'], current[1], type_poc):

                # Add the child to the graph, and to the array of node names
                G.add_node(nodeNames[child-1])
                if nodeNames[child - 1] in names_of_nodes:
                    x = 14
                else:
                    names_of_nodes.append(nodeNames[child - 1])

                # Determine if the child is an effect or mode and add it to the correct array
                '''if child < 27: effects.append(nodeNames[child - 1])
                else: modes.append(nodeNames[child - 1])'''
                if child < 49: modes.append(nodeNames[child - 1])
                else: effects.append(nodeNames[child - 1])

                # Add an edge between the current node and child in the tree we are building
                if type_poc == "child":
                    G.add_edge(nodeNames[current[0]-1], nodeNames[child-1])
                else:
                    G.add_edge(nodeNames[child-1], nodeNames[current[0]-1])

                queue.append([child, layer_fcn(current[1], type_poc)]) # Append the child to the queue
                nodeList[child - 1] = True # Change the status of the child to say we have visited it
                gens.update({nodeNames[child-1]: {"layer": layer_fcn(current[1], type_poc)}}) # Add child to dictionary of nodes

        # Remove the current node from the queue
        queue = queue[1:]

    # Return the tree we made, its adjacency matrix, the dictionary of nodes, the effects array, and the modes arrat
    return G, nx.to_numpy_array(G), gens, effects, modes, names_of_nodes


'''# **** Below code is still being worked on *****

arr, nodeNames = excel2Matrix("failureData.xlsx", "bigMatrix")
source = 11
G = nx.DiGraph()

effects = []
modes = []

G.add_node(nodeNames[source - 1]) 

if source < 49: modes.append(nodeNames[source - 1])
else: effects.append(nodeNames[source - 1])

adj = make_binary(arr).astype(int)

nodeList = np.reshape(np.repeat(False, arr.shape[0]), (arr.shape[0], 1))
nodes = diagonal_nodes(adj)

nodeList[source-1] = True
queue = [[source, 0]]
gens = {nodeNames[source-1]: {"layer": 0}}

while len(queue) > 0:
    current = queue[0]

    children_bool = adj[current[0]-1] @ nodes # vector of zeros and child names (numerical names)
    children = children_bool[np.nonzero(children_bool)] #list of just the child names (numerical names)

    for child in children:
        if nodeList[child - 1] == False or gens[nodeNames[child - 1]]['layer'] > current[1]:
            G.add_node(nodeNames[child-1])
            
            if child < 27: effects.append(nodeNames[child - 1])
            else: modes.append(nodeNames[child - 1])

            G.add_edge(nodeNames[current[0]-1], nodeNames[child-1])

            queue.append([child, current[1]+1]) # Append the child to the queue
            nodeList[child - 1] = True # Change the status of the child to say we have visited it
            gens.update({nodeNames[child-1]: {"layer": current[1]+1}}) # Add child to dictionary of nodes

    queue = queue[1:]'''