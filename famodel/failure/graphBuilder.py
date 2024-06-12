import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

''' graphBuilder---------------------------------------------------------------------------------------------------

            ****************************** Last Updated: 19 February 2024 ******************************

 Methods:
 1) excel2Matrix: inputs fileName, sheetName, boolean of nodeLabels --> output adjacency matrix and name array

 2) matrix2Graph: inputs adjacency matrix, name of nodes, boolean of nodeLabels --> output graph and adjacency 
 matrix

 3) export2graphml: inputs G, string (file name) --> writes graphML file with inputted file name (no return)

 4) get_node_dict: inputs adjacency matrix and name of nodes --> outputs array containing dictionaries with each
 entry containing the node's name, number, and whether it is an effect or mode.

 5) get_node_array: inputs adjacency matrix and name of nodes --> outputs array with dictionary entries of the 
 names of the nodes

 6) get_graph_edges: inputs adjacency matrix, string indicating type of output desired, and threshold --> outputs 
 an array of edges

 7) make_binary:  inputs adjacency matrix and float type --> outputs binarized adjacency matrix

 8) plot_graph: inputs graph and type (str) --> no output (prints graph)

 9) max_degrees: inputs adjacency matrix, names of nodes, threshold for binarization of adjacency matrix, boolean
 for showing names --> outputs tuples of max out degree, max in degree, and max overall degree with the nodes that
 have these max degrees

10) getProbabilities: inputs fileName, sheetName, boolean of nodeLabels --> output probabilities and name array

11) diagonal_nodes: inputs adjacency matrix --> outputs nodes diagonal matrix

12) cosine_similarity: inputs vector one and vector two --> outputs cosine similarity of vector one and two

13) make_undirected_unweighted: inputs adjacency matrix and (optional) threshold --> outputs undirected, unweighted
 version of the adjacency matrix

------------------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------------------


   excel2Matrix Documentation
 -------------------------------
 This code is meant to take in an excel file (already properly formatted by theh user) and output the corresponding
 graph. This assumes that the contents of the excel file is an adjacency matrix. The sheetName input identifies which
 sheet in the excel file the adjacency matrix is in. The nodeLabels input lets the user choose if they would like 
 words (input TRUE) or numeric (input FALSE) labels on the graph. The plot input lets the user toggle on/off plotting
 of the graph. Finally, the graph is returned.'''

def excel2Matrix(fileName, sheetName = None):
    # Read in excel file as a dataframe
    df = pd.read_excel(fileName, sheet_name=sheetName)

    # Convert the data frame to a numpy array
    arr = df.to_numpy()

    # Create an array of the names of the nodes for the graph.
    nodeNames = arr[:, 0].flatten()
    
    # return the adjacency matrix without labels and an array of the names of the nodes
    return arr[:, 1:], nodeNames



''' matrix2Graph Documentation
 -------------------------------
 This method takes in an adjacency matrix, the names of the nodes, whether or not to display the node names, and 
 whether or not to display the plot of the graph once constructed. The output is the graph and adjacency matrix.'''

def matrix2Graph(arr, nodeNames, nodeLabels = False):
    # Create the graph (DiGraph because the graph is directed and weighted) and add the nodes selected above
    G = nx.DiGraph()
    if not nodeLabels:
        nodeNames = range(arr.shape[0])
    G.add_nodes_from(nodeNames)

    # The following nested for loops find the edges of the graph and add them to the graph using the 'add_weighted_edges_from()'
    # command. 'i' indexes through the rows and 'j' indexes through the columns. If there is no edge (i.e. 0 in the adjacency
    # matrix), then move on to the next set of nodes. Otherwise, add the edge to the graph.
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if arr[i,j] == 0:
                continue
            else:
                if nodeLabels: G.add_weighted_edges_from([(nodeNames[i], nodeNames[j], arr[i,j])])
                else: G.add_weighted_edges_from([(i, j, arr[i,j])])

    # Return the graph and adjacency matrix
    return G, arr



''' export2graphML Documentation
 -------------------------------
 This method takes in a networkx graph and ouputs a graphML file.'''

def export2graphML(G, name):

    # create a name for the graphML file that you are going to write
    filename = name + "_FOWT_graphML.graphml"

    # wrtie the graphML file using networkx methods
    nx.write_graphml(G, filename, encoding='utf-8', prettyprint=True, infer_numeric_types=True, named_key_ids=True, edge_id_from_attribute=None)
    return



''' get_node_dict Documentation
 -------------------------------
 This method takes in the adjacency matrix and the names of all the nodes. It outputs an array of node dictionaries.
 Each key/value pair contains the node's number, name, and 'class_id' (which refers to if the node represents a 
 failure effect or failure mode).'''

def get_node_dict(arr, nodeNames):
    # Initialize array for the nodes
    nodes = []

    # Iterate through each node number
    for index in range(arr.shape[0]):

        # If the number is less than 26 (value chosen based on excel sheets), then the node is a failure effect.
        if index < 26:
            class_id = 'effect'

        # Otherwise, the node must be a failure mode.
        else:
            class_id = 'mode'

        # Add a dictionary with the node's number, name, selection criteria (in this case they can all be selected),
        # and class_id (whether they are an effect or mode).
        nodes.append({
                'data': {'id': str(index), 'label': nodeNames[index]},
                'selectable': True,
                'classes': class_id
            })
        
    # Return the array of dictionaries for the nodes
    return nodes



''' get_node_array Documentation
 -------------------------------
 This method takes in an adjacency matrix and list of node names. It outputs an array with each entry a
 key/value pair. Each key is the string 'label', but the values are the string names of the nodes.'''

def get_node_array(arr, nodeNames):
    # Initialize an array for the nodes
    nodes = []

    # Iterate through each node number
    for index in range(arr.shape[0]):

        # For every node, locate its name and append a dictionary with 'label' as the only key and
        # the name just acquired as the only value.
        nodes.append({'label': nodeNames[index]})

    # Return the array of nodes
    return nodes


''' get_grpah_edges Documentation
 -------------------------------
 This method takes in an adjacency matrix (array type) and a string indicating the type of edge list
 desired. This method outputs an array of all the edges. If the 'dictionary' type is indicated, then each
 entry in the array is formatted as
        {'data': {'id': name_of_source_and_target, 'source' name_of_source, 'target', name_of_target}}.
 If the 'array' type is indicated, then each entry in the array is formatted as
        [source, target]).
 Lastly, if the 'tuple' type is indicated, then each entry in the array is formatted as
        (source, target).'''

def get_graph_edges(arr, type, threshold):
    # Initialize an array for the edges
    edges = []

    # Iterate through each row of the adjacency matrix
    for i in range(arr.shape[0]):

        # For each row, iterate through each column of the matrix
        for j in range(arr.shape[1]):

            # If the entry in the adjacency matrix is greater than 0,5 (can change this threshold), then
            # add an edge to the array of edges. Each edge is added by appending an array containing the
            # source and target.
            if arr[i, j] > threshold:
                if type == 'dict':
                    edges.append({'data': {'id': str(i)+str(j), 'source': str(i), 'target': str(j)},
                              'classes': 'red'})
                elif type == 'array':
                    edges.append([i,j])
                elif type == 'tuple':
                    edges.append((i,j))

            # Otherwise, if the input does not indicate an edge, move on to the nect entry in the
            # adjacency matrix
            else:
                continue

    # Return the array of edges
    return edges



''' make_binary Documentation
 -------------------------------
 This method takes in an adjacency matrix (array typle) and a threshold value (float type). It outputs a
 new adjacency matrix such that all the entries are either 0 or 1.'''

def make_binary(arr, threshold = 0):
    # Initialize the new adjacency matrix by setting it equal to the original adjacency matrix
    arr_altered = arr

    # Iterate through each row of the matrix
    for i in range(0, arr_altered.shape[0]):

        # In each row, iterate through each column of the matrix
        for j in range(arr_altered.shape[0]):

            # If the entry is greater than the indicated threshold, set the entry of the new adjacency
            # matrix equal to 1
            if arr[i,j] > threshold:
                arr_altered[i,j] = 1

            # Otherwise, set the entry equal to 0
            else:
                arr_altered[i,j] = 0

    # Return the new adjacency matrix
    return  arr_altered



''' plot_graph Documentation
 -------------------------------
 This method takes in a plot and type of plot desired, then prints a plot of the graph.'''

def plot_graph(G, type, effects_mark, nodeNames):

    # This type places the failure effects on one side and the failure modes on the other in vertical columns
    if type == "bipartite":
        pos = nx.bipartite_layout(G, nodeNames[:effects_mark])  # positions for all nodes
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[effects_mark:], node_color="#fabc98", edgecolors="#c89679", node_shape="s")
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[:effects_mark], node_color="#98c5ed", edgecolors="#799dbd", node_shape="s")

        # This for loop places labels outside of the nodes for easier visualization
        for i in pos:
            x, y = pos[i]
            if x <= 0:
                plt.text(x-0.075,y,s=i, horizontalalignment='right')

            else:
                plt.text(x+0.075,y,s=i, horizontalalignment='left')
        
        #nx.draw_networkx_labels(G, pos, horizontalalignment='center')
        nx.draw_networkx_edges(G, pos)
        plt.box(False)
        plt.show()
        return

    # This type draws all the nodes in a circle
    elif type == "circular":
        pos = nx.circular_layout(G)  # positions for all nodes
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[:effects_mark], node_color="#98c5ed")
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[effects_mark:], node_color="#fabc98")
        nx.draw_networkx_labels(G, pos)
        nx.draw_networkx_edges(G, pos)
        plt.box(False)
        plt.show()
        return
    
    # This type attempts to have all the edges the same length
    elif type == "kamada_kawai":
        pos = nx.kamada_kawai_layout(G)  # positions for all nodes
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[:effects_mark], node_color="#98c5ed")
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[effects_mark:], node_color="#fabc98")
        nx.draw_networkx_labels(G, pos)
        nx.draw_networkx_edges(G, pos)
        plt.box(False)
        plt.show()
        return
    
    # This type attempts to have non-intersecting edges. This only works for planar graphs.
    elif type == "planar":
        pos = nx.planar_layout(G)  # positions for all nodes
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[:effects_mark], node_color="#98c5ed")
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[effects_mark:], node_color="#fabc98")
        nx.draw_networkx_labels(G, pos)
        nx.draw_networkx_edges(G, pos)
        plt.box(False)
        plt.show()
        return
    
    # This type plots the graph in a random configuration
    elif type == "random":
        pos = nx.random_layout(G)  # positions for all nodes
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[:effects_mark], node_color="#98c5ed")
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[effects_mark:], node_color="#fabc98")
        nx.draw_networkx_labels(G, pos)
        nx.draw_networkx_edges(G, pos)
        plt.box(False)
        plt.show()
        return
    
    # This type uses the eigenvectrs of the graph Laplacian in order to show possible clusters of nodes
    elif type == "spectral":
        pos = nx.spectral_layout(G)  # positions for all nodes
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[:effects_mark], node_color="#98c5ed")
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[effects_mark:], node_color="#fabc98")
        nx.draw_networkx_labels(G, pos)
        nx.draw_networkx_edges(G, pos)
        plt.box(False)
        plt.show()
        return
    
    # This type treats nodes as repelling objects and edges as springs (causing them to moving in simulation)
    elif type == "spring":
        pos = nx.spring_layout(G)  # positions for all nodes
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[:effects_mark], node_color="#98c5ed")
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[effects_mark:], node_color="#fabc98")
        nx.draw_networkx_labels(G, pos)
        nx.draw_networkx_edges(G, pos)
        plt.box(False)
        plt.show()
        return
    
    # This type places nodes in concentric circles
    elif type == "shell":
        pos = nx.shell_layout(G)  # positions for all nodes
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[:effects_mark], node_color="#98c5ed")
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[effects_mark:], node_color="#fabc98")
        nx.draw_networkx_labels(G, pos)
        nx.draw_networkx_edges(G, pos)
        plt.box(False)
        plt.show()
        return
    
    # This type places nodes in spiral layout
    elif type == "spiral":
        pos = nx.spiral_layout(G)  # positions for all nodes
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[:effects_mark], node_color="#98c5ed")
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[effects_mark:], node_color="#fabc98")
        nx.draw_networkx_labels(G, pos)
        nx.draw_networkx_edges(G, pos)
        plt.box(False)
        plt.show()
        return
    
    # This type places nodes in spiral layout
    elif type == "multipartite":
        
        pos = nx.multipartite_layout(G)  # positions for all nodes
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[:effects_mark], node_color="#98c5ed")
        nx.draw_networkx_nodes(G, pos, nodelist=nodeNames[effects_mark:], node_color="#fabc98")
        nx.draw_networkx_labels(G, pos)
        nx.draw_networkx_edges(G, pos)
        plt.box(False)
        plt.show()
        return

    # Plot the grpah with networkx and matplotlib using regular algorithm
    nx.draw(G, with_labels=True, font_weight='bold')
    plt.show()
    return



''' max_degrees Documentation
 --------------------------------
 This method takes in an adjacency matrix, the names of the nodes, a threshold for binarization
 of the adjacency matrix, and a boolean for labeling the nodes with names (True) or numbers (False).
 The output consists of three tuples: 
        1. (the nodes with the highest out degree, value of the highest out degree)
        1. (the nodes with the highest in degree, value of the highest in degree)
        1. (the nodes with the highest overall degree, value of the highest overall degree)'''

def max_degrees(arr, nodeNames, threshold = 0, name = False):
    
    # Create copy of the adjacency matrix so that we can alter it without losing information about
    # our original adjacency matrix.
    arr_altered = arr

    # Binarization of adjacency matrix. The threshold determines the cutoff for what will be labeled
    # as a 1 versus a 0.  Anything above the threshold will be a 1, and anything below the threshold
    # will be set to 0.
    for i in range(0, arr_altered.shape[0]):
        for j in range(arr_altered.shape[0]):
            if arr[i,j] > threshold:
                arr_altered[i,j] = 1
            else:
                arr_altered[i,j] = 0

    # Calculating out degrees and the maximum of the out degrees
    out_degrees = np.sum(arr_altered, axis=1)
    max_out = np.where(out_degrees == max(out_degrees))

    # Calculating in degrees and the maximum of the in degrees
    in_degrees = np.sum(arr_altered, axis=0)
    max_in = np.where(in_degrees == max(in_degrees))

    # Calculating overall degrees and the maximum of the overall degrees
    degrees = out_degrees + in_degrees
    max_deg = np.where(degrees == max(degrees))

    # If the user chooses to list the nodes by their proper name (rather than with numbers), then
    # we index into the array of node names, nodeNames, to find the names of these nodes. We then
    # return the three tuples about maximum out degree, in degree, and overall degree.
    if name: 
        out_name = nodeNames[max_out[0].tolist()]
        in_name = nodeNames[max_in[0].tolist()]
        deg_name = nodeNames[max_deg[0].tolist()]
        return  (out_name, max(out_degrees)), (in_name, max(in_degrees)), (deg_name, max(degrees))
    
    # Otherwise, if the user chooses not to label nodes with their proper names, then we return the
    # three tuples about maximum out degree, in degree, and overall degree. Rather than listing the
    # names of the nodes, their corresponding numbers are listed.
    return  (max_out, max(out_degrees)), (max_in, max(in_degrees)), (max_deg, max(degrees))



''' getProbabilities Documentation
 -------------------------------
 This method inputs an Excel file and reads the probabilities from the file. We return an array of these
 probabilities and an array of the node names (indexed the same as the probability array).'''

def getProbabilities(fileName, sheetName = None):
    df = pd.read_excel(fileName, sheet_name=sheetName)

    # Convert the data frame to a numpy array
    arr = df.to_numpy()

    # Create an array of the names of the nodes for the graph.
    nodeNames = arr[:, 0].flatten()
    
    # return the adjacency matrix without labels and an array of the names of the nodes
    return arr[:, 1], nodeNames


'''diagonal_nodes Documentation
 --------------------------------
 This method takes in an adjacency matrix and outputs a diagonal matrix. For an adjacency matrix of length k, the
 outputted diagonal matrix places values 1-k on the diagonal.'''

def diagonal_nodes(arr):
    # Return the diagonal matrix with i+1 in the i,i spot and zero elsewhere
    return np.diag([i+1 for i in range(arr.shape[0])])

'''cosine_similarity Documentation
 -----------------------------------
 This method takes in two vectors and outputs the cosine similarity between them. In linear algebra, cosine
 similarity is defined as the measure of how much two vectors are pointing in the same direction. It can also
 be thought of cosine of the angle between the two vectors.'''

def cosine_similarity(v1, v2):
    # Find the dot product between the two vectors inputted
    dot_prod = np.transpose(v1) @ v2

    # Calculate and return the cosine similarity of the two vectors
    return (dot_prod) / ((np.sqrt(np.sum(np.square(v1))) * np.sqrt(np.sum(np.square(v2)))))



''' make_undirected_unweighted Documentation
 ---------------------------------------------
 This method takes in an adjacency matrix and threshold. It outputs an adjusted adjacency matrix that
 corresponds to the undirected and unweighted graph of the one inputted.'''

def make_undirected_unweighted(arr, threshold = 0):
    # Make sure that the array is a numpy array
    arr = np.array(arr)

    # Add the array to the transpose of itself. This makes the graph undirected.
    arr = arr + arr.T

    # Use the method from graphBuilder.py to binarize the matrix. This makes the graph unweighted.
    make_binary(arr, threshold)

    # Since the array now only consists of zeros and ones, make sure that it is of integer type
    arr = arr.astype(int)

    # Return the adjusted matrix
    return arr