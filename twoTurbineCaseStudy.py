import numpy as np
from graphBuilder import *
from searchMethods import *

''' twoTurbineCaseStudy.py -----------------------------------------------------------------------------------------------

            ****************************** Last Updated: 18 April 2024 ******************************

 Methods:
 1) twoTurbine_bayesian_table: input adjusted adjacency matrix, adjacency matrix, current node, list of node names,
 adjusted list of node names --> outputs list of parents, probability table

-----------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------------------


   twoTurbine_bayesian_table Documentation
 ------------------------------------------
 This method inputs an adjusted adjacency matrix, regular adjacency matrix, current node (integer), array of node names (strings), 
 and adjusted node names (same indexing as adjusted adjacency matrix). We then calculate the probability table for the current node
 given its parents in the adjusted adjacency matrix. We use the weights from the regular adjacency matrix to determine the transitional
 probabilities between nodes. We output an array of the current node's parents and its probability table.'''


def twoTurbine_bayesian_table(a, arr, current, nodeNames, non):
    # arr, nodeNames = excel2Matrix("ExcelFiles/failureData.xlsx", "twoTurbines_simplified") # For debugging, feel free to uncomment
    adj = a.copy()
    nodes = diagonal_nodes(adj) # Diagonal matrix of node's numerical names +1
    adj = make_binary(adj, 0) # Binarize adjacency matrix
    parents_bool = nodes @ adj[:, current-1] # vector of zeros and parent names (numerical names)
    parents = list(parents_bool[np.nonzero(parents_bool)]) # list of just the parent names (numerical names)
    prob_table = np.zeros((2 ** len(parents), len(parents) + 2)) # Initialize the array for the node's conditional probability distribution

    current_index = np.where(nodeNames == non[int(current - 1)])[0][0]

    for iteration in range(2 ** len(parents)): # For each combination of parents having either failed or not
        true_false_array = [] # Iniitalize to record which parents have failed

        for p in range(len(parents)):
            parent_index = np.where(nodeNames == non[int(parents[p] - 1)])[0][0]
            prob_table[iteration][p] = int(iteration / (2 ** p)) % 2 # Append whether or not the parent node has failed
            prob_table[iteration][-2] += (1- int(iteration / (2 ** p)) % 2) * arr[parent_index][current_index] # Determine parent probability (given if the parent has failed or not)
            true_false_array.append((int(iteration / (2 ** p)) % 2)) # Append whether or not the parent node has failed
        
        if prob_table[iteration][-2] > 1: # If the value is greater than 1, set the probability to 1
            prob_table[iteration][-2]  = 1
        prob_table[iteration][-1] = 1 - prob_table[iteration][-2] # Add the probability of the node not failing to the array
    # print(prob_table)
    return parents, prob_table



'''
# The following code is for writing excel file with all pair-wise probabilities. To run, copy and paste the code below into main.py and hit 'Run.'

adjacency_matrix, nodeNames = excel2Matrix("ExcelFiles/failureData.xlsx", "twoTurbines_simplified")
current = 1
midpoint = True
num_iterations = 1
probabilities, nodeNames = getProbabilities("ExcelFiles/failureProbabilities.xlsx", sheetName = "Conditional Probabilities (2)")
all_probs = np.zeros(adjacency_matrix.shape)
with pd.ExcelWriter("bayesian_simulations3.xlsx") as writer:
    for i in range(1, len(nodeNames) + 1):
        print("----------------------------------", i, "----------------------------------")
        array_of_probs = np.zeros((len(nodeNames), 2))
        for j in range(1,len(nodeNames)+1):
            adjacency_matrix2 = adjacency_matrix.copy()
            probabilities = np.array(probabilities).copy()
            A, B = backward_bayesian_inference(adjacency_matrix2, nodeNames, [i], [i], [j], probabilities, start_bool = True, multi = False, poc="child", twoTurbine = True)
            array_of_probs[j-1] = A
        df = pd.DataFrame(np.array(array_of_probs))
        df.to_excel(writer, sheet_name="Probabilities given "+str(i))
        all_probs[:,i-1] = array_of_probs[:,0]
    df2 = pd.DataFrame(np.array(all_probs))
    df2.to_excel(writer, sheet_name="allProbs")
'''