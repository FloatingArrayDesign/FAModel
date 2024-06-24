import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from famodel.failure.failureProbabilities import *
from famodel.failure.twoTurbineCaseStudy import *
from famodel.project import Project
        

class failureGraph():
    def __init__(self, project_file):
        # Create project class instance from yaml file
        self.Array = Project(file=project_file)
        self.G = nx.DiGraph()
        self.ms = self.Array.ms # MoorPy object
        self.ms_initialized = False

        self.critical_failures = None
        self.criticality_type = None
        self.imput_and_susceptibility_table = None
        self.mean_is = None
        self.iteration = 0


    def create_failureGraph(self, matrix_file, matrix_sheet, probabilities_file, probability_sheet):
        '''Create a graph of failures based on the FAModel Project object
        Parameters
        ----------
        matrix_file : string
            Failure matrix file that encodes interaction between failure modes and effects
        matrix_sheet : string
            Name of sheet in Excel file to pull the failure matrix from
        probabilities_file : string
            File name for list of failure rates (iniital probabilities) of all failures
        probability_sheet : string
            Name of sheet in Excel file to pull the failure probabilities from
        '''

        print("\nBegin generating failure graph...")
        # Create adjacency matrix from failure matrix
        df = pd.read_excel(matrix_file, sheet_name=matrix_sheet)
        arr = df.to_numpy()[:,1:]
        nodeNames = df.to_numpy()[:, 0].flatten()
        self.nodeNames = nodeNames

        # Get initial failure probabilities for each failure mode and effect
        probabilities, array_of_probs = getProbabilities(probabilities_file, probability_sheet)
        init_prob_dict = {}
        for prob_index in range(len(probabilities)):
            init_prob_dict.update({nodeNames[prob_index]: probabilities[prob_index]})
        self.init_prob_dict = init_prob_dict

        # Determine angle of clashing we are interested in
        angle_degree = input("What angle do you want to use? (in degrees) ")
        if angle_degree == '': self.angle_radians = 0.0
        else: self.angle_radians = float(angle_degree)/360 * math.pi * 2

        # Initialize and create the dictionaries of the children and parents
        self.failures_c = {}
        self.failures_p = {}
        for i in range(arr.shape[0]):
            node_children = []
            node_parents = []
            for j in range(arr.shape[1]):
                if arr[i,j] > 0:
                    node_children.append(nodeNames[j])
                if arr[j,i] > 0:
                    node_parents.append(nodeNames[j])
            self.failures_c.update({nodeNames[i]: node_children})
            self.failures_p.update({nodeNames[i]: node_parents})

        # Get the systems (groups of failures by component type) and list of nodes that could impact the FAModel
        systems = self.get_systems(nodeNames)
        self.systems = systems

        # Initialize graph, boolean for plotting, and list of probabilities
        self.G = nx.DiGraph()

        # FIRST DEGREE NODES -------------------------------------------------------------------------------------------
        for platform in self.Array.platformList:
            attachments = self.Array.platformList[platform].attachments
            failure_probabilities = self.Array.platformList[platform].failure_probability
            platform_obj = self.Array.platformList[platform]
            nearby_platforms = []
            mooring_clashes = []
            cable_clashes = []
            num_cables = []

            # Create platform failure nodes
            for platform_failure in systems['platform']:
                if platform_failure in failure_probabilities.keys(): fail_prob = failure_probabilities[platform_failure]
                else: fail_prob = init_prob_dict[platform_failure]
                self.G.add_node(platform_failure + "\n" + str(platform), probability=fail_prob, obj=[platform_obj], failure=platform_failure, m_or_e=self.mode_effect_dict[platform_failure])
                self.G = self.addMoreEdges(platform_failure, platform, [platform])

            # FIRST DEGREE EDGES -------------------------------------------------------------------------------------------
            for attach1 in attachments.keys():
                attach1_name = str(attachments[attach1]['id'])
                attach1_type = ''
                failure_probabilities = attachments[attach1]['obj'].failure_probability
                if 'mooring' in str(type(attachments[attach1]['obj'])): 
                    if attachments[attach1]['obj'].shared: attach1_type = 'sharedmooring'
                    else: attach1_type = 'mooring'
                elif 'cable' in str(type(attachments[attach1]['obj'])): 
                    attach1_type = 'cable'
                    num_cables.append(attachments[attach1]['obj'])
                elif 'turbine' in str(type(attachments[attach1]['obj'])):
                    attach1_type = 'turbine'
                    # Create turbine failure nodes
                    # for turbine_failure in systems['turbine']:
                    #     if turbine_failure in failure_probabilities.keys(): fail_prob = failure_probabilities[turbine_failure]
                    #     else: fail_prob = init_prob_dict[turbine_failure]
                    #     self.G.add_node(turbine_failure + "\n" + str(platform),  probability=fail_prob, obj=[platform_obj], failure=turbine_failure, m_or_e=self.mode_effect_dict[turbine_failure])
                    #     self.G = self.addMoreEdges(turbine_failure, platform, [platform])

                # Create moroing/cable failure nodes
                for attach1_failure in systems[attach1_type]:
                    original_name = attach1_name
                    if 'connect' in attach1_failure: attach1_name = platform + attach1_name
                    if attach1_failure in failure_probabilities.keys(): fail_prob = failure_probabilities[attach1_failure]
                    else: fail_prob = init_prob_dict[attach1_failure]
                    self.G.add_node(attach1_failure + "\n" + attach1_name,  probability=fail_prob, obj=[attachments[attach1]['obj']], failure=attach1_failure, m_or_e=self.mode_effect_dict[attach1_failure])
                    self.G = self.addMoreEdges(attach1_failure, attach1_name, [platform, attach1_name])
                    attach1_name = original_name
                
                if attach1_type == 'turbine': continue

                # Create clashing failure nodes
                for attach2 in attachments.keys():
                    attach2_name = str(attachments[attach2]['id'])
                    attach2_type = ''
                    clash_name = str(attach1_name)+str(attach2_name)
                    if 'mooring' in str(type(attachments[attach2]['obj'])): attach2_type = 'mooring'
                    elif 'cable' in str(type(attachments[attach2]['obj'])): attach2_type = 'cable'
                    elif 'turbine' in str(type(attachments[attach2]['obj'])): 
                        attach2_type = 'turbine'
                        continue
                    for clash_failure in systems[(str(attach1_type)+str(attach2_type))]:
                        if 'shared' in attach1_type and all(abs(np.array(attachments[attach1]['obj'].rB[:2]) - np.array(attachments[attach2]['obj'].rB[:2])) < 100): reverse = True
                        else: reverse = False
                        if self.couldClash(clash_failure, attachments[attach1]['obj'], attachments[attach2]['obj'], reverse):
                            if clash_failure in failure_probabilities.keys(): fail_prob = failure_probabilities[clash_failure]
                            else: fail_prob = init_prob_dict[clash_failure]
                            self.G.add_node(clash_failure + "\n" + clash_name,  probability=fail_prob, obj=[attachments[attach1]['obj'], attachments[attach2]['obj']], failure=clash_failure, m_or_e=self.mode_effect_dict[clash_failure])
                            self.G = self.addMoreEdges(clash_failure, clash_name, [platform, attach1_name, attach2_name, clash_name])
                            if attach1_type == 'mooring' and attach2_type == attach1_type: mooring_clashes.append(clash_failure + "\n" + clash_name)
                            elif ('shared' not in attach1_type) and ('shared' not in attach2_type): cable_clashes.append(clash_failure + "\n" + clash_name)

                # SUBNODES AND SUBEDGES ------------------------------------------------------------------------------------
                subcomponents = attachments[attach1]['obj'].subcomponents
                component_num = 0
                for component in subcomponents:
                    if 'section' in str(type(component)).lower(): failure_probabilities = {}
                    else: failure_probabilities = component.failure_probability
                    component_num += 1
                    if 'mooring' in attach1_type:
                        if 'type' in component.keys():
                            # Create clump weight failure nodes
                            if 'str' in str(type(component['type'])) and 'weight' in component['type']: 
                                component_type = 'weight'
                                component_name = str(attach1_name + ' ' + component['type'] + ' ' + str(component_num))

                            # Create mooring material failure nodes
                            if 'dict' in str(type(component['type'])):
                                if 'polyester' in component['type']['material']: component_type = 'polyester'
                                elif 'chain' in component['type']['material']: component_type = 'chain'
                                elif 'rope' in component['type']['material']: component_type = 'rope'
                                component_name = attach1_name + ' ' + str(component['type']['name'])

                        # Create connector failure nodes
                        else:
                            component_type = 'connector'
                            component_name = attach1_name + ' connector'
                    
                    elif 'cable' in attach1_type:
                        # Create dynamic cable section failure nodes
                        if 'dynamic' in str(type(component)):
                            component_type = 'dynamic'
                            component_name = str(component.id)
                        
                        # Create static cable section failure nodes
                        elif 'static' in str(type(component)):
                            component_type = 'static'
                            component_name = str(component.id)

                        # Create offshore Joints failure nodes
                        elif 'joint' in str(type(component)).lower():
                            component_type = 'joints'
                            component_name = attach1_name + ' ' + str(component.id)

                    for component_failure in systems[component_type]:
                        if component_failure in failure_probabilities: fail_prob = failure_probabilities[component_failure]
                        else: fail_prob = init_prob_dict[component_failure]
                        self.G.add_node(component_failure + "\n" + component_name,  probability=fail_prob, obj=[component], failure=component_failure, m_or_e=self.mode_effect_dict[component_failure])
                        self.G = self.addMoreEdges(component_failure, component_name, [platform, attach1])


                # SECOND ORDER NODES -------------------------------------------------------------------------------------------
                attached_to = attachments[attach1]['obj'].attached_to
                for attach1A in attached_to:
                    attach1A_name = str(attach1A.id)
                    failure_probabilities = attach1A.failure_probability

                    # Create anchor failure nodes
                    if 'anchor' in str(type(attach1A)):
                        attach1A_type = 'anchor'
                        if len(attach1A.attachments) > 1: attach1A_type = 'sharedanchor'
                        for anchor_failure in systems[attach1A_type]:
                            if anchor_failure in attach1A.failure_probability.keys(): fail_prob = failure_probabilities[anchor_failure]
                            else: fail_prob = init_prob_dict[anchor_failure]
                            self.G.add_node(anchor_failure + "\n" + attach1A_name,  probability=fail_prob, obj=[attach1A], failure=anchor_failure, m_or_e=self.mode_effect_dict[anchor_failure])
                            self.G = self.addMoreEdges(anchor_failure, attach1A_name, [platform, attach1_name, attach1A_name])
                    
                    # Create edges between platforms
                    elif 'platform' in str(type(attach1A)): 
                        attach1A_type = 'platform'
                        attach1A_name = attach1A.id
                        for platform_failure in systems['platform']:
                            self.G = self.addMoreEdges(platform_failure, platform, [platform, attach1A_name])
                    
                    # Create substation/grid failure nodes
                    elif 'substation' in str(type(attach1A)):
                        attach1A_type = 'substation'
                        for grid_failure in systems['grid']:
                            if grid_failure in failure_probabilities.keys(): fail_prob = failure_probabilities[grid_failure]
                            else: fail_prob = init_prob_dict[grid_failure]
                            self.G.add_node(grid_failure + "\n" + attach1A_name,  probability=fail_prob, obj=[attach1A], failure=grid_failure, m_or_e=self.mode_effect_dict[grid_failure])
                            self.G = self.addMoreEdges(grid_failure, attach1A_name, [platform, attach1_name, attach1A_name])

            # Create mooring-mooring clashing failure node if no two mooring lines likely to clash
            failure_probabilities = self.Array.platformList[platform].failure_probability
            if len(mooring_clashes) < 1 and not(attach1_type == 'turbine'):
                if systems['mooringmooring'][0] in failure_probabilities.keys(): fail_prob = failure_probabilities[systems['mooringmooring'][0]]
                else: fail_prob = init_prob_dict[systems['mooringmooring'][0]]
                self.G.add_node(systems['mooringmooring'][0] + "\n" + str(platform),  probability=fail_prob, obj=[platform_obj], failure=systems['mooringmooring'][0], m_or_e=self.mode_effect_dict[systems['mooringmooring'][0]])
                self.G = self.addMoreEdges(systems['mooringmooring'][0], str(platform), [platform])

            # Create cable-mooring clashing failure nodes if no cable and mooring pairing likely to clash
            if len(cable_clashes) < 1 and not(attach1_type == 'turbine'):
                for cable_num_obj in num_cables:
                    cable_num = cable_num_obj.id
                    for clashing_failure in systems['cablemooring']:
                        if clashing_failure in failure_probabilities.keys(): fail_prob = failure_probabilities[clashing_failure]
                        else: fail_prob = init_prob_dict[clashing_failure]
                        self.G.add_node(clashing_failure + "\n" + str(platform) + ' ' + str(cable_num),  probability=fail_prob, obj=[platform_obj, cable_num], failure=clashing_failure, m_or_e=self.mode_effect_dict[clashing_failure])
                        self.G = self.addMoreEdges(clashing_failure, str(platform) + ' ' + str(cable_num), [platform])

        # If the user wants to input probabilities for specific edges, reweight edges based on the user's inputs
        user_inputs = input("Would you like to input probabilities into adjacency matrix? ")
        if (user_inputs == 'y' or user_inputs == 'yes') or user_inputs == 'True':
            for i in range(len(self.G.edges)):
                edge = list(self.G.edges)[i]
                self.G = self.edgeReweight(edge)

        # Initialize MoorPy array
        self.Array.getMoorPyArray(cables=1, pristineLines=1)
        self.Array.ms.initialize()
        self.Array.ms.solveEquilibrium()
        print()



    def get_systems(self, nodeNames):
        '''Create dictionary for each subsystem of failures
        Parameters
        ----------
        nodeNames : list
            List of all the failure names to use to create dictionary of subsystems
        '''
        # Create dictionary of which failures are modes and which are effects
        self.mode_effect_dict = {}
        for i in range(len(nodeNames)):
            if i < 26: self.mode_effect_dict.update({nodeNames[i]: 'effect'})
            else: self.mode_effect_dict.update({nodeNames[i]: 'mode'})

        # Systems and indices of corresponding failures in nodeNames
        turbine = [0,1,2,3,4,5,26,27,28,29]
        platform = [6,7,8,9,10,11,30,31,32]
        mooringmooring = [12]
        mooringcable = [13,14]
        rope = [35]
        polyester = [34]
        chain = [33]
        mooring = [15,16,17]
        connector = [36]
        clump_weight = [37]
        anchor = [19,20,38]
        cable = [21,22,23,40, 41,42, 45]
        dynamic = [43]
        static = [44]
        grid = [24,25]
        joints = [46]
        buoyancy = [40]
        sharedmooring = [15,16,18]
        sharedanchor = [19,20,39]

        # Dictionary of systems and their failures
        systems = {'turbine':nodeNames[turbine], 'platform':nodeNames[platform], 'mooringmooring':nodeNames[mooringmooring], 'rope':nodeNames[rope], 
                   'chain':nodeNames[chain],'polyester':nodeNames[polyester],'mooringcable':nodeNames[mooringcable], 'cablemooring':nodeNames[mooringcable], 
                   'mooring':nodeNames[mooring], 'connector':nodeNames[connector], 'weight':nodeNames[clump_weight], 'anchor':nodeNames[anchor], 'cable':nodeNames[cable], 
                   'grid':nodeNames[grid], 'dynamic':nodeNames[dynamic], 'static':nodeNames[static], 'buoyancy':nodeNames[buoyancy], 'cablecable': [], 
                   'sharedmooring':nodeNames[sharedmooring], 'sharedmooringcable':nodeNames[mooringcable], 'joints': nodeNames[joints], 'cablesharedmooring':nodeNames[mooringcable], 
                   'sharedmooringmooring':nodeNames[mooringmooring], 'mooringsharedmooring':nodeNames[mooringmooring], 'sharedanchor': nodeNames[sharedanchor]}
        return systems
    


    # def get_impact_nodes(self, nodeNames):
    #     '''Create dictionary for each node that tells us if the node could impact the FAModel
    #     Parameters
    #     ----------
    #     nodeNames : list
    #         List of all the failure names to use to create dictionary of subsystems
    #     '''
    #     # List of nodes that could impact hte FAModel
    #     could_impact = [2, 6, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 20, 21, 23]
    #     impactful_nodes = nodeNames[could_impact]

    #     # Create dictionary of nodes that tells us if the node impacts the FAModel or not
    #     impact_dict = {}
    #     for node in nodeNames:
    #         if node in impactful_nodes: impact_dict.update({node: True})
    #         else: impact_dict.update({node: False})
    #     return impact_dict


    def get_critical_node(self, param, mode = False):
        '''Identify and return the critical failure(s) of the failure graph
        Parameters
        ----------
        param : string
            Measurement for criticality (either initial probability, degree [in, out, or total], susceptibility, or impact)
        '''
        nodalNames = np.array(list(self.G.nodes))
        self.criticality_type = param
        critical_prob = [0, []]

        # Find critical failure for critical indicating the maximum initial probability
        if 'prob' in param:
            nodal_probabilities = nx.get_node_attributes(self.G, "probability")
            for failure_node in nodal_probabilities:
                if nodal_probabilities[failure_node] > critical_prob[0]:
                    if not mode or (mode and self.G.nodes[failure_node]['m_or_e'] == 'mode'):
                        critical_prob[0] = nodal_probabilities[failure_node]
                        critical_prob[1] = [failure_node]
                elif nodal_probabilities[failure_node] == critical_prob[0]:
                    critical_prob[1].append(failure_node)
        
        # Find the critical failure for crtitical refering to the maximum degree (either in-degree, out-degree, or total-degree)
        elif 'deg' in param:
            # Binaraize adjacency matrix
            arr_altered = nx.to_numpy_array(self.G)
            for i in range(0, arr_altered.shape[0]):
                for j in range(arr_altered.shape[0]):
                    if nx.to_numpy_array(self.G)[i,j] > 0: arr_altered[i,j] = 1
                    else: arr_altered[i,j] = 0

            # Get list of degrees
            if 'in' in param: degrees = np.sum(arr_altered, axis=0)
            elif 'out' in param: degrees = np.sum(arr_altered, axis=1)
            else: degrees = np.sum(arr_altered, axis=0) + np.sum(arr_altered, axis=1)
            
            # Determine max degree
            for node in range(len(degrees)):
                degree = degrees[node]
                if not mode or (mode and self.G.nodes[list(self.G.nodes)[node]]['m_or_e'] == 'mode'):
                    if degrees[node] > critical_prob[0]:
                        critical_prob[0] = degrees[node]
                        critical_prob[1] = [list(self.G.nodes)[node]]
                    elif degrees[node] == critical_prob[0]: critical_prob[1].append(list(self.G.nodes)[node])

        # Find the crtitical failure for ctitical refering to the susceptibility or impact of a node
        elif 'sus' in param or 'impact' in param:
            self.get_susceptibility_and_impact()
            if 'impact' in param: average_lst = np.mean(self.imput_and_susceptibility_table, axis=1)
            elif 'sus' in param: average_lst= np.mean(self.imput_and_susceptibility_table, axis=0)
            for node in range(len(average_lst)):
                if not mode or (mode and self.G.nodes[list(self.G.nodes)[node]]['m_or_e'] == 'mode'):
                    if average_lst[node] > critical_prob[0]:
                        critical_prob[0] = average_lst[node]
                        critical_prob[1] = [list(self.G.nodes)[node]]
                    elif average_lst[node] == critical_prob[0]: critical_prob[1].append(list(self.G.nodes)[node])
        
        self.critical_failures = critical_prob
        return critical_prob
        

    def get_susceptibility_and_impact(self):
        '''Run Bayesian inference over the graph to find all conditional probabilities (probability of A given B for all failures A and B), then average to determine
        the failure with the highest susceptibility and that with the highest impact

        Parameters
        ----------
        None
        '''
        print("\nStarting susceptibility and impact calculation... ")
        # Create list of node names and adjacency matrix from the failure graph
        nodeNamesArray = np.array(list(self.G.nodes))
        arr = nx.to_numpy_array(self.G)

        # If the user wants to input probabilities for specific edges, reweight edges based on the user's inputs
        user_inputs = input("Would you like to input probabilities into adjacency matrix? ")
        twoTurbine_calculationType = False
        if (user_inputs == 'y' or user_inputs == 'yes') or user_inputs == 'True':
            twoTurbine_calculationType = True
            for i in range(len(self.G.edges)):
                edge = list(self.G.edges)[i]
                self.G = self.edgeReweight(edge)
                # -------------------------> If using the two-Turbine test case, uncomment below for quicker reweighting: <-------------------------
                # if ('rift off' in edge[0].replace("\n", " ") or 'ncreased' in edge[0].replace("\n", " ")) or 'ynamics' in edge[0].replace("\n", " "): self.G = self.edgeReweight(edge)
                # elif ('apsize' in edge[0].replace("\n", " ") or '-cable' in edge[0].replace("\n", " ")) or 'ing line non' in edge[0].replace("\n", " "): self.G = self.edgeReweight(edge)
                # elif ('ragging' in edge[0].replace("\n", " ") or 'hain' in edge[0].replace("\n", " ")) or 'ire rope' in edge[0].replace("\n", " "): self.G = self.edgeReweight(edge)
                # elif ('ynthetic' in edge[0].replace("\n", " ") or 'able profile' in edge[0].replace("\n", " ")) or 'ared line' in edge[0].replace("\n", " "): self.G = self.edgeReweight(edge)
                # elif ('load on cable' in edge[0].replace("\n", " ") or 'eight' in edge[0].replace("\n", " ")): self.G = self.edgeReweight(edge)

        # Ask user if they are ready to continue to Bayesian network calculations (if not, quit)
        continue_input = input("Ready to continue? ")
        if 'n' in continue_input.lower(): quit()

        # Bayesian network calculation
        arr = nx.to_numpy_array(self.G)
        nodeNames = np.reshape(np.array(list(self.G.nodes)), (len(list(self.G.nodes)), ))
        all_probabilities = np.zeros(arr.shape) # Initialize a large array to put all the pairwise probabilities in
        for start_component in range(1,len(list(self.G.nodes))+1): # Iterate through each failure mode/effect in turbine
            a = arr.copy()
            non = nodeNames
            K, a, g, e, m, non = breadth_first_multi(a, nodeNames, [start_component], "child") # Generate tree for Bayesian network
            prblts = [] # Initialize array of node probabilities (in order of appearance in graph)
            for node in non:
                prblts.append(self.G.nodes[node]['probability']) # Add nodes to array of node probabilities
            prblts = np.array(prblts)

            probabilitiy_table = np.zeros((2, a.shape[0])) # Initialize table of inference probabilities
            nodes = diagonal_nodes(a) # Diagonal matrix of node names (numerical +1)
            a = make_binary(a, 0.5) # Binarize adjacency table

            nodeNamesArray = np.array(list(K.nodes)) # Initialize array of names of nodes
            nodeNamesArray = np.reshape(nodeNamesArray, (len(nodeNamesArray),)) # Make numpy array

            # Interence-----------------------------------------------------------------------
            for node in range(a.shape[0]):
                pts_bool = nodes @ a[:, node] # vector of zeros and child names (numerical names)
                pts = pts_bool[np.nonzero(pts_bool)] #list of just the child names (numerical names)

                if len(pts) < 1: # If no parents, add probability of failure happening to the probability table
                    probabilitiy_table[0][node] = self.G.nodes[list(self.G.nodes)[node]]['probability']
                    probabilitiy_table[1][node] = 1 - self.G.nodes[list(self.G.nodes)[node]]['probability']
                    continue
                
                if twoTurbine_calculationType:
                    parents, our_table = twoTurbine_bayesian_table(a, arr, node + 1, nodeNames, non) # Calculate the probability distribution table
                else:
                    parents, our_table = bayesian_table(a, node+1, True, nodeNames, True, prblts)
                mlt_table = np.ones((our_table.shape[0],2)) # Initialize table for multiplying across rows of probability distribution table

                # Calculate Probability Table ------------------------------------------------------------
                for i in range(our_table.shape[0]):
                    for j in range(our_table.shape[1] - 2):
                        parent = int(parents[j])
                        if our_table[i,j] == 0:
                            our_table[i,j] = probabilitiy_table[0][parent - 1]
                            if probabilitiy_table[0][parent - 1] == 0:
                                break
                        else:
                            our_table[i,j] = probabilitiy_table[1][parent - 1]
                            if (parent-1 == 0): # If the node's parent is the evidence, zero out the non-failing possibility
                                our_table[i,j] = 0
                        mlt_table[i,0] *= our_table[i,j] # Multiply the probabilities across the probability distribution table
                    mlt_table[i,1] = mlt_table[i,0] * our_table[i, -1] # Multiple by the probability of event not failing given combination of parent failure
                    mlt_table[i,0] *= our_table[i, -2] # Multiple by the probability of event failing given combination of parent failure
                sm_table = np.sum(mlt_table, axis = 0) #/np.sum(mlt_table) # Sum the products of probabilities across the columns
                probabilitiy_table[0][node] = sm_table[0] # Update the inference probability table with the probabilites just calculated
                probabilitiy_table[1][node] = sm_table[1]

                # Print and add probability of node to table
                print(start_component, node, " --> Probability of ", nodeNamesArray[node].replace("\n", " "), "=", sm_table)
                index2 = np.where(nodeNames == nodeNamesArray[node])[0][0]
                all_probabilities[0 * arr.shape[0] + start_component - 1][0 * arr.shape[0] + index2] = sm_table[0]/np.sum(sm_table)

        # Save the probabilities calcuated
        self.imput_and_susceptibility_table = all_probabilities

        # Calculate and return highest impact and susceptibility
        mean_impact = np.mean(all_probabilities, axis=1)
        max_impact = np.max(mean_impact)
        max_impact_index = np.where(mean_impact == max_impact)[0]

        mean_susceptibility = np.mean(all_probabilities, axis=0)
        max_susceptibility = np.max(mean_susceptibility, axis=0)
        max_impact_susceptibility = np.where(mean_impact == max_impact)[0]

        self.mean_is = np.hstack((np.reshape(mean_impact, (len(list(self.G.nodes)), 1)), np.reshape(mean_susceptibility, (len(list(self.G.nodes)), 1))))
        return [max_impact, nodeNames[max_impact_index]], [max_susceptibility, nodeNames[max_impact_susceptibility]]



    def couldClash(self, failure, a1, a2, reverse):
        '''Determine if two lines (either mooring or cable) could clash with each other
        Parameters
        ----------
        failure : string
            Name of failure (so that we can check if failure node already exists)
        a1      : object
            Either cable or mooring object for the first line
        a2      : object
            Either cable or mooring object for the second line
        reverse : boolean
            Determines orientation of the first line (True if orientation of vector needs to be reversed, False otherwise)
        '''
        # If the failure node already exists (perhaps by a different name), return that the lines cannot clash
        if a1 == a2 or (failure + "\n" + str(a1.id) + str(a2.id) in list(self.G.nodes) or failure + "\n" + str(a2.id) + str(a1.id) in list(self.G.nodes)): return False

        # Obtain the (x,y) coordinates of the start and end points of the lines
        a1_pnt1 = np.array(a1.rA[:2])
        a1_pnt2 = np.array(a1.rB[:2])
        a2_pnt1 = np.array(a2.rA[:2])
        a2_pnt2 = np.array(a2.rB[:2])

        # Determine the boundaries of movement of the lines (boundaries form a rectangle)
        a1MaxX, a1MinX, a1MaxY, a1MinY = self.get_min_max_vals(a1_pnt1, a1_pnt2, self.angle_radians, reverse)
        a2MaxX, a2MinX, a2MaxY, a2MinY = self.get_min_max_vals(a2_pnt1, a2_pnt2, self.angle_radians, False)

        # If the rectangles overlap, return TRUE (the lines can clash). Else, return FALSE (the lines cannot clash)
        overlap = False
        for corner_x in [a1MaxX,a1MinX,]:
                for corner_y in [a1MaxY, a1MinY]:
                    if (corner_x <= a2MaxX and corner_x >= a2MinX) and (corner_y <= a2MaxY and corner_y >= a2MinY): overlap = True
        return overlap



    def addMoreEdges(self, node, node_id, ids):
        '''Adds edges between current node and any already created nodes
        Parameters
        ----------
        node    : string
            Name of current failure
        node_id : string
            Name of current component that we are creating the failure for
        ids     : list of strings
            List of other component ids for components whose failure nodes have already been created
        '''
        # Create edges from current node to nodes already created (if created node is child of current node in reference matrix)
        for child in self.failures_c[node]:
            for id in ids:
                if child + "\n" + str(id) in self.G.nodes: self.G.add_edge(node + "\n" + str(node_id), child + "\n" + str(id), weight=0.001)

        # Create edges from nodes already created to current node (if created node is parent of current node in reference matrix)
        for parent in self.failures_p[node]:
            for id in ids:
                if parent + "\n" + str(id) in self.G.nodes: self.G.add_edge(parent + "\n" + str(id), node + "\n" + str(node_id), weight=0.001)
        return self.G



    def edgeReweight(self, edge):
        '''Ask user for new weight and set edge weight to user's input
        Parameters
        ----------
        edge : list of node names
            Edge that we want to reweight, comprised of the two nodes that it connects
        '''
        new_weight = input("Enter weight for (" + str(edge[0].replace("\n", " ")) + ", "  + str(edge[1].replace("\n", " ")) + ") (press \'Enter\' for default value)")
        if new_weight=='': new_weight=0.01
        self.G[edge[0]][edge[1]]['weight']=float(new_weight)
        return self.G



    def get_min_max_vals(self, pnt2, pnt1, angle_radians, reverse):
        '''Determine boundary of movement for a line (either cable or mooring)
        Parameters
        ----------
        pnt2          : list of floats
            Point (x,y) of beginning of line
        pn1           : list of floats
            Point (x,y) of end of line
        angle_radians : float
            Angle of movement (in radians)
        reverse : boolean
            Determines orientation of the first line (True if orientation of vector needs to be reversed, False otherwise)
        '''
        # If the vector's orientation needs to be reversed, reverse the order of the points
        if reverse:
            pnt_hold = pnt1
            pnt1 = pnt2
            pnt2 = pnt_hold

        # Find the vector that represents the line and the vector's length
        vector = pnt2 - pnt1
        length = math.sqrt(vector[0]**2 + vector[1]**2)

        # Find the angle (in polar coordinates) of the vector
        if vector[0] == 0.0 and vector[1] > 0: angle = math.pi/2
        elif vector[0] == 0.0 and vector[1] < 0: angle = 3*math.pi/2
        elif vector[0] == 0.0 and vector[1] == 0: angle = 0.0
        else: angle = math.atan(vector[1]/vector[0])
        if angle == 0 and vector[0] < 0: angle = math.pi
        if (angle > -math.pi*0.5 and angle < math.pi*0.5) and vector[0] < 0: angle += math.pi

        # Add and subtraact the angle of motion (while vector is in polar coordinates) to create two new vectors and then convert them back to rectangular coordinates
        new_vector1 = np.array([length * math.cos(angle - angle_radians), length * math.sin(angle - angle_radians)]) + np.array(pnt1)
        new_vector2 = np.array([length * math.cos(angle + angle_radians), length * math.sin(angle + angle_radians)]) + np.array(pnt1)
        if angle == math.pi and angle_radians==0: # Handle the issue of sin(0) not equal to 0 in python
            new_vector1 = np.array([length * math.cos(angle - angle_radians), 0]) + np.array(pnt1)
            new_vector2 = np.array([length * math.cos(angle + angle_radians), 0]) + np.array(pnt1)

        # Determine the bounds of the smallest rectangle which contains the original vector and the two new vectors
        max_x = max([new_vector1[0], new_vector2[0], pnt1[0], pnt2[0]])
        min_x = min([new_vector1[0], new_vector2[0], pnt1[0], pnt2[0]])
        max_y = max([new_vector1[1], new_vector2[1], pnt1[1], pnt2[1]])
        min_y = min([new_vector1[1], new_vector2[1], pnt1[1], pnt2[1]])
        if max_x == min_x:
            if max_x < 0: max_x = 0
            elif min_x > 0: min_x = 0
        return max_x, min_x, max_y, min_y
    


    def enact_failures(self, failure):
        '''Update the FAModel based on failure(s) occurring
        Parameters
        ----------
        failure : string
            Sting (name) of failure we would like to enact
        '''
        if failure not in list(self.G.nodes):
            return
        
        # --- Sections to be finished (but not by Emma) ---
        # Flood a ballast section in RAFT
        if 'watertightness' in self.G.nodes[failure]['failure'].lower():
            return
        # Remove a ballast section in RAFT
        if 'ballast system' in self.G.nodes[failure]['failure'].lower():
            return
        # Update RNA mass/inertia for blade or nacelle falling off/add load for hitting platform
        if 'RNA' in self.G.nodes[failure]['failure'].lower():
            tower_obj = self.G.nodes[failure]['obj'][0]
            turbine_num = int(failure[-1])
            raft_turbine = self.Array.array.fowtList[turbine_num]
            raft_turbine.mRNA += 0  # Update RNA mass --> need help on
            raft_turbine.IxRNA += 0 # Update RNA inerta --> need help on
            raft_turbine.IrRNA += 0 # Update RNA inerta --> need help on

            platform_obj = self.G.nodes[failure]['obj'][0].attachments[list(self.G.nodes[failure]['obj'][0].attachmentskeys())[0]]
            platform_obj.loads.update({'z': 0}) # Add load for hitting platform --> need help on

        # Update mass/inertia of tower and RNA to simulate a tower buckling and falling to the side
        if 'tower structur' in self.G.nodes[failure]['failure'].lower():
            tower_obj = self.G.nodes[failure]['obj'][0]
            turbine_num = int(failure[-1])
            raft_turbine = self.Array.array.fowtList[turbine_num]
            raft_turbine.calcStatics()
            raft_turbine.mtower += 0  # Update tower mass --> need help on
            raft_turbine.rCG_tow += 0 # Update tower center of mass --> need help on
            raft_turbine.mRNA += 0    # Update RNA mass --> need help on
            raft_turbine.IxRNA += 0   # Update RNA inerta --> need help on
            raft_turbine.IrRNA += 0   # Update RNA inerta --> need help on
        

        # --- Working sections ---
        # Detach line from platform
        if self.G.nodes[failure]['failure'].lower() in ['chain', 'wire rope', 'synthetic rope', 'clump weights or floats', 'connectors', 'dynamic cable', 'terminations']:
            failure_obj = self.G.nodes[failure]['obj'][0].part_of
            if 'termination' in failure.lower(): failure_obj = self.G.nodes[failure]['obj'][0]
            failure_obj.detachFrom('b')
            if 'termination' in failure.lower(): failure_obj = self.G.nodes[failure]['obj'][0].subcomponents[0]
            elif 'cable' in failure.lower(): failure_obj = self.G.nodes[failure]['obj'][0]
            # if not failure_obj.attached_to[-1]:
            #     failure_obj.attached_to = failure_obj.attached_to[:-1]
            # if 'moor' in str(type(failure_obj)): self.Array.ms.disconnectLineEnd(lineID=failure_obj.ss.number, endB=1)
            self.Array.ms.disconnectLineEnd(lineID=failure_obj.ss.number, endB=1)
            # print(failure_obj.attached_to)
            # self.Array.getMoorPyArray(cables = 1, pristineLines=1)

        # Detach cable at joint
        if self.G.nodes[failure]['failure'].lower() in ['offshore joints']:
            failure_obj = self.G.nodes[failure]['obj'][0]
            failure_obj.detach(failure_obj.attachments[failure_obj.attachments.keys[0]], 'a')
        
        # Anchor becomes a free point
        if 'anchor' in self.G.nodes[failure]['failure'].lower() and not('tether' in self.G.nodes[failure]['failure'].lower()):
            # anchor_pnt = self.G.nodes[failure]['obj'][0].r
            # for point in self.Array.ms.pointList:
            #         print(anchor_pnt[:2], point.r[:2])
            #         if all(abs(anchor_pnt[:2] - point.r[:2]) < 10):
            #                 point.type = 0
            self.G.nodes[failure]['obj'][0].mpAnchor.type = 0

        # Turbine not parked above cut-out wind speed (may want to add blades pitched at wrong angles)
        if self.G.nodes[failure]['failure'].lower() in ['turbine controls']:
            cut_out_wind_speed = 30 # in m/s --> need a way to determine this value
            new_case = self.Array.array.design['cases']['data'][-1]
            new_case[0] = cut_out_wind_speed + 5
            new_case[3] = 'operating'
            self.Array.array.design['cases']['data'].append(new_case)
            self.Array.array.analyzeCases()
        
        # Replace buoyancy section with regular cable section in MoorPy
        if 'buoyancy' in self.G.nodes[failure]['failure'].lower():
            all_buoys = True #bool('y' in input('Would you like to change all buoy sections? (y/n) '))
            names = []
            # Find buoyancy sections of cable
            for i in range(len(self.G.nodes[failure]['obj'][0].subcomponents[0].dd['sections'])):
                if 'buoy' in self.G.nodes[failure]['obj'][0].subcomponents[0].dd['sections'][i]['type']['name'].lower(): 
                    names.append(i)
            if not all_buoys: names = [names[0]]
            # Replace buoyancy section(s) with regular section
            for i in names:
                new_type = self.G.nodes[failure]['obj'][0].subcomponents[0].dd['sections'][i-1]['type']
                self.G.nodes[failure]['obj'][0].subcomponents[0].dd['sections'][i]['type'] = new_type
        # self.Array.getMoorPyArray(cables=1, pristineLines=1)
        self.Array.ms.initialize()
        self.Array.ms.solveEquilibrium()
        return

    

    def get_descendants(self, failure, threshold = 0.0):
        '''Find the children of a specific failure node
        Parameters
        ----------
        failure : string
            Name of the failure whose children you want to find
        threshold : float
            Value for which we will base our binarization of the adjacency matrix off of
            (anything > threshold will equal 1, anything < threshold will equal 0)
        '''
        # Get adjacency matrix and list of node names from graph
        adj = nx.to_numpy_array(self.G)
        arr = np.zeros(adj.shape)
        nodeNames = np.array(list(self.G.nodes))
        
        # Binarize adjacency matrix
        for i in range(len(nodeNames)):
            for j in range(len(nodeNames)):
                if adj[i][j] > threshold: arr[i][j] = 1
                else: arr[i][j] = 0

        # Find index of failure
        failure_index = np.where(nodeNames == failure)[0]

        # Create diagonal matrix of row numbers (index starts at 1)
        nodes = diagonal_nodes(arr)

        # Find children of failure
        child_bool = arr[failure_index] @ nodes # vector of zeros and child names (numerical names)
        children_ints = child_bool[np.nonzero(child_bool)] #list of just the child names (numerical names)
        children = [nodeNames[int(child_index - 1)] for child_index in children_ints]

        # Return list of children of failure
        return children

    def get_effects_identifiers(self, child):
        '''Check if the children (failures) of the failure enacted were reached
        Parameters
        ----------
        failure : string
            Name of failure we would like to enact
        '''        
        print("Starting effects check...'", child.replace('\n',' '))
        results = 'no results'

        # Check platform angles
        if 'capsize' in child.lower() or 'excess dynamics' in child.lower():
            platform_obj = self.G.nodes[child]['obj'][0]
            case_num = 0
            platform_num = int(platform_obj.id[-1])
            results = []
            self.Array.array.analyzeCases()
            results.append(self.Array.array.results['case_metrics'][case_num][platform_num - 1]['surge_max'])
            results.append(self.Array.array.results['case_metrics'][case_num][platform_num - 1]['sway_max'])
            results.append(self.Array.array.results['case_metrics'][case_num][platform_num - 1]['heave_max'])
            results.append(self.Array.array.results['case_metrics'][case_num][platform_num - 1]['roll_max'])
            results.append(self.Array.array.results['case_metrics'][case_num][platform_num - 1]['pitch_max'])
            results.append(self.Array.array.results['case_metrics'][case_num][platform_num - 1]['yaw_max'])
            # ---- Uncomment below if you want to include these measures ----
            # results.append(self.Array.array.results['case_metrics'][case_num][platform_num - 1]['AxRNA_max'])
            # results.append(self.Array.array.results['case_metrics'][case_num][platform_num - 1]['Mbase_max'])
            # results.append(self.Array.array.results['case_metrics'][case_num][platform_num - 1]['omega_max'])

        # Check z-location of the platform
        elif ('stability' in child.lower() or 'sink' in child.lower()) or 'hydrostatic' in child.lower():
            results = self.G.nodes[child]['obj'][0].body.r6[2]

        # Check for rope drag (taut/semitaut lines), line lay length (catenary/semitaut), shared line deep enough
        elif 'change in mooring profile' in child.lower():
            results = []
            # results = self.G.nodes[child]['obj'][0].ss.lineList[1]

        # Check watch circle
        elif 'drift' in child.lower() or 'clashing' in child.lower():
            output = self.G.nodes[child]['obj'][0].getWatchCircle()
            results = np.hstack((np.array(output[0]), np.array(output[1])))
            # results = 'WATCH CIRCLE STILL NOT WORKING'

        # Check loads on turbine
        elif 'turbine loads' in child.lower():
            results = self.G.nodes[child]['obj'][0].loads

        # Check if tensions exceed MBL or connector MBL?
        elif ('moor' in child.lower() and 'loads' in child.lower()) or 'nonfunctional' in child.lower():
            self.G.nodes[child]['obj'][0].updateTensions()
            results = [self.G.nodes[child]['obj'][0].loads['TAmax'], self.G.nodes[child]['obj'][0].loads['TBmax']]

            # ---------------- or ----------------
            # results = self.G.nodes[child]['obj'][0].ss.getTen(-1)
            # iLine = np.where(np.array(self.Array.ms.lineList) == cable.subcomponents[subcomponent].ss)[0]

        # Check if vertical tension for DEA, check if tensions exceed MBL
        elif 'anchor drag' in child.lower():
            anchor_obj = self.G.nodes[child]['obj'][0]
            anchor_obj.attachments[list(self.G.nodes[child]['obj'][0].attachments.keys())[0]]['obj'].updateTensions()
            results = self.G.nodes[child]['obj'][0].loads

            # ---------------- or ----------------
            # iLine = np.where(np.array(self.Array.ms.lineList) == cable.subcomponents[subcomponent].ss)[0]
            # results = self.G.nodes[child]['obj'][0].ss.getTen(-1)

        # Check if cable sag,hog,etc are within limits
        elif ('profile' in child.lower() or 'load' in child.lower()) and 'cable' in child.lower():
            results = []
            # Get platform(s) connected to cable
            for platform in self.G.nodes[child]['obj'][0].attached_to:
                    if 'platform' in str(type(platform)):
                        cables = []
                        # Determine which index the cable is in
                        for attachment in platform.attachments:
                                if 'cable' in str(type(platform.attachments[attachment]['obj'])).lower(): cables.append(platform.attachments[attachment]['obj'])
                        cable_index = np.where(np.array(cables) == self.G.nodes[child]['obj'][0])[0][0]
                        # Get the watch circle from the platform, which gives us the min curve and min sag of the array cable
                        a,b,c = platform.getWatchCircle(self.Array.ms)
                        results.append(c['minCurvSF'][cable_index])
                        results.append(c['minSag'][cable_index])
        # Plotting code - feel free to uncomment
        # self.Array.plot2d()
        # plt.show()
        return results
    


    def choose_new_failure_mode(self, failure):
        '''Choose new failures to enact based off of enacted failure
        Parameters
        ----------
        failure : string
            Name of failure we would like to enact/have enacted
        '''
        possible_failures = self.get_descendants(failure) # Finding children of failure effect
        new_failures = [] # List of new fialure modes

        # Going through all the children of the failure effect, add the failure modes in the list to the new_failures list
        for possible_failure in possible_failures:
            if self.G.nodes[possible_failure]["m_or_e"] == 'mode': new_failures.append(possible_failure)
        return new_failures



    def run_failure_simulation_iteration(self, param='prob', plot = False):
        '''Enact failure mode and determine if failure effects have occurred
        Parameters
        ----------
        param : string
            Measurement for criticality (either initial probability, degree [in, out, or total], susceptibility, or impact)
        '''
        # Determine failure modes to enact
        if not self.critical_failures: self.critical_failures = self.get_critical_node(param, mode=True)
        critical_node = self.critical_failures

        observed_effects = [] # List of observed failure effects
        new_critical_failures = []

        for node in critical_node[1]:
            print('\n----------- Enacting', node.replace('\n', ' '), '-----------')
            children = self.get_descendants(node) # Children of failure mode to check
            results_before = {}
            results_after = {}

            # Check starting state of failure effects
            for child in children: 
                results_before.update({child: self.get_effects_identifiers(child)})

            # Plot before failure enactment
            self.Array.updatePositions()
            if plot:
                self.Array.plot2d()
                plt.savefig("famodel/failure/Demos/" + param + "/" + str(self.iteration) + node.replace("\n", "-").replace("/","-") + "_demo_before.png")

            # Enact failure
            self.enact_failures(node)

            # Reevaluate state of FOWT
            self.Array.array.solveStatics(case=0) # Solve statics in RAFT
            self.Array.ms.solveEquilibrium() # Solve equilibruim in MoorPy
            for platform in self.Array.platformList:
                self.Array.platformList[platform].getWatchCircle() # Update/get watch circles for each platform

            # Check ending state of failure effects
            for child in children: 
                results_after.update({child: self.get_effects_identifiers(child)})
                
                # Identify if tensions are above MBL
                if (('moor' in child.lower() and 'loads' in child.lower()) or 'nonfunctional' in child.lower()) or 'anchor drag' in child.lower():
                    if ('moor' in child.lower() and 'loads' in child.lower()) or 'nonfunctional' in child.lower():
                        MBL = self.G.nodes[child]['obj'][0].dd['sections'][0]['type']['MBL']
                    elif 'anchor drag' in child.lower():
                        moor_line = self.G.nodes[child]['obj'][0].attachments[list(self.G.nodes[child]['obj'][0].attachments.keys())[0]]
                        MBL = moor_line.dd['sections'][0]['type']['MBL']
                    if 'e' in MBL:
                        MBL_info = MBL.split('e')
                        MBL = float(MBL_info[0]) ** int(MBL_info[1])
                    else: MBL = float(MBL)

                    if results_after[child][0] > MBL or results_after[child][1] > MBL:
                        observed_effects.append(child)
                        new_failure_modes = self.choose_new_failure_mode(child)
                        for nfm in new_failure_modes:
                            if not nfm in new_critical_failures and not nfm == node: new_critical_failures.append(nfm)
                        if len(new_failure_modes) > 1: print("Effect found -", child.replace('\n', ' '))
                    if child in list(self.G.nodes):
                        self.G.remove_node(child)

                # If the angle of roll is greater than 180 degrees (or pi in radians), then the turbine capsized
                elif 'capsize' in child.lower():
                    if results_after[child][3] > math.pi:
                        observed_effects.append(child)
                        new_failure_modes = self.choose_new_failure_mode(child)
                        for nfm in new_failure_modes:
                            if not nfm in new_critical_failures and not nfm == node: new_critical_failures.append(nfm)
                        if len(new_failure_modes) > 1: print("Effect found -", child.replace('\n', ' '))
                    if child in list(self.G.nodes):
                        self.G.remove_node(child)

                # Identify if state of failure effect has changed
                elif not(results_before[child] == results_after[child]):
                    observed_effects.append(child)
                    new_failure_modes = self.choose_new_failure_mode(child)
                    for nfm in new_failure_modes:
                        if not nfm in new_critical_failures and not nfm == node: new_critical_failures.append(nfm)

                    if child in list(self.G.nodes):
                        self.G.remove_node(child)
                    if len(new_failure_modes) > 1: print("Effect found -", child.replace('\n', ' '))
            
            # Check for new possibilities of clashing ndoes
            systems = self.systems
            connected_components = []
            for failure in list(self.G.nodes):
                # Find first line
                if ('cable' in str(self.G.nodes[failure]['obj'][0]) or 'moor' in str(self.G.nodes[failure]['obj'][0]))and not ('anchor' in str(type(self.G.nodes[failure]['obj'][0]))):
                    line1 = self.G.nodes[failure]['obj'][0]
                    # Find type of line
                    if 'cable'in str(self.G.nodes[failure]['obj'][0]): line1_type = 'cable'
                    elif line1.shared: line1_type = 'sharedmooring'
                    else: line1_type = 'mooring'
                    # Identify connected failures to attach edges to
                    if line1.attached_to:
                        for item in line1.attached_to:
                            if item: connected_components.append(item)
                        for subitem in line1.subcomponents:
                            if subitem: connected_components.append(subitem)
                        
                        # Find second line
                        for failure2 in list(self.G.nodes):
                            if (('cable' in str(self.G.nodes[failure2]['obj'][0]) or 'moor' in str(self.G.nodes[failure2]['obj'][0])) and not(line1 == self.G.nodes[failure2]['obj'][0])) and not ('anchor' in str(type(self.G.nodes[failure2]['obj'][0]))):
                                line2 = self.G.nodes[failure2]['obj'][0]
                                # Find type of line
                                if 'cable'in str(self.G.nodes[failure2]['obj'][0]): line2_type = 'cable'
                                elif 'moor' in str(self.G.nodes[failure2]['obj'][0]) and line2.shared: line2_type = 'sharedmooring'
                                else: line2_type = 'mooring'
                                # Identify connected failrues to attach edges to
                                if line2.attached_to:
                                    for item in line2.attached_to:
                                        if item: connected_components.append(item)
                                    for subitem in line2.subcomponents:
                                        if subitem: connected_components.append(subitem)
                                    connected_components.append(str(line1.id)+str(line2.id))

                                    # Check if the cables could clash. If so, create clashing node and constructed edges between appropriate failure nodes
                                    for clash_failure in systems[(str(line1_type)+str(line2_type))]:
                                        if 'shared' in line1_type and all(abs(np.array(line1.rB[:2]) - np.array(line2.rB[:2])) < 100): reverse = True
                                        else: reverse = False
                                        if self.couldClash(clash_failure, line1, line2, reverse) and not(clash_failure + "\n" + str(line1.id)+str(line2.id) in list(self.G.nodes)):
                                            if clash_failure in line1.failure_probability.keys(): fail_prob = line1.failure_probability[clash_failure]
                                            else: fail_prob = self.init_prob_dict[clash_failure]
                                            if not(('shared' in line1_type or 'shared' in line2_type) and ('anchor' in clash_failure.lower())):
                                                print('NEW CLASHING NODE -', clash_failure.replace('\n', ' ') + "\n" + str(line1.id)+str(line2.id))
                                                self.G.add_node(clash_failure + "\n" + str(line1.id)+str(line2.id),  probability=fail_prob, obj=[line1, line2], failure=clash_failure, m_or_e=self.mode_effect_dict[clash_failure])
                                                self.G = self.addMoreEdges(clash_failure, str(line1.id)+str(line2.id), connected_components)
            # Plot after failure enactment
            self.Array.updatePositions()
            if plot:
                self.Array.plot2d()
                plt.savefig("famodel/failure/Demos/" + param + "/" + str(self.iteration) + node.replace("\n", "-").replace("/","-") + "_demo_after.png")
            self.iteration += 1
        print()
        return new_critical_failures
    


    def run_failure_simulation(self, auto = True, plot = False):
        '''Run one iteration of enacting/checking failures, then decide which new failure mode to enact
        Parameters
        ----------
        None
        '''
        # Determine criticality measure or input critical node
        param = input('How would like to select the first failure mode to enact? (probability, degree, impact, susceptibility, or type-in) ')
        if 'type' in param.lower(): self.get_user_failure_node()
        else: self.criticality_type = param
        continue_on = True

        # While the user wants to continue and there are nodes left to enact, continue the following
        while len(self.G.nodes) > 0 and continue_on:

            # Run an interation of enacting a failure mode and checking for failure effects
            new_critical_failures = self.run_failure_simulation_iteration(param, plot)

            # If there are no new failures, determine which failures to continue with
            if len(new_critical_failures) < 1:
                print('There are no observed effects from the FAModel!')
                continue_bool = 'cr'
                if not auto: continue_bool = input('Would you like to choose a new critical node (type \'cr\'), enter a new critical node (type \'tcr\'), or quit (type \'q\')? ')
                
                # If the user wants us to determine the new critical node, find the new critical nodes based on previous criteria
                if continue_bool.lower() == 'cr':
                    stipulation = 'global' #input('Would you like child, global, or combined critical failures? (type \'child\', \'global\', or \'both\') ')
                    self.update_critical_node(stipulation, mode = True)

                # If the user wants to input their own failure, get their failure and set it as the new critical failure
                elif continue_bool == 'tcr': self.get_user_failure_node()
                
                # If the user wants to quit, quit the program
                elif continue_bool == 'q':
                    quit()

            # If there are new failure modes identified in the iteration, set them as the critical failures and remove old critical failures
            else:
                self.G.remove_nodes_from(self.critical_failures[1])
                self.critical_failures[1] = new_critical_failures
            # print("NEW CFs:", self.critical_failures)
            if len(self.critical_failures[1]) < 1: continue_on = False
            if not auto: continue_on = bool('y' in input('Are you ready to continue to next failure simulation? (y/n) '))



    def get_user_failure_node(self):
        '''Get the failure node from the user
        Parameters
        ----------
        None
        '''
        # Ask for failure node
        new_failure = input('Please type failure here: ')

        # If failure not in graph, keep asking
        while not(new_failure in list(self.G.nodes)) and not('q' in new_failure): 
            print('Failure not in graph.')

            # If the user wants to see the list of failure nodes, list them below
            show_nodes = bool('y' in input('Show list of failure nodes? (y/n) '))
            if show_nodes:
                for node in range(len(list(self.G.nodes))):
                    print([node, list(self.G.nodes)[node]])
            
            # Ask for number with associated failure (numbers shown in list produced above)
            new_failure_num = input('Please type number associated with failure (or \'q\' to quit): ')
            new_failure = list(self.G.nodes)[int(new_failure_num)]

        # If the user wants to quit, quit
        if 'q' in new_failure: quit()

        # Set new critical failures
        self.G.remove_nodes_from(self.critical_failures[1])
        self.critical_failures[1] = [new_failure]
        self.critical_failures[0] = 0

    

    def update_critical_node(self, criticality_stipulation, critical_failures = None, mode = False):
        '''Determine and return new critical failures (presumably after the previous critical failure(s) occur)
        Parameters
        ----------
        criticality_stipulation : string
            Specifies how to determine the new critical failures
        '''
        # If multiple critical nodes allowed, give the option to use both of the above methods to find critical nodes
        if ('child' in criticality_stipulation.lower() and 'global' in criticality_stipulation.lower()) or 'both' in criticality_stipulation.lower():
            clist1 = self.get_child_critical_node(critical_failures=critical_failures, mode = mode)
            clist2 = self.get_critical_node(self.criticality_type, mode)
            self.critical_failures[0] = (clist1[0] + clist2[0])/2
            self.critical_failures[1] = clist1[1] + clist2[1]

        # Find new critical failure (such that the new critical failure must be a direct effect of the previous one)
        elif 'child' in criticality_stipulation.lower(): self.get_child_critical_node(critical_failures=critical_failures, mode=mode)

        # Find critical node such that it doesn't have to be a direct effect of the previous critical failure
        elif 'global' in criticality_stipulation.lower(): 
            for critical_failure in self.critical_failures[1]:
                if critical_failure in list(self.G.nodes): self.G.remove_node(critical_failure)
            self.critical_failures = self.get_critical_node(self.criticality_type, mode)

        return self.critical_failures
    

    def get_child_critical_node(self, critical_failures = None, mode = False):
        '''Get the critical nodes when user desires future critical nodes be children of original critical failure(s)
        Parameters
        ----------
        None
        '''
        # Create list of direct effects from critical failure(s)
        nodalNames = np.array(list(self.G.nodes))
        critical_prob = [0, []]
        critical_successors = []
        if not critical_failures: critical_failures = self.critical_failures[1]
        for critical_failure in critical_failures:
            node_index = np.where(nodalNames == critical_failure)[0][0]
            nodes = diagonal_nodes(nx.to_numpy_array(self.G))
            a = make_binary(nx.to_numpy_array(self.G), 0)
            child_bool = nodes @ a[:, node_index]
            children = child_bool[np.nonzero(child_bool)]
            for child in children:
                if not(nodalNames[int(child)] in critical_successors) and not(nodalNames[int(child)] in self.critical_failures[1]):
                    if not mode or (mode and self.G.nodes[nodalNames[int(child)]]['m_or_e'] == 'mode'):
                        critical_successors.append(nodalNames[int(child)])
            if critical_failure in list(self.G.nodes):
                self.G.remove_node(critical_failure)
        
        # Find critical failure for critical indicating the maximum initial probability
        if 'prob' in self.criticality_type:
            nodal_probabilities = nx.get_node_attributes(self.G.subgraph(critical_successors), "probability")
            for failure_node in nodal_probabilities:
                if nodal_probabilities[failure_node] > critical_prob[0]: critical_prob = [nodal_probabilities[failure_node], [failure_node]]
                elif nodal_probabilities[failure_node] == critical_prob[0]: critical_prob[1].append(failure_node)
        
        # Find the critical failure for crtitical refering to the maximum degree (either in-degree, out-degree, or total-degree)
        elif 'deg' in self.criticality_type:
            out_deg, in_deg, deg = max_degrees(nx.to_numpy_array(self.G), nodalNames, threshold = 0, name = True)
            if 'in' in self.criticality_type: critical_prob = [in_deg[1], list(in_deg[0])]
            elif 'out' in self.criticality_type: critical_prob = [out_deg[1], list(out_deg[0])]
            else: critical_prob = [deg[1], list(deg[0])]
        
        # Find the crtitical failure for ctitical refering to the susceptibility or impact of a node
        elif 'sus' in self.criticality_type or 'impact' in self.criticality_type:
            if 'sus' in self.criticality_type: row_num = 1
            elif 'impact' in self.criticality_type: row_num = 0
            vals = np.zeros(len(critical_successors), 2)
            for node in critical_successors: vals[0] = self.mean_is[np.where(nodalNames == node)[0][0]]
            critical_prob = [np.max(vals), critical_successors[np.where(critical_successors == np.max(vals))]]

        self.critical_failures = critical_prob
        return critical_prob