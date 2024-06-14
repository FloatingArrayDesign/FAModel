import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from failure.failureProbabilities import *
from failure.twoTurbineCaseStudy import *
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

        # Get initial failure probabilities for each failure mode and effect
        probabilities, array_of_probs = getProbabilities(probabilities_file, probability_sheet)
        init_prob_dict = {}
        for prob_index in range(len(probabilities)):
            init_prob_dict.update({nodeNames[prob_index]: probabilities[prob_index]})

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

            # Create failure nodes
            for turbine_failure in systems['turbine']:
                if turbine_failure in failure_probabilities.keys(): fail_prob = failure_probabilities[turbine_failure]
                else: fail_prob = init_prob_dict[turbine_failure]
                self.G.add_node(turbine_failure + "\n" + str(platform),  probability=fail_prob, obj=[platform_obj], failure=turbine_failure, m_or_e=self.mode_effect_dict[turbine_failure])
                self.G = self.addMoreEdges(turbine_failure, platform, [platform])

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

                # Create moroing/cable failure nodes
                for attach1_failure in systems[attach1_type]:
                    original_name = attach1_name
                    if 'connect' in attach1_failure: attach1_name = platform + attach1_name
                    if attach1_failure in failure_probabilities.keys(): fail_prob = failure_probabilities[attach1_failure]
                    else: fail_prob = init_prob_dict[attach1_failure]
                    self.G.add_node(attach1_failure + "\n" + attach1_name,  probability=fail_prob, obj=[attachments[attach1]['obj']], failure=attach1_failure, m_or_e=self.mode_effect_dict[attach1_failure])
                    self.G = self.addMoreEdges(attach1_failure, attach1_name, [platform, attach1_name])
                    attach1_name = original_name

                # Create clashing failure nodes
                for attach2 in attachments.keys():
                    attach2_name = str(attachments[attach2]['id'])
                    attach2_type = ''
                    clash_name = str(attach1_name)+str(attach2_name)
                    if 'mooring' in str(type(attachments[attach2]['obj'])): attach2_type = 'mooring'
                    elif 'cable' in str(type(attachments[attach2]['obj'])): attach2_type = 'cable'
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
                    failure_probabilities = component.failure_probability
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
            if len(mooring_clashes) < 1:
                if systems['mooringmooring'][0] in failure_probabilities.keys(): fail_prob = failure_probabilities[systems['mooringmooring'][0]]
                else: fail_prob = init_prob_dict[systems['mooringmooring'][0]]
                self.G.add_node(systems['mooringmooring'][0] + "\n" + str(platform),  probability=fail_prob, obj=[platform_obj], failure=systems['mooringmooring'][0], m_or_e=self.mode_effect_dict[systems['mooringmooring'][0]])
                self.G = self.addMoreEdges(systems['mooringmooring'][0], str(platform), [platform])

            # Create cable-mooring clashing failure nodes if no cable and mooring pairing likely to clash
            if len(cable_clashes) < 1:
                for cable_num_obj in num_cables:
                    cable_num = cable_num_obj.id
                    for clashing_failure in systems['cablemooring']:
                        if clashing_failure in failure_probabilities.keys(): fail_prob = failure_probabilities[clashing_failure]
                        else: fail_prob = init_prob_dict[clashing_failure]
                        self.G.add_node(clashing_failure + "\n" + str(platform) + ' ' + str(cable_num),  probability=fail_prob, obj=[platform_obj, cable_num], failure=clashing_failure, m_or_e=self.mode_effect_dict[clashing_failure])
                        self.G = self.addMoreEdges(clashing_failure, str(platform) + ' ' + str(cable_num), [platform])



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


    def get_critical_node(self, param):
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
                    critical_prob[0] = nodal_probabilities[failure_node]
                    critical_prob[1] = [failure_node]
                elif nodal_probabilities[failure_node] == critical_prob[0]:
                    critical_prob[1].append(failure_node)
        
        # Find the critical failure for crtitical refering to the maximum degree (either in-degree, out-degree, or total-degree)
        elif 'deg' in param:
            out_deg, in_deg, deg = max_degrees(nx.to_numpy_array(self.G), nodalNames, threshold = 0, name = True)
            if 'in' in param:
                critical_prob[0] = in_deg[1]
                critical_prob[1] = list(in_deg[0])
            elif 'out' in param:
                critical_prob[0] = out_deg[1]
                critical_prob[1] = list(out_deg[0])
            else:
                critical_prob[0] = deg[1]
                critical_prob[1] = list(deg[0])
        
        # Find the crtitical failure for ctitical refering to the susceptibility or impact of a node
        elif 'sus' in param or 'impact' in param:
            max_impact, max_sus = self.get_susceptibility_and_impact()
            if 'impact' in param: critical_prob = max_impact
            elif 'sus' in param: critical_prob = max_sus
            else: return
        
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
                if ('rift off' in edge[0].replace("\n", " ") or 'ncreased' in edge[0].replace("\n", " ")) or 'ynamics' in edge[0].replace("\n", " "): self.G = self.edgeReweight(edge)
                elif ('apsize' in edge[0].replace("\n", " ") or '-cable' in edge[0].replace("\n", " ")) or 'ing line non' in edge[0].replace("\n", " "): self.G = self.edgeReweight(edge)
                elif ('ragging' in edge[0].replace("\n", " ") or 'hain' in edge[0].replace("\n", " ")) or 'ire rope' in edge[0].replace("\n", " "): self.G = self.edgeReweight(edge)
                elif ('ynthetic' in edge[0].replace("\n", " ") or 'able profile' in edge[0].replace("\n", " ")) or 'ared line' in edge[0].replace("\n", " "): self.G = self.edgeReweight(edge)
                elif ('load on cable' in edge[0].replace("\n", " ") or 'eight' in edge[0].replace("\n", " ")): self.G = self.edgeReweight(edge)

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
        max_impact_index = np.where(mean_impact == max_impact)[0][0]

        mean_susceptibility = np.mean(all_probabilities, axis=0)
        max_susceptibility = np.max(mean_susceptibility, axis=0)
        max_impact_susceptibility = np.where(mean_impact == max_impact)[0][0]

        self.mean_is = np.hstack(np.reshape(mean_impact, (len(list(self.G.nodes)), 1)), np.reshape(mean_susceptibility, (len(list(self.G.nodes)), 1)))
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
    

    def run_moorpy_simulation(self):
        print('HDIGH')
        self.Array.getMoorPyArray
        fig, ax = self.ms.plot(ax=ax, color='red')
        print(self.ms.bodyList[0].r6)
        self.ms.bodyList[0].f6Ext = np.array([3e6, 0, 0, 0, 0, 0])
        print(self.ms.bodyList[0].r6)
    

    def run_raft_simulation(self):
        return


    def enact_failures(self, failure, upper = False):
        '''Update the FAModel based on failure(s) occurring
        Parameters
        ----------
        failure : object
            Object of failure you would like to enact
        '''
        # Detach from platform
        if self.G.nodes[failure]['failure'].lower() in ['chain', 'wire rope', 'synthetic rope', 'clump weights or floats', 'connectors', 'dynamic cable', 'terminations']:
            failure_obj = self.G.nodes[failure]['obj'][0].part_of
            failure_obj.detachFrom('b')

        # Detach cable at joint
        if self.G.nodes[failure]['failure'].lower() in ['offshore joints']:
            failure_obj = self.G.nodes[failure]['obj'][0]
            failure_obj.detach(failure_obj.attachments[failure_obj.attachments.keys[0]], 'a')
        # Update MoorPy Array
        # self.Array.getMoorPyArray
    
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
        children = [nodeNames[int(child_index)] for child_index in children_ints]

        # Return list of children of failure
        return children



    def check_for_effects(self, failure):
        '''Check if the children (failures) of the failure enacted were reached
        Parameters
        ----------
        failure : string
            Name of failure we would like to enact
        '''
        x=0
        # Initialize list of observed effects (those that have occurred)
        observed_effects = []

        # Enact the failure and get its children
        without_failure, with_failure = self.enact_failures(failure)
        children_of_failure = self.get_descendants(failure)

        # Check if effect has occurred. If so, remove the effect from the graph and add it to the list of observed effects
        for child in children_of_failure:
            if 'capsize' in child.lower() or 'excess dynamics' in child.lower():
                platform_obj = self.G.nodes[child]['obj'][0]
                case_num = 0
                platform_num = int(platform_obj.id[-1])
                results = self.Array.array.results['case_metrics'][case_num][platform_num]['<degree of freedom>_max']

            if ('stability' in child.lower() or 'sink' in child.lower()) or 'hydrostatic' in child.lower():
                platform_obj = self.G.nodes[child]['obj'][0]
                z_location = platform_obj.body.r6[2]

            if 'change in mooring profile' in child.lower():
                mooring_obj = self.G.nodes[child]['obj'][0]
                z_location = mooring_obj.ss.lineList[1]

            if 'drift' in child.lower() or 'clashing' in child.lower():
                platform_obj = self.G.nodes[child]['obj'][0]
                x,y = platform_obj.getWatchCircle()


            elif without_failure == with_failure:
                self.G.remove_node(child)
                observed_effects.append(child)
        return observed_effects
    


    def choose_new_failure(self, failure):
        '''Choose new failures to enact based off of enacted failure
        Parameters
        ----------
        failure : string
            Name of failure we would like to enact/have enacted
        '''
        observed_effects = self.check_for_effects(failure)
        new_failures = []

        # If there are no observed effects, ask if user would like to use the new critical failure. If so, return the new critical failure modes
        if len(observed_effects) < 1:
            pick_cf = input("\nThere are no observed effects from " + str(failure).replace('\n', ' ') + ".\nWould you like to update critical failure? (y/n) ")
            if 'y' in pick_cf: 
                critical_failures = self.update_critical_node(self.criticality_type)[1]

                # Going through all the critical failures, add the failure modes in the list to the new_failures list
                for critical_failure in critical_failures:
                    if self.G.nodes[critical_failure]["m_or_e"] == 'mode': new_failures.append(critical_failure)
        
        else:
            # Find the descendants of each observed effect
            for observed_effect in observed_effects:
                possible_failures = self.get_descendants(observed_effect)

                # Going through all the descendants from the observed effect, add the failure modes in the list to the new_failures list
                for possible_failure in possible_failures:
                    if self.G.nodes[possible_failure]["m_or_e"] == 'mode': new_failures.append(possible_failure)
        return new_failures



    def update_critical_node(self, criticality_stipulation):
        '''Determine and return new critical failures (presumably after the previous critical failure(s) occur)
        Parameters
        ----------
        criticality_stipulation : string
            Specifies how to determine the new critical failures
        '''
        # If multiple critical nodes allowed, give the option to use both of the above methods to find critical nodes
        if ('child' in criticality_stipulation.lower() and 'global' in criticality_stipulation.lower()) or 'both' in criticality_stipulation.lower():
            clist1 = self.get_child_critical_node()
            clist2 = self.get_critical_node(self.criticality_type)
            self.critical_failures[0] = (clist1[0] + clist2[0])/2
            self.critical_failures[1] = clist1[1] + clist2[1]

        # Find new critical failure (such that the new critical failure must be a direct effect of the previous one)
        elif 'child' in criticality_stipulation.lower(): self.critical_failures = self.get_child_critical_node()

        # Find critical node such that it doesn't have to be a direct effect of the previous critical failure
        elif 'global' in criticality_stipulation.lower(): 
            for critical_failure in self.critical_failures[1]:
                self.G.remove_node(critical_failure)
            self.critical_failures = self.get_critical_node(self.criticality_type)

        return self.critical_failures
    

    def get_child_critical_node(self):
        '''Get the critical nodes when user desires future critical nodes be children of original critical failure(s)
        Parameters
        ----------
        None
        '''
        # Create list of direct effects from critical failure(s)
        nodalNames = np.array(list(self.G.nodes))
        critical_prob = [0, []]
        critical_successors = []
        for critical_failure in self.critical_failures[1]:
            node_index = np.where(nodalNames == critical_failure)[0][0]
            nodes = diagonal_nodes(nx.to_numpy_array(self.G))
            a = make_binary(nx.to_numpy_array(self.G), 0)
            child_bool = nodes @ a[:, node_index]
            children = child_bool[np.nonzero(child_bool)]
            for child in children:
                if not(nodalNames[int(child)] in critical_successors): critical_successors.append(nodalNames[int(child)])
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