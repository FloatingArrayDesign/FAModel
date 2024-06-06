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
        angle_degree = float(input("What angle do you want to use? (in degrees) "))
        self.angle_radians = angle_degree/360 * math.pi * 2

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
                self.G.add_node(platform_failure + "\n" + str(platform), probability=fail_prob, obj = [platform_obj])
                self.G = self.addMoreEdges(platform_failure, platform, [platform])

            # Create failure nodes
            for turbine_failure in systems['turbine']:
                if turbine_failure in failure_probabilities.keys(): fail_prob = failure_probabilities[turbine_failure]
                else: fail_prob = init_prob_dict[turbine_failure]
                self.G.add_node(turbine_failure + "\n" + str(platform),  probability=fail_prob, obj = [platform_obj])
                self.G = self.addMoreEdges(turbine_failure, platform, [platform])

            # FIRST DEGREE EDGES -------------------------------------------------------------------------------------------
            for attach1 in attachments.keys():
                attach1_name = str(attachments[attach1]['id'])
                attach1_type = ''
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
                    if attach1_failure in attachments[attach1]['obj'].failure_probability.keys(): fail_prob = failure_probabilities[attach1_failure]
                    else: fail_prob = init_prob_dict[attach1_failure]
                    self.G.add_node(attach1_failure + "\n" + attach1_name,  probability=fail_prob, obj = [attachments[attach1]['obj']])
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
                            if clash_failure in attachments[attach1]['obj'].failure_probability.keys(): fail_prob = failure_probabilities[clash_failure]
                            else: fail_prob = init_prob_dict[clash_failure]
                            self.G.add_node(clash_failure + "\n" + clash_name,  probability=fail_prob, obj = [attachments[attach1]['obj'], attachments[attach2]['obj']])
                            self.G = self.addMoreEdges(clash_failure, clash_name, [platform, attach1_name, attach2_name, clash_name])
                            if attach1_type == 'mooring' and attach2_type == attach1_type: mooring_clashes.append(clash_failure + "\n" + clash_name)
                            elif ('shared' not in attach1_type) and ('shared' not in attach2_type): cable_clashes.append(clash_failure + "\n" + clash_name)

                # SUBNODES AND SUBEDGES ------------------------------------------------------------------------------------
                subcomponents = attachments[attach1]['obj'].subcomponents
                component_num = 0
                for component in subcomponents:
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
                        if component_failure in component.failure_probability.keys(): fail_prob = failure_probabilities[component_failure]
                        else: fail_prob = init_prob_dict[component_failure]
                        self.G.add_node(component_failure + "\n" + component_name,  probability=fail_prob, obj = [component])
                        self.G = self.addMoreEdges(component_failure, component_name, [platform, attach1])


                # SECOND ORDER NODES -------------------------------------------------------------------------------------------
                attached_to = attachments[attach1]['obj'].attached_to
                for attach1A in attached_to:
                    attach1A_name = str(attach1A.id)

                    # Create anchor failure nodes
                    if 'anchor' in str(type(attach1A)):
                        attach1A_type = 'anchor'
                        if len(attach1A.attachments) > 1: attach1A_type = 'sharedanchor'
                        for anchor_failure in systems[attach1A_type]:
                            if anchor_failure in attach1A.failure_probability.keys(): fail_prob = failure_probabilities[anchor_failure]
                            else: fail_prob = init_prob_dict[anchor_failure]
                            self.G.add_node(anchor_failure + "\n" + attach1A_name,  probability=fail_prob, obj = [attach1A])
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
                            if grid_failure in attach1A.failure_probability.keys(): fail_prob = failure_probabilities[grid_failure]
                            else: fail_prob = init_prob_dict[grid_failure]
                            self.G.add_node(grid_failure + "\n" + attach1A_name,  probability=fail_prob, obj = [attach1A])
                            self.G = self.addMoreEdges(grid_failure, attach1A_name, [platform, attach1_name, attach1A_name])

            # Create mooring-mooring clashing failure node if no two mooring lines likely to clash
            if len(mooring_clashes) < 1:
                if systems['mooringmooring'][0] in failure_probabilities.keys(): fail_prob = failure_probabilities[systems['mooringmooring'][0]]
                else: fail_prob = init_prob_dict[systems['mooringmooring'][0]]
                self.G.add_node(systems['mooringmooring'][0] + "\n" + str(platform),  probability=fail_prob, obj= [platform_obj])
                self.G = self.addMoreEdges(systems['mooringmooring'][0], str(platform), [platform])

            # Create cable-mooring clashing failure nodes if no cable and mooring pairing likely to clash
            if len(cable_clashes) < 1:
                for cable_num_obj in num_cables:
                    cable_num = cable_num_obj.id
                    for clashing_failure in systems['cablemooring']:
                        if clashing_failure in failure_probabilities.keys(): fail_prob = failure_probabilities[clashing_failure]
                        else: fail_prob = init_prob_dict[clashing_failure]
                        self.G.add_node(clashing_failure + "\n" + str(platform) + ' ' + str(cable_num),  probability=fail_prob, obj = [platform_obj, cable_num])
                        self.G = self.addMoreEdges(clashing_failure, str(platform) + ' ' + str(cable_num), [platform])
        return self.G



    def get_systems(self, nodeNames):
        '''Create dictionary for each subsystem of failures
        Parameters
        ----------
        nodeNames : list
            List of all the failure names to use to create dictionary of subsystems
        '''
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


    def get_critical_node(self, param):
        '''Identify and return the critical failure(s) of the failure graph
        Parameters
        ----------
        param : string
            Measurement for criticality (either initial probability, degree [in, out, or total], susceptibility, or impact)
        '''
        nodalNames = np.array(list(self.G.nodes))
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
        for start_component in range(1,arr.shape[0]+1): # Iterate through each failure mode/effect in turbine
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

        # Calculate and return highest impact and susceptibility
        mean_impact = np.mean(all_probabilities, axis=1)
        max_impact = np.max(mean_impact)
        max_impact_index = np.where(mean_impact == max_impact)[0][0]

        mean_susceptibility = np.mean(all_probabilities, axis=0)
        max_susceptibility = np.max(mean_susceptibility, axis=0)
        max_impact_susceptibility = np.where(mean_impact == max_impact)[0][0]
        return [max_impact, nodeNamesArray[max_impact_index]], [max_susceptibility, nodeNamesArray[max_impact_susceptibility]]



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