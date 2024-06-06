import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from failure.failureProbabilities import *
from failure.twoTurbineCaseStudy import *
from famodel.project import Project

def couldClash(G, failure, a1, a2, reverse):
    if a1 == a2 or (failure + "\n" + str(a1.id) + str(a2.id) in list(G.nodes) or failure + "\n" + str(a2.id) + str(a1.id) in list(G.nodes)):
            return False
    a1_pnt1 = np.array(a1.rA[:2])
    a1_pnt2 = np.array(a1.rB[:2])
    a2_pnt1 = np.array(a2.rA[:2])
    a2_pnt2 = np.array(a2.rB[:2])
    # print('\n', a1_pnt1, a1_pnt2)
    # print(a2_pnt1, a2_pnt2)
    a1MaxX, a1MinX, a1MaxY, a1MinY = get_min_max_vals(a1_pnt1, a1_pnt2, angle_radians, reverse)
    a2MaxX, a2MinX, a2MaxY, a2MinY = get_min_max_vals(a2_pnt1, a2_pnt2, angle_radians, False)

    overlap = False
    for corner_x in [a1MaxX,a1MinX,]:
            for corner_y in [a1MaxY, a1MinY]:
                if (corner_x <= a2MaxX and corner_x >= a2MinX) and (corner_y <= a2MaxY and corner_y >= a2MinY):
                    overlap = True
                    # print("TRUE!!!!! ------------->")
                    # print(corner_x < a2MaxX, corner_x > a2MinX)
                    # print([corner_x, corner_y], a2MaxX, a2MinX, a2MaxY, a2MinY)
                    # print(a1MaxX, a1MinX, a1MaxY, a1MinY, '\n')
    return overlap

def addMoreEdges(nearby_platforms, G, Array, node, node_id, failures_c, failures_p, ids):
    for child in failures_c[node]:
        for id in ids:
            if child + "\n" + str(id) in G.nodes: 
                G.add_edge(node + "\n" + str(node_id), child + "\n" + str(id), weight=0.001)
    for parent in failures_p[node]:
        for id in ids:
            if parent + "\n" + str(id) in G.nodes: 
                G.add_edge(parent + "\n" + str(id), node + "\n" + str(node_id), weight=0.001)
    return G

def edgeReweight(G, edge):
    new_weight = input("Enter weight for (" + str(edge[0].replace("\n", " ")) + ", "  + str(edge[1].replace("\n", " ")) + ") (press \'Enter\' for default value)")
    if new_weight=='':
        new_weight=0.01
    G[edge[0]][edge[1]]['weight']=float(new_weight)
    return G

def get_y_Value(point1, point2, x):
    return ((point2[1]-point1[1])/(point2[0]-point1[0]))*(x-point1[0])+point1[1]

def get_min_max_vals(pnt2, pnt1, angle_radians, reverse):
    if reverse:
        pnt_hold = pnt1
        pnt1 = pnt2
        pnt2 = pnt_hold
    vector = pnt2 - pnt1
    # print('vect', vector)
    length = math.sqrt(vector[0]**2 + vector[1]**2)
    if vector[0] == 0.0 and vector[1] > 0: angle = math.pi/2
    elif vector[0] == 0.0 and vector[1] < 0: angle = 3*math.pi/2
    elif vector[0] == 0.0 and vector[1] == 0: angle = 0.0
    else: angle = math.atan(vector[1]/vector[0])
    if angle == 0 and vector[0] < 0:
        angle = math.pi
    if (angle > -math.pi*0.5 and angle < math.pi*0.5) and vector[0] < 0:
        angle += math.pi
    # print('angle', angle)
    new_vector1 = np.array([length * math.cos(angle - angle_radians), length * math.sin(angle - angle_radians)]) + np.array(pnt1)
    new_vector2 = np.array([length * math.cos(angle + angle_radians), length * math.sin(angle + angle_radians)]) + np.array(pnt1)
    if angle == math.pi and angle_radians==0:
        new_vector1 = np.array([length * math.cos(angle - angle_radians), 0]) + np.array(pnt1)
        new_vector2 = np.array([length * math.cos(angle + angle_radians), 0]) + np.array(pnt1)
    max_x = max([new_vector1[0], new_vector2[0], pnt1[0], pnt2[0]])
    min_x = min([new_vector1[0], new_vector2[0], pnt1[0], pnt2[0]])
    max_y = max([new_vector1[1], new_vector2[1], pnt1[1], pnt2[1]])
    min_y = min([new_vector1[1], new_vector2[1], pnt1[1], pnt2[1]])
    if max_x == min_x:
        if max_x < 0: max_x = 0
        elif min_x > 0: min_x = 0
    return max_x, min_x, max_y, min_y

# Create project class instance from yaml file
Array = Project(file='famodel/OntologySample600m.yaml')

# Create adjacency matrix from failure matrix
df = pd.read_excel("failure/failureData.xlsx", sheet_name="bMMatch")
arr = df.to_numpy()[:,1:]
nodeNames = df.to_numpy()[:, 0].flatten()

# Get initial failure probabilities for each failure mode and effect
probabilities, array_of_probs = getProbabilities("failure/failureProbabilities.xlsx", "Sheet3")


# Determine angle of clashing we are interested in
# angle_degree = float(input("What angle do you want to use? (in degrees) "))
# angle_radians = angle_degree/360 * math.pi * 2
angle_radians = 30/360 * math.pi * 2

# Initialize and create the dictionaries of the children, parents, and probabilities for each failure
failures_c = {}
failures_p = {}
probability_dict = {}

for i in range(arr.shape[0]):
    node_children = []
    node_parents = []
    for j in range(arr.shape[1]):
        if arr[i,j] > 0:
            node_children.append(nodeNames[j])
        if arr[j,i] > 0:
            node_parents.append(nodeNames[j])
    failures_c.update({nodeNames[i]: node_children})
    failures_p.update({nodeNames[i]: node_parents})
    probability_dict.update({nodeNames[i]: probabilities[i]})

# Create dictionary for each subsystem of failures
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
systems = {'turbine':nodeNames[turbine], 'platform':nodeNames[platform], 'mooringmooring':nodeNames[mooringmooring], 
           'rope':nodeNames[rope], 'chain':nodeNames[chain],'polyester':nodeNames[polyester],'mooringcable':nodeNames[mooringcable], 
           'cablemooring':nodeNames[mooringcable], 'mooring':nodeNames[mooring], 'connector':nodeNames[connector], 
           'weight':nodeNames[clump_weight], 'anchor':nodeNames[anchor], 'cable':nodeNames[cable], 'grid':nodeNames[grid], 
           'dynamic':nodeNames[dynamic], 'static':nodeNames[static], 'buoyancy':nodeNames[buoyancy], 'cablecable': [], 
           'sharedmooring':nodeNames[sharedmooring], 'sharedmooringcable':nodeNames[mooringcable], 'joints': nodeNames[joints],
           'cablesharedmooring':nodeNames[mooringcable], 'sharedmooringmooring':nodeNames[mooringmooring], 
           'mooringsharedmooring':nodeNames[mooringmooring], 'sharedanchor': nodeNames[sharedanchor]}


# Initialize graph, boolean for plotting, and list of probabilities
G = nx.DiGraph()
plot = False

# FIRST DEGREE NODES -------------------------------------------------------------------------------------------
for platform in Array.platformList:
    # print(platform, Array.platformList[platform].r)
    attachments = Array.platformList[platform].attachments
    nearby_platforms = []
    mooring_clashes = []
    cable_clashes = []
    num_cables = []

    # Create platform failure nodes
    for platform_failure in systems['platform']:
        G.add_node(platform_failure + "\n" + str(platform))
        G = addMoreEdges(nearby_platforms, G, Array, platform_failure, platform, failures_c, failures_p, [platform])

    # Create failure nodes
    for turbine_failure in systems['turbine']:
        G.add_node(turbine_failure + "\n" + str(platform))
        G = addMoreEdges(nearby_platforms, G, Array, turbine_failure, platform, failures_c, failures_p, [platform])

    # FIRST DEGREE EDGES -------------------------------------------------------------------------------------------
    for attach1 in attachments.keys():
        attach1_name = str(attachments[attach1]['id'])
        attach1_type = ''
        if 'mooring' in str(type(attachments[attach1]['obj'])): 
            if attachments[attach1]['obj'].shared: attach1_type = 'sharedmooring'
            else: attach1_type = 'mooring'
        elif 'cable' in str(type(attachments[attach1]['obj'])): 
            attach1_type = 'cable'
            num_cables.append(attachments[attach1]['obj'].id)

        # Create moroing/cable failure nodes
        for attach1_failure in systems[attach1_type]:
            original_name = attach1_name
            if 'connect' in attach1_failure: attach1_name = platform + attach1_name
            G.add_node(attach1_failure + "\n" + attach1_name)
            G = addMoreEdges(nearby_platforms, G, Array, attach1_failure, attach1_name, failures_c, failures_p, [platform, attach1_name])
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
                # print('could clash --', couldClash(G, clash_failure, attachments[attach1]['obj'], attachments[attach2]['obj']), '--',  attach1_name, attach2_name)
                if couldClash(G, clash_failure, attachments[attach1]['obj'], attachments[attach2]['obj'], reverse):
                    G.add_node(clash_failure + "\n" + clash_name)
                    G = addMoreEdges([attach1_name, attach2_name], G, Array, clash_failure, clash_name, failures_c, failures_p, [platform, attach1_name, attach2_name, clash_name])

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
                G.add_node(component_failure + "\n" + component_name)
                G = addMoreEdges(nearby_platforms, G, Array, component_failure, component_name, failures_c, failures_p, [platform, attach1])


        # SECOND ORDER NODES -------------------------------------------------------------------------------------------
        attached_to = attachments[attach1]['obj'].attached_to
        for attach1A in attached_to:
            attach1A_name = str(attach1A.id)

            # Create anchor failure nodes
            if 'anchor' in str(type(attach1A)):
                attach1A_type = 'anchor'
                if len(attach1A.attachments) > 1: attach1A_type = 'sharedanchor'
                for anchor_failure in systems[attach1A_type]:
                    G.add_node(anchor_failure + "\n" + attach1A_name)
                    G = addMoreEdges(nearby_platforms, G, Array, anchor_failure, attach1A_name, failures_c, failures_p, [platform, attach1_name, attach1A_name])
            
            # Create edges between platforms
            elif 'platform' in str(type(attach1A)): 
                attach1A_type = 'platform'
                attach1A_name = attach1A.id
                for platform_failure in systems['platform']:
                    G = addMoreEdges(nearby_platforms, G, Array, platform_failure, platform, failures_c, failures_p, [platform, attach1A_name])
            
            # Create substation/grid failure nodes
            elif 'substation' in str(type(attach1A)):
                attach1A_type = 'substation'
                for grid_failure in systems['grid']:
                    G.add_node(grid_failure + "\n" + attach1A_name)
                    G = addMoreEdges(nearby_platforms, G, Array, anchor_failure, attach1A_name, failures_c, failures_p, [platform, attach1_name, attach1A_name])

    if len(mooring_clashes) < 1:
        G.add_node(systems['mooringmooring'][0] + "\n" + str(platform))
        G = addMoreEdges(nearby_platforms, G, Array, systems['mooringmooring'][0], str(platform), failures_c, failures_p, [platform])

    if len(cable_clashes) < 1:
        for cable_num in num_cables:
            for clashing_failure in systems['cablemooring']:
                G.add_node(clashing_failure + "\n" + str(platform) + ' ' + str(cable_num))
                G = addMoreEdges(nearby_platforms, G, Array, clashing_failure, str(platform) + ' ' + str(cable_num), failures_c, failures_p, [platform])


print('\n')
for node in G.nodes:
    print(node.replace('\n', ' '))