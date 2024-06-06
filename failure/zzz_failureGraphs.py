import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from failure.failureProbabilities import *
from failure.twoTurbineCaseStudy import *
from famodel.project import Project


def addMoreEdges(nearby_platforms, G, Array, node, node_id, failures_c, failures_p, ids):
    node_index = np.where(nodeNames == node)[0][0]
    platforms = list(Array.platformList.keys())
    if len(nearby_platforms) > 1:
        platforms = nearby_platforms
    cables = list(Array.cableList.keys())
    all_ids = []
    combos = []
    for p in platforms:
        for c in cables:
            if c in p:
                combos.append(str(p) + ',' + str(c))
                all_ids.append(c)
                for m in Array.platformList[p]:
                    combos.append(str(m) + ',' + str(c))
                    all_ids.append(str(m))
        print("ATTACHMENTS: ", Array.platformList[p].attachments)
        for a in Array.platformList[p].anchorList:
            all_ids.append(a)
    all_ids = platforms+all_ids+list(combos)
    for child in failures_c[node]:
        child_index = np.where(nodeNames == child)[0][0]
        original_ids = ids
        if arr[node_index, child_index] > 1.5:
            ids = all_ids
        for id in ids:
            if child + "\n" + str(id) in G.nodes: 
                G.add_edge(node + "\n" + str(node_id), child + "\n" + str(id), weight=0.001)
        ids = original_ids
    for parent in failures_p[node]:
        parent_index = np.where(nodeNames == parent)[0][0]
        original_ids = ids
        if arr[parent_index, node_index] > 1.5:
            ids = all_ids
        for id in ids:
            if parent + "\n" + str(id) in G.nodes: 
                G.add_edge(parent + "\n" + str(id), node + "\n" + str(node_id), weight=0.001)
        ids = original_ids
    return G

def edgeReweight(G, edge):
    new_weight = input("Enter weight for (" + str(edge[0].replace("\n", " ")) + ", "  + str(edge[1].replace("\n", " ")) + ") (press \'Enter\' for default value)")
    if new_weight=='':
        new_weight=0.01
    G[edge[0]][edge[1]]['weight']=float(new_weight)
    return G

def get_y_Value(point1, point2, x):
    return ((point2[1]-point1[1])/(point2[0]-point1[0]))*(x-point1[0])+point1[1]

def get_min_max_vals(pnt1, pnt2, angle_radians):
    vector = pnt2 - pnt1
    length = math.sqrt(vector[0]**2 + vector[1]**2)
    if vector[0] == 0.0 and vector[1] > 0: angle = math.pi/2
    elif vector[0] == 0.0 and vector[1] < 0: angle = 3*math.pi/2
    elif vector[0] == 0.0 and vector[1] == 0: angle = 0.0
    else: angle = math.atan(vector[1]/vector[0])
    if angle == 0 and vector[0] < 0:
        angle = math.pi
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
df = pd.read_excel("/Users/eslack/Documents/Code/ExcelFiles/failureData.xlsx", sheet_name="bMMatch")
arr = df.to_numpy()[:,1:]
nodeNames = df.to_numpy()[:, 0].flatten()

# Get initial failure probabilities for each failure mode and effect
probabilities, array_of_probs = getProbabilities("failureProbabilities.xlsx", "Sheet3")



# Determine if using the twoTurbine case study model
tT_input = input("Would you like to calculate for the twoTurbine case study? ")
if 'y' in tT_input: twoTurbine = True
else: twoTurbine = False

# Determine angle of clashing we are interested in
angle_degree = float(input("What angle do you want to use? (in degrees) "))
angle_radians = angle_degree/360 * math.pi * 2

# Determine what the nearby platforms are (True ==> all other platforms, False ==> physically close turbines)
nba_input = input("Would you like one turbine to directly affect all turbines? ")
nba = False
if 'y' in nba_input.lower():
    nba = True



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
clashing = [12,13,14]
mooring_materials = [33,34,35]
mooring = [15,16,17, 18]
connector = [36]
clump_weight = [37]
anchor = [19,20,38,39]
cable = [21,22,23,41,42, 45,46]
cable_mat = [43,44]
grid = [24,25]
buoyancy = [40]

systems = {'turbine':nodeNames[turbine], 'platform':nodeNames[platform], 'clashing':nodeNames[clashing], 'moor_mat':nodeNames[mooring_materials], 
           'mooring':nodeNames[mooring], 'connector':nodeNames[connector], 'weight':nodeNames[clump_weight], 
           'anchor':nodeNames[anchor], 'cable':nodeNames[cable], 'grid':nodeNames[grid], 'cable_mat':nodeNames[cable_mat],
           'buoyancy':nodeNames[buoyancy]}


# Initialize graph, boolean for plotting, and list of probabilities
G = nx.DiGraph()
plot = False
probabilities = []

for platform in Array.platformList:
    # Initialize list of mooring-mooring clashing nodea and list of mooring-cable and anchor-cable clashing nodes
    mooring_clashes = []
    cable_clashes = []

    # Determine the nearby platforms
    nearby_platforms = []
    for platform2 in Array.platformList:
        platform_xy = Array.platformList[platform].r
        platform2_xy = Array.platformList[platform2].r
        dist_btwn_turbines = math.sqrt((platform_xy[0] - platform2_xy[0])**2 + (platform_xy[1] - platform2_xy[1])**2)
        if dist_btwn_turbines <= (1600 * math.sqrt(2)):
            nearby_platforms.append(platform2)
    if 'y' in nba_input:
        nearby_platforms = []

    # Create platform nodes
    for platform_failure in systems['platform']:
        if not(platform_failure + "\n" + str(platform) in list(G.nodes)): 
            probabilities.append(probability_dict[platform_failure])
        G.add_node(platform_failure + "\n" + str(platform))
        G = addMoreEdges(nearby_platforms, G, Array, platform_failure, platform, failures_c, failures_p, [platform])

    # Create turbine nodes
    for turbine_failure in systems['turbine']:
        if not(turbine_failure + "\n" + str(platform) in list(G.nodes)): probabilities.append(probability_dict[turbine_failure])
        G.add_node(turbine_failure + "\n" + str(platform))
        G = addMoreEdges(nearby_platforms, G, Array, turbine_failure, platform, failures_c, failures_p, [platform])
    
    print(Array.platformList[platform].getTopLevelEdge)
    for anchor in Array.platformList[platform].anchorList:
        # Create anchor nodes
        for anchor_failure in systems['anchor']:
            if len(Array.platformList[platform].anchorList[anchor].mooringList) > 1 and "ingle" in anchor_failure: continue
            elif len(Array.platformList[platform].anchorList[anchor].mooringList) <= 1 and "hared" in anchor_failure: continue
            if not(anchor_failure + "\n" + str(anchor) in list(G.nodes)): probabilities.append(probability_dict[anchor_failure])
            G.add_node(anchor_failure + "\n" + str(anchor))
            G = addMoreEdges(nearby_platforms, G, Array, anchor_failure, anchor, failures_c, failures_p, [platform, anchor])

        for mooring in Array.platformList[platform].anchorList[anchor].mooringList:
            # Create mooring nodes based on specific materials
            for section in Array.platformList[platform].anchorList[anchor].mooringList[mooring].dd["sections"]:
                for material in systems['moor_mat']:
                    if ((("ope" in section["type"]["material"]) and ("ire" in material)) or (("oly" in section["type"]["material"]) and ("ynth" in material)) )or (("hain" in section["type"]["material"]) and ("hain" in material)):
                        if not(material + "\n" + str(mooring) in list(G.nodes)): probabilities.append(probability_dict[material])
                        G.add_node(material + "\n" + str(mooring))
                        G = addMoreEdges(nearby_platforms, G, Array, material, mooring, failures_c, failures_p, [platform, anchor, mooring])

            # Create other mooring nodes
            for mooring_failure in systems['mooring']:
                if (Array.platformList[platform].anchorList[anchor].mooringList[mooring].shared) and ("ing line non" in mooring_failure.replace("\n", " ")): continue
                elif (not Array.platformList[platform].anchorList[anchor].mooringList[mooring].shared) and ("ared line" in mooring_failure.replace("\n", " ")): continue
                else:
                    if not(mooring_failure + "\n" + str(mooring) in list(G.nodes)): probabilities.append(probability_dict[mooring_failure])
                    G.add_node(mooring_failure + "\n" + str(mooring))
                    G = addMoreEdges(nearby_platforms, G, Array, mooring_failure, mooring, failures_c, failures_p, [platform, anchor, mooring])

            # Create mooring-clashing nodes
            for mooring2 in Array.platformList[platform].anchorList[anchor].mooringList:
                if not(mooring == mooring2):
                    print('not shared', mooring, mooring2)
                    mooring_pnt1 = Array.platformList[platform].anchorList[anchor].mooringList[mooring].rA[:2]
                    mooring_pnt2 = Array.platformList[platform].anchorList[anchor].mooringList[mooring].rB[:2]
                    mooring2_pnt1 = Array.platformList[platform].anchorList[anchor].mooringList[mooring2].rA[:2]
                    mooring2_pnt2 = Array.platformList[platform].anchorList[anchor].mooringList[mooring2].rB[:2]
                    mMaxX, mMinX, mMaxY, mMinY = get_min_max_vals(mooring_pnt1, mooring_pnt2, angle_radians)
                    cMaxX, cMinX, cMaxY, cMinY = get_min_max_vals(mooring2_pnt1, mooring2_pnt2, angle_radians)
                    x_overlap = (cMaxX < mMaxX and cMaxX > mMinX) or (cMinX > mMinX and cMinX < mMaxX)
                    y_overlap = (cMaxY < mMaxY and cMaxY > mMinY) or (cMinY > mMinY and cMinY < mMaxY)
                    x_on_top = (cMaxX <= mMaxX and cMaxX >= mMinX) and (cMinX >= mMinX and cMinX <= mMaxX)
                    y_on_top = (cMaxY <= mMaxY and cMaxY >= mMinY) and (cMinY >= mMinY and cMinY <= mMaxY)
                    if (x_overlap and y_overlap) or (x_on_top and y_on_top):
                        if not(systems['clashing'][0] + "\n" + str(mooring) + "," + str(mooring2) in list(G.nodes)): 
                            probabilities.append(probability_dict[systems['clashing'][0]])
                        G.add_node(systems['clashing'][0] + "\n" + str(mooring) + "," + str(mooring2))
                        G = addMoreEdges(nearby_platforms, G, Array, systems['clashing'][0], str(mooring) + str(mooring2), failures_c, failures_p, [platform, anchor, mooring, mooring2])
                        mooring_clashes.append(str(mooring) + str(mooring2))
            
            # Assign connector failures
            for connector in Array.platformList[platform].anchorList[anchor].mooringList[mooring].dd['connectors']:
                for connector_failure in systems['connector']:
                    if not(connector_failure + "\n" + str(mooring) in list(G.nodes)): probabilities.append(probability_dict[connector_failure])
                    G.add_node(connector_failure + "\n" + str(mooring))
                    G = addMoreEdges(nearby_platforms, G, Array, connector_failure, mooring, failures_c, failures_p, [platform, anchor, mooring, connector])

            # Create cable-clashing nodes
            for cable in Array.cableList:
                if platform in Array.cableList[cable].dd['platforms']:
                    
                    # Cable failures
                    for cable_failure in systems['cable']:
                        if 'disconnect' in cable_failure.replace('\n',' ').lower():
                            if not(cable_failure + "\n" + str(cable) + ' ' + str(platform) in list(G.nodes)): 
                                probabilities.append(probability_dict[cable_failure])
                            G.add_node(cable_failure + "\n" + str(platform)+','+str(cable))
                            G = addMoreEdges(nearby_platforms, G, Array, cable_failure, str(platform)+','+str(cable), failures_c, failures_p, [platform, anchor, mooring, cable])
                        else:
                            if not(cable_failure + "\n" + str(cable) in list(G.nodes)): 
                                probabilities.append(probability_dict[cable_failure])
                            G.add_node(cable_failure + "\n" + str(cable))
                            G = addMoreEdges(nearby_platforms, G, Array, cable_failure, cable, failures_c, failures_p, [platform, anchor, mooring, cable])

                    # Cable buoyancy modules and materials
                    for cable_section in Array.cableList[cable].dd['cables']:
                        # Create buoyancy module nodes (one node per cable)
                        if 'buoyancy_sections' in cable_section.dd:
                            n_modules = cable_section.dd['buoyancy_sections'][0]['N_modules']
                            for buoy_failure in systems['buoyancy']:
                                if not(buoy_failure + "\n" + str(cable) in list(G.nodes)): 
                                    probabilities.append(probability_dict[buoy_failure])
                                G.add_node(buoy_failure + "\n" + str(cable))
                                G = addMoreEdges(nearby_platforms, G, Array, buoy_failure, cable, failures_c, failures_p, [platform, anchor, mooring, cable])

                        # Create either dynamic or static cable failure node (one of each per cable at a maximum)
                        if cable_section.dd['cable_type']['dynamic']:
                            if not(systems['cable_mat'][0] + "\n" + str(cable) in list(G.nodes)): 
                                probabilities.append(probability_dict[systems['cable_mat'][0]])
                            G.add_node(systems['cable_mat'][0] + "\n" + str(cable))
                            G = addMoreEdges(nearby_platforms, G, Array, systems['cable_mat'][0], cable, failures_c, failures_p, [platform, anchor, mooring, cable])
                        else:
                            if not(systems['cable_mat'][1] + "\n" + str(cable) in list(G.nodes)): 
                                probabilities.append(probability_dict[systems['cable_mat'][1]])
                            G.add_node(systems['cable_mat'][1] + "\n" + str(cable))
                            G = addMoreEdges(nearby_platforms, G, Array, systems['cable_mat'][1], cable, failures_c, failures_p, [platform, anchor, mooring, cable])

                    # Cable clashing failures
                    cable_pnt1 = np.array(Array.platformList[Array.cableList[cable].dd['platforms'][0]].r)
                    if not ('substation' in Array.cableList[cable].dd['platforms'][1].lower()):
                        cable_pnt2 = np.array(Array.platformList[Array.cableList[cable].dd['platforms'][1]].r)
                    else: cable_pnt2 = np.array(Array.substationList[Array.cableList[cable].dd['platforms'][1]].r)
                    mooring_pnt1 = Array.platformList[platform].anchorList[anchor].mooringList[mooring].rA[:2]
                    mooring_pnt2 = Array.platformList[platform].anchorList[anchor].mooringList[mooring].rB[:2]
                    cMaxX, cMinX, cMaxY, cMinY = get_min_max_vals(cable_pnt1, cable_pnt2, angle_radians)
                    mMaxX, mMinX, mMaxY, mMinY = get_min_max_vals(mooring_pnt1, mooring_pnt2, angle_radians)
                    x_overlap = (cMaxX < mMaxX and cMaxX > mMinX) or (cMinX > mMinX and cMinX < mMaxX)
                    y_overlap = (cMaxY < mMaxY and cMaxY > mMinY) or (cMinY > mMinY and cMinY < mMaxY)
                    x_on_top = (cMaxX <= mMaxX and cMaxX >= mMinX) and (cMinX >= mMinX and cMinX <= mMaxX)
                    y_on_top = (cMaxY <= mMaxY and cMaxY >= mMinY) and (cMinY >= mMinY and cMinY <= mMaxY)
                    if (x_overlap and y_overlap) or ((cMaxX == mMaxX and cMinX == mMinX)and(cMaxY == mMaxY and cMinY == mMinY)):
                        for cable_clashing_failure in systems['clashing'][1:]:
                            if not(cable_clashing_failure + "\n" + str(cable)+str(mooring) in list(G.nodes)): 
                                probabilities.append(probability_dict[cable_clashing_failure])
                            G.add_node(cable_clashing_failure + "\n" + str(cable)+str(mooring))
                            G = addMoreEdges(nearby_platforms, G, Array, cable_clashing_failure, str(cable)+str(mooring), failures_c, failures_p, [platform, cable, anchor, mooring, str(cable)+str(mooring)])
                            cable_clashes.append(str(cable)+str(mooring))
                    if len(cable_clashes) < 1:
                        if not(systems['clashing'][1] + "\n" + str(platform) + "," + str(cable) in list(G.nodes)): 
                            probabilities.append(probability_dict[systems['clashing'][1]])
                        G.add_node(systems['clashing'][1] + "\n" + str(platform) + "," + str(cable))
                        G = addMoreEdges(nearby_platforms, G, Array, systems['clashing'][1], str(platform) + ',' + str(cable), failures_c, failures_p, [platform, anchor, cable, mooring, str(cable) + ' ' + str(platform)])
                        if not(systems['clashing'][2] + "\n" + str(platform) + ',' + str(cable) in list(G.nodes)):
                            probabilities.append(probability_dict[systems['clashing'][2]])
                        G.add_node(systems['clashing'][2] + "\n" + str(platform) + ',' + str(cable))
                        G = addMoreEdges(nearby_platforms, G, Array, systems['clashing'][2], str(platform) + ',' + str(cable), failures_c, failures_p, [platform, anchor, cable, mooring, str(cable) + ' ' + str(platform)])
    
    # If there are no specific mooring-mooring clashing, create a generic mooring-mooring clashing node for the platform
    if len(mooring_clashes) < 1:
        fail_list = [platform, anchor]
        for a_variable in Array.platformList[platform].mooringList.keys():
            fail_list.append(a_variable)
        if not(systems['clashing'][0] + "\n" + str(platform) in list(G.nodes)): 
            probabilities.append(probability_dict[systems['clashing'][0]])
        G.add_node(systems['clashing'][0] + "\n" + str(platform))
        G = addMoreEdges(nearby_platforms, G, Array, systems['clashing'][0], platform, failures_c, failures_p, fail_list)

    # Shared mooring line failures
    for mooring in Array.platformList[platform].mooringList:
        if Array.platformList[platform].mooringList[mooring].shared:
            for mooring_failure in systems['mooring']:
                if not twoTurbine or 'shared line' in mooring_failure.replace('\n', ' ').lower():
                    if not(mooring_failure + "\n" + str(mooring) in list(G.nodes)): 
                        probabilities.append(probability_dict[mooring_failure])
                    G.add_node(mooring_failure + "\n" + str(mooring))
                    G = addMoreEdges(nearby_platforms, G, Array, mooring_failure, mooring, failures_c, failures_p, [platform, mooring])

            # Create mooring nodes based on specific materials
            for section in Array.platformList[platform].mooringList[mooring].dd["sections"]:
                for material in systems['moor_mat']:
                    if ((("ope" in section["type"]["material"]) and ("ire" in material)) or (("oly" in section["type"]["material"]) and ("ynth" in material)) )or (("hain" in section["type"]["material"]) and ("hain" in material)):
                        if not(material + "\n" + str(mooring) in list(G.nodes)): probabilities.append(probability_dict[material])
                        G.add_node(material + "\n" + str(mooring))
                        G = addMoreEdges(nearby_platforms, G, Array, material, mooring, failures_c, failures_p, [platform, mooring])

            # Create mooring clashing nodes
            for mooring2 in Array.platformList[platform].anchorList[anchor].mooringList:
                if not(mooring == mooring2):
                    mooring_pnt1 = np.array(Array.platformList[platform].mooringList[mooring].rA[:2])
                    mooring_pnt2 = np.array(Array.platformList[platform].mooringList[mooring].rB[:2])
                    mooring2_pnt1 = np.array(Array.platformList[platform].mooringList[mooring2].rA[:2])
                    mooring2_pnt2 = np.array(Array.platformList[platform].mooringList[mooring2].rB[:2])
                    mMaxX, mMinX, mMaxY, mMinY = get_min_max_vals(mooring_pnt1, mooring_pnt2, angle_radians)
                    cMaxX, cMinX, cMaxY, cMinY = get_min_max_vals(mooring2_pnt1, mooring2_pnt2, angle_radians)
                    x_overlap = (cMaxX < mMaxX and cMaxX > mMinX) or (cMinX > mMinX and cMinX < mMaxX)
                    y_overlap = (cMaxY < mMaxY and cMaxY > mMinY) or (cMinY > mMinY and cMinY < mMaxY)
                    x_on_top = (cMaxX <= mMaxX and cMaxX >= mMinX) and (cMinX >= mMinX and cMinX <= mMaxX)
                    y_on_top = (cMaxY <= mMaxY and cMaxY >= mMinY) and (cMinY >= mMinY and cMinY <= mMaxY)
                    if (x_overlap and y_overlap) or (x_on_top or y_on_top):
                        if not(systems['clashing'][0] in list(G.nodes)): probabilities.append(probability_dict[systems['clashing'][0]])
                        G.add_node(systems['clashing'][0] + "\n" + str(mooring) + "," + str(mooring2))
                        G = addMoreEdges(nearby_platforms, G, Array, systems['clashing'][0], str(mooring) + "," + str(mooring2), failures_c, failures_p, [platform, mooring, mooring2])
                        mooring_clashes.append(str(mooring) + str(mooring2))

            # Cable clashing failures
            for cable in Array.cableList:
                if platform in Array.cableList[cable].dd['platforms']:
                    cable_pnt1 = np.array(Array.platformList[Array.cableList[cable].dd['platforms'][0]].r)
                    if not ('substation' in Array.cableList[cable].dd['platforms'][1].lower()):
                        cable_pnt2 = np.array(Array.platformList[Array.cableList[cable].dd['platforms'][1]].r)
                    else: cable_pnt2 = np.array(Array.substationList[Array.cableList[cable].dd['platforms'][1]].r)
                    
                    if not ('substation' in mooring[0].lower()):mooring_pnt1 = np.array(Array.platformList[mooring[0]].r)
                    else: mooring_pnt1 = np.array(Array.substationList[mooring[0]].r)
                    if not ('substation' in mooring[1].lower()): mooring_pnt2 = np.array(Array.platformList[mooring[1]].r)
                    else: mooring_pnt2 = np.array(Array.substationList[mooring[1]].r)

                    x_overlap = (cMaxX < mMaxX and cMaxX > mMinX) or (cMinX > mMinX and cMinX < mMaxX)
                    y_overlap = (cMaxY < mMaxY and cMaxY > mMinY) or (cMinY > mMinY and cMinY < mMaxY)
                    cMaxX, cMinX, cMaxY, cMinY = get_min_max_vals(cable_pnt1, cable_pnt2, angle_radians)
                    mMaxX, mMinX, mMaxY, mMinY = get_min_max_vals(mooring_pnt1, mooring_pnt2, angle_radians)
                    if (x_overlap and y_overlap) or ((cMaxX == mMaxX and cMinX == mMinX)and(cMaxY == mMaxY and cMinY == mMinY)):
                        for cable_clashing_failure in systems['clashing'][1:]:
                            if not (systems['clashing'][1] in list(G.nodes)): probabilities.append(probability_dict[systems['clashing'][1]])
                            G.add_node(systems['clashing'][1] + "\n" + str(mooring) + "," + str(cable))
                            G = addMoreEdges(nearby_platforms, G, Array, systems['clashing'][1], str(mooring) + "," + str(cable), failures_c, failures_p, [platform, cable, mooring, str(mooring) + "," + str(cable)])
                            G = addMoreEdges(nearby_platforms, G, Array, systems['clashing'][1], str(platform) + "," + str(cable), failures_c, failures_p, [platform, cable, mooring, str(platform) + "," + str(cable)])
                            cable_clashes.append(str(cable)+str(mooring))

    # Grid/substation failures
    for grid_failure in systems['grid']:
        if not(grid_failure in list(G.nodes)):
            probabilities.append(probability_dict[grid_failure])
        G.add_node(grid_failure + "\n" + "")
        G = addMoreEdges(nearby_platforms, G, Array, grid_failure, "", failures_c, failures_p, [platform])



# If interested in the twoTurbine case study, alter a few of the failure nodes
if twoTurbine:
    list_of_failures_to_rid = ['Tether & anchor systems array_cable10', 'Cable protection system array_cable10', 
                           'Terminations array_cable10', 'Offshore joints array_cable10']
    for failure in list(G.nodes):
        if 'connector' in failure.replace('\n', ' ').lower():
            G.remove_node(failure)
        if 'buoyancy modules' in failure.replace('\n', ' ').lower():
            G.remove_node(failure)
        if failure.replace('\n', ' ') in list_of_failures_to_rid:
            G.remove_node(failure)



# Print the number of nodes and (if the user wants) the list of nodes
print('\nNumber of nodes -',len(G.nodes),'\n')
user_input3 = input("Would you like to see the list of failures? ")
if 'y' in user_input3.lower() or 'rue' in user_input3.lower():
    for edge in list(G.nodes):
        print(edge.replace("\n", " "))

# Print the number of nodes and (if the user wants) the list of edges
print('\nNumber of edges -',len(G.edges),'\n')
user_input4 = input("Would you like to see the list of edges? ")
if 'y' in user_input4.lower() or 'rue' in user_input4.lower():
    itervar = 0
    for edge in list(G.edges):
        print(edge)
        itervar += 1
        if (itervar + 1) % 1000 == 0: user_input45 = input("Continue? ")



# If the user wants to input probabilities for specific edges, reweight edges based on the user's inputs
user_inputs = input("Would you like to input probabilities into adjacency matrix? ")
twoTurbine_calculationType = False
if (user_inputs == 'y' or user_inputs == 'yes') or user_inputs == 'True':
    twoTurbine_calculationType = True
    for i in range(len(G.edges)):
        edge = list(G.edges)[i]
        if ('rift off' in edge[0].replace("\n", " ") or 'ncreased' in edge[0].replace("\n", " ")) or 'ynamics' in edge[0].replace("\n", " "):
            G = edgeReweight(G, edge)
        elif ('apsize' in edge[0].replace("\n", " ") or '-cable' in edge[0].replace("\n", " ")) or 'ing line non' in edge[0].replace("\n", " "):
            G = edgeReweight(G, edge)
        elif ('ragging' in edge[0].replace("\n", " ") or 'hain' in edge[0].replace("\n", " ")) or 'ire rope' in edge[0].replace("\n", " "):
            G = edgeReweight(G, edge)
        elif ('ynthetic' in edge[0].replace("\n", " ") or 'able profile' in edge[0].replace("\n", " ")) or 'ared line' in edge[0].replace("\n", " "):
            G = edgeReweight(G, edge)
        elif ('load on cable' in edge[0].replace("\n", " ") or 'eight' in edge[0].replace("\n", " ")):
            G = edgeReweight(G, edge)



# Ask user if they are ready to continue to Bayesian network calculations (if not, quit)
continue_input = input("Ready to continue? ")
if 'n' in continue_input.lower():
    quit()


# Bayesian network calculation
arr = nx.to_numpy_array(G)
nodeNames = np.reshape(np.array(list(G.nodes)), (len(list(G.nodes)), ))

poc = "child"
filename = "turbineInference_" + poc + "FAModelTwoTurbine" + "_" + str(nba_input) + ".xlsx"
for start_component in range(1,arr.shape[0]+1): # Iterate through each failure mode/effect in turbine
    print(start_component)

with pd.ExcelWriter(filename) as writer:
    all_probabilities = np.zeros(arr.shape) # Initialize a large array to put all the pairwise probabilities in
    for start_component in range(1,arr.shape[0]+1): # Iterate through each failure mode/effect in turbine
        a = arr.copy()
        non = nodeNames
        K, a, g, e, m, non = breadth_first_multi(a, nodeNames, [start_component], poc) # Generate tree for Bayesian network
        prblts = [] # Initialize array of node probabilities (in order of appearance in graph)
        for node in non:
            node_index = np.where(nodeNames == node)[0][0]
            prblts.append(probabilities[node_index]) # Add nodes to array of node probabilities
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
                probabilitiy_table[0][node] = probabilities[np.where(nodeNames == nodeNamesArray[int(node)])[0][0]]
                probabilitiy_table[1][node] = 1 - probabilities[np.where(nodeNames == nodeNamesArray[int(node)])[0][0]]
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

    # Write array to dataframe
    df3 = pd.DataFrame(all_probabilities)
    df3.to_excel(writer, sheet_name="allProbs")


# If we want to plot the network, plot it now
if plot:
    nx.draw_networkx(G)
    plt.show()