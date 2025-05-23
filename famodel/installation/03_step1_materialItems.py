"""
Script to create material packages for installation.

Steps
-----
1. Load the Project object.
2. Group components into material packages for transport and installation.

Dependencies
-----------
famodel.project.Project
"""

'''Loads an ontology file and performs material itemization.'''
import matplotlib.pyplot as plt
import os
from fadesign.conceptual.conceptDesign import ConceptDesign
import numpy as np
from famodel.project import Project
from famodel.mooring.mooring import Mooring
from fadesign.conceptual import metrics as mtr
import pickle
import networkx as nx


def create_mtrlPkg(moor):
    '''Clusters components that needs to be mobilized/installed simultaneously. '''
    
    pkg = {}
    if moor.shared:
        # installation itemization for shared line
        for i, sec in enumerate(moor.dd['sections']):
            pkg[f"sec_{i}"] = {
                "obj": sec,
                "mass": sec['type']['m']*sec['L']/1e3,  # mass [t]
                "length": sec['L'],                     # length [m]
                "dependencies": [],
            }
            if i>1:
                pkg[f"sec_{i}"]['dependencies'].append(f"sec_{0}")
                
        for i, conn in enumerate(moor.dd['connectors']):
            if conn['m'] > 0:
                pkg[f"conn_{i}"] = {
                    "obj": conn,
                    "mass": conn['m']/1e3, # mass [t]
                    # "load": ?,  # pressure [t/m^2]    # NOTE: I am not checking load because clump weights could be easily sized in a way to reduce load if needed.
                    "dependencies": [f"sec_{0}"], 
                }
    else:
        for att in moor.attached_to:
            if type(att).__name__ == "Anchor":
                M = att.mass / 1e3        # mass [t]
                D = att.dd['design']['D'] # diameter [m]
                L = att.dd['design']['L'] # length [m]
                area_xy = D * L
                pkg = {
                    f"anchor_{att.id}": {
                        "obj": att,
                        "mass": M,
                        "load": M/area_xy,
                        "space": area_xy,
                        "dependencies": []
                    },
                }
                break

        for i, sec in enumerate(moor.dd['sections']):
            pkg[f"sec_{i}"] = {
                "obj": sec,
                "mass": sec['type']['m']*sec['L']/1e3,
                "length": sec['L'],
                "dependencies": [f"anchor_{att.id}"],
            }
        for i, conn in enumerate(moor.dd['connectors']):
            if conn['m'] > 0:
                pkg[f"conn_{i}"] = {
                    "obj": conn,
                    "mass": conn['m']/1e3,
                    "load": conn['m']/conn['v'],
                    "dependencies": [f"anchor_{att.id}"],
                }         
    return pkg


# CONCEPT CATEGORY & NAME
category    = "ILIA"  # ILIA | ILSA | SLIA | SLSA | hybrid 
conceptName = "3_ILIA_al[cr]_grd"

# FILE LOCATIONS
filePath     = os.path.dirname(os.path.abspath(__file__))
inputFile    = os.path.join(filePath, category, f"proj_{conceptName}.yaml")

# Load project
proj = Project(file=inputFile)

# Step1: create Material Items
mtrlPkgs = []
for moor in proj.mooringList.values():
    mtrlPkgs.append(create_mtrlPkg(moor))

output = os.path.join(filePath, "mtrlPkgs.pkl")
with open(output, 'wb') as f:
    pickle.dump(mtrlPkgs, f)