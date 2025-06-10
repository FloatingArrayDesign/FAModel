"""
loads material package saved in step1, loads vessel description, and creates action items for transport, mobilization, and installation.

Steps
-----
1. Define tasks for each stage of the installation process.
2. Set dependencies between tasks.

Dependencies
-----------
None
"""


import matplotlib.pyplot as plt
import os
from fadesign.conceptual.conceptDesign import ConceptDesign
import numpy as np
from famodel.project import Project
from famodel.mooring.mooring import Mooring
from fadesign.conceptual import metrics as mtr
import yaml
import networkx as nx
import pickle 
from . import install_helpers as inst


# FILE LOCATIONS
filePath     = os.path.dirname(os.path.abspath(__file__))
itemFile     = os.path.join(filePath, "temp_files", "mtrlPkgs.pkl")

with open(itemFile, 'rb') as f:
    pkgs = pickle.load(f)

vesselFile   = os.path.join(filePath, "input_files/agent_yamls", "vesselDesc.yaml")
vesselNames = ['AHTS']

with open(vesselFile) as file:
    vesselDisc = yaml.load(file, Loader=yaml.FullLoader)


for vesselName in vesselNames:
    vessel = vesselDisc[vesselName]
    vessel['state'] = {
        'remaining_cargo': vessel['storage_specs']['max_cargo'],
        'remaining_deck_space': vessel['storage_specs']['max_deck_space'],
        'remaining_spool_capacity': vessel['storage_specs']['max_spool_capacity'],
        'assigned_materials': []
        }
    
# Action list:
# Transport
distance2port = 250  # km [not sure where this information would be handled]
transport_V1 = inst.tranportTo_actionItem(vessel, distance2port)
# Mobilization
pkg = pkgs[0] # example
mobilize_V1, vessel = inst.mobilizeM_actionItem(vessel, pkg)

# TODO: This seems like the same as vessel.mobilize? Are these two different approaches that are duplicates? How do we map this out? 

print(vessel['state'])
# Installation
pkg = pkgs[0] # example
install_V1, vessel = inst.install_actionItem(vessel, pkg)
print(vessel['state'])

inst.visualizeAction(transport_V1)
plt.show() 

inst.visualizeAction(mobilize_V1)
plt.show() 

inst.visualizeAction(install_V1)
plt.show()
