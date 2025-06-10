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
import yaml
import pickle 
from . import vessel as V

# FILE LOCATIONS
filePath     = os.path.dirname(os.path.abspath(__file__))
itemFile     = os.path.join(filePath, "temp_files", "mtrlPkgs.pkl")

# Load package info
with open(itemFile, 'rb') as f:
    pkgs = pickle.load(f)

# Load vessel info
vesselsFile   = os.path.join(filePath, "input_files/agent_yamls", "vesselDesc.yaml")
vesselNames = ['AHTS']

with open(vesselsFile) as file:
    vesselDisc = yaml.load(file, Loader=yaml.FullLoader)

# Create a dictionary to hold vessel objects with the vessel names as keys
vesselDict = {name: V.Vessel(vesselDisc = vesselDisc[name]) for name in vesselNames}  # Filter vessel descriptions to only include specified names


# Action list: for each vessel
for name in vesselNames:
    vessel = vesselDict[name] # Get the vessel object

    # Mobilization (needs to initialize before transport because transport depends on mobilization ActionItems)
    pkg = pkgs[0] # example
    mobilize = vessel.get_mobilize_action(pkg = pkg)
    print("Vessel State: \n", vessel.state)

    # Transport
    distance2port = 250  # km [not sure where this information would be handled]
    transit_to = vessel.get_transit_to_action(distance2port = distance2port)
    print("Vessel State: \n", vessel.state)

    # Installation
    pkg = pkgs[0] # example
    install = vessel.get_install_action(pkg = pkg)
    print("Vessel State: \n", vessel.state)

    mobilize.visualize()
    plt.show() 

    transit_to.visualize()
    plt.show() 

    install.visualize()
    plt.show()

# Pickle the vessels for loading later
output = os.path.join(filePath, "temp_files", "vessels.pkl")
with open(output, 'wb') as f:
    pickle.dump(vesselDict, f)