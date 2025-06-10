"""
an example of how install manager is used to register a vessel and port, and schedule an event and run (not finished yet).
This is loads in the vessels from 04 and their respective actions.

Steps
-----
1. Load material packages and action items.
2. Simulate the installation process using the InstallManager.

Dependencies
-----------
None
"""

# from fadesign.conceptual.installation.vessel import Vessel
# from fadesign.conceptual.installation.port import Port
from .install_manager import InstallManager as IM
import os
import pickle 

# FILE LOCATIONS
filePath     = os.path.join(os.path.dirname(os.path.abspath(__file__)))
port1_file   = os.path.join(filePath, "input_files/agent_yamls", "humboldt_bay.yaml")
vesselsFile = os.path.join(filePath, "temp_files", "vessels.pkl")

# Load vessel info
with open(vesselsFile, 'rb') as f:
    vesselDict = pickle.load(f)

# Initialize
im = IM()

# Register Port
im.registerPort(port1_file)

# Register Vessels (set up for a list of vessels, currently only one)
for keys in vesselDict.keys():
    vessel = vesselDict[keys]
    # Register the vessel with the InstallManager
    im.registerVessel(vessel)

# Register vessel action (INCOMPLETE) - have not finished how to call the actions in a time loop, how time is handled, and stnadard action input/outputs. Have not considered port locations. 
# im.scheduleEvent(im.now, im.vesselList['AHTS'], action='mob', params={"portLocation": [0, 0]}) # Incomplete. This calls vessel.mob(), which is a function that's decoupled from the other actions in the class. 
im.scheduleEvent(im.now, im.vesselList['AHTS'], action=vessel.mobilize, params={"portLocation": [0, 0]}) # Incomplete. This calls vessel.mob(), which is a function that's decoupled from the other actions in the class. 
# Register 
im.run()