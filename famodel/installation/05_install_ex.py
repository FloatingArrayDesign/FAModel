"""
an example of how install manager is used to register a vessel and port, and schedule an event and run (not finished yet).
This is independent of 04_step2_actionItems.py and 03_step1_materialItems.py

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
from pyproj import Proj, Transformer
import pandas as pd
import os

# FILE LOCATIONS
filePath     = os.path.dirname(os.path.abspath(__file__))
vessel1_file = os.path.join(filePath, "ahts.yaml")
port1_file   = os.path.join(filePath, "humboldt_bay.yaml")

# Initialize
im = IM()

# Register Port
im.registerPort(port1_file)
# Register Vessels
im.registerVessel(vessel1_file)

# Register vessel mob activity
im.scheduleEvent(im.now, im.vesselList['3gs'], action='mob', params={"portLocation": [0, 0]}) # TODO: how does this relate to 04? In 04 you create actions. Here you schedule events and it seems like the action is defined as a method of the class (not a separate object).
# Register 
im.run()