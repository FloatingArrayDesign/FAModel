"""
Script to integrate all steps into a complete installation simulation.

Steps
-----
1. Load material packages and action items.
2. Simulate the installation process using the InstallManager.

Dependencies
-----------
None
"""

'''an example of how install manager is used to register a vessel and port, and schedule an event and run (not finished yet).'''
# from fadesign.conceptual.installation.vessel import Vessel
# from fadesign.conceptual.installation.port import Port
from fadesign.installation.install_manager import InstallManager as IM
from pyproj import Proj, Transformer
import pandas as pd
import os

# FILE LOCATIONS
filePath     = os.path.dirname(os.path.abspath(__file__))
vessel1_file = os.path.join(filePath, "ahts.yaml")
port1_file   = os.path.join(filePath, "port_of_humboldt.yaml")

# Initialize
im = IM()

# Register Port
im.registerPort()
# Register Vessels
im.registerVessel(vessel1_file)

# Register vessel mob activity
im.scheduleEvent(im.now, im.vessels['ahts1'], action='mob', params={"portLocation": [0, 0]})
# Register 
im.run()