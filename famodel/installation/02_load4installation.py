"""
This example loads a project class and conducts some analysis and performance metrics. 

NOTE: If the project yaml file does not exist, please run 01_unload_driver.py before this 
so you have a project that you can load here in this example. The project yaml file generated 
by the unload driver is named proj_{<concept name>}. Insert that in the inputFile in this script.
"""

import matplotlib.pyplot as plt
import os
from fadesign.conceptual.conceptDesign import ConceptDesign
import numpy as np
from famodel.project import Project
from fadesign.conceptual import metrics as mtr

# CONCEPT CATEGORY & NAME
category    = "ILIA"  # ILIA | ILSA | SLIA | SLSA | hybrid 
conceptName = "3_ILIA_1P2C_grd_mini"

# FILE LOCATIONS
filePath     = os.path.dirname(os.path.abspath(__file__))
inputFile    = os.path.join(filePath, category, f"proj_{conceptName}.yaml")
mtricFile    = os.path.join(filePath, category, f"{conceptName}_mtr.xlsx")
windFile     = os.path.join(filePath, "sites/wind_data/humboldt_rose_1.csv")
currentFile  = os.path.join(filePath, "sites/surface_currents/HumboldtBay_currentRose.csv")

# Load project
proj = Project(file=inputFile)
model = proj.array
model.mooring_currentMod = 0     
model.ms.moorMod = 0
proj.trimGrids()

proj.plot3d(draw_boundary=True, boundary_on_bath=True, draw_bathymetry=True)
plt.show()
# Compute metrics
proj = mtr.metrics(proj, mtricFile, windFile=windFile, currentFile=currentFile)

# Plots
proj.plot3d(fowt=True, draw_boundary=False, boundary_on_bath=False, draw_bathymetry=False)
proj.plot2d(plot_boundary=False, plot_bathymetry=True)
plt.show()

# unload to yaml
proj.unload(file=inputFile)