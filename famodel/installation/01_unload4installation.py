"""
Script to unload a project class from a ConceptDesign class.

This just sets up a project, does not test the install model. 

Steps
-----
1. Load the ConceptDesign object.
2. Extract the project class for further processing.

Dependencies
-----------
fadesign.conceptual.conceptDesign.ConceptDesign
"""
import matplotlib.pyplot as plt
import os
from fadesign.conceptual.conceptDesign import ConceptDesign
from fadesign.conceptual import metrics as metr
import numpy as np
from famodel.project import Project
import pandas as pd

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']

# CONCEPT CATEGORY & NAME
# LIST OF CONCEPTS for installation (keep it updated)

# category | conceptName
# ILIA | 3_ILIA_1P2C_grd_mini (caseI for installation)

category    = "ILIA"
conceptName = "3_ILIA_al[cr]_grd"

# FILE LOCATIONS
filePath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "input_files")
inputFile = os.path.join(filePath, category, f"{conceptName}.yaml")

# Load the design
concept = ConceptDesign(baseDir=filePath, filename=inputFile, plot=False)
concept.design(plot=False)

# unload to yaml
outputFilename = f"proj_{conceptName}.yaml"
outputFile = os.path.join(filePath, category, outputFilename)
print(f"Unloading project to {outputFilename}")
concept.project.unload(file=outputFile)