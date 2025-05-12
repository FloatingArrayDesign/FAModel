"""
This example unloads a project class from a conceptDesign class. 
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
conceptName = "3_ILIA_1P2C_grd_mini"

# FILE LOCATIONS
filePath = os.path.dirname(os.path.abspath(__file__))
inputFile = os.path.join(filePath, category, f"{conceptName}.yaml")

concept = ConceptDesign(baseDir=filePath, filename=inputFile, plot=False)
concept.design(plot=False)

# unload to yaml
outputFilename = f"proj_{conceptName}.yaml"
outputFile = os.path.join(filePath, category, outputFilename)
print(f"Unloading project to {outputFilename}")
concept.project.unload(file=outputFile)