"""
Simple driver file to create an FAModel project and make a moorpy system model
of the array, then run a watch circle analysis.
The input file only contains the bare minimum information to build a moorpy
array with moorings and dynamic cables (static cables are not modeled in MoorPy)
"""

from famodel import Project
import matplotlib.pyplot as plt

# define name of ontology input file
input_file = '04_array-watch-circle.yaml'

# initialize Project class with input file, we don't need RAFT for this so mark False
project = Project(file=input_file,raft=False)

# create moorpy array 
project.getMoorPyArray() # Platform hydrostatics information from ontology now used instead of default UMaine VolturnUS-S

# - - - Let's run a watch circle of the array - - - 
project.arrayWatchCircle() # run a watch circle using default thrust force (no wave/current/dynamic forces)

# - - - Now analyze the results
project.plot2d() # plot the motion envelopes

# Print some outputs for a random mooring line, anchor, and cable
print(f"Max load on a mooring line of platform fowt1: {project.mooringList['fowt1a'].loads}\n")
print(f"Safety factor for max load on mooring line fowt1a: {project.mooringList['fowt1a'].safety_factors}\n")
print(f"Max load on an anchor of platform fowt1: {project.anchorList['fowt1a'].loads}\n")
print(f"Max load on a dynamic cable attached to platform fowt2: {project.cableList['cable1'].subcomponents[-1].loads}\n")
print(f"Safety factors for dynamic cable attached to platform fowt2: {project.cableList['cable1'].subcomponents[-1].safety_factors}")

plt.show()

