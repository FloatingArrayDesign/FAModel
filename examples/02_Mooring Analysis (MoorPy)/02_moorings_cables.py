# -*- coding: utf-8 -*-
"""
Simple driver file to create an FAModel project and make a moorpy system model
of the array including moorings and dynamic cables.
The input file only contains the bare minimum information to build a moorpy
array with moorings and dynamic cables (static cables are not modeled in MoorPy)
"""

from famodel import Project
import matplotlib.pyplot as plt
import numpy as np

# define name of ontology input file
input_file = '02_moorings_cables.yaml'

# initialize Project class with input file, we don't need RAFT for this so mark False
project = Project(file=input_file,raft=False)

# create moorpy array
project.getMoorPyArray()

# - - - Let's do a quick simulation of force on the moorpy array - - - 
ms = project.ms                                             # create a variable shortcut to moorpy system model
print(f"Body initial position is {ms.bodyList[0].r6}")      # print initial position of a platform 
fig, ax = ms.plot()                                         # plot the system in original configuration
ms.bodyList[0].f6Ext = np.array([3e6, 0, 0, 0, 0, 0])       # apply an external force on the body [N]
ms.solveEquilibrium3(DOFtype='both')                        # equilibrate
fig, ax = ms.plot(ax=ax, color='red')                       # plot the system in displaced configuration (on the same plot, in red)

print(f"Body offset position is {ms.bodyList[0].r6}")       # print offset position of a platform for comparison


plt.show()

