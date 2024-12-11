# -*- coding: utf-8 -*-
"""
Create a platform object and its associated moorings and anchors 
from a moorpy system. Relocate system to desired location

moorpy system must have only one platform, and no shared anchors or moorings allowed.

If using a moordyn file, it must have a body, or you will need to edit the ms to
include a body after loading the moordyn file.
"""

from famodel.project import Project
import moorpy as mp
import os
import matplotlib.pyplot as plt

#### INPUTS ####
input_directory = 'Inputs/' # relative location of directory for input files (yaml, bath files, etc)
filename = 'MoorDyn_semitaut200m.dat' # moordyn file to create a moorpy system
rep_pf_name = 'FOWT1' # platform to replicate (look at yaml file array data table to get platform names)
new_pf_loc = [100,100]


# change to input directory
os.chdir(input_directory)

# create moorpy system
ms = mp.System(file=filename)
ms.initialize()
ms.solveEquilibrium()
ms.plot()

# create empty project class
project = Project(depth=ms.depth, raft=0)

# add bathymetry if you want 
project.loadBathymetry('bathymetry200m_sample.txt')

# add soil info as needed
project.loadSoil('soil_sample.txt')

# add platform, mooring, and anchor objects from ms
project.addPlatformMS(ms, r=new_pf_loc)

# get new moorpy array and plot
project.getMoorPyArray(plt=1)
plt.show()
