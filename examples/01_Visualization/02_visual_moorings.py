# -*- coding: utf-8 -*-
"""
Simple driver file to create a 2d plot of an platform locations with
mooring lines in an array.
The input file only contains the bare minimum information to build a 2d plot 
of the turbine locations and moorings (no cables, platform design, turbines, 
                                       site condition information, etc.)
"""

from famodel import Project
import matplotlib.pyplot as plt

# define name of ontology input file
input_file = '02_visual_moorings.yaml'

# initialize Project class with input file, we don't need RAFT for this so mark False
project = Project(file=input_file,raft=False)

# plot
project.plot2d()

# to moorings plot in 3d, we'll need to add depth and create a moorpy model of the system
project.depth = 200 # depth added because we did not include the site conditions section of the yaml
project.getMoorPyArray()
project.plot3d()

plt.show()

