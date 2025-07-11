"""
Simple driver file to create a 2d plot of locations of miscellaneous platform types in an array.
The input file only contains the bare minimum information to build a 2d plot 
of the platform locations (no moorings, cables, anchors, platform design, topsides, 
                          site condition information, etc.)
"""

from famodel import Project
import matplotlib.pyplot as plt

# define name of ontology input file
input_file = '06_2D-visual_misc_platforms.yaml'

# initialize Project class with input file, we don't need RAFT for this so mark False
project = Project(file=input_file,raft=False)

# plot
project.plot2d()

plt.show()