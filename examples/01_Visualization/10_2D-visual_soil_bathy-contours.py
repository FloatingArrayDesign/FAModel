# -*- coding: utf-8 -*-
"""
Simple driver file to create a 2d plot of soil types.
The input file only contains the bare minimum information to build a 2d plot 
of the soil and bathymetry contour lines (no platforms, moorings, cables, platform design, turbines, 
                                       site condition information, etc.)
"""

from famodel import Project
import matplotlib.pyplot as plt
import os

# define name of ontology input file
dir = os.path.dirname(os.path.realpath(__file__))
input_file = os.path.join(dir,'10_2D-visual_soil_bathy-contours.yaml')

# initialize Project class with input file, we don't need RAFT for this so mark False
project = Project(file=input_file,raft=False)

# plot the soil types with the bathymetry contour lines (must turn off bathymetry colorplot)
project.plot2d(plot_soil=True, plot_bathymetry=False, 
               plot_bathymetry_contours=True,
               bathymetry_levels=10) # choose # of contour lines to show

# Let's change the number of contour levels
project.plot2d(plot_soil=True, plot_bathymetry=False, 
               plot_bathymetry_contours=True,
               bathymetry_levels=5) # choose # of contour lines to show


plt.show()
