# -*- coding: utf-8 -*-
"""
Simple driver file to create an array with moorings and rotate the platforms and array.

This allows you to rotate platforms including all of their moorings, anchors, and fairleads
The input file only contains the bare minimum information to build a 2d plot 
of the turbine locations and moorings with fairleads (no cables, platform design, turbines, 
                                       site condition information, etc.)
"""

from famodel import Project
import matplotlib.pyplot as plt

# define name of ontology input file
input_file = '01_Fairleads.yaml'

# initialize Project class with input file, we don't need RAFT for this so mark False
project = Project(file=input_file,raft=False)

# plot
project.plot2d()

# let's rotate a platform but keep in the same x,y position
pf_loc = project.platformList['fowt0'].r
new_heading = 143 # [deg]
project.platformList['fowt0'].setPosition(r=pf_loc, heading=new_heading,
                                          degrees=True, project=project)

# plot again to see the difference
project.plot2d()

# let's now change the platform's position
new_r = [-2000, -2200]
project.platformList['fowt0'].setPosition(r=new_r, project=project)

# plot again
project.plot2d()


