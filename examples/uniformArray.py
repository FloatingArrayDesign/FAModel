# -*- coding: utf-8 -*-
"""
Example that shows how to create a uniform array 
from a yaml with NO Array table

Second portion of the code updates the uniform array with project.updateUniformArray() 
to add skew, different platform heading patterns, offsets, etc
"""
from famodel.project import Project
import os
import numpy as np
import matplotlib.pyplot as plt

# point to location of yaml file with uniform array info
dir = os.path.dirname(os.path.realpath(__file__))
filename = '\OntologySample200m_uniformArray.yaml' # yaml file to make initial platform(s)
# This yaml file does not contain explicit locations of each platform in the array table,
# but rather has a 'uniform_array' section that describes # of rows, cols, spacing, etc.
# This info is then used to automatically make a uniform array when the yaml file is loaded


# load in yaml
project = Project(file=dir+filename,raft=True)
project.plot2d()
# plot the system
project.plot3d(fowt=True)
project.trimGrids()
project.unload('checkyaml.yaml')

proj2 = Project(file='checkyaml.yaml')

#%% Update uniform array - edits in progress
print('\nUpdating uniform array shape\n')
# inputs for updateUniformArray
nrows = 5 # number of rows (total # of platforms should be same as input file)
ncols = 5 # number of columns
spacing_x = 1700 # spacing between columns
spacing_y = 1600 # spacing between rows
rotation_angle = 10 # angle to rotate the entire grid
skew_x = 2 # skew angle in x direction
skew_y = 2 # skew angle in y-direction
heading_pattern = [[0,180,0,180,0],[180,0,180,0,180]] # heading pattern. Length of outer list is # of rows it takes to repeat 
# the pattern, length of inner list is either 1 (for all same heading in the row) or the # of platforms in a row
offset_x = 600 # offset from boundary in x direction
offset_y = 1500 # offset from boundary in y direction
center = [0,0] # center point of array

# call function to update array
project.updateUniformArray(nrows, ncols, [spacing_x,spacing_y], rotation_angle, skew_x,
                            skew_y, offset_x, offset_y, heading_pattern, center=center)

model = project.array # this is the RAFT model!!
# update platform locations and headings in RAFT
for i,body in enumerate(model.fowtList):
    # set position
    pf_heading = project.platformList['fowt'+str(i)].phi
    body.setPosition([project.platformList['fowt'+str(i)].r[0],project.platformList['fowt'+str(i)].r[1],0,0,0,pf_heading])
    body.heading_adjust = np.degrees(project.platformList['fowt'+str(i)].phi)
# plot again
project.plot3d(fowt=True,draw_boundary=False,boundary_on_bath=False)

#%% Run FLORIS
print('\nRunning FLORIS\n')
config_file = 'gch.yaml'
turb_file = 'iea_15MW.yaml'
wr = 'maine_rose.csv'

project.getFLORISArray(config_file,[turb_file],[0,10.59,25],[0,1.95e6,1.9E6])
project.getFLORISMPequilibrium(10.59,0,.06,3,150,plotting=True)

#%% Watch circles for full array and plot envelopes
print('\nObtaining watch circle and motion envelopes\n')
project.arrayWatchCircle(ang_spacing=15)
project.plot2d()

plt.show()