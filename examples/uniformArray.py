# -*- coding: utf-8 -*-
"""
Example that shows how to create a uniform array 
from a yaml with NO Array table

Second portion of the code updates the uniform array with project.updateUniformArray() 
to add skew, different platform heading patterns, offsets, etc
"""
from famodel import Project
import os

# point to location of yaml file with uniform array info
input_directory = 'Inputs/' # relative location of directory for input files (yaml, bath files, etc)
filename = 'OntologySample200m_uniformArray.yaml' # yaml file to make initial platform(s)
# This yaml file does not contain explicit locations of each platform in the array table,
# but rather has a 'uniform_array' section that describes # of rows, cols, spacing, etc.
# This info is then used to automatically make a uniform array when the yaml file is loaded

# switch to directory of input files
os.chdir(input_directory)

# load in yaml
project = Project(file=filename)

# plot the system
project.getMoorPyArray()
project.plot3d()

#%% Update uniform array
# inputs for updateUniformArray
nrows = 5 # number of rows (total # of platforms should be same as input file)
ncols = 5 # number of columns
spacing_x = 1700 # spacing between columns
spacing_y = 1600 # spacing between rows
rotation_angle = 10 # angle to rotate the entire grid
skew_x = 5 # skew angle in x direction
skew_y = 5 # skew angle in y-direction
heading_pattern = [[0,180,0,180,0],[180,0,180,0,180]] # heading pattern. Length of outer list is # of rows it takes to repeat 
# the pattern, length of inner list is either 1 (for all same heading in the row) or the # of platforms in a row
offset_x = 600 # offset from boundary in x direction
offset_y = 600 # offset from boundary in y direction
center = [0,0] # center point of array

# call function to update array
project.updateUniformArray(nrows, ncols, [spacing_x,spacing_y], rotation_angle, skew_x,
                           skew_y, offset_x, offset_y, heading_pattern, center=center)
# plot again
project.plot3d()
