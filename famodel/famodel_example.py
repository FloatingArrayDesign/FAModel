# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 09:53:14 2024

@author: lsirkis
"""

###### FAModel sample farm driver file ######
from famodel.project import Project
import matplotlib.pyplot as plt
from copy import deepcopy
# Create the famodel farm and save it under the variable name project
print('Creating farm project')
project = Project(file='OntologySample.yaml',raft=0)
# create MoorPy array
print('Creating MoorPy array')
project.getMoorPyArray(pristineLines=1,plt=1) # use unmodified mooring lines

# plot the famodel farm
print('Plotting farm')
project.plot3d()

####### create lists of x-y coordinates of different objects ########
print('Obtaining coordinates of objects in farm')
# anchors
anchor_coords = []                              # initialize an empty list to store anchor locations in
for anch in project.anchorList.values():      # go through each anchor object in the farm
    anchor_coords.append(anch.r)                # add the x-y location of this anchor object to the anchor_xy list
    
# platforms
turb_coords = []                            # initialize an empty list to store platform locations in
for PF in project.platformList.values():    # go through each platform object in the farm
    turb_coords.append(deepcopy(PF.r))                # add the x-y location of this substation object to the anchor_xy list
    # add buffer zone around anchors etc
    PF.getbufferzones()
    
# substations
substation_coords = []
for sub in project.substationList.values():    # go through each substation object in the farm
    substation_coords.append(deepcopy(sub.r))                # add the x-y location of this substation object to the anchor_xy list
    
# get site boundary x-y coordinates
boundary_coords = project.boundary
    
    
    
# # ########## get watch circles and envelopes for platforms and mooring lines #############
# print('Creating watch circles and envelopes')
# # This step will take a long time...
# for moor in project.mooringList.values():
#     moor.getEnvelope()

project.plot2d()  # this should also plot the watch circles/envelopes!

plt.show()