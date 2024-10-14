# -*- coding: utf-8 -*-
"""
This example shows how to duplicate a single platform object and relocate it.
Any cables will not be replicated. This function does not work with shared moorings or 
shared anchors
"""
from famodel.project import Project
import os

#### INPUTS ####
input_directory = 'Inputs/' # relative location of directory for input files (yaml, bath files, etc)
filename = 'OntologySample200m_noshared.yaml' # yaml file to make initial platform(s)
rep_pf_name = 'FOWT1' # platform to replicate (look at yaml file array data table to get platform names)
new_pf_loc = [0,0]

# switch to directory of input files
os.chdir(input_directory)

# first load in single platform from yaml
project = Project(file=filename)

# plot the system for comparison later
project.getMoorPyArray(plt=1)

# now choose a platform to replicate
rep_pf = project.platformList[rep_pf_name]

# call duplicate function
project.duplicate(rep_pf,r=new_pf_loc)

# make new moorpy array
project.getMoorPyArray()

# plot the new system
project.plot3d()