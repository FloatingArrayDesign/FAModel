"""
Simple driver file to create a RAFT model of a 2-platform array without turbines.
This file then adds a wave spectrum to run, analyzes that case, and plots the results.

For more information on using RAFT, please see RAFT documentation at https://github.com/WISDEM/RAFT
"""

from famodel import Project
import matplotlib.pyplot as plt
import os

# define name of ontology input file
dir = os.path.dirname(os.path.realpath(__file__))
input_file = os.path.join(dir,'01_platform.yaml')

# initialize Project class with input file, we don't need RAFT for this so mark False
project = Project(file=input_file,raft=True)

# pull out RAFT object, this will only have a raft platform (no turbine), because we didn't specify a RAFT turbine definition in the ontology file
raft_model = project.array # store short cut to raft model 

# - - - Let's try running a case - - - 
## First, let's add a case in to the model's design dictionary (since we didn't add this in the ontology)
raft_model.design['cases'] = {}
raft_model.design['cases']['keys'] = ['wind_speed', 'wind_heading', 'turbulence', 'turbine_status', 'yaw_misalign', 'wave_spectrum', 'wave_period', 'wave_height', 'wave_heading']
raft_model.design['cases']['data'] = [[     0,         0,             0,             'off',        0,             'JONSWAP',         12,           6,              0        ]]

# analyze our case
raft_model.analyzeCases(display=True) # display what's happening for fun

# plot RAFT results
raft_model.plotResponses()
raft_model.plot()

plt.show()