# -*- coding: utf-8 -*-
"""
Create moorpy system of one platform in the array and save as a MoorDyn file
"""

from famodel import Project
import os
import moorpy as mp
from moorpy.helpers import subsystem2Line

# get locations of files
cwd = os.getcwd()
input_dir = cwd+'/Inputs/'
yaml_file = 'OntologySample200m_uniformArray.yaml'

os.chdir(input_dir)
# create project
project = Project(file=yaml_file)

# create a moorpy system for a single platform
project.platformList['fowt0'].mooringSystem(project=project)

# save location
ms = project.platformList['fowt0'].ms

# convert subsystems in ms to lines
for line in range(len(ms.lineList)):
    subsystem2Line(ms,0)

# unload
ms.unload('output_MD.dat')
