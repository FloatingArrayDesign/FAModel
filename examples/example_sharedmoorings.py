# -*- coding: utf-8 -*-
"""
This example shows shared mooring systems, including a shared mooring line and 
a shared anchor. The purpose of this example is to show the Ontology inputs for 
shared systems and how these systems might appear in FAModel
"""
import famodel
from famodel import Project
import matplotlib.pyplot as plt
import os
from moorpy.helpers import subsystem2Line
import moorpy as mp

# point to location of yaml file with uniform array info
dir = os.path.dirname(os.path.realpath(__file__))
filename = os.path.join(dir,'OntologySample600m_shared.yaml') # yaml file for project


# load in yaml
project = Project(file=filename, raft=True)


# plot in 2d and 3d
project.plot2d()
project.plot3d(plot_fowt=True)

# convert any subsystems to moorpy lines so we can unload to a moordyn file

i=0
for line in project.ms.lineList:
    if hasattr(project.ms.lineList[i],'lineList'):
        # this is a subsystem
        subsystem2Line(project.ms,i)
    else:
        i+=1

# unload to a moordyn file
project.ms.unload('testMD.dat')

# load in again to check it looks correct
ms = mp.System(file='testMD.dat')
ms.initialize()
ms.solveEquilibrium()
ms.plot()

plt.show()