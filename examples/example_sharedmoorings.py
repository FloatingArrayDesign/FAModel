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

# point to location of yaml file with uniform array info
dir = os.path.dirname(os.path.realpath(__file__))
filename = '\OntologySample600m_shared.yaml' # yaml file for project


# load in yaml
project = Project(file=dir+filename, raft=True)


# plot in 2d and 3d
project.plot2d()
project.plot3d(fowt=True)

plt.show()