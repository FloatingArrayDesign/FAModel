# tests anchor capacity and load functionality
from famodel.project import Project
import numpy as np


def test_capacities():
    # load in famodel project
    project = Project(file='testOntology.yaml',raft=False)
    for anch in project.anchorList.values():
        anch.getAnchorCapacity()
        