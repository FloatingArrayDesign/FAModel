
"""
Test mooring loading, configurations, methods
"""
import pytest

import numpy as np

from famodel.project import Project

from famodel.platform.fairlead import Fairlead
from famodel.platform.platform import Platform

import os


@pytest.fixture
def setup_project():
    dir = os.path.dirname(os.path.realpath(__file__))
    return(Project(file=os.path.join(dir,'mooring_ontology.yaml'), raft=False))
    
    
def test_num_moorings(setup_project):
    
    assert(len(setup_project.mooringList)==11)
    
def test_moor_heading(setup_project):
    
    moor = setup_project.mooringList['FOWT1a']
    dists = moor.rA[:2]-moor.rB[:2]
    heading = np.pi/2 - np.arctan2(dists[1], dists[0])
    pf = setup_project.platformList['FOWT1']
    assert(heading == np.radians(45+180))
    assert(heading == pf.mooring_headings[0]+pf.phi)
    assert(heading == np.radians(moor.heading))
    
def test_platform_connection(setup_project):
    
    moor = setup_project.mooringList['FOWT1a']
    assert(moor.attached_to[0]==setup_project.anchorList['FOWT1a'])
    assert(moor.attached_to[1]==setup_project.platformList['FOWT1'])
    
def test_fairlead_connection(setup_project):
    
    end_sub = setup_project.mooringList['FOWT1a'].subcomponents[-1]
    assert(len(end_sub.attachments)==2)
    assert(np.any([isinstance(att['obj'], Fairlead) for att in end_sub.attachments.values()]))
    
def test_rA_depth(setup_project):
    moor = setup_project.mooringList['FOWT1a']
    loc = moor.rA
    true_depth = setup_project.getDepthAtLocation(loc[0],loc[1])
    assert(moor.rA[2] == -true_depth)
    assert(moor.dd['zAnchor'] == -true_depth)
    assert(moor.z_anch == -true_depth)
    setup_project.getMoorPyArray()
    assert(moor.ss.rA[2] == -true_depth)
    
'''
    
def test_end_locs(self):
    moor = self.project.mooringList['fowt1a']
    assert(moor.rB == )
    assert(moor.rA == )
'''
    
def test_num_sections(setup_project):
    moor = setup_project.mooringList['FOWT1a']
    setup_project.getMoorPyArray()
    assert(len(moor.dd['sections'])==len(moor.ss.lineList))
    assert(len(moor.dd['sections'])==2)
    
def test_num_connectors(setup_project):
    moor = setup_project.mooringList['FOWT1a']
    assert(len(moor.dd['connectors'])==3)
    
def test_shared_connections(setup_project):
    
    moor = setup_project.mooringList['FOWT1-FOWT2']
    assert(len(moor.subcomponents[0].attachments)==2)
    assert(np.any([isinstance(att['obj'], Fairlead) for att in moor.subcomponents[0].attachments.values()]))
    assert(isinstance(moor.attached_to[0], Platform))
    assert(isinstance(moor.attached_to[1], Platform))
    
def test_shared_flag(setup_project):
    
    moor = setup_project.mooringList['FOWT1-FOWT2']
    assert(moor.shared == 1)
 
'''
def test_shared_depth(self):
    
    moor = self.project.mooringList['fowt1b']
    self.project.getMoorPyArray()
    assert(moor.ss)
'''
    
