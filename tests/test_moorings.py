
"""
Test mooring loading, configurations, methods
"""
import pytest

import numpy as np

from famodel.project import Project

from famodel.platform.fairlead import Fairlead
from famodel.platform.platform import Platform
from famodel.mooring.connector import Section, Connector
from famodel.famodel_base import Node, Edge

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
    
def test_fairlead_position(setup_project):
    moor = setup_project.mooringList['FOWT1a']
    fl = moor.subcomponents[-1].attachments['FOWT1_F1']
    assert(fl.r==setup_project.mooringList['FOWT1a'].rB)
    pf = setup_project.platformList['FOWT1']
    head_fl = np.radians(90-30)
    head_pf = np.radians(90)-pf.phi
    assert(fl.r==[58*cos(head_fl+head_pf),58*sin(head_fl+head_pf),-14])
    
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
    
# - - - -tests in progress- - - -
    

def bridle_project():
    dir = os.path.dirname(os.path.realpath(__file__))
    return(Project(file=os.path.join(dir,'mooring_ontology_parallels.yaml'), raft=False))

def test_bridle_setup(bridle_project):
    moor = bridle_project.mooringList['FOWT2a']
    # check subcons_B is a list of length 2
    assert(len(moor.subcons_B)==2)
    # check each item in subcons_B is attached to 2 things (fairlead and another subcomponent)
    for sub in moor.subcons_B:
        assert(len(sub.attachments)==2)
        for att in sub.attachments.values():
            assert(isinstance(att['obj'],[Fairlead,Section]))
    pf = moor.attached_to[1]
    fl_attachment = [False, False]
    for i,sub in enumerate(moor.subcons_B):
        for att in pf.attachments.values():
            if sub.id in att['obj'].attachments:
                fl_attachment[i] = True
            
    assert(all(fl_attachment))
    
def test_bridle_end_locs(bridle_project):
    moor = bridle_project.mooringList['FOWT1a']
    # check rB is at midpoint of fairlead locs
    fl_locs = []
    for sub in moor.subcons_B:
        att = [att['obj'] for att in sub.attachments.values() if isinstance(att['obj'],Fairlead)]
        fl_locs.append(att[0].r)
    from famodel.helpers import calculate_midpoint
    midpoint = calculate_midpoint(fl_locs)
    assert(midpoint==moor.rB)
    # check 
    # check location of anchor is correct
    u = np.array([np.cos(np.radians(moor.heading)),np.sin(np.radians(moor.heading))])
    anch_loc = np.hstack((np.array(midpoint[:2])+moor.span*u,-bridle_project.depth))
    assert(anch_loc==moor.rA)
    
    

 
'''
def test_shared_depth(self):
    
    moor = self.project.mooringList['fowt1b']
    self.project.getMoorPyArray()
    assert(moor.ss)
'''
    
