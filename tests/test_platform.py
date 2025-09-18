# -*- coding: utf-8 -*-
"""
Created on Mon Jul 14 13:33:12 2025

@author: lsirkis
"""
import pytest
import os
from famodel import Project
import numpy as np
from copy import deepcopy
from numpy.testing import assert_allclose

@pytest.fixture
def project():
    dir = os.path.dirname(os.path.realpath(__file__))
    return(Project(file=os.path.join(dir,'platform_ontology.yaml'), raft=False))

def test_platform_location(project):
    assert_allclose(project.platformList['FOWT2'].r, [1600, 0, 0])
    
def test_basic_pf_relocation(project):
    project.platformList['FOWT1'].setPosition(r=[20, 20, -10], heading=30, degrees=True)
    assert_allclose(project.platformList['FOWT1'].r,[20, 20, -10])
    assert pytest.approx(project.platformList['FOWT1'].phi) == np.radians(30)
    
def test_fl_relloc(project):
    fl = project.platformList['FOWT2'].attachments['FOWT2_F2']
    fl_loc = [np.cos(np.radians(60))*40.5,
              np.sin(np.radians(60))*40.5,
              -20]
    assert_allclose(fl['r_rel'],fl_loc)

def test_pf_relocation(project):
    new_r = [1500, 1500]
    new_head = 15
    moor = project.mooringList['FOWT2a']
    moor_head_start = deepcopy(project.mooringList['FOWT2a'].heading)
    fl = project.platformList['FOWT2'].attachments['FOWT2_F2']
    project.platformList['FOWT2'].setPosition(r = new_r, 
                                              heading=new_head, 
                                              degrees=True, 
                                              project=project)
    assert pytest.approx(moor.heading) == new_head+moor_head_start
    fl_loc_new = [1500 + np.cos(np.radians(60-new_head))*40.5,
                  1500 + np.sin(np.radians(60-new_head))*40.5,
                  -20]
    fl = project.platformList['FOWT2'].attachments['FOWT2_F2']
    assert_allclose(fl['obj'].r, fl_loc_new)
    moor_head_new = np.radians(90-(new_head+moor_head_start))
    new_x = fl['obj'].r[0] + np.cos(moor_head_new)*moor.span
    new_y = fl['obj'].r[1] + np.sin(moor_head_new)*moor.span
    new_anch_r = [new_x,
                  new_y,
                  -project.getDepthAtLocation(new_x,new_y)]
    project.plot2d()
    assert_allclose(project.anchorList['FOWT2a'].r, new_anch_r)
    

    

