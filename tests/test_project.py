# tests Project class functionality
# (very early start)

import pytest

import numpy as np

from numpy.testing import assert_allclose

import matplotlib.pyplot as plt

from famodel import Project

from famodel.mooring.mooring import Mooring
from famodel.platform.platform import Platform
from famodel.anchors.anchor import Anchor
from famodel.mooring.connector import Connector
from famodel.substation.substation import Substation
from famodel.cables.cable import Cable
from famodel.cables.dynamic_cable import DynamicCable
from famodel.cables.static_cable import StaticCable
from famodel.cables.cable_properties import getCableProps, getBuoyProps
from famodel.cables.components import Joint
from famodel.turbine.turbine import Turbine

"""

def test_tensions_swap():
    '''Compares two equivalent catenary mooring lines that are defined in opposite directions.'''
       
    ms = mp.System(depth=60)

    #ms.lineTypes['chain'] = getLineProps(120, name='chain')    # add a line type
    ms.setLineType(120, 'chain', name='chain')    # add a line type
    
    ms.addPoint(1, [  0, 0, -60])
    ms.addPoint(1, [100, 10, -30])

    # line sloping up from A to B, and another in the opposite order
    ms.addLine(120, 'chain', pointA=1, pointB=2)
    ms.addLine(120, 'chain', pointA=2, pointB=1)

    ms.initialize()
    
    # compare tensions
    assert_allclose(np.hstack([ms.lineList[0].fA, ms.lineList[0].fB]),
                    np.hstack([ms.lineList[1].fB, ms.lineList[1].fA]), rtol=0, atol=10.0, verbose=True)
"""

def test_bathymetry():
    # initialize a blank project
    project = Project()
    
    # load bathymetry file (MoorDyn style)
    project.loadBathymetry('tests/bathymetry_sample.txt')
    
    # sample the depth somewhere
    x=20
    y=50
    z = project.getDepthAtLocation(x, y)
    print(f"The depth at {x}, {y} is {z:.3f} m.")
    
    # plot to see if it worked
    #project.plot3d()
    plt.show()

def test_create_components():
    project = Project(file='tests/testOntology.yaml',raft=0)
       
    # check number of mooring lines
    assert len(project.mooringList) == 8
    
    # check number of anchors
    assert len(project.anchorList) == 6
    
    # check number of cables
    assert len(project.cableList) == 3
    
    # check number of turbines
    assert len(project.turbineList) == 3
    
    # check number of substations
    assert len(project.substationList) == 1
    
    # check number of platforms
    assert len(project.platformList) == 3
    
def test_check_connections():
    project = Project(file='tests/testOntology.yaml',raft=0)
    # check connections
    for i,anch in enumerate(project.anchorList.values()):
        
        # check last anchor is a shared anchor:
        if i == 5:
            assert len(anch.attachments) == 2
            
        # check anchor is attached to a mooring line and mooring line is attached to the anchor
        for j,att in enumerate(anch.attachments):
            # get mooring line attached to anchor
            moor = anch.attachments[att]['obj']
            # check that it is a mooring line
            assert isinstance(moor,Mooring)
            # check that mooring line is attached to anchor
            assert moor.attached_to[0] == anch
            
    for i, pf in enumerate(project.platformList.values()):
        
        # check number of things connected to platform
        if i == 1:
            assert len(pf.attachments) == 5 # 3 lines, 1 turbine, 1 cable 
        else:
            assert len(pf.attachments) == 6 # 3 lines, 1 turbine, 2 cables
            
def test_headings_repositioning():
    project = Project(file='tests/testOntology.yaml',raft=0)
    # check angles and repositioning for regular mooring line and shared mooring line, reg cable and suspended cable
    assert_allclose(np.hstack((project.mooringList['FOWT1a'].rA,project.mooringList['FOWT1-FOWT2'].rA)),
                    np.hstack(([-828.637,-828.637,-600],[40.5,0,-20])),rtol=0,atol=0.5)
    assert_allclose(np.hstack((project.cableList['array_cable10'].subcomponents[0].rB,project.cableList['suspended_cable11'].subcomponents[0].rB)),
                    np.hstack(([640.5,0,-600],[0,1615.5,-20])),rtol=0,atol=0.5)
    
def test_marine_growth():
    project = Project(file='tests/testOntology.yaml',raft=0)
    # check correct mg gets added to specified mooring lines and cables for ss_mod
    project.getMarineGrowth(lines=['FOWT1a',['suspended_cable11',0]])
    # pull out a mooring line and a cable to check
    Moor = project.mooringList['FOWT1a'].ss.lineList[0].type['d_vol']
    mgMoor = project.mooringList['FOWT1a'].ss_mod.lineList[1].type['d_vol']
    
    Cab = project.cableList['suspended_cable11'].subcomponents[0].ss.lineList[0].type['d_vol']
    mgCab = project.cableList['suspended_cable11'].subcomponents[0].ss_mod.lineList[0].type['d_vol']
    
    assert_allclose(np.hstack((mgMoor,mgCab)),np.hstack((Moor+0.1156,Cab+0.4)),rtol=0,atol=0.05)
    
def test_seabed():
    '''test seabed properties are properly loaded and included in anchor design dictionaries'''
    # check soil at a location
    project = Project(file='tests/testOntology.yaml',raft=0)
    soilInfo = project.getSoilAtLocation(-828.637,-828.637)
    assert soilInfo[0] == 'rock'
    assert soilInfo[1]['UCS'] == 7
    assert soilInfo[1]['Em'] == 50
    anch = project.anchorList['FOWT1a']
    assert anch.dd['soil_properties']['UCS'] == 7
    assert anch.dd['soil_type'] == 'rock'
    
    
    

    

if __name__=='__main__':
    '''
    '''
