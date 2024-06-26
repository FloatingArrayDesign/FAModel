# tests Project class functionality
# (very early start)

import pytest

import numpy as np

from numpy.testing import assert_allclose

import matplotlib.pyplot as plt

from famodel import Project

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

if __name__=='__main__':
    
    # initialize a blank project
    project = Project()
    
    # load bathymetry file (MoorDyn style)
    project.loadBathymetry('bathymetry_sample.txt')
    
    # sample the depth somewhere
    x=20
    y=50
    z = project.getDepthAtLocation(x, y)
    print(f"The depth at {x}, {y} is {z:.3f} m.")
    
    # plot to see if it worked
    #project.plot3d()
    plt.show()
    
    
    ######## test load design features ################
    project.load('testOntology.yaml')
    
    from famodel.mooring.mooring import Mooring
    from famodel.platform.platform import Platform
    from famodel.mooring.anchor import Anchor
    from famodel.mooring.connector import Connector
    from famodel.substation.substation import Substation
    from famodel.cables.cable import Cable
    from famodel.cables.dynamic_cable import DynamicCable
    from famodel.cables.static_cable import StaticCable
    from famodel.cables.cable_properties import getCableProps, getBuoyProps
    from famodel.cables.components import Joint
    from famodel.turbine.turbine import Turbine
    
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
    
    
    
    
    # check connections
    for i,anch in enumerate(project.anchorList.values()):
        
        # check last anchor is a shared anchor:
        if i == 5:
            assert len(anch.attachments) == 2
            
        # check anchor is attached to a mooring line and mooring line is attached to the anchor
        for j in range(0,len(anch.attachments)):
            # get mooring line attached to anchor
            moor = anch.attachments[j]['obj']
            # check that it is a mooring line
            assert isinstance(moor,Mooring)
            # check that mooring line is attached to anchor
            assert moor.attached_to[0] == anch
            
    for i, pf in enumerate(project.platformList.values()):
        
        # check number of things connected to platform
        if i == 0:
            assert len(pf.attachments) == 5 # 3 lines, 1 turbine, 1 cable 
        else:
            assert len(pf.attachments) == 6 # 3 lines, 1 turbine, 2 cables
        
    # check angles and repositioning for regular mooring line and shared mooring line
    assert_allclose(project.mooringList[0].rA,[-828.637,-828.637,-600],rtol=0,atol=0.1)
    assert_allclose(project.mooringList[-1].rA,)

    
    
        
        
        