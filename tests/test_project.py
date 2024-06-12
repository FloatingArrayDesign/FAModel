# tests Project class functionality
# (very early start)

import pytest

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
    
    # check number of mooring lines
    assert len(project.mooringList) == 8