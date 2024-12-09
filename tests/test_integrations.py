# tests Project class functionality
# (very early start)

import pytest

import numpy as np

from numpy.testing import assert_allclose

import matplotlib.pyplot as plt

from famodel.project import Project

from famodel.mooring.mooring import Mooring


def test_MoorPy_integration():
    project = Project(file='tests/testOntology.yaml',raft=0)
    project.getMoorPyArray(cables=1,plt=1)
    # check a random mooring line for ss
    assert project.mooringList['FOWT1a'].ss is not None
    
def test_RAFT_integration():
    project = Project(file='tests/testOntology.yaml')
    assert project.array is not None
    
'''def test_FLORIS_integration():'''

def test_lineDesign_integration():
    
    # make a dummy design dictionary for Mooring to make a Subsystem with
    dd = dict(sections={}, connectors={})
    dd['sections']   = [{} for i in range(1)]
    dd['connectors'] = [{} for i in range(2)]
    import moorpy as mp
    ms = mp.System(depth=200)
    # the sizing function coefficients to use in the design
    lineProps = [ms.setLineType(185,material='chain',name='chain')]
    lengths = [850]
    
    # Assign section properties for use in Mooring's Subsystem.makeGeneric call
    for i in range(1):
        dd['sections'][i]['type'] = lineProps[i]
        dd['sections'][i]['L'] = lengths[i]
    
    # # Assign props for intermediate points/connectors
    # for i in range(self.nLines-1):
    #     # if this is an intermediate line
    #     pointDict = getPointProps(Ws[ i + 1*(self.shared==1)]) 
        
    #     dd['connectors'][i+1]['m'] = pointDict['m']
    #     dd['connectors'][i+1]['v'] = pointDict['v']
    #     # CdA?
    
    # General mooring dimension info
    dd['span'    ] = 779.6
    dd['zAnchor' ] = -200
    dd['rad_fair'] = 58
    dd['z_fair'  ] = -14
    
    # create a mooring object
    moor = Mooring(dd=dd)
    
    # now make Subsystem, self.ss ??
    ss = moor.createSubsystem(case=0)
    assert ss is not None
    
    

if __name__=='__main__':
      
    
    ######## create Sample ontology and make raft, moorpy models
    project = Project(file='testOntology.yaml',raft=1)
    
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
    
    # test integration of floris model
    
    
    # test integration with failure modes and effects model
    
    # test integration with cable layout optimization model
    
    # test integration with lineDesign and cableDesign
    
    
    

    
    
        
        
        