# tests anchor capacity and load functionality
from famodel.project import Project
import numpy as np
import os


def test_anchor_loads():
    # load in famodel project 
    dir = os.path.dirname(os.path.realpath(__file__))
    project = Project(file=os.path.join(dir,'testOntology.yaml'), raft=False)
    project.getMoorPyArray(cables=1)
    anch = project.anchorList['FOWT1a']
    
    # get lug loads on anchor
    anch.getLugForces(plot=False)
    
    assert('Ha' in anch.loads)
    assert('Hm' in anch.loads)
    assert(anch.loads['Ha'] != anch.loads['Hm'])

def test_anchor_capacities():
    # load in famodel project (suction pile anchor)
    dir = os.path.dirname(os.path.realpath(__file__))
    project = Project(file=os.path.join(dir,'testOntology.yaml'), raft=False)
    project.getMoorPyArray(cables=1)
    anch = project.anchorList['FOWT1a']
    
    # fill in load dictionary to skip watch circle run
    loads = {'Ha':4522222,'Va':3948278}
    # get capacity and safety factor
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    # try suction pile with sand
    soil = anch.soilProps
    soil['sand'] = soil.pop(next(iter(soil.keys())))
    soil['sand']['phi'] = 33
    soil['sand']['Dr'] = 50
    # get capacity and safety factor
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    # check plate anchor type   
    newdd = anch.dd
    newdd['type'] = 'plate'
    newdd['design'] = {'type':'plate','A':20,'zlug':10,'beta':10}
    anch.soilProps['clay'] = anch.soilProps.pop('sand')
    # new loads
    loads['Ha'] = 1000000
    loads['Va'] = 0
    # get capacity and safety factor
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    # check drilled and grouted anchor type (need to change material to rock)
    loads = {'Ha':4522222,'Va':3948278} # go back to original loads    
    newdd['type'] = 'dandg_pile'
    newdd['design'] = {'type':'dandg_pile','L':50,'D':3,'zlug':0}
    soil['rock'] = soil.pop('clay') # soil_properties has default rock info in there already, just change name
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    # check driven pile anchor in rock and clay
    newdd['type'] = 'driven'
    soil['weak_rock'] = soil.pop('rock')
    newdd['design'] = {'type':'driven','L':20,'D':1.5,'zlug':-3}
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    soil['clay'] = soil.pop('weak_rock')
    newdd['design'] = {'type':'driven','L':30,'D':2,'zlug':3}
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    soil['sand'] = soil.pop('clay')
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    soil['sand']['Dr'] = 50
    
    # check helical pile anchor with sand
    newdd['type'] = 'helical_pile'
    newdd['design'] = {'type':'helical_pile','L':25.1,'d':1,'D':5.01,'zlug':5}
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    # check helical pile anchor with clay
    soil['clay'] = soil.pop('sand')
    newdd['type'] = 'helical_pile'
    newdd['design'] = {'type':'helical_pile','L':25.1,'d':1,'D':5.01,'zlug':5}
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    # check torpedo anchor
    newdd['type'] = 'torpedo_pile'
    newdd['design'] = {'type':'torpedo_pile','D1':3,'D2':1.1,'L1':10,'L2':4,'zlug':16}
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    
    
    

    
        