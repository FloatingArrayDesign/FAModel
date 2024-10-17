# tests anchor capacity and load functionality
from famodel.project import Project
import numpy as np


def test_anchor_loads():
    # load in famodel project 
    project = Project(file='tests/testOntology.yaml', raft=False)
    project.getMoorPyArray()
    anch = project.anchorList['FOWT1a']
    
    # get lug loads on anchor
    anch.getLugForces(plot=False)
    
    assert('Ha' in anch.loads)
    assert('Hm' in anch.loads)
    assert(anch.loads['Ha'] != anch.loads['Hm'])

def test_anchor_capacities():
    # load in famodel project (suction pile anchor)
    project = Project(file='tests/testOntology.yaml', raft=False)
    project.getMoorPyArray(cables=1)
    anch = project.anchorList['FOWT1a']
    
    # fill in load dictionary to skip watch circle run
    loads = {'Ha':4522222,'Va':3948278}
    # get capacity and safety factor
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    # try suction pile with sand
    newdd = anch.dd
    newdd['soil_type'] = 'sand'
    newdd['soil_properties']['phi'] =33
    newdd['soil_properties']['beta'] = 5
    # get capacity and safety factor
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    # check plate anchor type   
    newdd['type'] = 'plate'
    newdd['design'] = {'type':'plate','A':20,'zlug':10,'beta':10}
    newdd['soil_type'] = 'clay'
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
    newdd['soil_type'] = 'rock' # soil_properties has default rock info in there already, just change name
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    # check driven pile anchor in rock and clay
    newdd['type'] = 'driven'
    newdd['soil_type'] = 'weak_rock'
    newdd['design'] = {'type':'driven','L':20,'D':1.5,'zlug':-3}
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    newdd['soil_type'] = 'clay'
    newdd['design']['zlug'] = 3
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    
    # check torpedo anchor
    newdd['type'] = 'torpedo_pile'
    newdd['design'] = {'type':'torpedo_pile','D1':3,'D2':1.1,'L1':10,'L2':4,'zlug':16}
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    # check helical pile anchor
    newdd['type'] = 'helical_pile'
    newdd['design'] = {'type':'helical_pile','L':25.1,'d':1,'D':5.01}
    anch.getAnchorCapacity(loads=loads, plot=False)
    anch.getFS(loads=loads)
    
    

    
        