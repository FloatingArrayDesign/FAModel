# tests anchor capacity and load functionality
from famodel.project import Project
import numpy as np
import os
import matplotlib.pyplot as plt
import pytest

# --- Helper goes at module level ---
def assign_soil(anchor, soil_label, project):
    soil_def = project.soilProps[soil_label]
    layers = soil_def['layers']
    print('[DEBUG] assign_soil: soil_label =', soil_label)
    print('[DEBUG] assign_soil: soil_def =', soil_def)
    profile_map = [{
        'name': 'CPT_Assigned',
        'x': 0, 'y': 0,
        'layers': layers}]
    
    anchor.setSoilProfile(profile_map)
    anchor.profile_name = 'CPT_Assigned'

@pytest.fixture
def project():
    dir = os.path.dirname(os.path.realpath(__file__))
    return(Project(file=os.path.join(dir,'testOntology.yaml'), raft=False))

def test_anchor_loads(project):
    # load in famodel project 
    project.getMoorPyArray(cables=1)
    anch = project.anchorList['FOWT1a']
   
    assign_soil(anch, 'mud_soft', project)
    
    # Force calculation
    anch.getMudlineForces(max_force=False)
    
    # Extract mudline loads
    Hm = anch.loads.get('Hm')
    Vm = anch.loads.get('Vm')
    zlug = anch.dd['design']['zlug']

    # Compute lug loads
    _, Ha, Va = anch.getLugForces(Hm, Vm, zlug, plot=False)
    anch.loads['Ha'] = Ha
    anch.loads['Va'] = Va

    # Assertions
    assert 'Ha' in anch.loads
    assert 'Hm' in anch.loads
    assert anch.loads['Ha'] != anch.loads['Hm']
    
def test_anchor_capacities(project):
    # load in famodel project (suction pile anchor)
    project.getMoorPyArray(cables=1)
    anch = project.anchorList['FOWT1a']

    assign_soil(anch, 'mud_firm', project)
    
    # fill in load dictionary to skip watch circle run
    loads = {'Ha':4.5e6, 'Va':1.9e6}
    # get capacity and safety factor
    Ha = loads['Ha']
    Va = loads['Va']
    zlug = anch.dd['design']['zlug']
    anch.getCapacityAnchor(Ha, Va, zlug=zlug, plot=False)
    anch.getSafetyFactor()

    # try SUCTION PILE with sand
    assign_soil(anch, 'sand', project)
    # get capacity and safety factor
    Ha = loads['Ha']
    Va = loads['Va']
    zlug = anch.dd['design']['zlug']
    anch.getCapacityAnchor(Ha, Va, zlug=zlug, plot=False)
    anch.getSafetyFactor()
    
    # check PLATE ANCHOR type   
    anch.dd['type'] = 'plate'
    anch.dd['design'] = {'B':5, 'L':2, 'zlug':10, 'beta':10}
    assign_soil(anch, 'mud_soft', project)
    # new horizontal load
    loads['Ha'] = 2e6
    loads['Va'] = 0
    Ha = loads['Ha']
    Va = loads['Va']
    zlug = anch.dd['design']['zlug']
    # get capacity and safety factor
    anch.getCapacityAnchor(Ha, Va, zlug=zlug, plot=False)
    anch.getSafetyFactor()
       
    # check DRILLED & GROUTED PILE (need to change material to rock)
    loads = {'Ha':4.5e5, 'Va':1.9e5} # again assign new loads    
    anch.dd['type'] = 'dandg'
    anch.dd['design'] = {'L':10, 'D':3, 'zlug':0}
    assign_soil(anch, 'weak_rock', project)
    Ha = loads['Ha']
    Va = loads['Va']
    zlug = anch.dd['design']['zlug']
    anch.getCapacityAnchor(Ha, Va, zlug=zlug, plot=False)
    anch.getSafetyFactor()
    
    # check DRIVEN PILES in soils and rock
    anch.dd['type'] = 'driven'
    anch.dd['design'] = {'L':20, 'D':1.5, 'zlug': 3}
    assign_soil(anch, 'mud_firm', project)
    Ha = loads['Ha']
    Va = loads['Va']
    zlug = anch.dd['design']['zlug']
    anch.getCapacityAnchor(Ha, Va, zlug=zlug, plot=False)
    anch.getSafetyFactor()
    
    anch.dd['design'] = {'L':30, 'D':2, 'zlug':3}
    assign_soil(anch, 'sand', project)
    Ha = loads['Ha']
    Va = loads['Va']
    zlug = anch.dd['design']['zlug']
    anch.getCapacityAnchor(Ha, Va, zlug=zlug, plot=False)
    anch.getSafetyFactor()
    
    assign_soil(anch, 'weak_rock', project)
    Ha = loads['Ha']
    Va = loads['Va']
    # change the padeye back to mudline elevation in rock
    anch.dd['design']['zlug'] = 0
    zlug = anch.dd['design']['zlug'] 
    anch.getCapacityAnchor(Ha, Va, zlug=zlug, plot=False)
    anch.getSafetyFactor()
       
    # check HELICAL PILE with sand
    anch.dd['type'] = 'helical'
    anch.dd['design'] = {'L':15, 'd':1.25, 'D':2.00, 'zlug':3}
    assign_soil(anch, 'sand', project)
    Ha = loads['Ha']
    Va = loads['Va']
    zlug = anch.dd['design']['zlug']
    anch.getCapacityAnchor(Ha, Va, zlug=zlug, plot=False)
    anch.getSafetyFactor()
    
    # check HELICAL PILE with clay
    anch.dd['type'] = 'helical'
    anch.dd['design'] = {'L':12, 'd':0.5, 'D':1.5, 'zlug':3}
    assign_soil(anch, 'mud_firm', project)
    Ha = loads['Ha']
    Va = loads['Va']
    zlug = anch.dd['design']['zlug']
    anch.getCapacityAnchor(Ha, Va, zlug=zlug, plot=False)
    anch.getSafetyFactor()
    
    # check TORPEDO PILE
    anch.dd['type'] = 'torpedo'
    anch.dd['design'] = {'D1':3, 'D2':1.1, 'L1':10, 'L2':4, 'zlug':16}
    assign_soil(anch, 'mud_soft', project)
    Ha = loads['Ha']
    Va = loads['Va']
    zlug = anch.dd['design']['zlug']
    anch.getCapacityAnchor(Ha, Va, zlug=zlug, plot=False)
    anch.getSafetyFactor()
    

    
    

    
        