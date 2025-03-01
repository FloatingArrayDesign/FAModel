# -*- coding: utf-8 -*-
"""
Example showing how to call forces and 
anchor capacity functions, along with safety factors and material costs.
"""
# import necessary packages
from famodel.project import Project
import os

os.chdir('./Inputs/')

# set yaml file location and name
ontology_file = 'OntologySample200m_noshared.yaml'

# create project class
project = Project(file=ontology_file)

# let's choose a single anchor from the array to look at
anch = project.anchorList['FOWT1a']

# now let's get the mudline and lug forces on this anchor
anch.getLugForces() # getLugForces calls getMudlineForces() to get the anchor forces at both locations

# get anchor capacity for one anchor (this case is for suction pile in clay)
anch.getAnchorCapacity()
# get anchor cost
anch.getCost()
print('\nClay suction pile capacity is: ',anch.anchorCapacity)
print('Clay suction pile safety factor is: ',anch.getFS())
print('Clay suction pile cost is: ', anch.cost,'\n')

# try suction pile with sand
newdd = anch.dd
newdd['soil_type'] = 'sand'
newdd['soil_properties']['phi'] = 33
newdd['soil_properties']['Dr'] = 70
newdd['soil_properties']['delta'] = 25

# update anchor loads at lug point (mudline load should be constant), then get anchor capacity
anch.getLugForces()
anch.getSize(startGeom=anch.dd['design'])
# anch.getAnchorCapacity()
anch.getCost()
print('\nSand suction pile capacity is: ',anch.anchorCapacity,' N')
print('Sand suction pile safety factor is: ',anch.getFS())
print('Sand suction pile cost is: ', anch.cost,' USD\n')

# check plate anchor type   
newdd['type'] = 'DEA'
newdd['design'] = {'type':'DEA','A':20,'zlug':10,'beta':10}
newdd['soil_type'] = 'clay'

anch.getLugForces()
anch.getAnchorCapacity()
anch.getCost()
print('\nClay plate capacity is: ',anch.anchorCapacity,' N')
print('Clay plate safety factor is: ',anch.getFS())
print('Clay plate cost is: ', anch.cost,' USD\n')

# check drilled and grouted pile anchor type
newdd['type'] = 'dandg_pile'
newdd['design'] = {'type':'dandg_pile','L':50,'D':3,'zlug':0}
newdd['soil_type'] = 'rock' # soil_properties has default rock info in there already, just change name

anch.getLugForces()
anch.getAnchorCapacity()
print('\nRock drilled and grouted pile capacity is: ',anch.anchorCapacity,' N')
print('Rock drilled and grouted pile safety factor is: ',anch.getFS())

# check driven pile anchor in rock
newdd['type'] = 'driven'
newdd['soil_type'] = 'weak_rock'
newdd['design'] = {'type':'driven','L':20,'D':1.5,'zlug':-3} # zlug should be negative (above mudline) for rock!

anch.getLugForces()
anch.getAnchorCapacity()
print('\nWeak rock driven pile capacity is: ',anch.anchorCapacity,' N')
print('Weak rock driven pile safety factor is: ',anch.getFS())

# check driven pile anchor in clay
newdd['soil_type'] = 'clay'
newdd['design'] = {'type':'driven','L':40,'D':4,'zlug':10}

anch.getLugForces()
anch.getAnchorCapacity()
print('\nClay driven pile capacity is: ',anch.anchorCapacity,' N')
print('Clay driven pile safety factor is: ',anch.getFS())

# check driven pile anchor in sand
newdd['soil_type'] = 'sand'
newdd['soil_properties']['Dr'] = 50

anch.getLugForces()
anch.getAnchorCapacity()
print('\nSand driven pile capacity is: ',anch.anchorCapacity,' N')
print('Sand driven pile safety factor is: ',anch.getFS())

# check helical pile anchor with sand
newdd['type'] = 'helical_pile'
newdd['design'] = {'type':'helical_pile','L':25.1,'d':1,'D':5.01, 'zlug':5}

anch.getLugForces()
anch.getAnchorCapacity()
print('\nSand helical pile capacity is: ',anch.anchorCapacity,' N')
print('Sand helical pile safety factor is: ',anch.getFS())

# check helical pile anchor with clay
newdd['soil_type'] = 'clay'
newdd['type'] = 'helical_pile'
newdd['design'] = {'type':'helical_pile','L':25.1,'d':1,'D':5.01,'zlug':5}

anch.getLugForces()
anch.getAnchorCapacity()
print('\nClay helical pile capacity is: ',anch.anchorCapacity,' N')
print('Clay helical pile safety factor is: ',anch.getFS())

# check torpedo anchor in clay
newdd['type'] = 'torpedo_pile'
newdd['design'] = {'type':'torpedo_pile','D1':3,'D2':1.1,'L1':10,'L2':4,'zlug':16}

anch.getLugForces()
anch.getAnchorCapacity()
print('\nClay torpedo pile capacity is: ',anch.anchorCapacity,' N')
print('Clay torpedo pile safety factor is: ',anch.getFS())





