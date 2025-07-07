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
ontology_file = 'OntologySample200m_1turb.yaml'

# create project class
project = Project(file=ontology_file)
project.getMoorPyArray()

# let's choose a single anchor from the array to look at
anch = project.anchorList['fowt0a']

# now let's get the mudline and lug forces on this anchor
anch.getLugForces() # getLugForces calls getMudlineForces() to get the anchor forces at both locations

# establish a factor of safety in horizontal (Ha) and vertical (Va) directions
minfs = {'Ha': 1.8, 'Va': 2}

# let's get the loads with the factor of safety included 
loads_with_FS = {'Ha':anch.loads['Ha']*minfs['Ha'],'Va':anch.loads['Va']*minfs['Va']}

# get anchor capacity for one anchor (this case is for suction pile in clay)
anch.getAnchorCapacity(loads=loads_with_FS) # loads are used in capacity calculation, so let's send in the loads with factor of safety applied

# get anchor cost
startGeom = [10,2,6.6]
geomKeys = ['L','D','zlug']
geomBounds = [(5, 50), (1, 7), (3.3,16.7)]
FSDiff_max = {'Ha':5,'Va':5}
anch.getSize(startGeom,geomKeys,geomBounds,minfs=minfs,FSdiff_max=FSDiff_max, plot=True)
anch.getCost()
print('\nClay suction pile capacity is: ',anch.anchorCapacity)
print('Clay suction pile safety factor is: ',anch.getFS())
print('Clay suction pile cost is: ', anch.cost,'\n')
# try suction pile with sand
newdd = anch.dd
anch.soilProps['sand'] = anch.soilProps.pop('mud_firm')
anch.soilProps['sand']['phi'] = 33
anch.soilProps['sand']['Dr'] = 70
anch.soilProps['sand']['delta'] = 25
# update anchor loads at lug point (mudline load should be constant), then get anchor capacity
anch.getLugForces()
anch.getSize(startGeom,geomKeys,geomBounds,plot=True)
anch.getAnchorCapacity(loads=loads_with_FS)
anch.getCost()
print('\nSand suction pile capacity is: ',anch.anchorCapacity,' N')
print('Sand suction pile safety factor is: ',anch.getFS())
print('Sand suction pile cost is: ', anch.cost,' USD\n')

# check plate anchor type   
newdd['type'] = 'DEA'
newdd['design'] = {'type':'DEA','A':20,'zlug':20,'beta':10}
anch.soilProps['clay'] = anch.soilProps.pop('sand')

startGeom = [10,20]
geomKeys = ['A','zlug']

anch.getLugForces()
anch.getAnchorCapacity(loads=loads_with_FS)
# let's fix the zlug for the plate anchor - set fix_zlug=True to prevent it being changed
anch.getSize(startGeom,geomKeys,minfs={'Ha':2,'Va':0}, fix_zlug = True)
anch.getCost()
print('\nClay plate capacity is: ',anch.anchorCapacity,' N')
print('Clay plate safety factor is: ',anch.getFS())
print('Clay plate cost is: ', anch.cost,' USD\n')

# check drilled and grouted pile anchor type
newdd['type'] = 'dandg_pile'
newdd['design'] = {'type':'dandg_pile','L':50,'D':3,'zlug':0}
anch.soilProps['rock'] = anch.soilProps.pop('clay') # soil_properties has default rock info in there already, just change name

# startGeom = [5,50]
# geomKeys = ['L','D']
anch.getLugForces()
# anch.getSize(startGeom,geomKeys,minfs={'Ha':2,'Va':2})
anch.getAnchorCapacity(loads=loads_with_FS)
print('\nRock drilled and grouted pile capacity is: ',anch.anchorCapacity,' N')
print('Rock drilled and grouted pile safety factor is: ',anch.getFS())

# check driven pile anchor in rock
newdd['type'] = 'driven'
anch.soilProps['weak_rock'] = anch.soilProps.pop('rock')
newdd['design'] = {'type':'driven','L':20,'D':1.5,'zlug':-3} # zlug should be negative (above mudline) for rock!

anch.getLugForces()
anch.getAnchorCapacity(loads=loads_with_FS)
print('\nWeak rock driven pile capacity is: ',anch.anchorCapacity,' N')
print('Weak rock driven pile safety factor is: ',anch.getFS())

# check driven pile anchor in clay
anch.soilProps['clay'] = anch.soilProps.pop('weak_rock')
newdd['design'] = {'type':'driven','L':40,'D':4,'zlug':10}

anch.getLugForces()
anch.getAnchorCapacity(loads=loads_with_FS)
print('\nClay driven pile capacity is: ',anch.anchorCapacity,' N')
print('Clay driven pile safety factor is: ',anch.getFS())

# check driven pile anchor in sand
anch.soilProps['sand'] = anch.soilProps.pop('clay')
anch.soilProps['sand']['Dr'] = 50

anch.getLugForces()
anch.getAnchorCapacity(loads=loads_with_FS)
print('\nSand driven pile capacity is: ',anch.anchorCapacity,' N')
print('Sand driven pile safety factor is: ',anch.getFS())

# check helical pile anchor with sand
newdd['type'] = 'helical_pile'
newdd['design'] = {'type':'helical_pile','L':25.1,'d':1,'D':5.01, 'zlug':5}

anch.getLugForces()
anch.getAnchorCapacity(loads=loads_with_FS)
print('\nSand helical pile capacity is: ',anch.anchorCapacity,' N')
print('Sand helical pile safety factor is: ',anch.getFS())

# check helical pile anchor with clay
anch.soilProps['clay'] = anch.soilProps.pop('sand')
newdd['type'] = 'helical_pile'
newdd['design'] = {'type':'helical_pile','L':25.1,'d':1,'D':5.01,'zlug':5}

anch.getLugForces()
anch.getAnchorCapacity(loads=loads_with_FS)
print('\nClay helical pile capacity is: ',anch.anchorCapacity,' N')
print('Clay helical pile safety factor is: ',anch.getFS())

# check torpedo anchor in clay
newdd['type'] = 'torpedo_pile'
newdd['design'] = {'type':'torpedo_pile','D1':3,'D2':1.1,'L1':10,'L2':4,'zlug':16}

anch.getLugForces()
anch.getAnchorCapacity(loads=loads_with_FS)
print('\nClay torpedo pile capacity is: ',anch.anchorCapacity,' N')
print('Clay torpedo pile safety factor is: ',anch.getFS())





