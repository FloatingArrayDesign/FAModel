''' Example driver file for creating an FAModel project from a YAML file.

This particular example uses the OntologySample200m.yaml file as input.

The output of this file should be a 2-turbine array with a shared mooring line 
going between the turbines, and the bathymetry should be 200m. 

To run without RAFT installed, skip Section 2. To create a Project 
that will automatically create a RAFT model, run Section 2. 

Section 3 shows various modeling capabilities of FAModel
'''

# import necessary packages
from famodel.project import Project
import os

os.chdir('./Inputs/')

# set yaml file location and name
ontology_file = 'OntologySample200m.yaml'

#%% Section 1: Project without RAFT
print('Creating project without RAFT\n')
# create project object
project = Project(file=ontology_file,raft=0)
# create moorpy system of the array, include cables in the system
project.getMoorPyArray(cables=1)
# plot in 3d, using moorpy system for the mooring and cable plots
project.plot3d()



#%% Section 2: Project with RAFT
print('\nCreating project with RAFT \n')
#create project object, automatically create RAFT object (and automatically create moorpy system in the process!)
project = Project(file=ontology_file,raft=1)
# plot in 3d, use moorpy system for mooring and cables, use RAFT for platform, tower, and turbine visuals
project.plot3d(fowt=True)

# get location of RAFT model (stored as array property in project class)
model = project.array
print('Running RAFT case')
# run cases
model.analyzeCases()
# plot results
model.plotResponses()
model.plot()


#%% Section 3: Other capabilities

#### get motion envelopes of platforms and moorings ####
print('\nGetting motion envelopes of platforms and moorings\n')
# loop over each platform object in the farm
for platform in project.platformList.values():
    # get watch circle
    platform.getWatchCircle()
# plot motion envelopes with 2d plot
project.plot2d()

#### get anchor capacities, loads, and safety factors ####
print('\nGetting anchor capacities, loads, and safety factors\n')
# loop over each anchor in the farm
capacities = []
loads = []
sfs = []
for anchor in project.anchorList.values():
    # get anchor capacity
    anchor.getAnchorCapacity()
    capacities.append(anchor.anchorCapacity)
    # get anchor loads at mudline and anchor lug depth (if applicable)
    loads.append(anchor.getLugForces())
    # get safety factor
    sfs.append(anchor.getFS())
    
print('\nAnchor safety factors: ',sfs)
    
#### add marine growth to the mooring lines and cables ####
print('\nAdding marine growth\n')
# marine growth dictionary is read in from YAML, see Ontology ReadMe for description
project.getMarineGrowth()
# moorpy system lines with marine growth are stored in the respective objects under ss_mod (pristine lines are stored under ss)
# check the difference in nominal diameter for a given line:
reg_line_d = project.mooringList['FOWT1a'].ss.lineList[1].type['d_nom']
mg_line_d = project.mooringList['FOWT1a'].ss_mod.lineList[-1].type['d_nom']
print('\nPristine line polyester nominal diameter just below surface: ',reg_line_d,' Marine growth line polyester nominal diameter just below surface: ',mg_line_d)

#### 

