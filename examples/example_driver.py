''' Example driver file for creating an FAModel project from a YAML file.

This particular example uses the OntologySample200m.yaml file as input.

The output of this file should be a 2-turbine array with a shared mooring line 
going between the turbines, and the bathymetry should be 200m. 

To run without RAFT installed, skip Section 2. To create a Project 
that will automatically create a RAFT model, run Section 2. 

Section 3 runs FLORIS/FAModel interface

Section 4 shows various modeling capabilities of FAModel
    - watch circle and motion envelopes of mooring lines
    - calculating anchor capacities and safety factors
    - resizing an anchor for a desired safety factor
    - adding marine growth to cables and mooring lines
'''

# import necessary packages
from famodel.project import Project
import os
import matplotlib.pyplot as plt

os.chdir('./Inputs/')

# set yaml file location and name
ontology_file = 'OntologySample200m.yaml'

#%% Section 1: Project without RAFT
print('Creating project without RAFT\n')
print(os.getcwd())
# create project object
project = Project(file=ontology_file,raft=False)
# create moorpy system of the array, include cables in the system
project.getMoorPyArray(cables=1)
# plot in 3d, using moorpy system for the mooring and cable plots
project.plot2d()
project.plot3d()

#%% Section 2: Project with RAFT
print('\nCreating project with RAFT \n')
#create project object, automatically create RAFT object (and automatically create moorpy system in the process!)
project = Project(file=ontology_file,raft=True)
# plot in 3d, use moorpy system for mooring and cables, use RAFT for platform, tower, and turbine visuals
project.plot3d(fowt=True,draw_boundary=False,boundary_on_bath=False,save=True)

# get location of RAFT model (stored as array property in project class)
model = project.array
model.mooring_currentMod = 0 # temp requirement to work with changes in RAFT
model.ms.moorMod = 0 # temp requirement to work with changes in RAFT
print('Running RAFT case')
# run cases
model.analyzeCases()
# plot results
model.plotResponses()

#%% Section 3: FLORIS
print('Running FLORIS')
config_file = 'gch.yaml' # configuration for running floris
turb_file = 'iea_15MW.yaml' # turbine file 

project.getFLORISArray(config_file,[turb_file],[0,10.59,25],[0,1.95e6,1.9E6])
project.getFLORISMPequilibrium(10.59,0,.06,3,150,plotting=True)


#%% Section 4: Other capabilities

#### get motion envelopes of platforms and moorings ####
print('\nGetting motion envelopes of platforms and moorings\n')
# get watch circle for all platforms in the farm
project.arrayWatchCircle(ang_spacing=20)
# save envelopes from watch circle information for each mooring line
for moor in project.mooringList.values():
    moor.getEnvelope()

# plot motion envelopes with 2d plot
project.plot2d(save=True,plot_bathymetry=False)


#### get anchor capacities, loads, and safety factors ####
print('\nGetting anchor capacities, loads, and safety factors\n')
# let's look at one anchor in the farm

# define anchor to analyze
anchor = project.anchorList['FOWT1a']
# get anchor capacity
anchor.getAnchorCapacity()
capacities = anchor.anchorCapacity
# get anchor loads at mudline and anchor lug depth (if applicable)
loads = anchor.getLugForces()
# size an anchor
starting_geometry = [15,20] # geometry values
starting_geom_labels = ['A','zlug'] # corresponding labels for the geometry list
min_safety_factors = {'Ha':2,'Va':2} # minimum safety factors
FSdiff_max = {'Ha':.1,'Va':.1} # allowable difference between actual and desired FS for final result
anchor.getSize(starting_geometry, starting_geom_labels, minfs=min_safety_factors,
               FSdiff_max=FSdiff_max)
# get safety factor
sfs = anchor.getFS()
    
print('\nAnchor safety factors: ',sfs) # NOTE that Va will show as 'inf' because there is no vertical force on the anchor.
    
#### add marine growth to the mooring lines and cables ####
print('\nAdding marine growth\n')
# marine growth dictionary is read in from YAML, see Ontology ReadMe for description
project.getMarineGrowth(display=False)
# moorpy system lines with marine growth are stored in the respective objects under ss_mod (pristine lines are stored under ss)
# check the difference in nominal diameter for a given line:
reg_line_d = project.mooringList['FOWT1a'].ss.lineList[1].type['d_nom']
mg_line_d = project.mooringList['FOWT1a'].ss_mod.lineList[-1].type['d_nom']
print('\nPristine line polyester nominal diameter just below surface: ',reg_line_d)
print('Marine growth line polyester nominal diameter just below surface: ',mg_line_d)






