# -*- coding: utf-8 -*-
"""
Example on how to make a project class and fill in manually rather than using a yaml file.
"""

from famodel import Project
import os
import moorpy as mp
import matplotlib.pyplot as plt

### INPUTS ###
dir = os.path.dirname(os.path.realpath(__file__))
bathfile = os.path.join(dir,'bathymetry200m_Array.txt')
soilfile = os.path.join(dir,'soil_sample.txt')
boundfile = os.path.join(dir,'boundary_sample.csv')
moordynfile = os.path.join(dir,'Moordyn_semitaut200m.dat')
n_pfs = 4
pf_spacing = 1600


# create empty project
project = Project()


### populate the Project class ###

# easiest way to do this is with a moorpy system 
ms = mp.System(file=moordynfile)
# create platform, moorings, anchors from moorpy system
pflocs = [-800,-800,0]
project.addPlatformMS(ms,pflocs)

# lets add a few more platforms, this time by duplicating the platform we just made
for i in range(n_pfs-1):
    # update location to place platform (this logic set will make a z-formation)
    if i%2 == 0:
        pflocs[0] += pf_spacing
    else:
        pflocs[1] += pf_spacing
    # copy platform we made from ms
    project.duplicate(project.platformList['fowt0'],r=pflocs)
    
# add bathymetry from file (alternatively, set depth in project initialization input in line above)
project.loadBathymetry(filename=bathfile)
# add soil info from file (only required if doing anchor capacity calc's)
project.loadSoil(filename=soilfile)
# add boundary from file
project.loadBoundary(filename=boundfile)
    

# create moorpy system and plot
project.getMoorPyArray()
project.plot3d()
plt.show()




