import os
from project import Project
from shapely.geometry import Point
import matplotlib.pyplot as plt


#### Model of the StandMoor baseline design array ####

lease_area_name = 'Humboldt_SW'
bathymetry_moorpy = 'bathymetry_humboldt-sw_gebco.txt'
bathymetry_gebco = 'bathymetry/gebco_2023_n41.3196_s40.3857_w-125.2881_e-123.9642.asc'

# initialize a new Project object
humboldt_sw = Project(file=os.path.dirname(os.path.realpath(__file__)) + '\Project_Input_File.yaml')         # could have also set the 'centroid' parameter here

# set the lease area boundary of the project 
#humboldt_sw.setLeaseArea(lease_area_name)

# set the bathymetry of the surrounding project
#humboldt_sw.loadBathymetry(bathymetry_gebco, moorpy_bathymetry_filename=bathymetry_moorpy)

# plot the project in 3D
humboldt_sw.plot3d(area=True, bathymetry=bathymetry_moorpy, boundary=False, area_on_bath=True)





#fig, ax = humboldt_sw.plot({'boundary':{'color':'r'}})

# set any maps to include in the plots
#deepfarm.addMap2GDF(filename=os.path.dirname(os.path.realpath(__file__))+'/cb_2018_us_state_20m.shp', states=['California'])

# set the farm layout of the project
#deepfarm.setFarmLayout(style='shared', nrows=10, ncols=10, turbine_spacing=2000, nOSS=2)

# plot whatever you want of the Project using geopandas
#fig, ax = deepfarm.plot({'centroid': {'color':'black', 'label':True}, 'map': {'color':'tab:blue', 'boundary':True}, 'farm': {'turbine': {'color':'r'}, 'oss': {'color':'b'}}})
#fig, ax = deepfarm.plot({'centroid': {'color':'black', 'label':True}, 'farm': {'turbine': {'color':'r'}, 'oss': {'color':'b'}}})
#fig, ax = deepfarm.plot({'map': {'color':'tab:blue', 'boundary':True}})


plt.show()

a = 2



