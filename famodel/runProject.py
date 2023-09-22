
from project import Project
from shapely.geometry import Point
import matplotlib.pyplot as plt


#### Test script for plotting/analysis of the DeepFarm project ####

# specify the project's area (either a centroid or a lease area)
farm_centroid = (-124.73094, 40.133304)

# initialize a new Project object using the project's area of interest
deepfarm = Project(centroid=farm_centroid)

# set any maps to include in the plots
deepfarm.addMap2GDF(filename='cb_2018_us_state_20m.shp', states=['California'])

# set the farm layout of the project
deepfarm.setFarmLayout(style='grid', nrows=10, ncols=10, turbine_spacing=2000, nOSS=2)

# plot whatever you want of the Project using geopandas
fig, ax = deepfarm.plot({'centroid': {'color':'black'}, 'map': {'color':'tab:blue', 'boundary':True}, 'farm': {'color':'r'}})

# add other things that you want to the plot, if desired
nrel_channel = Point(-120.66088, 34.188565)
nrel_humboldt = Point(-124.73094, 40.133304)
nrel_crescent_city = Point(-124.76659, 41.699739)
hawaii = Point(-157.83, 21)
pointlist = [ nrel_channel, nrel_humboldt, nrel_crescent_city, hawaii ]

deepfarm.addPoints(ax=ax, pointlist=pointlist, kwargs={'pointlist': {'color':'g', 'marker':'x'}})
deepfarm.addState(ax=ax, states=['Hawaii'], kwargs={'Hawaii': {'color':'tab:blue', 'boundary':True}})


plt.show()

a = 2




#### future things ####
# set the bathymetry of the area
# set the soil conditions of the area
# set the metocean conditions of the area
# create different layout options for the turbines
# calculate the closest distance to shore
# populate the model with ports, vessels, etc.
# calculate the depth at any point in the project (lat/long or relative to centroid)

