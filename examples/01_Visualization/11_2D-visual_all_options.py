# -*- coding: utf-8 -*-
"""
Simple driver file to create 2D plots exploring the different settings options
for plot2d().
"""

from famodel import Project
import matplotlib.pyplot as plt
import os

# define name of ontology input file
dir = os.path.dirname(os.path.realpath(__file__))
input_file = os.path.join(dir,'11_2D-visual_all_options.yaml')

# initialize Project class with input file, we don't need RAFT for this so mark False
project = Project(file=input_file,raft=False)

# create a dictionary of plot settings. We list here all the options showing default settings - comment/change different ones to see what they do!
plot_2d_settings = {
    'ax' : None, # matplotlib.pyplot axis to plot on. If None, a new figure & axis is created.
    'plot_soil' : False, # bool, if True plots soil
    'plot_bathymetry' : True, # bool, if True plots bathymetry. plot_soil & plot_bathymetry can not both be True, but you can plot bathy contours with soil.
    'plot_boundary' : True,  # bool, if True plots lease boundary
    'color_lineDepth' : False, # If True, color mooring lines based on depth. Only works if plot_bathymetry=False.
    'plot_bathymetry_contours' : False, # bool, if True plots bathymetry contour lines. Can be used with plot_soil or plot_bathymetry
    'bare' : False, # bool, if True does not plot extra things like colorbars
    'axis_equal' : True, # bool, if True plots x and y axes at same scale to prevent distortions
    'save' : False, # bool, if True saves the figure
    'figsize' : (8,8),  # the dimensions of the figure to be plotted
    'env_color' : [.5,0,0,.8], # platform motion envelope color (need to run project.arrayWatchCircle() to get this)
    'fenv_color' : [.6,.3,.3,.6], # mooring motion envelope color (need to run project.arrayWatchCircle() to get this)
    'alpha' : 0.5, # the opacity of the plot
    'return_contour' : False, # bool, if True returns the bathymetry or soil filled contour (only 1 of these can be plotted at a time)
    'cmap_cables' : None, # matplotlib colormap string for cable conductors sizes
    'cmap_soil' : None, # matplotlib colormap string for soil types
    'plot_platforms' : True, # bool, if True plots the platform locations
    'plot_anchors' : True, # bool, if True plots the anchor locations
    'plot_moorings' : True, # bool, if True plots the mooring lines
    'plot_cables' : True, # bool, if True plots the cables
    'cable_labels' : False, # bool, if True adds labels of cable names to the plot
    'depth_vmin' : None, # minimum depth to plot. If None, defaults to minimum depth in the bathymetry information
    'depth_vmax' : None, # maximum depth to plot. If None, defaults to maximum depths in the bathymetry information
    'bathymetry_levels' : 50, # number of bathymetry contour levels, for either contour lines or the filled contour
    'plot_legend' : True, # bool, if True plots the legend
    'legend_x' : 0.5, # x location of the legend
    'legend_y' : -0.1, # y location of the legend
    'plot_landmask': False, # bool, if True plots a gray mask over land areas     
    'max_line_depth': None,  # max depth for line coloring if color_lineDepth is True
    'only_shared' : False,   # if color_lineDepth is True, only color shared lines
    'linewidth_multiplier' : 2,  # multiplier for line widths if color_lineDepth is True
    'soil_alpha' : 0.5 # opacity of soil colormap
    }


project.plot2d(**plot_2d_settings) # unpack settings dictionary and plot

# Let's change some settings!
plot_2d_settings['cable_labels'] = True # add cable labels
plot_2d_settings['plot_bathymetry'] = False # turn off bathymetry
plot_2d_settings['plot_soil'] = True # turn on soil 
plot_2d_settings['bathymetry_levels'] = 6 # use 6 bathymetry levels for contour lines
plot_2d_settings['plot_bathymetry_contours'] = True # turn on contour lines
plot_2d_settings['cmap_soil'] = 'Greens' # change soil colormap to shades of green
plot_2d_settings['soil_alpha'] = 0.35 # change soil opacity to 35%

project.plot2d(**plot_2d_settings)



plt.show()
