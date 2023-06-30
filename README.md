# FAModel

The FAModel (or Floating Array Model) package serves as a high-level library for
modeling a floating wind array. It combines site condition information and a 
description of the floating array design, and contains functions for evaluating
the array's behavior considering the site conditions. For example, it combines
information about site soil conditions and an array's anchor characteristics to
estimate the holding capacity of each anchor.

The library works in conjunction with the tools RAFT and MoorPy to model floating
wind turbines and mooring systems, respectively.

An example of use of these tools to model three mooring lines over the bathymetry 
of the Humboldt lease area is shown below.

![Humboldt](famodel/seabed/images/slopeview4.PNG)


## Subpackages

The library has a core Project class for organizing information and an evolving
collection of subpackages for specific functions. The two current subpackges are:

- anchors: contains modules for anchor capacity calculations
- seabed: contains modules for seabed bathymetry and boundary information

Please navigate into the subfolders above for additional information.


## Authors

The NREL Floating Wind Array Design team.