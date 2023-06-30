# Project class and floating array subpackages

This in-progress package currently includes the start of an
overall floating array project class and a series of data
structures and functions for bathymetry, soil conditions, 
and anchor holding capacity.

## Core Class

### Project

The project class provides a standard data structure that combines
design information (through RAFT) and site information. 
Currently the site information focuses on the seabed and lease area
boundary. In the future, metocean data will also be included.

The project class includes a grid over which seabed bathymetry and
soil information is stored. The soil information follows a property
parameterization with the following fields:
- soil_class - soil classification name (TBD)
- soil_gamma - soil effective unit weight [kPa] (all soils)
- soil_Su0   - undrained shear strength at mudline [kPa] (clay soils)
- soil_K     - undrained shear strength gradient [kPa/m] (clay soils)
- soil_alpha - soil skin friction coefficient [-] (clay soils)
- soil_phi   - angle of internal friction [deg] (sand soils)


## Subpackages

### anchors

The anchors subpackage contains modules for anchor capacity calculations
as a function of soil type. It has the following:
- level-1 capacity curves for general anchor types and soil classes
- level-2 capacity functions for select anchor types and quantitative
  soil properties
- a general anchor capacity calculation 'switchboard' that routes to the
  appropriate functions for specified anchor and soil types.
  
### seabed

The seabed subpackage provides a set of functions for inputting and 
processing "seabed" information, including bathymetry, soil properties,
and other spatial properties of a lease area such as the lease area
boundary. 



