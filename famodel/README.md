# Project class and floating array subpackages

This in-progress package currently includes the start of an
overall floating array project class and a series of data
structures and functions for bathymetry, soil conditions, 
and anchor holding capacity.

## The Project Class

The Project class provides a standard data structure that combines
design information (through RAFT) and site information. 
Currently the site information focuses on the seabed and lease area
boundary. In the future, metocean data will also be included.

The project class includes a grid over which seabed bathymetry and
soil information is stored. The soil information follows a property
parameterization with the following fields:
- soil_class - soil classification name:
  - soft clay
  - medium clay
  - hard clay
  - sand
- soil_gamma - soil effective unit weight [kPa] (all soils)
- soil_Su0   - undrained shear strength at mudline [kPa] (clay soils)
- soil_K     - undrained shear strength gradient [kPa/m] (clay soils)
- soil_alpha - soil skin friction coefficient [-] (clay soils)
- soil_phi   - angle of internal friction [deg] (sand soils)

## Subpackages

### Anchors

The anchors subpackage contains modules for anchor capacity calculations
as a function of soil type. It has the following:
- level-1 capacity curves for the following general anchor types and soil classes:
  - DEAs in soft clay, medium clay, hard clay, or sand
  - VLAs in soft clay or medium clay
  - Suction anchors in soft clay or medium clay
  - SEPLAs in soft clay
- level-2 capacity functions for select anchor types and quantitative
  soil properties
- a general anchor capacity calculation 'switchboard' that routes to the
  appropriate functions for specified anchor and soil types.
  
### Seabed

The seabed subpackage provides a set of functions for inputting and 
processing "seabed" information, including bathymetry, soil properties,
and other spatial properties of a lease area such as the lease area
boundary. 

## Project Methods

### setGrid
        
Set up the rectangular grid over which site or seabed
data will be saved and worked with. Directions x and y are 
generally assumed to be aligned with the East and North 
directions, respectively, at the array reference point.
    
### loadBoundary
Load a lease area boundary for the project from an input file.
        
### loadBathymetry

Load bathymetry information from an input file (format TBD), convert to
a rectangular grid, and save the grid to the floating array object (TBD).
        
### loadSoil

Load geoetechnical information from an input file (format TBD), convert to
a rectangular grid, and save the grid to the floating array object (TBD).

The input file should provide rows with the following entries:
- x coordinate
- y coordinate
- class  - soil classification name ('clay', 'sand', or 'rock' with optional modifiers)
- gamma* - soil effective unit weight [kPa] (all soils)
- Su0*   - undrained shear strength at mudline [kPa] (clay 
- K*     - undrained shear strength gradient [kPa/m] (clay 
- alpha* - soil skin friction coefficient [-] (clay soils)
- phi*   - angle of internal friction [deg] (sand soils)

Some (*) parameters are optional depending on the soil class and mode.   

Irregular sampling points will be supported and interpolated to a 
rectangular grid.

### getSoilAtLocation

Interpolate soil properties at specified location from the soil
properties grid and return a dictionary of soil properties that
can be used in anchor capacity calculations.

### calcAnchorCapacity

Compute holding capacity of a given anchor based on the soil
info at its position. The anchor object's anchor properties and
location will be used to determine the holding capacity, which
will be saved to the anchor object.
        

