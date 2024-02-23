# Project class and floating array subpackages

This in-progress package currently includes the start of an
overall floating array project class and a series of data
structures and functions for bathymetry, soil conditions, 
and anchor holding capacity.

## The Project Class

The Project class provides a standard data structure that combines
design information (through an ontology yaml file or RAFT) and site information. 
Currently the site information focuses on the seabed and lease area
boundary. In the future, metocean data will also be included. The design information 
includes anchors, platforms, connectors, cables, and mooring lines.

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

### loadDesign

Load in the design information from a YAML file or dictionary. Loads in cables,
mooring lines, anchors, and connectors. Creates objects for each mooring line, anchor,
connector, and cable. Lists of mooring line objects, anchor objects, and connector objects 
are created and stored in the project (array level) object.

### loadSite

Load site information from a YAML file or dictionary. Loads in bathymetry from a file using
the loadBathymetry method, or from xyz data. Loads in boundaries from a file using the loadBoundary 
method, or loads in the boundaries from xy data. Additional site information such as general depth,
water density, air density, and air dynamic viscosity are also loaded here.

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
        
### getDepthAtLocation

Computes the depth at a specified x-y location based on the 
bathymetry grid stored in the project. The seabed normal vector
can also be obtained.

### seabedIntersect

Computes the location at which a ray crosses the seabed

### projectAlongSeabed

Obtain the z-coordinates of a projection along the seabed surface for a set of x-y coordinates.

### makeDistanceMatrix

Computes the horizontal distance between every turbine's undisplaced location 
and returns a matrix of these distances for the array.

### calcCableLength

Calculates the cable's length based on its routing

### checkCableExclusions

Checks whether a cable crosses over any exclusions or other out of bounds areas

### plot3d

Plots aspects of the Project object in matplotlib 3D

### plot2d

Plots aspects of the Project object into matplotlib in 2D

### getMoorPyArray

Creates an array in MoorPy based on the mooring, anchor, connector, and platform objects listed 
in the project object. Also allows the option to plot the resulting MoorPy array.

### getFromDict

Streamlines getting values from the design dictionary from YAML file, including error checking.

## Ontology YAML

Contains all information on the design of the floating array.

Key information for building the floating array design yaml file by section:

### Site
This section contains information about the site, either through listing of data file names, or site information listed in table format.
Subsections include general information, boundaries, exclusions, bathymetry, seabed, metocean, and resource.

### Array
This section contains a table of information about the turbine array layout. Each row 
in the table represents one turbine. The turbineID refers to the turbine type in the turbine 
section, the platformID refers to the platform type from the platform section, the mooringID refers 
to the mooring system number in the mooring_systems section. The x_location and y_location refer to the 
xy coordinates for the center of the platform, and the heading adjust describes any rotation of the platform's 
mean position.

If 0 is provided for the mooring system ID, this designates that the mooring details for that 
turbine are in the array_mooring section, which is generally in use for shared moorings or anchors.

### Array Mooring
This section contains array-level data on anchors and lines. The anchors are listed in the anchor_data table,
where the xy location, anchor type, and embedment information are provided. The line_data table provides 
information on each line, including which platform row number (from the array table) and/or anchor number (from the
anchor_data table) is attached to the line, what the headings for each end are, and whether there is any length adjustment.

headingA refers to the mooring heading of the turbine the line is connected to for end A of the line, if it is connected to 
a platform (for all non-shared moorings, end A is connected to the anchor, and in this instance NA can be used for headingA). 
Similarly, headingB refers to the mooring heading of the platform the line is connected to for end B of the line. End B must be 
connected to a platform.

### Array Cables

### Turbines

### Platforms

### Mooring Systems
This section lists the mooring system for each turbine, not including any shared lines or lines with shared anchors.
The list of mooring systems is numbered by msX where X is the mooring system ID number provided in the array table.
Each mooring system provides the ID of the configuration, the heading, the anchor type, and whether there is any length adjustment.

Platforms that have shared mooring or shared anchors for some of their lines may have a listing in the mooring system (if the MooringID 
is not 0 for that platform in the array section), but this mooring system will only contain mooring line configurations that are not shared 
moorings or anchors.

### Mooring Line Configurations
This section lists the mooring configurations that are called out i nthe mooring_systems section or the array_mooring section. For each 
configuration, a descriptive name, fairlead and anchor radius, fairlead depth, and list of sections is provided. Items in the sections list 
may be mooring lines (in which case type and length are given as keys) or connectors (in which case the key connectorType must be used). Connectors are 
optional, and some, all, or none of the connectors for a line may be provided. However, only one connector may be listed between each line type listed.

**Please note that the order of the sections list is very important. Entries in the sections list must be provided in order from end A 
(anchor end for non-shared lines) to end B (fairlead end).**

For shared line configurations, a key 'symmetric' may be used to describe whether the provided line configuration is half of a symmetric line 
(symmetric: True) or a full line configuration (symmetric: False) in the same level as the name and radii information. If a connector is given 
as the last item in the section list for a shared symmetric line, it is assumed that the provided connector is located in the center of the line.

### Mooring Line Types
This section provides design details for each line type listed in the line configurations section. Necessary information to be listed for each
line type includes nominal diameter, volume-equivalent diameter, mass per unit length, quasi-static stiffness, minimum breaking load, and material. 
Optional information includes cost per unit length, dynamic stiffness, dynamic stiffness mean-load multiplier, and material details.

### Mooring Connector Types
This section provides design details for each mooring connector type listed in the line configurations section. There are no requirements on what
must be listed for each connector, but optional inputs include mass, volume, and CdA (product of cross sectional area and drag coefficient).

### Anchor Types
This section provides design detials for each anchor type listed in the mooring_systems section or the array_mooring anchor_data table. There are no 
requirements on what must be listed for each anchor, but optional inputs include diameter, length, mass, area, thickness, and embedment length.

### Cables

### Dynamic Cable Configurations

### Cable Types

### Cable Appendages

### Cable Joints