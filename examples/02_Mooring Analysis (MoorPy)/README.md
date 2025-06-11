# Mooring Analysis
Mooring line and cable motions and forces can be determined in FAModel through the FAModel MoorPy integration.

A MoorPy model of the array can automatically be developed from FAModel, and it is then saved as a property of the project object called 'ms'. 

Users can perform their own analyses with the moorpy model, or use one of the built-in analysis tools such as:
- arrayWatchCircle() - array watch circle analysis from turbine thrust forces
    - this will save motion envelopes of moorings and platforms, as well as provide maximum forces for each mooring, cable, and anchor in the array.


## Build MoorPy System Model of Array
A MoorPy model of the array can be built with the method getMoorPyArray(). This can include information on moorings, dynamic cables, platform hydrostatics, connectors, and anchors, but the getMoorPyArray() method is flexible and will build only what you have included.

            
               
| Options | Min. Required Information | Ontology Sections to Include |
| ---------------- | ---------------|---------------|
| Mooring lines    | Platform x,y coordinates<br>Platform entity type (i.e. 'FOWT')<br>Platform fairlead radius<br>Platform fairlead depth<br>Mooring line span<br>Anchor type<br>Water depth/bathymetry | Array data table<br>Platforms<br>mooring_systems AND/OR array_mooring table<br>mooring_line_configs<br>mooring_line_types<br>anchor_types<br>site->general->water depth OR site->bathymetry|
| Cables           | Platform x,y coordinates<br>Platform entity type (i.e. 'FOWT')<br>Cable attachments<br>Dynamic cable span<br>Dynamic cable properties<br>Static cable properties (if applicable)<br>Water depth/bathymetry               | Array data table<br>Platforms<br>dynamic_cable_configs<br>\*\*cable_types<br>array_cables table AND/OR cables<br>\*\*cable_appendages <br>site->general->water depth OR site->bathymetry|

It is important to note that each FAModel component object (except static cables) has a MoorPy equivalent (body, point, subsystem, or line) and that equivalent is created automatically based on the information stored in FAModel component objects. The moorpy component objects are stored in the FAModel component objects, and the main moorpy system (which connects all of the MoorPy objects) is stored in the project class under the property ms.

## Array Watch Circle

The built-in FAModel arrayWatchCircle() method can be called to run a watch circle analysis on the array. This will apply a thrust force (either provided or default to IEA 15 MW reference turbine) at a set of angles and record the motion envelopes of the platforms and moorings, as well as record the maximum anchor, mooring, and cable forces. It will also compute and record the safety factors of the moorings and cables (including curvature safety factor for cables).
These results can be found in the loads property of the mooring, anchor, and cable objects.
The motion envelopes can be visualized by running plot2d() or the coordinates can be accessed in the envelopes property of the object.