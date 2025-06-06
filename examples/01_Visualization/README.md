# Mooring Analysis
Mooring line and cable motions and forces can be determined in FAModel through the FAModel MoorPy integration.

A MoorPy model of the array can automatically be developed from FAModel, and it is then saved as a property of the project object called 'ms'. 

Users can perform their own analyses with the moorpy model, or use one of the built-in analysis tools such as:
- array watch circle analysis from turbine thrust forces and save motion envelopes of moorings, as well as provide maximum forces.


## Build MoorPy System Model of Array
A MoorPy model of the array can be built with the method getMoorPyArray(). This can include information on moorings, dynamic cables, platform hydrostatics, connectors, and anchors, but the getMoorPyArray() method is flexible and will build only what you have included.

            
               
| Options | Min. Required Information | Ontology Sections to Include |
| ---------------- | ---------------|---------------|
| Platform locations | Platform x,y coordinates<br>Platform entity type (i.e. 'FOWT') | Array data table<br>Platforms |
| Mooring lines    | Platform x,y coordinates<br>Platform entity type (i.e. 'FOWT')<br>Platform fairlead radius<br>Platform fairlead depth<br>Mooring line span<br>Anchor type | Array data table<br>Platforms<br>mooring_systems AND/OR array_mooring table<br>mooring_line_configs<br>mooring_line_types<br>anchor_types|
| Cables           | Platform x,y coordinates<br>Platform entity type (i.e. 'FOWT')<br>Cable attachments<br>Dynamic cable span<br>Dynamic cable properties<br>Static cable properties (if applicable)               | Array data table<br>Platforms<br>dynamic_cable_configs<br>\*\*cable_types<br>array_cables table AND/OR cables<br>\*\*cable_appendages |
 Bathymetry       | MoorDyn bathymetry file                                                                                                                                                                        | site -> bathymetry                                                                                                                    |
| Soil             | MoorDyn soil file                                                                                                                                                                              | site -> seabed                                                                                                                        |
| Lease boundaries | Lease  coordinates                                                                                                                                                                             | site -> boundaries                                                                                                                    |

## 3D Plot

3D plots of the array are also flexible, with the required inputs depending on what the user wants to visualize.

The following table describes the required inputs for each optional use of the 3D plotting feature. See the ontology YAMLs in this folder to understand how to format this input information in an ontology file. The name of the ontology file matches the name of the option in the table below. To see all of these options together, check out the ontology YAML entitled 'plotting_3d_all.yaml'. The python script 'plot_3d_examples.py' shows example plots for each option, as well as a plot combining all of the options.

| Options          | Minimum Required Inputs                                                                                                                                                                        | Ontology Sections                                                                                                                     |
| ---------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| Mooring lines    | Platform x,y coordinates<br>Platform entity type (i.e. 'FOWT')<br>Platform fairlead radius<br>Platform fairlead depth<br>Mooring line span<br>Mooring line configuration<br>Mooring Line Types | Array data table<br>Platforms<br>mooring_systems AND/OR array_mooring table<br>mooring_line_configs<br>mooring_line_types             |
| Cables           | Platform x,y coordinates<br>Platform entity type (i.e. 'FOWT')<br>Cable attachments<br>Dynamic cable span<br>Dynamic cable properties<br>Static cable properties (if applicable)               | Array data table<br>Platforms<br>dynamic_cable_configs<br>\*\*cable_types<br>array_cables table AND/OR cables<br>\*\*cable_appendages |
| Platforms        | RAFT platform description<br>Platform x,y coordinates<br>Platform entity type (i.e.'FOWT')<br>Platform fairlead radius<br>Platform fairlead depth                                              | Platforms<br>Array data table                                                                                                         |
| Turbines         | RAFT platform description<br>Platform x,y coordinates<br>Platform entity type (i.e.'FOWT')<br>Platform fairlead radius<br>Platform fairlead depth<br>RAFT turbine description                  | Platforms<br>Array data table<br>Topsides                                                                                             |
 Bathymetry       | MoorDyn bathymetry file                                                                                                                                                                        | site -> bathymetry                                                                                                                    |
| Soil             | MoorDyn soil file                                                                                                                                                                              | site -> seabed                                                                                                                        |
| Lease boundaries | Lease  coordinates                                                                                                                                                                             | site -> boundaries                                                                                                                    |