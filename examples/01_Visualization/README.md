# Visualization 
FAModel arrays can be visualized in 2D or 3D from the project class.
The plotting tools are designed to be flexible and can plot varied levels of detail depending on what information is provided.


## Plotting Options and Input Requirements
Plots of an array can be as complex or simple as the user wants, with the required inputs depending on what needs to be visualized.

The following table describes the required inputs for each optional use of the 2D and 3D plotting features. See the ontology YAMLs in this folder to understand how to format this input information in an ontology file. The name of the ontology file matches the name of the driver file that will run the visualization. To see all of these options together, check out the ontology YAML entitled 'OntologySample200m.yaml' in the Inputs sub-folder of the examples folder. The python script 'example_driver.py' shows plots for all of these options together, as well as other functionality in FAModel.

            
               
| Options | Min. Required Information | Ontology Sections to Include |
| ---------------- | ---------------|---------------|
| Platform locations<br> (2D only) | Platform x,y coordinates<br>Platform entity type (i.e. 'FOWT') | Array data table<br>Platforms |
| Mooring lines    | Platform x,y coordinates<br>Platform entity type (i.e. 'FOWT')<br>Platform fairlead radius<br>Platform fairlead depth<br>Mooring line span | Array data table<br>Platforms<br>mooring_systems AND/OR array_mooring table<br>mooring_line_configs<br>mooring_line_types|
| Anchors          | Platform x,y coordinates<br>Platform entity type (i.e. 'FOWT')<br>Platform fairlead radius<br>Platform fairlead depth<br>Mooring line span |Array data table<br>Platforms<br>mooring_systems AND/OR array_mooring_table<br>mooring_line_configs<br>mooring_line_types<br>anchor_types|
| Cables           | Platform x,y coordinates<br>Platform entity type (i.e. 'FOWT')<br>Cable attachments<br>Dynamic cable span<br>Dynamic cable properties<br>Static cable properties (if applicable)               | Array data table<br>Platforms<br>dynamic_cable_configs<br>\*\*cable_types<br>array_cables table AND/OR cables<br>\*\*cable_appendages |
|Bathymetry       | MoorDyn bathymetry file                                                                                                                                                                         site -> bathymetry                                                                                                                    |
| Lease boundaries | Lease  coordinates                                                                                                                                                                             | site -> boundaries                                                                                                                    |
| Substations/<br>Various platform types  | Platform x,y coordinates<br>Platform entity type (i.e. 'Substation') | Array data table<br>Platforms |
|Platform shape| Platform x,y coordinates<br> Platform entity type<br> RAFT platform description<br>water depth| Array data table<br>Platforms (with full RAFT description)<br>Site->general->water depth|
|Turbine shape| Platform x,y coordinates<br>Platform entity type<br>RAFT platform description<br>RAFT turbine description<br>water depth|Array data table<br>Platforms (with full RAFT description)<br>Topsides (with full turbine RAFT description)<br>Site->general->water depth| 

