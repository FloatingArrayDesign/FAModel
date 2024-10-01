# FAModel Examples

This folder contains examples for various functionalities of FAModel.
These examples are designed to help users get acquainted with the package and how it can be used.
The sample .yaml files may also be used as templates for the user's own floating array yamls.
The examples are described below:

## FAModel project from a yaml file example
This example shows how to make an FAModel project from a yaml file. The sample shows various possibilities for mooring and anchoring including shared moorings and shared anchors. It also shows some functionalities of FAModel including:
- integration with MoorPy
- integration with RAFT
- adding marine growth to moorings and cables
- calculating capacity, loads, and safety factors on an anchor
- determining the watch circle of a platform, with maximum tensions and minimum safety factors for all mooring lines and cables attached to the platform
- plotting the 2d plan view of a farm, including: 
    - bathymetry
    - boundaries 
    - mooring lines, cables (including routing), platforms, and substations
    - motion envelopes for mooring lines and platforms
- plotting the 3d plan view of a farm, including: 
    - bathymetry
    - boundaries
    - mooring lines, platforms, and substations
    - dynamic cables
    - static cables, including burial depth and routing
## FAModel project manual example
This example shows how to make an FAModel project with platforms, anchors, and mooring lines without a yaml file. 

## Anchor capacity examples
This example calculates the anchor capacity for each different type of anchor, determines the load on the anchor, and then determines the safety factor. 
This is provided to show an example of what information is required for each anchor type, and how results can differ for different anchor types.