# Moorings, Sections, and Connectors

This file contains information on Mooring class and the subcomponent classes, [Section](#the-section-class) and [Connector](#the-connector-class).
The layout of this file is as follows:
## Layout
* [The Mooring Class](#the-mooring-class)
	* [Mooring Properties](#mooring-properties)
	* [Mooring Methods](#mooring-methods)
* [The Section Class](#the-section-class)
* [The Connector Class](#the-connector-class)


## The Mooring Class

The Mooring class provides a data structure for design information of a mooring line. Design information 
includes a design dictionary with the following details:
- Anchoring radius
- Anchor depth
- Fairlead radius
- Fairlead depth 
- Detail on each line section
	- Line section type
		- diameter (nominal and volume-equivalent)
		- material
		- cost
		- mass
		- MBL
		- EA
	- Line section length
	
The Mooring object contains subcomponent objects that represent each component of the full mooring line. Line segments are Section objects, while connectors between segments and at the ends of the lines are Connector objects. These segments alternate.

## Mooring Properties
- dd
: design description dictionary
- n_sec
: number of sections
- i_con
: indices of connectors in the subcomponents list
- i_sec
: indices of sections in the subcomponents list
- rad_anch
: anchoring radius
- rad_fair
: fairlead radius
- z_anch
: anchoring depth
- z_fair
: fairlead depth
- rA : end A absolute coordinates
- rB : end B absolute coordinates
- heading : compass heading from B to A
- adjuster : custom function that can adjust mooring
- shared : int for anchored line (0), shared line (1) or half of a shared line (2)
- symmetric : boolean for if the mooring line is symmetric shared line
- rho : water density
- g : acceleration due to gravity
- envelopes : 2D motion envelopes, buffers, etc.
- loads : dictionary of loads on the mooring line
- reliability : dictionary of reliability information on the line
- cost : dictionary of line costs
- failure_probability : dictionary of failure probabilities


## Mooring methods

### update
Update the Mooring object based on the current state of the design dictionary, or, a new design dictionary.

### setSectionLength
Sets length of the section, including in the subsystem if there is one

### setSectionType
Sets lineType of section, including in the subsystem if there is one

### reposition
Adjusts mooring position based on changed platform location or heading, It can cll a custom "adjuster" function if one is provided. Otherwise it iwll jsut update the end positions.

### setEndPosition
Set the position of an end of the mooring

### getCost
Finds the cost based on the MoorPy subsystem cost estimates

### updateTensions
Gets tensions from subsystem and updates the max tensions dictionary if it is larger than a previous tension

### createSubsystem

Create a MoorPy subsystem for a line configuration from the design dictionary. Regular or suspended lines 
may be used.

### addMarineGrowth
Re-creates sections part of design dictionary to account for marine growth on the subsystem, then calls createSubsystem() to recreate the line

### addCorrosion
Calculates MBL of chain line with corrosion included

### getEnvelope
Computes the motion envelope of the Mooring based on the watch 
circle(s) of what it's attached to. If those aren't already 
calculated, this method will call the relevant getWatchCircle method

## The Connector Class

The Connector class provides a data structure for design information of a connector. 
The Connector class inherits from dict and Node.
The design dictionary includes
the following details: 
- connector type
	- mass
	- volume 
	- CdA
	
The connector class also contains an xyz location of the connector, and a connector object in MoorPy.

## Connector methods

### makeMoorPyConnector

Create a MoorPy connector object in a MoorPy system. Mass, volume, and CdA are added as available in the design dictionary.

## The Section Class

The Section class provides a data structure for the mooring line section material and length. The Section class inherits from dict and Edge.

The line material properties (linear mass, material, MBL, Cd, etc) are stored in the type dictionary of the Section class.