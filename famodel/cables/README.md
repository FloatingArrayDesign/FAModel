# The Cable Class

The [Cable class](./cable.py) contains properties, methods, and subcomponents aimed at describing and modeling a power cable between two Nodes (i.e., between two platforms or a platform and a substation). The Cable class inherits from Edge.

There are three classes that can make up the subcomponents of a Cable object:
- [Dynamic Cable](#dynamic-cable-class)
- [Static Cable](#static-cable-class)
- [Joint](#joint-class)

The dynamic portion of the cable is a separate subcomponent from the static portion to differentiate these sections in modeling marine growth, cable motions, routing, and more.

These subcomponents alternate between Edges (either a DynamicCable or StaticCable object) and Nodes (Joints). A typical cable may have dynamic cables connected to the platforms/substations on either end, and a static cable is connected to the dynamic sections with a joint on the seabed. This would be characterized with the following subcomponents list: 
[DynamicCable, Joint, StaticCable, Joint, DynamicCable]. However, the sections of a cable object are customizeable as needed. They are defined by the user in the cables section of the Ontology file.

 A suspended cable may only have one subcomponent (a DynamicCable object) as there is not static portion of the cable, and there may not be a need for a joint.

### Cable Properties

- dd
: design description dictionary. This is composed of the following sections:
    - Joints : joint objects  in this cable
    - Cables : static and dynamic cable objects in this cable
    - routing : any cable routing information
- n_sec
: number of subcomponent cables (DynamicCable or StaticCable objects) the Cable is composed of
- upstream_turb_count
: number of turbines upstream from this cable
- i_con
: list of Joint indices in the subcomponents list
- i_sec
: list of StaticCable and DynamicCable indices in the subcomponents list
- rA 
: xyz location of end A of cable
- rB
: xyz location of end B of the cable
- x
: cable routing locations in x-direction
- y
: cable routing locations in y-direction
- r
: cable routing radii
- L
: total length of the cable
- failure_probability
: dictionary of failure probabilities for the cable
- cost
: total cost of the cable components (USD)




### Cable Methods

[Back to Top](#the-cable-class)
## Dynamic Cable Class

The [DynamicCable class](./dynamic_cable.py) contains properties and methods for a dynamic section of a power cable. The DynamicCable class inherits from Edge.

The DynamicCable class includes information on the bare cable properties as well as any buoyancy module sections.

### Dynamic Cable Properties

- dd
: design description dictionary. This is composed of the following sections:
    - cable_type : the bare cable properties
    - buoyancy_sections : list of buoyant sections of the cable including
        - N_modules: number of modules
        - spacing: spacing between module centerpoints
        - L_mid: location along the cable of the midpoint of the buoyancy section
        - module_props: buoyancy module properties
- ss
: pristine MoorPy subsystem associated with this dynamic cable
- ss_mod
: modified MoorPy subsystem associated with this dynamic cable. 
Often this is used for adding marine growth
- rA
: xyz location for end A of the dynamic cable
- rB
: xyz location for end B of the dynamic cable
- headingA
: heading from the Node connected to end A of the dynamic cable (or overall Cable object)
- headingB
: heading from the Node connected to end B of the dynamic cable (or overall Cable object) []
- L
: total length of the dynamic cable [m]
- rad_fair
: radius of fairlead connection [m]
- rad_anch
: distance from center of connected platform or substation to center of node connection at other end (either a Joint or another platform/substation). There are no anchors for a cable, this property name is solely for compatibility with certain functions and MoorPy [m]
- z_anch
: depth of end connection point. There are no anchors for a cable, this property name is solely for compatibility with certain functions and MoorPy [m]
- z_fair
: depth of fairlead connection point [m]
- adjuster
: custom function to adjust the cable for different depths
- shared
: int for if a regular (0), suspended (1), or symmetric half of suspended (2) cable
- rho
: density of fluid this cable is submerged in [kg/m^3]
- g
: acceleration due to gravity [m/s^2]
- loads
: dictionary of loads on the dynamic cable
- reliability
: dictionary of reliability information for the dynamic cable
- cost
: dictionary of costs for the dynamic cable
- failure_probability
: dictionary of failure probabilities for the dynamic cable



### Dynamic Cable Methods

 - makeCableType()
 : Processes dictionary info to make cableType dictionary
 - updateSubsystem()
 : Adjusts the subsystem properties when the buoyancy section info changes.
 The contents of self.dd['buoyancy_sections'] should already be changed before
 calling this method. The length (L) could also change
 - getBuoyancyOverlaps()
 : Calculates the margin for any buoyancy sections overlapping with
 each other or with the cable ends. The returned values could be used
 as constraints to avoid overlap in an optimization
 - calcEquivBuoyancy()
 : Calculates new cable information that includes equivalent buoyancy of buoys
 - reposition()
 : Adjusts dynamic cable position based on changed platform location or
 heading. It can call a custom "adjuster" function if one is
 provided. Otherwise it will just update the end positions
 - setEndPosition()
 : Sets the position of an end of the dynamic cable.
 - getCost()
 : Gets tensions from subsystem and updates the max tensions dictionary if it is larger than a previous tension
 - createSubsystem()
 : Creates a subsystem for cable and buoyancy section(s) configuration from the design dictionary
 - addMarineGrowth()
 : Creates a new design dictionary (adds a section to old one) to account for marine growth on the subystem, then calls createSubsystem() to recreate the cable. Marine growth on a cable buoyancy section is handled differently than that of a bare cable section. See [Marine Growth on DynamicCables](#marine-growth-on-dynamiccables) for details on the modeling of marine growth on cables in FAModel.
 - symmetricalMirror()
 : Mirrors a half symmetrical cable design dictionary to show the entire cable, creates subsystem if asked.

### Marine Growth on DynamicCables
Marine growth is modeled on bare dynamic cables similarly to the way it is modeled on mooring lines - using thicknesses and densities defined in a dictionary for different depth ranges. 

For buoyant sections of dynamic cables, the marine growth on the buoyancy modules significantly alters the method of calculating marine growth mass and weight. A thickness of marine growth is added to the entire surface area of the buoy. In this instance, buoys are assumed cylinders with a hole through the axis that is the same diameter as the bare cable. See the figure below:


The mass of marine growth is determined for each buoy, as well as for the bare cable length between buoys. The mass and weight of the marine growth for entire buoyant section length is then determined, and added to the pristine cable and buoy mass and weight. The calculations below describe this process in detail.

The volume of marine growth on the buoy $V_{mg_{buoy}}$ can be described by the following equation:

```math
V_{mg_{buoy}} = \frac{\pi}{4}((d_{buoy}+2t)^2-d_{buoy}^2)L_{buoy} + \frac{\pi}{2}((d_{buoy}+2t)^2-d_{cable}^2)t

```
Where $d_{buoy}$ is the diameter of the buoy cylinder, $t$ is the thickness of marine growth, $L_{buoy}$ is the length of the buoy and $d_{cable}$ is the diameter of the bare cable.

With the volume of marine growth, the mass of marine growth can easily be determined using the density of marine growth ($\rho_{mg}$)
```math
m_{mg_{buoy}} = \rho_{mg}V_{mg_{buoy}}

```
The weight of the marine growth on the buoy is calculated as follows:
```math
w_{mg_{buoy}} = (\rho_{mg}-\rho_{water})gV_{mg_{buoy}}
```
The volume of marine growth on the bare cable is calculated as :
```math
V_{mg_{cable}} = \frac{\pi}{4}((d_{cable}+2t)^2-d_{cable}^2)(L_{bare}-2t)*(N_{modules}-1)
```
Where 

[Back to Top](#the-cable-class)


## Static Cable Class

The StaticCable class contains properties and methods for a static section of a power cable. The StaticCable class inherits from Edge.

### Static Cable Properties

### Static Cable Methods

[Back to Top](#the-cable-class)
## Joint Class

[Back to Top](#the-cable-class)
