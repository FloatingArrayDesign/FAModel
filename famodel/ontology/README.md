# Floating Array Ontology

This subpackage of FAModel contains information about the floating array 
ontology--a way of recording information that describes a floating wind 
farm project, including both site condition information and design information. 
This ontology is in a draft form and will continue to be revised based
on feedback from prospective users and collaborating projects.

The goal of the ontology is to provide a standardized format for recording
and exchanging a description of a floating wind farm design. This capability
is aligned with the work of IEA Wind Task 49, which focuses on integrated 
design of floating wind arrays. The ontology proposed here draws on elements
from two established ontologies developed under a previous IEA Wind Task.
Task 37 developed [plant-level and turbine-level ontologies](https://windio.readthedocs.io).
The current floating array ontology has a number of additions and differences that 
better suit the scope and emphasis of floating wind arrays. The sections are as follows:

* [Site](#site)
  * [General                         ](#general)
  * [Boundaries                      ](#boundaries)
  * [Exclusions                      ](#exclusions)
  * [Bathymetry                      ](#bathymetry)
  * [Seabed                          ](#seabed)
  * [Metocean                        ](#metocean)
  * [Resource                        ](#resource)
  * [RAFT Cases                      ](#raft-cases)
  * [RAFT Settings                   ](#raft-settings)
* [Array                             ](#array)
  * [Array Layout                    ](#array-layout)
  * [Array Mooring                   ](#array-mooring)
  * [Array Cables                    ](#array-cables)
* [Substation                        ](#substation)
* [Turbine(s)                        ](#turbines)
* [Platform(s)                       ](#platforms)
* [Mooring                           ](#mooring)
  * [Mooring Systems                 ](#mooring-systems)
  * [Mooring line configurations     ](#mooring-line-configurations)
  * [Mooring line section properties ](#mooring-line-section-properties)
  * [Mooring Connectors              ](#mooring-connectors)
  * [Anchor types                    ](#anchor-types)
* [Cables                            ](#cables)
  * [Top Level Cables                ](#top-level-cables)
  * [Cable Configurations            ](#cable-configurations)
  * [Cable Cross Sectional Properties](#cable-types)
  * [Cable Appendages                ](#cable-appendages)
  * [Cable Joints                    ](#cable-joints)

The following sections give an overview of the array ontology makeup with examples. 


## Site
The site section contains information on the site conditions, including the boundaries of the array and exclusion zones. It also contains seabed conditions,
metocean data, and wind resource information for the selected site. 

### General
The general section includes water depth, density of water, density of air, and viscosity of air.

```yaml
    general:
        water_depth : 200        # [m]      uniform water depth
        rho_water   : 1025.0     # [kg/m^3] water density
        rho_air     : 1.225      # [kg/m^3] air density
        mu_air      : 1.81e-05   #          air dynamic 
```

### Boundaries
The boundaries section contains the boundaries of the array. This information is provided with a list of polygon vertices in order, which are then connected linearly 
to define the boundaries of the array. This information can be used to check that all the floating wind turbines are fully contained within the array boundaries.

```yaml
    boundaries:  # project or lease area boundary, via file or vertex list
        file:   # filename of x-y vertex coordinates [m]
        x_y:    # list of polygon vertices in order [m]
           -[x1, y1]
           -[x2, y2]
           -[x3, y3]
           -[x4, y4]
```

### Exclusions
The exclusions section contains any information on exclusion zones, which can be used to check that the power cables and anchors are avoiding
those areas. For simplicity, the exclusion zones allow for two types: circle and polygon. The circle option allows the user to input the x,y coordinates
of the center of the circle and the radius of the circle. The polygon option allows the user to input a list of x,y coordinates which are then 
connected linearly to define the exclusion zone. The user can define as many exclusion zones as needed with a "name" to distinguish them. 

```yaml
    exclusions:
      - name:
        type: circle
        x_y_r:
          - [x1, y1, r1] 
      
      - name:
        type: polygon
        x_y_r:
            -[x1, y1]
            -[x2, y2]
```

### Bathymetry
The bathymetry section currently just provides a link to a MoorDyn-style
bathymetry grid file. There is potential redundancy with the seabed section
and we likely want to support other file formats as well. To be improved.
	
```yaml
    bathymetry:
        file:   # MoorDyn-style bathymetry file
```


### Seabed
The seabed section contains information on the depth and soil type throughout the array. 
The user provides a list of x, y, z and "soil type" points to define the depth
and soil type classification. Alternatively, a file can be specified containing fully
gridded data about depth and soil properties on a rectangular grid, with the option
for quantitative soil index properties.
In between provided points, the depth is linearly interpolated. 
For now, the soil type will be assumed to match that of the nearest point. 
	
```yaml
    seabed:
        keys: [x,  y, depth,  soil_type]
        data:
          - [ x1, y1, z1, "soft clay"]
          - [ x2, y2, z1, "medium clay"]    
          
        filename: 'seabed_grid.csv'  # gridded data of seabed depth and soil classification   
```

### Metocean
The metocean section contains wind, wave, and current parameters that inform the design load cases. This section is further subdivided into extremes, probabilistics bins, and time series.
The extreme subsection contains metocean parameters for specified return periods which are used to check the strength of the design in extreme conditions.
The probabilistic bins section contains a list of metocean parameters, including directionality and probability. This section is used to inform a fatigue analysis,
so the probability of the bins should total to 1. Finally, the timeseries section contains the metocean parameters in time series form, which 
is needed for logistics analysis to inform vessel availability. To reduce the number of lines in the file, the timeseries section inputs a 
csv filename. 

```yaml
    metocean:
        extremes:  # extreme values for specified return periods (in years)
            return periods :   [ 1, 5, 50 ]  # in years
            all :   # unconditional extreme values
                Hs:    [     ,     ,     ]
                Tp:    [     ,     ,     ]
                Gamma: [     ,     ,     ]
                WS:    [     ,     ,     ]
                TI:    [     ,     ,     ]
                Shear: [     ,     ,     ]
                CS:    [     ,     ,     ]
                8.5	9.8, 10.4, 11.8, 12.4, 13.7, 16.8, 18.1, 18.6, 19.8, 20.3, 21.4

                
            WS_2_4 :  # conditional values from wind speed range of 2 - 4 m/s
                Hs:    [     ,     ,     ]
                Tp:    [     ,     ,     ]
                Gamma: [     ,     ,     ]
                TI:    [     ,     ,     ]
                Shear: [     ,     ,     ]
                CS:    [     ,     ,     ]
            
            WS_4_6 :  # conditional values from wind speed range of 4 - 6 m/s
                ...
            
            
        joint_probabality_bins:  # a set of cases representing the joint metocean probability distribution
            keys : [ prob , Hs  , Tp, WindSpeed, TI, Shear, Gamma, CurrentSpeed, WindDir, WaveDir, CurrentDir  ]
            data :
                -  [ 0.010  ,   ,   ]
                -  [ 0.006  ,   ,   ]
                -  [ 0.005  ,   ,   ]
                
        time_series :
            filename: 'metocean_timeseries.csv'
```

### Resource
The resource section contains information on the wind resource for the site, 
which can be used to calculate AEP. To again reduce the number of lines
in the file, the resource section inputs a filename which contains the resource data.
This resource data will follow the [WindIO plant ontology](https://windio.readthedocs.io/en/latest/source/plant.html).

```yaml
    resource :
        filename: 'windresource'
```

### RAFT Cases
The RAFT cases section contains the parameters for any load cases that are intended to be run in [RAFT](https://openraft.readthedocs.io), and as such follows the format specified by RAFT. 
The section inputs a list where each entry corresponds to a load case. Note that turbulence can be input 
as a percent or as a string representing a turbulence model such as IIB_NTM.
```yaml
RAFT_cases:
        keys : [wind_speed, wind_heading, turbulence, turbine_status, yaw_misalign, wave_spectrum, wave_period, wave_height, wave_heading  ]
        data :  #   m/s        deg    % or e.g. IIB_NTM    string            deg         string          (s)         (m)         (deg)
            -  [    10.5,         0,            0.01,    operating,          0,        JONSWAP,         12,         6,         0       ]
```

### RAFT Settings
The RAFT settings section contains the general parameters used for RAFT simulations, such as cutoff frequencies, 
Initial amplitudes for each degree of freedom at all frequencies, and the number of iterations to solve the model dynamics. 
As with the previous section, the format follows that specified by [RAFT](https://openraft.readthedocs.io).

```yaml
RAFT_settings:   
        min_freq     :  0.001    #  [Hz]       lowest frequency to consider, also the frequency bin width     
        max_freq     :  0.10    #  [Hz]       highest frequency to consider
        XiStart      :   0      # sets initial amplitude of each DOF for all frequencies
        nIter        :   4      # sets how many iterations to perform in Model.solveDynamics()
``` 

## Array

This part of the ontology includes a section for the tubine layout, as
well as optional array-level descriptions of the mooring system and 
array cabling.

### Array Layout
The array section summarizes the floating wind turbines in the array. The 
section inputs a list where each entry corresponds to a wind turbine. The ID serves as a method to identify the specific turbine system. 
As such, each list entry should have a unique ID, but the ID type (string, int, etc) is up to the user. The turbineID and platformID are specified for each list entry,
connecting to details in the [Turbine](#turbines) and [Platform](#platforms) sections. This allows the user to easily 
specify different turbine or platform types throughout the array. 
Similarly, the mooringID is included and refers to the [mooring_systems](#mooring-systems) section.
This allows the user to set up general mooring systems to be used throughout the array. Additionally, the x and y locations are input and the heading adjustment.
The heading adjustment refers to a rotation of the mooring system, relative to how it is defined in the mooring_systems section. This allows the user to 
easily define a single mooring system for various rotations throughout the array.

Alternatively, the mooringID can be set to zero and the mooring system can be 
input in the [array_mooring](#array-mooring) section.

```yaml
array:
    keys : [ID, turbineID, platformID, mooringID,   x_location,     y_location,   heading_adjust]
    data : # ID#   ID#        ID#        ID#           [m]             [m]           [deg]
        -  [fowt1,  1,         1,         ms1,         0,             0,           180  ]    
        -  [f2,     1,         2,         ms2,         1600,          0,            0   ]  
```


### Array Mooring
The array mooring section allows the user to input array-level mooring system details, instead of the more generalized mooring systems in mooring_systems.
This section inputs a list of x,y anchor positions, anchor type, and embedment depth. The anchor type links to the list in the anchor_types section.
Additionally, a list of mooring lines can be input with specified attachments at FOWTs and anchors. If there is an anchor connected to this line, it must be listed 
in end A, not end B. All anchors listed in line_data end A must have a matching ID in the anchor_data table, and all FOWTs listed in line_data end A or end B 
must have a matching ID in the array_data table. The anchor and fowt IDs must all be unique. The mooring lines each have a mooring configuration ID which links to the mooring_line_configs section. 
There is also an option to adjust the length of the line, depending on the spacing. 

```yaml
array_mooring:
    anchor_keys : 
          [ID, type,  x,  y,  embedment ]
    anchor_data :
        - [  anch1,  suction1,   ,   ,     ]
        - [  anch2,  suction1,   ,   ,     ]
    
    line_keys : 
          [MooringConfigID  ,  end A,   end B,  lengthAdjust]
    line_data :
        - [ semitaut-poly_1      ,  anch1,  fowt1,   0]
        - [ semitaut-poly_1      ,  anch2,  fowt1,   0]
        - [ semitaut-poly_2      ,  fowt1,  f2,   0]
```

### Array Cables

This section provides a straightforward and compact way to define the power
cables in the array. For each end (A and B) of the cable, it specifies the
turbine (matching and ID in the array table) or substation (matching an ID in the substation section)
attached to, the [Top Level Cables](#top-level-cables) 
used, the heading of the cable at the attachment of each end, and length adjustment. A 0-degree heading aligns with 
North and rotates clockwise.
Additional detail related to subsections of the top level cable, joints, and any cable 
routing are in the [Cables](#cables) section. 
The CableID refers to an entry in the [Top Level Cables](#top-level-cables) section.

```yaml
array_cables:   
    keys:  [ CableID,       AttachA,   AttachB,      headingA, headingB, lengthAdjust]
    data:
      - [ array_cable1,     fowt1,     f2,        270,      270,      0] 
      - [ suspended_cable1, f2,     substation1,  90,        270,      0]  
```

## Substation(s)

The substation section defines the substations used in the array. The substation key (name) must be unique 
from the ID keys in the array layout table.

## Turbine(s)

The turbine section can contain either a single turbine or a list of turbines, 
depending on the scenario. Note that if multiple turbines are listed, the section title must be 'turbines' instead of 'turbine'.
By default, the format follows that of [RAFT](https://openraft.readthedocs.io)
However, support will be added for linking to turbine design descriptions that follow
the [WindIO](https://windio.readthedocs.io) ontology format, which is also used 
by [WEIS](https://weis.readthedocs.io).

## Platform(s)

This section defines the floating support structures used in the design. As in
the previous section, it can contain a single platform or a list of platforms. 
By default, the format here follows that used by 
[RAFT](https://openraft.readthedocs.io) input files, with the addition of 'rFair' and 'zFair' entries to the 
dictionary for each platform in the first level of each platform listed. 
However, support will be added for also linking to platform descriptions that follow
the [WindIO](https://windio.readthedocs.io) ontology format, which is also used 
by [WEIS](https://weis.readthedocs.io).

```yaml 
platform:

    potModMaster :   1      # [int] master switch for potMod variables; 0=keeps all member potMod vars the same, 1=turns all potMod vars to False (no HAMS), 2=turns all potMod vars to True (no strip)
    dlsMax       :  5.0     # maximum node splitting section amount for platform members; can't be 0
    rFair        :  58 
    zFair        :  -15
    
    members:   # list all members here
        
      - name      :  center_column             # [-]    an identifier (no longer has to be number)       
        type      :  2                         # [-]    
        rA        :  [ 0, 0, -20]              # [m]    end A coordinates
        rB        :  [ 0, 0,  15]              # [m]    and B coordinates
        shape     :  circ                      # [-]    circular or rectangular
        gamma     :  0.0                       # [deg]  twist angle about the member's z-axis
        potMod    :  True                      # [bool] Whether to model the member with potential flow (BEM model) plus viscous drag or purely strip theory
        # --- outer shell including hydro---
        stations  :  [0, 1]                    # [-]    location of stations along axis. Will be normalized such that start value maps to rA and end value to rB
        d         :  10.0                      # [m]    diameters if circular or side lengths if rectangular (can be pairs)
        t         :  0.05                      # [m]    wall thicknesses (scalar or list of same length as stations)
        Cd        :  0.6                       # [-]    transverse drag coefficient       (optional, scalar or list of same length as stations)
        Ca        :  0.93                      # [-]    transverse added mass coefficient (optional, scalar or list of same length as stations)
        CdEnd     :  0.6                       # [-]    end axial drag coefficient        (optional, scalar or list of same length as stations)
        CaEnd     :  1.0                       # [-]    end axial added mass coefficient  (optional, scalar or list of same length as stations)
        rho_shell :  7850                      # [kg/m3] 
        # --- handling of end caps or any internal structures if we need them ---
        cap_stations :  [ 0    ]               # [m]  location along member of any inner structures (in same scaling as set by 'stations')
        cap_t        :  [ 0.001  ]             # [m]  thickness of any internal structures
        cap_d_in     :  [ 0    ]               # [m]  inner diameter of internal structures (0 for full cap/bulkhead, >0 for a ring shape)

```


## Mooring

This part of the ontology contains a set of sections that define parts of the 
mooring system, down to the definition of line cross sectional properties and 
anchor characteristics.

### Mooring Systems

This section describes the mooring systems that could be used for individual turbines and repeated throughout the array. Each mooring system contains a 
list of mooring lines, which contains the mooring configuration ID, the heading, the anchor type, and a possible length adjustment. The 
mooring configuration ID links to the details about the segments lengths and types in the mooring line configurations section. The heading refers to the angle of the mooring line and it rotates 
clockwise from North. The anchor type links to details about the anchor 
size and dimensions in the [anchor types section](#anchor-types). The length adjustment
is an optional parameter that can adjust the mooring line length for a shallower or deeper depth, for example. 

```yaml
mooring_systems:
    ms1:
        name: 3-line taut polyester mooring system
        
        keys: [MooringConfigID,  heading, anchorType, lengthAdjust] 
        data:
          - [  semitaut-poly_1,   30 ,    suction 1,   0 ]
          - [  semitaut-poly_1,  150 ,    suction 1,   0 ]
          - [  semitaut-poly_1,  270 ,    suction 1,   0 ]
```

### Mooring Line Configurations

The mooring line configurations lists the segment lengths and line types that make up each mooring line. Each line has a name that can then be specified 
as the MooringConfigID in the mooring systems section. The span is specified for each configuration, which represents the distance in the x-y plane between 
the two connection points of the line - i.e. between fairlead and anchor, or for shared lines, fairlead and fairlead. Fairlead radius and fairlead depth are specified in the [Platform](#platforms) section.
 Each line contains a list of sections that details the line section type and length. The line type name
connects to information in the mooring [line section properties](#mooring-line-section-properties) if the keyword 'type' is used. If 'mooringFamily' is instead 
specified (as in the catenary_1 example below), the mooring line properties are determined from MoorPy MooringProps values. In the latter case, the nominal diameter in m 'd_nom' must also be provided for that line section. 
Additionally, before and after each line section has an optional input which can list the 
ID of a [connector type](#mooring-connectors), such as an H-link or buoy. 
This information allows the weight and buoyancy of the connections to be included 
in the model, and provides clarity on the location of the connector relative to different line sections. 
There is also a True/False options for whether the section length is adjustable. 

Shared or suspended lines may also have an optional 'symmetric' input which, if set to true, signifies that 
the line is symmetric and only the first half of the line is provided in the 'sections' list. When loaded in to the project class, the mooring object will automatically be fully filled out by mirroring 
the first half of the line. If a connector is provided as the last item in the sections list for a symmetric line, it is assumed that the middle line is two identical lines with the given connector between, otherwise 
the middle line (last line given in the list) is doubled in length in the mirroring process. For example, the 'rope_shared' config in the yaml below would produce a symmetric shared line with sections in the following order
a 150 m section of rope, a clump weight, a 1172 m section of rope (note the doubled length), a clump weight, and finally a 150 m section of rope.



```yaml
  mooring_line_configs:
    
    catenary_1:
        name: catenary configuration 1
        
        span: 779.6 
        
        sections:
          - mooringFamily: chain
            d_nom: 0.185 # m
            length: 850
            adjustable: True
			
	semitaut-poly_1:  # mooring line configuration identifier
    
        name: Semitaut polyester configuration 1  # descriptive name
        
        span: 642
        
        sections:                 #in order from anchor to fairlead
          - type: chain_155mm       # ID of a mooring line section type
            length: 497.7            # [m] usntretched length of line section
            adjustable: True      # flags that this section could be adjusted to accommodate different spacings...
          - connectorType: h_link   
          - type: polyester_182mm        # ID of a mooring line section type
            length: 199.8           # [m] length (unstretched)


    rope_shared:
        name: shared rope 
		symmetric: True
        
        span: 1484

        
        sections:
          - type: rope
            length: 150
          - connectorType: clump_weight_80
          - type: rope
            length: 586

```    
    
### Mooring line section properties

The mooring line section properties contains the properties needed to accurately model the mooring lines. Each mooring line type is listed with 
a name that can be referenced in the [mooring line configurations section](#mooring-line-configurations). 
For each line type, the nominal and volume equivalent diameter are listed, 
as well as the mass density, stiffness, cost, MBL, and material name. The ontology supports either a single stiffness value (like for chain)
or the static-dynamic stiffness of fiber lines. An example of this is shown below. 

Alternatively, the mooring line parameters can be provided in a table-based format to reduce the number of lines.

```yaml
mooring_line_types:

    polyester_226mm:
        d_nom:    0.262      # [m] nominal diameter
        d_vol:    0.2258     # [m] volume-equivalent diameter
        m:        55.0       # [kg/m] mass per unit length (linear density)
        EA:       164e6      # [N] quasi-static stiffness
        MBL:    11.75e6      # [N] minimum breaking load
        EAd:    164.6e6      # [N] dynamic stiffness
        EAd_Lm:    0.34      # [-] dynamic stiffness mean-load multiplier
        cost:      194       # [$/m] cost per unit length
        material: polyester  # [-] material composition descriptor
        
    chain_170mm::
        d_nom:    0.170      # [m] nominal diameter
        d_vol:    0.306      # [m] volume-equivalent diameter
        m:        575.0      # [kg/m] mass per unit length (linear density)
        EA:      2468e6      # [N] quasi-static stiffness
        MBL:     25.2e6      # [N] minimum breaking load
        cost:      1486      # [$/m] cost per unit length
        material:  R3 studless chain  # [-] material composition descriptor

    # alternative table-based format
    keys :  [ name,   EA ,  MBL,  ...]
    data :
        -   [ poly1  , 3e7, 10e6, ... ]
        -   [ chain27, 3e9, 10e6, ... ]
```

### Mooring Connectors

This section lists properties of the mooring connections that are referenced in the mooring line configurations section. 
Each connector has a name that is used to identify it, as well as a mass (m) and volume (v). Optionally, the CdA of the connector 
can be specified to model the drag on the component. 
```yaml
 mooring_connector_types:
    
    h_link:
        m : 140    # [kg]  component mass
        v : 0.13   # [m^3] component volumetric displacement
        
    clump_weight_20:
        m : 20000  # [kg]
        v :  0.8   # [m^3]
        
    buoy_10:
        m  :  560   # [kg]  component mass
        v  : 10.2   # [m^3] component volumetric displacement
        CdA :  3.5   # [m^2] product of cross sectional area and drag coefficient
```

## Anchor types

The anchor types section lists dimensions and embedment depth for each anchor type. The anchor types section
allows the user to input the diameter, length, area, thickness, and embedment depth. All parameters are optional,
because the applicable information depends on the anchor type. 
The parameters align with the FAModel 
[intermediate anchor model](https://github.com/FloatingArrayDesign/FAModel/tree/main/famodel/anchors#parameters-needed-for-level-2-anchor-capacity-models). 

```yaml        
anchor_types:
    suction1:
        name   : standard suction pile
        d      :    # [m] Diameter
        L      :    # [m] Length
        t      :    # [mm] Thickness
        h      :    # [m] Embedment depth
```


## Cables

This section describes the cables through the array including both static 
and dynamic portions of cables. At the top level, each array cable going
between a pair of turbines (or a turbine and a substation) is defined. Each 
subsection of cable is then defined in the [Cable Configurations](#cable-configs) 
section, including buoyancy module layout. 


### Top Level Cables

A top-level cableis defined as the full assembly of electrical connection equipment between
two turbines or a turbine and a substation. 'type' links to the cable configuration description 
of the subsections of the cable, which can be dynamic or static cables. The 'connectorType' refers to 
an entry in the [Cable Joints](#cable-joints) section. The first entry in the sections list is connected 
to end A, while the last entry is connected to end B.

```yaml
 cables:

    array_cable1:      
        name: array cable with lazy wave - static - lazy wave sections # descriptive cable name
        # type : static_cable_80   # cable section type ID
   
        sections: # refers to a configuration in cable_configs
        # first entry is at end A, last entry is at end B
          - type: dynamic_lazy_wave1
     
          - connectorType: joint_1
            
          - type: static_1
     
          - connectorType: joint_1
     
          - type: dynamic_lazy_wave1

      
    suspended_cable1:
        sections: 
          - type: dynamic_suspended_1
```

### Cable Configurations

This section lists the cable configurations used in the array design.
The 'type' listed in the entry is either 'static' or 'dynamic'.
The 'cableFamily' key is used when the cable cross-sectional information 
will be imported from the cableProps_default yaml (in which case, an area A 
must be provided in mm^2), and the value for cableFamily must match an entry in the cableProps_default yaml.
Alternatively, if the cable cross-sectional properties will be provided in the [Cable Cross Sectional Properties](#cable-cross-sectional-properties)
section, the key 'typeID' will be used in place of 'cableFamily', and will
refer to an entry in the Cable Cross Sectional Properties list. 

The sections list provides details on the layout of buoyancy modules and clump weights, including the 
distance of the buoyancy section midpoint from end A, the number of modules, the spacing between modules, 
and the volume. The volume is only needed if the buoyancy module properties will be imported 
from the cableProps_defaul yaml. As with the cable properties, the 'type' in the sections list must refer to 
an entry in either the [Cable Appendages](#cable-appendages) section or in the cableProps_default.yaml.

Static cables can have routing information listed as vertex points along the cable route, and the radius of curve. 
Static cable burial information can also be provided.

Similar to mooring lines, the span refers to the end to end distance of the line in the x-y plane.

```yaml
    basic1:
        name: basic cable configuration, essentially a straight line between platforms
        span: 1600 # [m]
        type: dynamic
        zJTube: -30 # depth out of J-tube that cable starts in m
  
        cableFamily: dynamic_cable_33
        length: 1700 # [m]
        A: 100 # cable conductor cross-sectional area [mm^2] (Required for types listed in cable props yaml)
          
          
          
        
    dynamic_lazy_wave1:
        name: Lazy wave configuration 1 (simpler approach)
        voltage: 66 # [kV]
        span : 600    # [m] horizontal distance to end of dynamic cable
        zJTube: -20 # depth the cable comes out of the j-tube
        type: dynamic # dynamic or static
        A: 300
        

        cableFamily: dynamic_cable_33        # ID of a cable section type1
        length: 900                   # [m] length (unstretched) 
       
        sections: 
          - type: Buoyancy_750m #_w_buoy # (section properties including averaged effect of buoyancy modules)
            L_mid: 450 # [m] from end A
            N_modules: 5 # number of modules in this buoyancy section
            spacing: 21 # [m]
            V: 2   # [m^2]         
       
    
    static_1:
        name: Static cable configuration 1 
        voltage: 66 # [kV]
        span: 350
        type: static
        
        typeID: static_cable_36
        length: 2200
            
        routing_x_y_r:  # optional vertex points along the cable route. Nonzero radius wraps around a point at that radius.
          - [1000, 1200, 20] 
          - [2000, 1500, 20] 
        
        burial:  # optional definition of cable burial depth over its length
            station: [0, 1]                # length along cable, normalized by first and last value
            depth  : [0.1, 0.2]            # [m] burial depth
```       
    
### Cable Cross Sectional Properties

This section details the cross-sectional properties of each cable type.

```yaml
cable_types:

    dynamic_cable_66 :     # cable type identifier        
        dynamic :   True   # Flag for dynamic cable (default static)
        DC   :     False   # Flag for DC (default AC)
        kV   :        66   # [kV] voltage rating
        A    :       300   # [mm^2] cross-sectional area of each conductor (3 conductors)
        D    :      0.20   # [m] outer diameter
        m    :     30.59   # [kg/m] mass per unit length
        EA   :    700e+3   # [kN] axial stiffness 
        EI   :      10.0   # [kN.m^2] bending stiffness
        MBL  :       100   # [kN] minimum breaking load
        MBR  :       2.0   # [m] minimum bending radius

    static_cable_36:  
        dynamic :  False
        DC   :     False	
        kV   :        36		
        A    :       300		
        D    :    0.1248		
        m    :     12.90		
        EA   :    245e+3		
        EI   :      5.10		
        MBL  :      54.0		
        MBR  :     1.875		
```

### Cable Appendages

This section lists any cable appendages that might be added to the cables,
such as buoyancy modules or cable protection system components. Each entry
is given an identifier and can have a variety of parameters that describe
its lumped properties, such as mass, volume, and drag coefficient-area
product. These appendages are used in the 
[cable_configs](#cable-configurations) section.

```yaml
  cable_appendages:

    buoyancy_module_1:
        mass:    2700   # [kg]  mass
        volume: 8.615   # [m^3] volumetric displacement 
        CdA:     3.8    # [m^2] product of cross-sectional area and drag coefficient
        length:  2.2    # [m]   length taked up along cable
```

### Cable Joints

This section lists any cable joints that might connect cable subsections. Each entry is given 
and identifier and parameters to describe the joint. These joints are used in the [Top Level Cables] 
(#top-level-cables) section.
```yaml
cable_joints:
    joint_1:
        mass: 100000 # TBD
```
