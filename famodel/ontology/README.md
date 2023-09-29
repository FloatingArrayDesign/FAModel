# Floating Array Ontology

This subpackage of FAModel contains information about the floating array ontology--a way of recording information that describes a floating wind farm project, including both site condition information and design information.
The following sections outline the array ontology makeup. 

## Site
The site section contains information on the site conditions, including the boundaries of the array and exclusion zones. It also contains seabed conditions,
metocean data, and wind resource information for the selected site. 

### General
The general section includes water depth, density of water, density of air, and viscosity of air.

```python
    general:
        water_depth : 200        # [m]      uniform water depth
        rho_water   : 1025.0     # [kg/m^3] water density
        rho_air     : 1.225      # [kg/m^3] air density
        mu_air      : 1.81e-05   #          air dynamic 
```

### Boundaries
The boundaries section contains the boundaries of the array. This information is provided with a list of polygon vertices in order, which are then connected linearly 
to define the boundaries of the array. This information can be used to check that all the floating wind turbines are fully contained within the array boundaries.

```python
    boundaries: #list of polygon vertices in order
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

```python
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

### Seabed
The seabed section contains information on the depth and soil type throughout the array. The user provides a list of x,y,z and "soil type" points to define the depth. 
In between provided points, the depth is linearly interpolated. For now, the soil type will be assumed to match that of the nearest point. 
	
```python
	seabed:
        -[x1, y1, z1, "soft clay"]
        -[x2, y2, z2, "medium clay"]    
```

### Metocean
The metocean section contains wind, wave, and current parameters that inform the design load cases. This section is further subdivided into extremes, probabilistics bins, and time series.
The extreme subsection contains metocean parameters for specified return periods which are used to check the strength of the design in extreme conditions.
The probabilistic bins section contains a list of metocean parameters, including directionality and probability. This section is used to inform a fatigue analysis,
so the probability of the bins should total to 1. Finally, the timeseries section contains the metocean parameters in time series form, which 
is needed for logistics analysis to inform vessel availability. To reduce the number of lines in the file, the timeseries section inputs a 
csv filename. 

```python
	    metocean:
        extremes:  # extreme values for specified return periods (in years)
            keys :   [ Hs  , Tp  , WindSpeed, TI, Shear, Gamma, CurrentSpeed ]
            data :
                1:   [     ,     ,     ]
                10:  [     ,     ,     ]
                50:  [     ,     ,     ]
                500: [     ,     ,     ]
                
        probabalistic_bins:
            keys : [ prob , Hs  , Tp, WindSpeed, TI, Shear, Gamma, CurrentSpeed, WindDir, WaveDir, CurrentDir  ]
            data :
                -  [ 0.010  ,   ,   ]
                -  [ 0.006  ,   ,   ]
                -  [ 0.005  ,   ,   ]
                
        time_series :
            filename: 'metocean_timeseries.csv'
```

### Resource
The resource section contains information on the wind resource for the site, which can be used to calculate AEP. To again reduce the number of lines
in the file, the resource section inputs a filename which contains the resource data.

```python
    resource :
        filename: 'windresource'
```

## Array
The array section summarizes the floating wind turbines in the array. The section inputs a list where each entry corresponds to a wind turbine.
The turbineID and platformID are specified for each, connecting to details in the Turbine and Platform sections. This allows the user to easily 
specify different turbine or platform types throughout the array. Similarly, the mooringID is included and refers to the mooring_systems section.
This allows the user to set up general mooring systems to be used throughout the array. Additionally, the x and y locations are input and the heading adjustment.
The heading adjustment refers to a rotation of the mooring system, relative to how it is defined in the mooring_systems section. This allows the user to 
easily define a single mooring system for various rotations throughout the array.

Alternatively, the mooringID can be set to zero and the mooring system can be input in the array_mooring section.

```python
    array:         # [copy from RAFT style for the moment]
    keys : [turbineID, platformID, mooringID,   x_location,     y_location,   heading_adjust]
    data : #    ID#        ID#        ID#          [m]             [m]           [deg]
        -  [     1,         1,         1,             0,              0,          180   ]    # 2 array, shared moorings
        -  [     2,         1,         2,          1600,              0,            0   ]  
```



### Array Mooring
The array mooring section allows the user to input array-level mooring system details, instead of the more generalized mooring systems in mooring_systems.
This section inputs a list of x,y anchor positions, anchor type, and embedment depth. The anchor type links to the list in the anchor_types section.
Additionally, a list of mooring lines can be input with specified attachements at numbered FOWTs and anchors. The mooring lines each have a mooring 
configuration ID which links to the mooring_line_configs section. There is also an option to adjust the length of the line, depending on the spacing. 

```python
    array:         # [copy from RAFT style for the moment]
    keys : [turbineID, platformID, mooringID,   x_location,     y_location,   heading_adjust]
    data : #    ID#        ID#        ID#          [m]             [m]           [deg]
        -  [     1,         1,         1,             0,              0,          180   ]    # 2 array, shared moorings
        -  [     2,         1,         2,          1600,              0,            0   ]  


array_mooring:  # this is an array-level list, in addition to any per-mooring-system ones
    anchor_keys : 
          [ID, type,  x,  y,  embedment ]
    anchor_data :
        - [  1,  suction1,   ,   ,     ]
        - [  2,  suction1,   ,   ,     ]
          
    line_keys : 
          [MooringConfigID  ,  end A,   end B,  lengthAdjust?]
    line_data :
        - [ Taut-Poly_1      ,  FOWT 1,  Anch 1,   0]
        - [ Taut-Poly_1      ,  FOWT 1,  Anch 2,   0]
        - [ Taut-Poly_2      ,  FOWT 1,  Anch 3,   0]

```

### Turbine(s)
The turbine section can contain either a single turbine or a list of turbines, depending on the scenario. The information and layout of this 
section was taken from RAFT. 

###Platform(s)


## Mooring

### Mooring Systems

This section describes the mooring systems that could be for any individual turbine.

```python
mooring_systems:  # this is where individual mooring systems can be listed
    ms1:
        name: a great mooring system
        
        keys: [MooringConfigID,  heading, anchorType, lengthAdjust?] 
        data:
          - [  taut-poly_1,   60 ,    suction 1,   0 ]
          - [  taut-poly_1,  180 ,    suction 1,   0 ]
          - [  taut-poly_1,  300 ,    suction 1,   0 ]
```

### Mooring line configurations

```python
  mooring_line_configs:
    
    taut-poly_1:  # mooring line configuration identifier
    
        name: Taut polyester configuration 1  # descriptive name
        
        sections:
          - type: chain_160       # ID of a mooring line section type
            length: 80            # [m] usntretched length of line section
            connector: h_link     # ID of a connector type at the end of the line section (optional)
            adjustable: True      # flags that this section could be adjusted to accommodate different spacings...
            
          - type: poly_180        # ID of a mooring line section type
            length: 762           # [m] length (unstretched)
            connector: shackle    # ID of a connector type (optional)
            
        attachment:
            type: ?  # fairlead/pivot/padeye/other? (optional)
            coordinate: [58,0,-14]?  # relative position on platform??


    Name: shared-2-clump
        name: Shared line with two clump weights
        symmetric: True
        
        sections:
          - type: poly_180   
            length: 80       
            connector: clump_weight_20
            
          - type: poly_180
            length: 762   
        
        attachment:
            type:
            coordinate:
```    
    
### Mooring line section properties

```python
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
    keys :   name,   EA ,  MBL ...]
    data :
        -    poly1  , 3232, 23
        -    chain27, 3232, 23
```


### Mooring Connectors

```python
 mooring_connector_types:
    
    h_link:
        mass   : 140    # [kg]  component mass
        volume : 0.13   # [m^3] component volumetric displacement
        
    clump_weight_20:
        mass   : 20000  # [kg]
        volume :  0.8   # [m^3]
        
    buoy_10:
        mass   :  560   # [kg]  component mass
        volume : 10.2   # [m^3] component volumetric displacement
        CdA    :  3.5   # [m^2] produce of cross sectional area and drag coefficient
```

## Anchor types

```python        
anchor_types:
    Name: suction1
        Diameter
        Length
        Embedment depth
        …
```


## Cables

This section describes the cables through the array including both static 
and dynamic portions of cables. At the top level, each array cable going
between a pair of turbines (or a turbine and a substation) is defined. 
This definition can either occur in the [array_cables]() section or the
[cables]() section. The latter provides additional options for defining
the cable routing and burial depth.

### Detailed Array Cable Descriptions

The cables section contains a list of every cable in the array. Here, a cable
is defined as the full assembly of electrical connection equipment between
two turbines or a turbine and a substation. Similar to the [array_cables]() 
section, 'type' links to the cross-section property description of the static
portion of the cable. endA and endB sections define what each end of the cable
is attached to, at what heading it comes off at, and what dynamic cable
profile it uses. Additional fields specify the routing of the static portion
of the cable and the burial depth as as function of cable length.

```python
 cables:

  - name : array_cable1      # descriptive cable name
    type : static_cable_80   # cable section type ID
    
    endA: 
        attachID: turbine_1            # FOWT/substation/junction ID
        heading:  180                  # [deg] heading of attachment at end A
        dynamicID: dynamic_lazy_wave1  # ID of dynamic cable configuration at this end
    
    endB:
        attachID: turbine_2            # FOWT/substation/junction ID
        heading:  30                   # [deg] heading of attachment at end B
        dynamicID: dynamic_lazy_wave1  # ID of dynamic cable configuration at this end
    
    routing_x_y_r:  # optional vertex points along the cable route. Nonzero radius wraps around a point at that radius.
      - [1000, 1200, 20] 
      - [2000, 1500, 20] 
    
    burial:  # optional definition of cable burial depth over its length
        station: [0, 1]                # length along cable, normalized by first and last value
        depth  : [0.1, 0.2]            # [m] burial depth

  - name : array_cable_2     # descriptive cable name
    type : static_cable_80   # cable section type ID
    ...
```

## Dynamic Cable Configurations

This section lists the dynamic cable configurations used in the array design.
Similar to the mooring_line_configs section, it details the assembly of 
cable section that make up a dynamic cable profile, with links to the cross
sectional cable properties. Dynamic cable configurations have some special
properties including specification of the voltage, and the option of 
specifying 'appendages' along the cable length, which can represent discrete
objects like buoyancy modules.

```python
 dynamic_cable_configs:

    dynamic_lazy_wave1
        name: Lazy wave configuration 1 (simpler approach)
        voltage: 66 # [kV]
        span :     # [m] horizontal distance to end of dynamic cable
        
        sections:
          - type: dynamic_cable_27        # ID of a cable section type1
            length: 200                   # [m] length (unstretched)
            
          - type: dynamic_cable_27_w_buoy # (section properties including averaged effect of buoyancy modules)
            length: 300                  
            
          - type: dynamic_cable_27 
            length: 200            
            
        attachment:
            type: j-tube
            coordinate:   # relative location
    

    dynamic_lazy_wave2
        name: Lazy wave configuration 1 (more detailed approach)
        voltage: # [kV]
        span :     # [m] horizontal distance to end of dynamic cable
        
        sections:
          - type: dynamic_cable_27        # ID of a cable section type1
            length: 200                   # [m] length (unstretched)
            appendages:
                type: buoyancy_module_1
                locations: [10,12,13.5,15,18]
                
        attachment:
            type: j-tube
            coordinate:   # relative location
```       
    
### Cable Cross Sectional Properties

This section details the cross-sectional properties of each cable type.

```python
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
product. These appendages are used in the dynamic_cable_configs section.

```python
  cable_appendages:

    buoyancy_module_1:
        mass:    2700   # [kg]  mass
        volume: 8.615   # [m^3] volumetric displacement 
        CdA:     3.8    # [m^2] product of cross-sectional area and drag coefficient
		length:  2.2    # [m]   length taked up along cable
```
