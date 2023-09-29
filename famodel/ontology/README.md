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



## Array Mooring
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
