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
The boundaries section contains the boundaries of the array. This information is provided with a list of polygon vertices in order. The vertices are then connected linearly 
to define the boundaries of the array. This information can be used to check that all the floating wind turbines are fully contained within the array boundaries.

```python
    boundaries: #list of polygon vertices in order
       -[x1, y1]
       -[x2, y2]
       -[x3, y3]
       -[x4, y4]
```