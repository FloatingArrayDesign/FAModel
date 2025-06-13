# Anchors Library

This subpackage of FAModel contains the anchor class and all modules for the  capacity under extreme loads and the installation assessments

## Seabed Conditions

Introduction to different soil types

Heterogenous soil (mixed layers). Map of soil properties for horizontal and vertical spatial-variability.

The reference elevation of the pile is the pile head (z = 0 m), from here all elevations are derived. Thus, Z0 (mudline elevation) is the distance between the pile head and the top of the first layer of soil. Main padeye locations depend on their relative elevation to z0, if zlug > z0 mooring line is embedded below the mudline elevation
### Soil properties
##### Input
- profile_map 
	- name: CPT or reference in the system (-)
  - x, y: coordinates of the anchor within the lease area (m), (m)
  - layers 
    - top, bottom: depth for top and bottom layers (m), (m)
    - soil_type: clay/mud, sand and (weak) rock (-)  
    - soil properties: clay/mud, sand and (weak) rock, see soil properties for further details
- location_name
- z: depth of the layer (m)
- soil types:
  - clay/mud:
    - gamma: submerged soil unit weight (kN/m³)
    - Su: undrained shear strength (kPa)
  - sand:
    - gamma: submerged soil unit weight: (kN/m³)
    - phi: internal friction angle (deg)
    - Dr: relative density (%)
  - rock:
    - UCS: unconfined compressive strength at failure (MPa)
    - Em: rock mass modulus (MPa)
>[!NOTE]
Driven piles are only possible on weak rock, defined here as up to UCS = 5 MPa

> [!IMPORTANT] 
Units within FAModel follow the SI exclusively. The input soil parameters units follow common industry convention. Soil parameters conversion units to Pa and N/m³ take place in the capacity_soils module exclusively. There is no need to convert units.

    profile_map = [
        {
            'name': 'CPT_1',
            'x': 498234, 'y': 5725141,
            'layers': [
                {
                    'top': 1.0, 'bottom': 6.0,
                    'soil_type': 'clay',
                    'gamma_top': 8.0, 'gamma_bot': 8.0,
                    'Su_top': 10, 'Su_bot': 50},
                {
                    'top': 6.0, 'bottom': 15.0,
                    'soil_type': 'sand',
                    'gamma_top': 8.0, 'gamma_bot': 8.0,
                    'phi_top': 32, 'phi_bot': 38,
                    'Dr_top': 70, 'Dr_bot': 75},
                {
                    'top': 15.0, 'bottom': 30.0,
                    'soil_type': 'clay',
                    'gamma_top': 8.0, 'gamma_bot': 9.0,
                    'Su_top': 100, 'Su_bot': 200}]
        }
    ]
    Note:
    - z0 = 1 m, meaning the pile head is 1 m above the mudline
    - soil_type: clay, sand, clay
    - this method allows different soil types and gaps in between soil layers
##### Output
- z0: depth of the mudline relative to the pile head (m)
- soil types:
  - clay/mud:
    - f_gamma: effective unit soil weigtht at depth (N/m³)
    - f_Su: undrained shear strength at depth (Pa)
    - f_sigma_v_eff: effective vertical stress at depth (Pa)
    - f_alpha: adhesion factor from API correlation (-)
  - sand:
    - f_gamma: effective unit soil weigtht at depth (N/m³)
    - f_phi: friction angle at depth (deg)
    - f_sigma_v_eff: effective vertical stress at depth (Pa)
    - f_Dr: relative density at depth (%)
    - f_delta: skin friction angle at depth (-)
  - rock:
    - f_UCS: unconfined compressive strength at depth (Pa)
    - f_Em: rock mass modulus at depth (deg)
------------------------------------------------------------------------------

Soil classification for clay, sand and rock can be found in [Soil Classification Parameters](#soil-classification-parameters).
>[!NOTE] 
>Some anchor capacity functions require input loads at the anchor lug point. These loads can be sent in to the getAnchorCapacity() method, or the getAnchorCapacity() method will calculate the loads by calling getLugLoads(). 
The input loads must be maximum or large loads on the anchor.
		
### Soil classification parameters

The soft, medium and hard clay soil classes are distinguished by the following parameter ranges: 
| clay/mud     | N-Value  |     Eff. unit weight, gamma (kN/m³)  | Void ratio, e (-) | Water content,  (%) | Undrained shear strength, Su (kN/m2) |
|:-----------------:|:--------:|:---------------------------------:|:----------:|:--------------------------------------:|:-------------------------------:|
|     Very Soft     |  0 - 2   |             5.5 - 8.5             |  0.9 - 1.4 |                30 - 50                 |            0 - 12.5             |
|       Soft        |  2 - 4   |             5.5 - 8.5             |  0.9 - 1.4 |                30 - 50                 |            12.5 - 25            |
|       Medium      |  4 - 8   |             5.5 - 8.5             |  0.9 - 1.4 |                30 - 50                 |             25 - 50             |
|       Stiff       |  8 - 15  |             8.5 - 12              |    ~ 0.6   |                20 - 30                 |            50 - 100             |
|     Very Stiff    | 15 - 30  |             8.5 - 12              |    ~ 0.6   |                20 - 30                 |            100 - 200            |
|        Hard       |   < 30   |             8.5 - 12              |    ~ 0.6   |                20 - 30                 |              > 200              |


Sand can also be classified ranging from soft to hard. and are chracterize by the following ranges:

| sand |  N-Value | Eff. unit weight, gamma (kN/m³) | Void ratio, e (-) | Water content,  (%)| Eff. friction angle, phi (deg) | Relative density, Dr (%) |
|:----------------:|:--------:|:----------------------------:|:----------:|:--------------------------------------:|:-------------------:|:--------------------:|
|   Very   Loose   |    > 4   |            7 - 9.5           |    ~ 0.8   |                 25 - 30                |        < 30         |         < 15         |
|       Loose      |  4 - 10  |            7 - 9.5           |    ~ 0.8   |                 25 - 30                |       30 - 35       |        15 - 35       |
|     Compact      | 10 - 30  |          9.5 - 11.5          |   ~ 0.45   |                 12 - 16                |       35 - 40       |        35 - 65       |
|      Dense       | 30 - 50  |          9.5 - 11.5          |   ~ 0.45   |                 12 - 16                |       40 - 45       |       65 - 85        |
|    Very Dense    |   < 50   |          9.5 - 11.5          |   ~ 0.45   |                 12 - 16                |        > 45         |         > 85         |

## Anchor Types
The supported anchor types are listed below with their associated FAModel names in italics. Anchors types have specific [anchor capacity](#anchor-capacity-modules) and [anchor installation](#anchor-installation-modules) application modules, these are shown for clarity below as well.

|                                                        | Capacity       | Installation |
|--------------------------------------------------------|----------------|--------------|
|*DEA*   (drag-embedment anchors)                        | Plate          | Drag         | 
|*SEPLA* (suction embedded plate anchors)                | Plate          | Suction      | 
|*DEPLA* (dynamically embedded plate anchors)            | Plate          | Dynamic      |
|*VLA* (vertically loaded anchors)                       | Plate          | Drag         |
|*suction* (suction caisson/suction bucket anchors) | Suction        | Suction      |
|*torpedo* (torpedo pile anchors)                   | Torpedo        | Dynamic      |
|*helical* (helical pile anchors)                   | Helical/Driven | Torque-crowd |
|*driven*  (driven pile anchors)                    | Driven         | Driven       |
|*dandg* (drilled and grouted pile anchors)         | Drilled&Grout  | Drilled      |

### Anchor geometrical properties
#### DEA/SEPLA/DEPLA/VLA/plate
##### Input
- soil condition: 
    - z, gamma, Su: clay soil parameters (m), (kN/m3), (kPa) 
- geometry:
   - B: width of plate - dimension contained in the vertical plane (m)
   - L: length of plate - dimension perpendicular to the vertical plane (m)
   - zlug: embedded depth of bridle/padeye below mudline (m)
   - beta: angle of plate with horizontal plane (deg)
- loads:
  - Ha, Va: horizontal and vertical loads on padeye of anchor (N), (N)
##### Output

  
#### suction_pile (suction caisson/ suction bucket anchors)
##### Input
- soil condition:
    - location_name:
    - x, y: CPT or reference name
    - layers:
      - z, gamma, Su: clay soil parameters (m), (kN/m3), (kPa) 
      - z, gamma, phi, Dr: sand soil parameters (m), (kN/m3), (deg), (%) 
    
- geometry:
    - D: diameter of pile (m)
    - L: length of pile (m)
    - zlug: embedded depth of padeye below mudline (m)
- loads:
    - Ha, Va: horizontal and vertical loads on padeye of anchor (N), (N)
##### Output
  
#### torpedo_pile (torpedo pile anchors)
##### Input
- soil condition: 
    - z, gamma, Su: clay soil parameters () 
- geometry
   - D1: wing diameter (m)
   - D2: shaft diameter (m)
   - L1: wing length (m)
   - L2: shaft length (m)
   - zlug: embedded depth of padeye below mudline (m)
- loads
  - Ha, Va: horizontal and vertical loads on padeye of anchor (N), (N)

#### helical_pile (helical pile anchors)
##### Input
- soil condition: 
    - z, gamma, Su: clay soil parameters (m), (kN/m3), (kPa) 
    - z, gamma, phi, Dr: sand soil parameters (m), (kN/m3), (deg), (%) 
- geometry
    - D: helix diameter (m)
    - L: shaft length (m)
    - d: shaft diameter (m)
    - zlug: embedded depth of padeye below mudline (m)
- loads
    - Ha, Va: horizontal and vertical loads on padeye of anchor (N), (N)
  
#### driven_pile (driven pile anchors)
##### Input
- soil condition: 
    - z, gamma, Su: clay soil parameters (m), (kN/m3), (kPa) 
    - z, gamma, phi, Dr: sand soil parameters (m), (kN/m3), (deg), (%) 
    - z, UCS, Em: (weak) rock parameters (m), (MPa), (MPa)  
- geometry
    - L: length of pile (m)
    - D: diameter of pile (m)
    - zlug: embedded depth of padeye below mudline (m)
- loads
    - Ha, Va: horizontal and vertical loads on padeye of anchor (N), (N)
##### Output

> [IMPORTANT!] The general output is a lateral and rotational displacement or bending moment. In getCapacityAnchor, the driven pile capacity function is called in a while loop with incremented horizontal 
 input forces until one of the displacements goes past set failure criteria, thus providing a horizontal force capacity output [N]. Vertical capacity [N] is already calculated within the driven pile capacity function.
  
 For non-rock soil, the hinge (bending moment) is also considered as a failure mode along with the lateral and rotational displacement 
 
#### dandg_pile (drilled and grouted pile anchors)
##### Input
- soil condition: 
    - z, UCS, Em: (weak) rock parameters (m), (MPa), (MPa)  
- geometry
   - L: length of pile (m)
   - D: diameter of pile (m)
   - zlug: lug location (m)
- loads
   - Ha, Va: horizontal and vertical loads on padeye of anchor
   
## Loads
Loads derived from MoorPy and DRAFT are considered at a fixed point at mudline elevation. These loads need to be transfered from mudline to lug penetration when the main padeye is below the mudline. Transfer function: soil properties (profile)  mooring line properties (line_type, d and w), loads and zlug 
> [!IMPORTANT] It is cautious to condiser as input load the tension load at mudline since the load will be equal or larger to the tension at lug penetration. Conversely, the angle of the load at lug penetration will equal or larger to the angle at mudline. Therefore, yielding to more vertical componenent. Therefore, Tm >= Ta and thetam <= thetaa 

##### Input
- profile_map: soil profile 
- Tm: tension of the load on mudline (N)
- thetam: angle of the load on mudline (deg)
- zlug: main padeye embeddment (m)
- line_type: type of mooring line ('chain' or 'wire') (-)
- d: mooring line diameter (m)
- w: mooring line unit weight (N/m)
	
> [NOTE] Load components: Hm, Vm: horizontal and vertical load components on mudline (N), (N) and Ha, Va: horizontal and vertical load components on padeye of anchor (N), (N)
##### Output
- Ta: tension of the load on padeye of anchor (N)
- thetaa: angle of the load on padeye of anchor (deg)
- length: length of the embedded line (m)


> [!NOTE] check getLugForces for more details on this transfer function from mudline to lug elevation (below the seabed)
The getTransferLoad function requires **maximum** mudline forces as an input. These forces can be sent in as a dictionary, or anchor.loads dictionary will be searched for 'Hm' and 'Vm' values with additional 
key-value pair 'mudline_force_type':'max' to indicate these mudline forces are maximums.
 
If there are no max mudline forces in the anchor.loads dictionary, getMudlineForces(max_force=True) will be called. Stores results in loads dictionary. 
If lug is at mudline or no lug provided, equates mudline forces with lug forces. 
>[!NOTE]
>The getTransferFunction function called by getLugForces() is tuned to work with maximum loads on the anchor. Some anchor configuration, load, and soil condition combinations may produce invalid results in getTransferFunction. 
For example, the output Va may show as negative. In that case, getLugForces() will warn the user of the invalidity of the result and assign the anchor lug forces = mudline forces.
>[!NOTE] 
>Some anchor capacity functions require input loads at the anchor lug point. These loads can be sent in to the getAnchorCapacity() method, or the getAnchorCapacity() method will calculate the loads by calling getLugLoads(). 
The input loads must be maximum or large loads on the anchor.

------------------------------------------------------------------------------

> [!NOTE] 
> Load inputs to the capacity functions (with the exception of driven & drilled and grouted anchors) are in kN, while the anchor loads dictionary is in N. This conversion is automatically completed in the getAnchorCapacity() 
function so no manual load conversion is required. Load outputs are automatically converted in the getAnchorCapacity function where necessary. 

## Equipment

### Installation vessel
#### Pullard force
Drag installation. Input
  - Tmax: maximum pullard force (N)

#### Crane capacity
Dynamic installation. Output
  - Wp: dynamically installed plate/pile weight (N)

### Installation device
#### Suction pump
Suction installation. Output
  - delta_u_suction: maximum underpressure given by the suction pump during installation (Pa)
  - delta_u_retrieve: maxumum overpressure given by the suction pump during retrieval/removal (Pa)

#### Hydraulic drive head
Torque-crowd installation. Output
  - Torque: torque (Nm)
  - Force: crowd compressive force (N)

#### Hammer
Driven installation. Input
  - hammer_params: 
    - r_m: ram mass (kg) 
    - h: strock height (m) 
    - efficiency: efficiency of the hammer (-)

#### Drill head
Drilled installation

-----------------------------------------------------------------------------
  
## Anchor Class
The anchor class contains properties and methods related to mooring anchors.

The anchor class stores properties and methods that enable a wide range of modeling. 
The [anchor capacity modules](#anchor-capacity-modules) and the [anchor installation modules](#anchor-installation-modules) are integrated with FAModel through the anchor class and its methods.
 
### Anchor modules
Introduction
 
Inspection of the folder: anchors_famodel

|                                              |                |              |
|----------------------------------------------|----------------|--------------|
|![Plate anchor](images/Plateanchors/Plate.png)|![Suction pile anchor](images/Suctionpiles/Suction.png)|![Torpedo pile anchor](images/Torpedopiles/Torpedo.png)| 
|![Helical pile anchor](images/Helicalpiles/Helical.png)|![Driven pile anchor](images/Drivenpiles/Driven.png)|![Drilled and grouted pile anchor](images/Drilledandgroutedpiles/Drilled.png)| 


#### Anchor capacity modules
Analytical static capacity models for extreme load conditions. These models include static soil-structure interaction but the cyclic loading conditions are not covered yet. They will need to follow from further research.

- **capacity_plate** : 
	- getCapacityPlate(profile_map, location_name, D, L, zlug, Ha, Va, plot)
  - capacityPlate dic

- **capacity_suction** : 
	- getCapacitySuction(profile_map, location_name, D, L, zlug, Ha, Va, thetalug, psilug, plot)
  - capacitySuction dic
  
- **capacity_torpedo** : 
	- getCapacityTorpedo(profile_map, location_name, D1, D2, L1, L2, zlug, ballast, Ha, Va, plot)
  - capacityTorpedo dic

- **capacity_helical** : 
	- getCapacityHelical(profile_map, location_name, D, L, d, zlug, Ha, Va, plot)
  - capacityHelical dic

- **capacity_driven** : 
	- getCapacityDriven(profile_map, location_name, L, D, zlug, Ha, Va, plot)
  - capacityDriven dic
 
- **capacity_dandg** : 
	- getCapacityDandG(profile_map, location_name, L, D, zlug, Ha, Va, plot)
  - capacityDandG dic

- **capacity_load** : 
	- getTransferFunction(profile_map, Tm, thetam, zlug, line_type, d, w, plot)
  - capacityLoads dic

#### Anchor installation modules
Analytical installation models for main anchor types.

- **installation_drag** : 
	- getInstallationPlate(profile_map, location_name, B, Lf, Ls, Lca, Lj, plot)
  - installationDrag dic

- **installation_suction** : 
	- getInstallationSuction(profile_map, location_name, D, L, plot)
  - installationSuction dic

- **buckling_suction** : 
	- getBucklingSuction(profile_map, location_name, D, L, plot)
  - installationBuckling dic

- **installation_dynamic** : 
	- getInstallationTorpedo(profile_map, location, D1, D2, L1, L2, ballast, drop_height, plot)
  - installationDynamic dic

- **installation_torque** : 
	- getInstallationHelical(profile_map, location_name, D1, D2, L1, L2, zlug, ballast, Ha, Va, plot)
  - installationTorque dic

- **installation_driven** : 
	- getInstallationDriven(profile_map, location_name, D, L, hammer_params, J_shaft, J_toe, plot)
  - installationDriven dic

- **installation_drill** : 
	- getInstallationDrill(profile_map, location_name, D, L, driller_params, plot)
  - installationDrill dic

#### Anchor support modules

- **anchor_soils** : 
	- clay_profile(profile) 
	- sand_profile(profile) 
	- rock_profile(profile)

- **anchor_solvers** : 
	- fd_solver(n, N, h, D, t, fy, EI, Ha, Va, zlug, z0, k_secant) 

- **anchor_pycurves** : 
	- py_Matlock(z, D, Su, sigma_v_eff, gamma, z0, return_curve)
	- py_API(z, D, phi, sigma_v_eff, Dr, z0, return_curve)
	- py_Reese(z, D, UCS, Em, z0, return_curve) 
	- py_Lovera(z, D, UCS, Em, z0, delta_grout, E_grout, delta_crushed, return_curve)

- **anchor_plots** : 
	- plot_driven()
	- plot_suction()
	- plot_helical()
	- plot_plate()
	- plot_load()
	- plot_pycurve()

### Anchor class methods

- **setSoilProfile()** : Assign a bilinearly interpolated soil profile from the 4 nearest CPTs
- **makeMoorPyAnchor()** : Creates a MoorPy point object representing the anchor in a moorpy system
- **getMudlineForces()** : Finds forces on anchor at mudline using MoorPy Point.getForces method. Use max_force=True to obtain the maximum forces on that anchor from the platform.getWatchCircle() method. 
For more information on the getWatchCircle() calculations, see the [Platform ReadMe](../platform/README.md). An additional anchor.loads dictionary entry is included to describe the mudline load type. 
'mudline_load_type'='max' if max_force=True, and 'mudline_load_type'='current_state' if max_force=False.
- **getLugForces()** : Finds forces at the anchor lug location with getTransferFunction function in capacity_loads.py
- **getCapacityAnchor()** : Calls anchor capacity functions for the correct anchor type using the getLugForces embedded in the method
- **getSizeAnchor()** : Calls sizing anchor functions for the correct anchor type using the getLugForces embedded in the method
- **getFS()** : Computes safety factor for loads on the anchor
- **getCostAnchor()** : Finds costs of anchor from MoorProps and stores in design dictionary
- **getCombinedPlot()** : Create a plot showing the anchor and the inverse catenary overlay in the same coordinate system.

### Anchor Object Properties

- **r** : anchor [x,y,z] position
- **dd** : anchor design dictionary, containing geometric properties, soil properties at the anchor location, cost
- **ms** : moorpy system associated with this anchor point
- **aNum** : anchor index in array mooring list (generally only used for shared moorings)
- **mpAnchor** : moorpy point object that models this anchor
- **anchorCapacity** : dictionary with horizontal and vertical capacity of the anchor. Generally these are loads in [N], but can also be displacements (generally for driven or drilled and grouted piles)
- **loads** : dictionary of loads on the anchor, and the method used to obtain these loads (static or dynamic modeling). 
Loads include mooring line tension (T) with the angle of the load (theta) as well as horizontal (H) and vertical (V) loads components. 
The keys for these loads will either include an m (for loads at the mudline) or a (for loads at the anchor lug / main padeye).
- **soilProps** : dictionary of soil property information at the location of the anchor
- **failure_probability** : dictionary of probabilities for failure of the anchor

### Anchor Type Requirements

Different geometric properties and soil conditions are needed for each anchor type. See the [Anchor Capacity Modules](#anchor-capacity-modules) section for details on the requirements of each anchor type.

## Anchor Capacity Modules
The following list shows the required soil conditions, load information and geometrical properties involved in the anchors' capacity calculations. 


#### Output notes
 The general output is a lateral and rotational displacement. In getAnchorCapacity, the drilled and grouted pile function is called in a while loop with incremented horizontal input forces 
 until one of the displacements goes past set failure criteria, thus providing a horizontal force capacity output [N]. Vertical capacity [N] is already calculated within the driven pile capacity function. 
 


