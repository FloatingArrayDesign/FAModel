
from famodel import Project
from famodel.anchors.anchor import Anchor

# Step 1: Initialize and load soil
proj = Project()
proj.loadSoil(
    filename='inputs/GulfOfMaine_soil_layered_100x100.txt',
    soil_mode='uniform')  # 'uniform' soil does not need from the profile_source yaml file

for label, props in proj.soilProps.items():
    print(f"{label}: {props}")
    
# Convert to profile_map format so anchor capacity models can use it
proj.convertUniformToLayered(default_layer=50.0)

for label, props in proj.profile_map.items():
    print(f"{label}: {props}")

# Step 2: Create and register an anchor at a known position in the grid
anchor = Anchor(
    dd = {'type': 'suction', 'design': {'D': 3.5, 'L': 12.0, 'zlug': 9.67}},
    r  = [54.0, -4450.0, 0.0])

# Step 3: Assign local soil profile from project (nearest neighbor lookup)
soil_id, _ = proj.getSoilAtLocation(anchor.r[0], anchor.r[1])
soil_profile = proj.profile_map[soil_id]  # get compatible layered format
anchor.soilProps = {soil_id: soil_profile}
anchor.setSoilProfile([{'name': soil_id, 'layers': soil_profile}])  # ensures `anchor.soil_profile` is set

# Step 4: Assign loads and line
anchor.loads = {'Hm': 1e6, 'Vm': 5e4}
anchor.line_type = 'chain'
anchor.d = 0.16
anchor.w = 5000.0

# Step 5: Run capacity check and optimization
anchor.getLugForces(
    Hm=anchor.loads['Hm'], Vm=anchor.loads['Vm'], 
    zlug = anchor.dd['design']['zlug'], 
    d=anchor.d, w=anchor.w,  
    plot=True)

anchor.getCapacityAnchor(
    Hm=anchor.loads['Hm'], Vm=anchor.loads['Vm'], 
    zlug = anchor.dd['design']['zlug'], 
    line_type=anchor.line_type, d=anchor.d, w=anchor.w,   
    mass_update=True,
    plot=True)
anchor.getCost()
print(f'Material cost: {anchor.cost["Material cost"]:.2f} USD [2024]')

results = anchor.getSizeAnchor(
    geom = [anchor.dd['design']['L'], anchor.dd['design']['D']],
    geomKeys = ['L', 'D'],
    geomBounds = [(8.0, 15.0), (2.0, 4.0)],
    lambdap_con = [3, 6],
    zlug_fix = False,
    safety_factor = {'SF_combined': 1},
    plot = True)

# Step 6: Report
print('\nFinal Optimized Anchor:')
print('Design:', anchor.dd['design'])
print('Capacity Results:', anchor.anchorCapacity)
anchor.getCost()
print(f'Material cost: {anchor.cost["Material cost"]:.2f} USD [2024]')


