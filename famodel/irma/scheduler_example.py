from famodel.irma.scheduler import Scheduler
import numpy as np



# weather
weather = [1, 1, 1, 1, 1]

# tasks
tasks = [
{
    'name': "install_mooring",
    'requirements': ['mooring_reel', 'positioning']
},
{
    'name': "install_anchor",
    'requirements': ['anchor_handling','positioning']
}
]

# assets
assets = [
{
    'name': 'AHTS', 
    'capabilities': ['anchor_handling', 'mooring_reel', 'positioning'],
    'daily_cost': 50000,
    'max_weather': 2
},
{
    'name': 'MPSV', 
    'capabilities': ['mooring_reel', 'positioning'],
    'daily_cost': 25000,
    'max_weather': 1
}
]

# task-asset matrix
task_asset_matrix = np.array([
                [(2000, 2), (1000, 3), (2500, 3)],
                [(1500, 3), (-1,  -1), (4000, 2)]
])

# asset groups
asset_groups = [
{
    'assets': ['AHTS'], 
},
{
    'assets': ['MPSV'], 
},
{
    'assets': ['AHTS','MPSV'], 
},
]

# task dependencies
task_dependencies = {
'install_mooring': ['install_anchor']  # Mooring installation depends on anchor installation
}

# dependency types
dependency_types = {
    'install_anchor->install_mooring': 'start_start'  # Anchor must finish before mooring starts
}

offsets = {
    #'install_anchor->install_mooring': 1    # Mooring installation to start 1 period after Anchor installation
    'install_anchor->install_mooring': (1, 'exact')    # Tuple format: (value, type)
}

# calculate the minimum duration
min_duration = np.min(task_asset_matrix[:, :, 1][task_asset_matrix[:, :, 1] > 0])  # minimum non-zero duration

# intialize and run the scheduler
scheduler = Scheduler(task_asset_matrix, tasks, assets, task_dependencies, dependency_types, offsets, weather, min_duration, asset_groups=asset_groups)
scheduler.optimize()

a = 2

