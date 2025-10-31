"""
Capability-Based Asset Group Generator for MILP Scheduler

This module provides intelligent asset group generation based on capability matching,
designed to work with offshore installation scheduling problems. It creates sparse
task-asset matrices by pre-filtering operationally feasible combinations.

Key Features:
- Capability-based matching between tasks and assets
- Smart pre-filtering to avoid 2^N explosion
- Strategic batch support for aggregated installation tasks
- Sparse matrix generation for computational efficiency
- Operational feasibility validation

Strategic batch support allows multiple task alternatives like:
- "install_1_anchor" vs "install_4_anchors" with mutual exclusion
- Pre-aggregated task strategies for large-scale projects
"""

import numpy as np
from itertools import combinations


class TaskAssetGroupGenerator:
    """
    Generator for asset group assignments to tasks based on shared capabilities.
    
    Creates sparse task-asset matrices by matching task capability requirements 
    to asset capabilities, with support for flexible batch sizes per task.
    
    Features:
    - Configurable batch sizes (anchors per task)
    - Capability-based asset matching
    - Smart pre-filtering for computational efficiency
    - Automatic handling of remainder tasks for uneven divisions
    """
    
    def __init__(self, max_group_size=3):
        """Initialize the generator with constraints."""
        self.max_group_size = max_group_size
        self.verbose = False
        
    def generate_asset_groups(self, task_definitions, asset_definitions):
        """Generate asset groups and return scheduler-ready inputs."""
        if self.verbose:
            print("=== Capability-Based Asset Group Generation ===")
            
        # === VALIDATION PHASE ===
        validated_tasks = self._validate_task_definitions(task_definitions)
        validated_assets = self._validate_asset_definitions(asset_definitions)
        
        if self.verbose:
            print(f"\n=== VALIDATION RESULTS ===")
            print(f"Validated tasks: {len(validated_tasks)}")
            print(f"Validated assets: {len(validated_assets)}")
        
        # === ASSET GROUP GENERATION ===
        asset_groups = self._generate_feasible_asset_groups(validated_assets)
        
        if self.verbose:
            print(f"\n=== ASSET GROUP RESULTS ===")
            print(f"Feasible asset groups: {len(asset_groups)}")
        
        # === TASK-ASSET MATCHING ===
        task_asset_matches = self._match_tasks_to_asset_groups(validated_tasks, asset_groups)
        task_asset_matrix = self._build_sparse_matrix(validated_tasks, asset_groups, task_asset_matches)
        
        # === SCHEDULER INPUT PREPARATION ===
        scheduler_inputs = {
            "task_asset_matrix": task_asset_matrix,
            "tasks": [task["name"] for task in validated_tasks],
            "assets": [group["name"] for group in asset_groups],
            "asset_groups": asset_groups
        }
        
        if self.verbose:
            self._print_efficiency_stats(scheduler_inputs)
        
        return scheduler_inputs
    
    def _validate_task_definitions(self, task_definitions):
        """Validate and standardize task capability requirements."""
        validated_tasks = []
        
        if self.verbose:
            print(f"\n=== TASK VALIDATION ===")
        
        for task_name, requirements in task_definitions.items():
            if "required_capabilities" not in requirements:
                raise ValueError(f"Task '{task_name}' missing required_capabilities")
            
            # === TASK CONFIGURATION ===
            validated_task = {
                "name": task_name,
                "required_capabilities": set(requirements["required_capabilities"]),
                "min_weather_rating": requirements.get("min_weather_rating", 1),
                "max_duration": requirements.get("max_duration", 24),
                "complexity_factor": requirements.get("complexity_factor", 1.0),
                "batch_size": requirements.get("batch_size", 1)
            }
            validated_tasks.append(validated_task)
            
            if self.verbose:
                batch_info = f" (batch_size={validated_task['batch_size']})" if validated_task['batch_size'] > 1 else ""
                print(f"  Task: {task_name} requires {validated_task['required_capabilities']}{batch_info}")
        
        return validated_tasks
    
    def _validate_asset_definitions(self, asset_definitions):
        """Validate and standardize asset capabilities."""
        validated_assets = []
        
        if self.verbose:
            print(f"\n=== ASSET VALIDATION ===")
        
        for asset_name, capabilities in asset_definitions.items():
            if "capabilities" not in capabilities:
                raise ValueError(f"Asset '{asset_name}' missing capabilities")
            
            # === ASSET CONFIGURATION ===
            validated_asset = {
                "name": asset_name,
                "capabilities": set(capabilities["capabilities"]),
                "max_weather": capabilities.get("max_weather", 1),
                "base_cost": capabilities.get("base_cost", 10000),
                "daily_rate": capabilities.get("daily_rate", 5000),
                "availability": capabilities.get("availability", 1.0)
            }
            validated_assets.append(validated_asset)
            
            if self.verbose:
                print(f"  Asset: {asset_name} provides {validated_asset['capabilities']}")
        
        return validated_assets
    
    def _generate_feasible_asset_groups(self, validated_assets):
        """Generate all operationally feasible asset group combinations."""
        asset_groups = []
        
        if self.verbose:
            print(f"\n=== ASSET GROUP GENERATION ===")
        
        # === INDIVIDUAL ASSET GROUPS ===
        for asset in validated_assets:
            group = {
                "name": asset["name"],
                "assets": [asset["name"]],
                "combined_capabilities": asset["capabilities"],
                "min_weather": asset["max_weather"],
                "total_cost": asset["base_cost"],
                "total_daily_rate": asset["daily_rate"],
                "group_size": 1
            }
            asset_groups.append(group)
        
        if self.verbose:
            print(f"  Individual asset groups: {len(asset_groups)}")
        
        # === MULTI-ASSET COMBINATIONS ===
        combination_count = 0
        for size in range(2, min(len(validated_assets) + 1, self.max_group_size + 1)):
            for combo in combinations(validated_assets, size):
                if self._is_operationally_feasible(combo):
                    group_props = self._calculate_group_properties(combo)
                    asset_groups.append(group_props)
                    combination_count += 1
        
        if self.verbose:
            print(f"  Multi-asset combinations: {combination_count}")
            print(f"  Total asset groups: {len(asset_groups)}")
        
        return asset_groups
    
    def _is_operationally_feasible(self, asset_combination):
        """Apply operational feasibility filters to asset combinations."""
        
        # === FEASIBILITY FILTER 1: Weather compatibility ===
        # All assets must handle similar weather conditions
        weather_ratings = [asset["max_weather"] for asset in asset_combination]
        if max(weather_ratings) - min(weather_ratings) > 2:
            return False
        
        # === FEASIBILITY FILTER 2: Capability overlap ===
        # Avoid redundant capabilities (inefficient combinations)
        all_capabilities = [asset["capabilities"] for asset in asset_combination]
        total_capabilities = set().union(*all_capabilities)
        individual_count = sum(len(caps) for caps in all_capabilities)
        
        # Reject if overlap is too high (more than 70% overlap)
        overlap_ratio = (individual_count - len(total_capabilities)) / individual_count
        if overlap_ratio > 0.7:
            return False
        
        # === FEASIBILITY FILTER 3: Cost efficiency ===
        # Combination shouldn't be extremely expensive
        total_cost = sum(asset["base_cost"] for asset in asset_combination)
        avg_individual_cost = total_cost / len(asset_combination)
        if avg_individual_cost > 100000:
            return False
        
        # === FEASIBILITY FILTER 4: Group size limits ===
        # Practical limits on group size
        if len(asset_combination) > self.max_group_size:
            return False
            
        return True
    
    def _calculate_group_properties(self, asset_combination):
        """Calculate combined properties for an asset group."""
        assets = list(asset_combination)
        asset_names = [asset["name"] for asset in assets]
        
        # === CAPABILITY COMBINATION ===
        combined_capabilities = set()
        for asset in assets:
            combined_capabilities.update(asset["capabilities"])
        
        # === GROUP PROPERTY CALCULATION ===
        min_weather = min(asset["max_weather"] for asset in assets)
        total_cost = sum(asset["base_cost"] for asset in assets)
        total_daily_rate = sum(asset["daily_rate"] for asset in assets)
        group_name = "+".join(asset_names)
        
        return {
            "name": group_name,
            "assets": asset_names,
            "combined_capabilities": combined_capabilities,
            "min_weather": min_weather,
            "total_cost": total_cost,
            "total_daily_rate": total_daily_rate,
            "group_size": len(assets)
        }
    
    def _match_tasks_to_asset_groups(self, validated_tasks, asset_groups):
        """Match tasks to feasible asset groups based on capability requirements."""
        task_asset_matches = {}
        
        for task in validated_tasks:
            task_name = task["name"]
            feasible_groups = []
            
            for group in asset_groups:
                if self._can_group_handle_task(group, task):
                    # Calculate cost and duration for this task-group combination
                    cost, duration = self._calculate_task_cost_duration(task, group)
                    feasible_groups.append({
                        "group_name": group["name"],
                        "cost": cost,
                        "duration": duration
                    })
            
            task_asset_matches[task_name] = feasible_groups
            
            if self.verbose:
                print(f"  Task '{task_name}' can be handled by {len(feasible_groups)} asset groups")
        
        return task_asset_matches
    
    def _can_group_handle_task(self, asset_group, task):
        """Check if an asset group can handle a specific task."""
        # Capability check
        required_caps = task["required_capabilities"]
        available_caps = asset_group["combined_capabilities"]
        
        if not required_caps.issubset(available_caps):
            return False
        
        # Weather rating check
        if asset_group["min_weather"] < task["min_weather_rating"]:
            return False
        
        return True
    
    def _calculate_task_cost_duration(self, task, asset_group):
        """Calculate cost and duration for a task-asset group combination."""
        
        # === DURATION CALCULATION ===
        base_duration = task["max_duration"]
        complexity_factor = task["complexity_factor"]
        batch_size = task.get("batch_size", 1)
        
        # Duration scales with complexity and batch size, but with efficiency gains
        batch_efficiency = 1.0 if batch_size == 1 else (batch_size * 0.8)  # 20% efficiency gain for batches
        duration = base_duration * complexity_factor * batch_efficiency
        
        # === COST CALCULATION ===
        setup_cost = asset_group["total_cost"] * 0.1  # 10% of asset cost as setup
        operational_cost = asset_group["total_daily_rate"] * (duration / 24)  # Daily rate prorated
        batch_cost_factor = batch_size * 0.9  # 10% cost efficiency for larger batches
        
        total_cost = (setup_cost + operational_cost) * batch_cost_factor
        
        return round(total_cost, 2), round(duration, 2)
    
    def _build_sparse_matrix(self, validated_tasks, asset_groups, task_asset_matches):
        """Build sparse task-asset matrix avoiding mostly (-1,-1) entries."""
        num_tasks = len(validated_tasks)
        num_groups = len(asset_groups)
        
        # Initialize with infeasible values using object array
        matrix = np.empty((num_tasks, num_groups), dtype=object)
        matrix.fill((-1, -1))
        
        # Fill in feasible combinations
        for task_idx, task in enumerate(validated_tasks):
            task_name = task["name"]
            feasible_groups = task_asset_matches.get(task_name, [])
            
            for match in feasible_groups:
                # Find asset group index
                group_idx = next(i for i, group in enumerate(asset_groups) 
                               if group["name"] == match["group_name"])
                matrix[task_idx, group_idx] = (match["cost"], match["duration"])
        
        return matrix
    
    def _print_efficiency_stats(self, scheduler_inputs):
        """Print efficiency statistics about the generated matrix."""
        matrix = scheduler_inputs["task_asset_matrix"]
        total_entries = matrix.size
        
        # Count feasible entries by checking each element
        feasible_entries = 0
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                if matrix[i, j] != (-1, -1):
                    feasible_entries += 1
        
        sparsity = 1 - (feasible_entries / total_entries)
        
        print(f"\n=== Efficiency Statistics ===")
        print(f"Task-Asset Matrix: {matrix.shape}")
        print(f"Total entries: {total_entries}")
        print(f"Feasible entries: {feasible_entries} ({feasible_entries/total_entries:.1%})")
        print(f"Sparsity: {sparsity:.1%} (reduced computational load)")
        print(f"Asset groups: {len(scheduler_inputs['assets'])}")

def generate_capability_based_groups(task_definitions, asset_definitions, max_group_size=3, verbose=True):
    """
    Convenience function to generate capability-based asset groups.
    
    Args:
        task_definitions (dict): Task capability requirements
        asset_definitions (dict): Asset capabilities  
        max_group_size (int): Maximum assets per group
        verbose (bool): Print detailed output
    
    Returns:
        dict: Scheduler inputs ready for use with the MILP scheduler
    
    Example:
        task_defs = {
            "install_anchor_task_1": {
                "required_capabilities": ["anchor_handling", "positioning"],
                "min_weather_rating": 1,
                "max_duration": 12,
                "batch_size": 1
            }
        }
        
        asset_defs = {
            "anchor_vessel": {
                "capabilities": ["anchor_handling", "positioning"],
                "max_weather": 2,
                "base_cost": 30000
            }
        }
        
        scheduler_inputs = generate_capability_based_groups(task_defs, asset_defs)
    """
    generator = TaskAssetGroupGenerator(max_group_size=max_group_size)
    generator.verbose = verbose
    return generator.generate_asset_groups(task_definitions, asset_definitions)


if __name__ == "__main__":
    # Configurable anchor installation demo
    # 
    # CONFIGURATION PARAMETERS:
    # - num_anchors: Total number of anchors to install for your project
    # - anchors_per_task: Batch size - how many anchors each task will install
    #
    # EXAMPLES:
    # - num_anchors=100, anchors_per_task=1   → 100 individual tasks
    # - num_anchors=100, anchors_per_task=4   → 25 batch tasks (4 anchors each)
    # - num_anchors=100, anchors_per_task=100 → 1 mega-batch task
    # - num_anchors=7, anchors_per_task=3     → 2 tasks (3 anchors) + 1 task (1 anchor)

    num_anchors = 4        # Total number of units to install
    anchors_per_task = 1   # Batch size: anchors installed per task
    
    # Calculate strategy based on batch size
    if anchors_per_task == 1:
        strategy = 1  # Individual tasks
    elif anchors_per_task > 1 and anchors_per_task < num_anchors:
        strategy = 2  # Intermediate batches
    elif anchors_per_task == num_anchors:
        strategy = 3
    else:
        raise ValueError("Input strategy is not yet supported")
    
    print(f"=== Anchor Installation Demo ===")
    print(f"Total anchors: {num_anchors}")
    print(f"Anchors per task: {anchors_per_task}")
    print(f"Strategy: {strategy}\n")

    task_definitions = {}      # initialize the dictionary of tasks
    
    if strategy == 1:
        # Strategy 1: Multiple individual anchor installation tasks
        num_tasks = num_anchors // anchors_per_task
        print(f"Strategy 1: {num_tasks} tasks, each installing {anchors_per_task} anchor(s)")
        
        for i in range(1, num_tasks + 1):
            task_name = f"install_anchor_{i}"
            task_definitions[task_name] = {
                "required_capabilities": ["anchor_handling", "positioning"],
                "max_duration": 12 * anchors_per_task,  # Scale duration with batch size
                "batch_size": anchors_per_task
            }
        
    elif strategy == 2:
        # Strategy 2: Intermediate batches
        num_batches = num_anchors // anchors_per_task
        print(f"Strategy 2: {num_batches} batch tasks, each installing {anchors_per_task} anchors")
        
        for i in range(1, num_batches + 1):
            task_name = f"install_batch_{i}"
            task_definitions[task_name] = {
                "required_capabilities": ["anchor_handling", "positioning"],
                "max_duration": 8 * anchors_per_task,  # Batch efficiency
                "batch_size": anchors_per_task
            }
        
        # Handle remainder anchors
        remainder = num_anchors % anchors_per_task
        if remainder > 0:
            task_name = f"install_batch_{num_batches + 1}"
            task_definitions[task_name] = {
                "required_capabilities": ["anchor_handling", "positioning"],
                "max_duration": remainder * 8,
                "batch_size": remainder
            }
    
    elif strategy == 3:
        # Single batch for all anchors
        print(f"Strategy 3: 1 batch task installing {num_anchors} anchors at once")
        task_definitions["install_anchors"] = {
            "required_capabilities": ["anchor_handling", "positioning"],
            "max_duration": num_anchors * 8,  # Batch efficiency: 8 hours per anchor
            "batch_size": num_anchors
        }
        
    else:
        raise ValueError("Strategy must be 1, 2, or 3")
    
    # Same asset definitions for both strategies
    asset_definitions = {
        "anchor_vessel": {
            "capabilities": ["anchor_handling", "positioning"],
            "max_weather": 2,
            "base_cost": 30000,
            "daily_rate": 15000
        },
        "positioning_vessel": {
            "capabilities": ["positioning"],
            "max_weather": 3,
            "base_cost": 15000,
            "daily_rate": 6000
        }
    }
    
    # Generate asset groups for anchor tasks
    scheduler_inputs = generate_capability_based_groups(
        task_definitions, 
        asset_definitions, 
        max_group_size=2
    )
    
    print(f"=== Results ===")
    print(f"Total anchors configured: {num_anchors}")
    print(f"Anchors per task: {anchors_per_task}")
    print(f"Generated tasks: {scheduler_inputs['tasks']}")
    print(f"Number of tasks: {len(scheduler_inputs['tasks'])}")
    print(f"Asset groups: {len(scheduler_inputs['assets'])}")
    print(f"Matrix shape: {scheduler_inputs['task_asset_matrix'].shape}")