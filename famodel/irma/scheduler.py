# author: @rdavies, 9-8-2025

# Scheduler class for managing actions and tasks
# WIP, to be merged into Irma later on

'''
--- TODO List ---
- [] How to expand this to multiple assets per task?
- [] Eventually enable parallel tasks and multiple assets per task
- [] Convert input tasks and assets from dicts to Task and Asset objects
 - [] When tasks and assets are converted from lists to objects, update the type hints for task and asset list at class initialization.
- [] Add a delay cost, i.e. a cost for each time period where X = 0
- [] Do we want to return any form of info dictionary?
- [] Figure out if this can be parallelized
'''

from famodel.irma.task import Task
from famodel.irma.assets import Asset
from scipy import optimize
from scipy.optimize import milp
import numpy as np
import os

wordy = 1  # level of verbosity for print statements

class Scheduler:

    # Inputs are strictly typed, as this is an integer programming problem (ignored by python at runtime, but helpful for readability and syntax checking).
    def __init__(self, task_asset_matrix : np.ndarray, tasks : list[str] = [], assets : list[dict] = [], task_dependencies = {}, weather : list[int] = [], period_duration : float = 1, **kwargs):
        '''
        Initializes the Scheduler with assets, tasks, and constraints.

        Inputs
        ------
        task_asset_matrix : array-like
            A 3D array of (cost, duration) tuples indicating the cost and duration for each asset to perform each task.
            Must be len(tasks) x len(assets) x 2.
        tasks : list
            A list of Task objects to be scheduled.
        assets : list
            A list of Asset objects to be scheduled.
        task_dependencies : dict
            A dictionary mapping each task to a list of its dependencies.
        weather : list
            A list of weather windows. The length of this list defines the number of discrete time periods available for scheduling.
        period_duration : float
            The duration of each scheduling period. Used for converting from periods to real time.
        kwargs : dict
            Additional keyword arguments for future extensions.

        Returns
        -------
        None
        '''

        if wordy > 0:
            print("Initializing Scheduler...")

        self.task_asset_matrix = task_asset_matrix
        self.tasks = tasks
        self.assets = assets
        self.weather = weather
        self.task_dependencies = task_dependencies

        self.num_tasks = len(self.tasks)
        self.num_assets = len(self.assets)

        self.period_duration = period_duration # duration of each scheduling period. Used for converting from periods to real time.
        self.num_periods = len(weather)  # number of scheduling periods

        # --- Check for valid inputs ---

        # check for valid task_asset_matrix dimensions (must be len(tasks) x len(assets) x 2)
        if self.task_asset_matrix.ndim != 3 or self.task_asset_matrix.shape[0] != len(self.tasks) or self.task_asset_matrix.shape[1] != len(self.assets) or self.task_asset_matrix.shape[2] != 2:
            raise ValueError("task_asset_matrix must be a 3D array with shape (len(tasks), len(assets), 2).")
        
        # check for integer matrix, try to correct
        if self.task_asset_matrix.dtype != np.dtype('int'):
            try:
                self.task_asset_matrix = self.task_asset_matrix.astype(int)
            except:
                raise ValueError("task_asset_matrix must be a 3D array of integers with shape (len(tasks), len(assets), 2).")
            else: 
                print("Input task_asset_matrix was not integer. Converted to integer type.")

        # check for valid tasks and assets
        if not all(isinstance(task, str) for task in self.tasks):
            raise ValueError("All elements in tasks must be strings.")
        if not all(isinstance(asset, dict) for asset in self.assets):
            raise ValueError("All elements in assets must be dictionaries.")

        # check for valid weather
        if not all(isinstance(w, int) and w >= 0 for w in self.weather):
            raise ValueError("All elements in weather must be non-negative integers representing weather severity levels.")
        
        # check period duration is valid
        if self.period_duration <= 0:
            raise ValueError("period_duration must be positive non-zero.")
        
        # --- Debug helpers ---
        # make a list of indicies to help with building constraints
        self.x_indices = []
        for p in range(self.num_periods):
            for t in range(self.num_tasks):
                for a in range(self.num_assets):
                    self.x_indices.append(f"x_[{p}][{t}][{a}]")

        if wordy > 0:
            print(f"Scheduler initialized with {self.num_periods} time periods, {self.num_tasks} tasks, and {self.num_assets} assets.")

    def set_up_optimizer(self, goal : str = "cost"):
        '''
        Workspace for building out an optimizer. Right now, assuming the goal is minimize cost. This could easily be reworked to minimize duration, or some other value.

        This is a binary-integer linear programming problem, which can be solved with scipy.optimize.milp.

        Inputs
        ------  
        goal : str
            The optimization goal, minimize either "cost" or "duration". Default is "cost".

        Returns
        -------
        values : np.ndarray
            The values vector for the optimization problem.
        constraints : list
            A list of constraints for the optimization problem.
        integrality : np.ndarray
            An array that sets decision variables as integers.
        bounds : scipy.optimize.Bounds
            The bounds for the decision variables (0-1).
        '''

        if wordy > 0:
            print("Setting up the optimizer...")

        # Solves a problem of the form minimize:     v^T * x
        #                        subject to:   A_ub * x <= b_ub
        #                                      A_eq * x == b_eq
        #                                      A_lb * x >= b_lb         
        #                                       lb <= x <= ub        # These are constrained as integers on range 0-1

        # --- Check and Process Inputs ---
        if goal == "cost":
            goal_index = 0
        elif goal == "duration":
            goal_index = 1
        else:
            raise ValueError("goal must be either 'cost' or 'duration'.")
        
        # --- Build the objective function --- 

        # v^T * x

        # shape of V = shape of X, these will be flattened to equal-length vectors for the solver
        V = np.zeros((self.num_periods, self.num_tasks, self.num_assets), dtype=int) # Values matrix: value of asset a assigned to task t in period p
        X = np.zeros(V.shape, dtype=int) # Decision variable matrix:  X[p, t, a] = 1 if asset a is assigned to task t in period p, else 0

        # Values vector: In every planning period, the value of assigning asset a to task t is the same. Constraints determine which periods are chosen.
            # Note: Intentionally using values here instead of "cost" to avoid confusion between the program 'cost' of a pairing (which could be financial cost, duration, or some other target metric for minimization) to the solver and the financial cost of a asset-task pairing.
        for p in range(self.num_periods):
            for t in range(self.num_tasks):
                for a in range(self.num_assets):
                    V[p, t, a] = self.task_asset_matrix[t, a, goal_index]  # cost

        
        if wordy > 1:
            print("Values matrix V (periods x tasks x assets) of length " + str(V.flatten().shape[0]) + " created")
            print("Decision variable matrix X (periods x tasks x assets) of same shape initialized to zeros.")

        # Decision variable start as 0's, nothing decided. Constrainted to 0 or 1 by bounds and integrality in flattening objective function section.

        # --- Flatten the objective function for the solver ---

        # values vector (horizontal vector)
        values = V.flatten()  # Flatten the values tensor (num_periods x num_tasks x num_assets) into a 1D array for the solver. Solver requires 1D problem

        # decision vars (vertical vector)
        decision_vars = X.flatten()  # Flatten the decision variable tensor (num_periods x num_tasks x num_assets) into a 1D array for the solver. Solver requires 1D problem.

        # lb <= x <= ub
            # Constrain decision variables to be 0 or 1
        bounds = optimize.Bounds(0, 1)  # 0 <= x_i <= 1
        integrality = np.full_like(decision_vars, True)  # x_i are integers

        if wordy > 0:
            print("Bounds and integrality for decision variables set. Begining to build constraints...")

        # --- build the constraints ---

        # A_ub * x <= b_ub
        # A_eq * x == b_eq
        # A_lb * x >= b_lb 

        '''
        A note on constraints: There are two constraint matrices, the equality constraints (A_eq, b_eq) and the upper bound constraints (A_ub, b_ub).
        Each row in the coefficient matrices corresponds to a constraint, and each column corresponds to a decision variable. Thus the number of columns 
        is equal to the number of decision variables (num_periods * num_tasks * num_assets), and the number of rows is equal to the number of constraints.
        Similarly, the length of the limits matrices (b_eq, b_ub) is equal to the number of constraints.

        The equality constraints are expressed in the form A_eq * x = b_eq. Where A_eq is the coefficient matrix and b_eq is the limits matrix. 
        For example, the constraints 5x+3y=15 and x-y=1 can be expressed as:
        A_eq = [[5, 3],
                [1, -1]]
        b_eq = [15, 1] 

        Similarly, the upper bound constraints are expressed in the form A_ub * x <= b_ub. Where A_ub is the coefficient matrix and b_ub is the limits matrix.
        For example, the constraints 2x+3y<=12 and x+y<=5 can be expressed as:
        A_ub = [[2, 3],
                [1, 1]]
        b_ub = [12, 5]

        The lower bound constraints (A_lb and b_lb) follow the same form as the upper bound constraints.

        The lower and upper bound constraints on the decision variables (lb <= x <= ub) is handled above, limiting them to integer values of 0 or 1.

        The shape of decision vars is: 
            period_0_task_0_asset_0, period_0_task_0_asset_1, ..., period_0_task_1_asset_0, ..., period_1_task_0_asset_0, ...
        or equivalently:
            period_0_task_0_asset_0, ..., period_0_task_0_asset_max, period_0_task_1_asset_0, ..., period_0_task_max_asset_max, period_1_task_0_asset_0, ..., period_max_task_max_asset_max,

        Constraints column (decision variable) indexing used in definitions below:
        x = x_000, ..., x_00a, x_001, ..., x_00a, ..., x_0t1, ..., x_0ta,  # task asset pairings in period 0
            x_100, ..., x_10a, x_101, ..., x_10a, ..., x_1t1, ..., x_1ta,  # task asset pairings in period 1
            ...,
            x_p00, ..., x_p0a, x_p21, ..., x_p2a, ..., x_pt1, ..., x_pta   # task asset pairings in period p

        where:
        - 0:p is the period index, 0:t is the task index, and 0:a is the asset index
        - x is the flattened decision variable tensor X[p, t, a]
        '''

        # Empty list of constraint coefficient matrices
        A_ub_list = []
        A_eq_list = []
        A_lb_list = []

        # Empty list of constraint limit vectors
        b_ub_list = []
        b_eq_list = []
        b_lb_list = []

        # 0) Total number of assignments needs to be less than or equal to the number of periods available
        '''
        the sum of the total amount of periods assigned to tasks cannot be more than the total number of periods available:
        (x_000 + ... + x_pta) <= N
        '''

        # 1 row
        A_ub_0 = np.ones((1, len(decision_vars)), dtype=int)  # Every period assigned to a task counts as 1 towards the total assigned periods. This assumes one pair per period
        b_ub_0 = np.array([self.num_periods])

        A_ub_list.append(A_ub_0)
        b_ub_list.append(b_ub_0)

        if wordy > 0:
            print("Constraint 0 built.")

        # 1) asset can only be assigned to a task if asset is capable of performing the task (value of pairing is non-negative)
        '''
        if task j cannot be performed by asset k, then x_pjk = 0 for all periods p

        (x_0jk + ... + x_pjk) = 0   # for all tasks j in range(0:t) and assets k in range(0:a) where task_asset_matrix[j, k, goal_index] < 0
        '''
        
        mask = np.zeros(X.shape, dtype=int)

        for p in range(self.num_periods):
            mask[p,:,:] = self.task_asset_matrix[:, :, goal_index] < 0  # Create a mask of invalid task-asset pairings where cost is negative (indicating invalid)

        # 1 row
        A_eq_1 = mask.flatten().reshape(1, -1)  # Flatten the mask to match the decision variable shape, and reshape to be a single row
        b_eq_1 = np.zeros(A_eq_1.shape[0], dtype=int)

        if wordy > 2:  # example debugging code for looking at indicies. Can be applied to any constraint matrix if row index is adjusted accordingly
            print(f"Task {t}:")
            for i in range(len(self.x_indices)):
                print(f"  {self.x_indices[i]}: {A_eq_1[0, i]}")

        A_eq_list.append(A_eq_1)
        b_eq_list.append(b_eq_1)

        # 2) task dependencies must be respected (i.e., a task cannot start until all its dependencies have been satisfied)
        # TODO: enforce task dependencies

        # 3) assets cannot be assigned in a time period where the weather is above the maximum capacity
        # TODO: weather is disabled until this is added

        # 4) assets cannot be assigned to multiple tasks in the same time period
        '''
        this is a simplification, eventually we may want to allow multiple assets per task, or parallel tasks

        Sum of all pairings in a period must be <= 1:
        (x_000 + ... + x_0ta) <= 1   # for period 0
        (x_100 + ... + x_1ta) <= 1   # for period 1
        ...
        (x_p00 + ... + x_pta) <= 1   # for period p
        '''
        
        # num_periods rows
        A_ub_4 = np.zeros((self.num_periods, len(decision_vars)), dtype=int)
        b_ub_4 = np.ones(self.num_periods, dtype=int)  # right-hand side is 1 for each period

        for p in range(self.num_periods):
            # Create a mask for all variables in period p
            mask = np.zeros(X.shape, dtype=int) 
            mask[p, :, :] = 1             # Set all task-asset pairings in period p to 1 (so they are included in the sum)
            A_ub_4[p, :] = mask.flatten()

        A_ub_list.append(A_ub_4)
        b_ub_list.append(b_ub_4)

        if wordy > 0:
            print("Constraint 4 built.")

        # 5) The total number of tasks assigned cannot be greater than the number of tasks available
        # TODO: enforce task limits

        # 6) The total number of assets assigned cannot be greater than the number of assets available
        # TODO: enforce asset limits

        # 7) There is a penalty associated with a time period with no assigned tasks
        # TODO: This delay costs enforces tasks are finished as soon as possible

        # 8) All tasks must be assigned to a time period once
        '''

        The sum of all decision variables for each task must equal 1, indicating that each task is assigned to exactly one time period. This does not support parallel tasks.

        (x_000 + ... + x_p0a) = 1   # for task 0
        (x_010 + ... + x_p1a) = 1   # for task 1
        ...
        (x_pta + ... + x_pta) = 1   # for task t
        '''

        # num_tasks rows
        A_eq_8 = np.zeros((self.num_tasks, len(decision_vars)), dtype=int)
        b_eq_8 = np.ones(self.num_tasks, dtype=int)

        for t in range(self.num_tasks):
            # Create a mask for all variables for task t
            mask = np.zeros(X.shape, dtype=int)
            mask[:, t, :] = 1             # Set all periods and assets for task t to 1 (so they are included in the sum)
            A_eq_8[t, :] = mask.flatten()

            if wordy > 0:
                print(f"Task {t}:")
                for i in range(len(self.x_indices)):
                    print(f"  {self.x_indices[i]}: {A_eq_8[t, i]}")

        A_eq_list.append(A_eq_8)
        b_eq_list.append(b_eq_8)

        if wordy > 0:
            print("Constraint 8 built.")

        # 9) A task must be assigned to the number of time periods equal to its duration for the asset assigned to it
        # TODO: this enforces task lengths for assets that are assigned



        if wordy > 0:
            print("All constraints built. Stacking and checking constraints...")

        # --- Build the constraints ---
        # A series of linear constraints required by the solver by stacking the constraint matrices and limits vectors
        # The number of rows in these matrices is equal to the number of constraints, so they can be vertically stacked

        # Check num columns of all constraint matrices matches number of decision variables before stacking
        for i, A in enumerate(A_ub_list):
            if A.size > 0 and A.shape[1] != decision_vars.shape[0]:
                raise ValueError(f"Upper bound constraint matrix {i} has incorrect number of columns. Expected {decision_vars.shape[0]}, got {A.shape[1]}.")
        for i, A in enumerate(A_eq_list):
            if A.size > 0 and A.shape[1] != decision_vars.shape[0]:
                raise ValueError(f"Equality constraint matrix {i} has incorrect number of columns. Expected {decision_vars.shape[0]}, got {A.shape[1]}.")
        for i, A in enumerate(A_lb_list):
            if A.size > 0 and A.shape[1] != decision_vars.shape[0]:
                raise ValueError(f"Lower bound constraint matrix {i} has incorrect number of columns. Expected {decision_vars.shape[0]}, got {A.shape[1]}.")

        # Stack, check shapes of final matrices and vectors, and save the number of constraints for later use
        if len(A_ub_list) > 0:
            A_ub = np.vstack(A_ub_list)  # upperbound coefficient matrix
            b_ub = np.concatenate(b_ub_list)  # upperbound limits vector
            if A_ub.shape[0] != b_ub.shape[0]:
                raise ValueError(f"A_ub and b_ub have inconsistent number of rows. A_ub has {A_ub.shape[0]}, b_ub has {b_ub.shape[0]}.")
            self.num_ub_constraints = A_ub.shape[0]
        else:
            self.num_ub_constraints = 0

        if len(A_eq_list) > 0:
            A_eq = np.vstack(A_eq_list)  # equality coefficient matrix
            b_eq = np.concatenate(b_eq_list)  # equality limits vector
            if A_eq.shape[0] != b_eq.shape[0]:
                raise ValueError(f"A_eq and b_eq have inconsistent number of rows. A_eq has {A_eq.shape[0]}, b_eq has {b_eq.shape[0]}.")
            self.num_eq_constraints = A_eq.shape[0]
        else:
            self.num_eq_constraints = 0
        
        if len(A_lb_list) > 0:
            A_lb = np.vstack(A_lb_list)  # lowerbound coefficient matrix
            b_lb = np.concatenate(b_lb_list)  # lowerbound limits vector
            if A_lb.shape[0] != b_lb.shape[0]:
                raise ValueError(f"A_lb and b_lb have inconsistent number of rows. A_lb has {A_lb.shape[0]}, b_lb has {b_lb.shape[0]}.")
            self.num_lb_constraints = A_lb.shape[0]
        else:
            self.num_lb_constraints = 0

        if wordy > 0:
            print(f"Final constraint matrices built with {self.num_ub_constraints} upperbound constraints, {self.num_eq_constraints} equality constraints, and {self.num_lb_constraints} lowerbound constraints.")

        # Build constraint objects if they exist
        constraints = []
        if self.num_ub_constraints > 0:
            constraints.append(optimize.LinearConstraint(A = A_ub, ub = b_ub))
        if self.num_eq_constraints > 0:
            constraints.append(optimize.LinearConstraint(A = A_eq, lb = b_eq, ub = b_eq)) # equality constraints have same lower and upper bounds (thuis equality)
        if self.num_lb_constraints > 0:
            constraints.append(optimize.LinearConstraint(A = A_lb, lb = b_lb))

        # --- Save the optimization problem parameters for later use ---
        self.values = values
        self.constraints = constraints
        self.integrality = integrality
        self.bounds = bounds

        if wordy > 0:
            print("Optimizer set up complete.")

    def optimize(self, threads = -1):
        '''
        Run the optimizer

        Inputs
        ------
        threads : int, None
            Number of threads to use (<0 or None to auto-detect).

        Returns
        -------
        None
        '''

        # --- set up the optimizer ---
        # if the optimizer has not been set up yet, set it up
        if not hasattr(self, 'values') or not hasattr(self, 'constraints') or not hasattr(self, 'integrality') or not hasattr(self, 'bounds'):
            self.set_up_optimizer()

        if wordy > 0:
            print("Starting optimization...")

        # --- Check for valid inputs ---
        if not isinstance(threads, int) and threads is not None:
            raise ValueError("threads must be an integer or None.")

        # detect max number of threads on system if requested for passing into solver
        if threads < 0 or threads is None:
            threads = os.cpu_count()
            if threads is None:
                raise ValueError("Could not detect number of CPU threads on system.")

        # --- call the solver ---
        res = milp(
            c=self.values,                # milp function doesnt not care about the shape of values, just that it is a 1D array 
            constraints=self.constraints,
            integrality=self.integrality, # milp function doesnt not care about the shape of values, just that it is a 1D array 
            bounds=self.bounds
        )

        if wordy > 0:
            print("Solver complete. Analyzing results...")
            print("Results: \n", res)

        # --- process the results ---
        if res.success:
            # Reshape the flat result back into the (num_periods, num_tasks, num_assets) shape

            if wordy > 0:
                print("Decision variable [periods][tasks][assets]:")
                for i in range(len(self.x_indices)):
                    print(f"  {self.x_indices[i]}: {res.x[i]}")

            X_optimal = res.x.reshape((self.num_periods, self.num_tasks, self.num_assets))
            self.schedule = X_optimal
            if wordy > 0:
                print("Optimization successful. The following schedule was generated:")
                for p in range(self.num_periods):
                    print_string = f"Period {p+1}:"
                    whitespace = " " * (3 - len(str(p+1)))  # adjust spacing for single vs double digit periods. Limited to 99 periods.
                    print_string += whitespace
                    displayed = False
                    for t in range(self.num_tasks):
                        for a in range(self.num_assets):
                            if X_optimal[p, t, a] == 1:
                                task_name = self.tasks[t]
                                asset_name = self.assets[a]['name'] if 'name' in self.assets[a] else f"Asset {a+1}"
                                cost = self.task_asset_matrix[t, a, 0]
                                duration = self.task_asset_matrix[t, a, 1]
                                print_string += f"Task '{task_name}' assigned to Asset '{asset_name}' (Cost: {cost}, Duration: {duration})"
                                displayed = True
                                break # break because we only support one asset per task per period for now
                            else:
                                print_string += "No assignment"
                                displayed = True
                                break # break because we only support one asset per task per period for now
                        if displayed:
                            break

                    print(print_string)

        if wordy > 0:
            print("Optimization complete.")
        

if __name__ == "__main__":

    # A dummy system to test the scheduler

    # 21 weather periods = 21 time periods
    weather = [3]*5 + [2]*1 + [1]*10 + [2]*3 + [3]*2 # Three weather types for now. Example weather windows. The length of each window is equal to min_duration
    
    # Example tasks, assets, dependencies, and task_asset_matrix. Eventually try with more tasks than assets, more assets than tasks, etc.
    tasks = [
        "task1",
        "task2",
        "task3"
    ]
    assets = [
        {"name": "asset1", "max_weather" : 3},
        {"name": "asset2", "max_weather" : 2},
        {"name": "asset3", "max_weather" : 1}
    ]

    # task dependencies
    task_dependencies = {
        "task1": [],
        "task2": ["task1_start"],
        "task3": ["task1_end", "task2_start"]
    }

    # cost and duration tuples for each task-asset pair. -1 indicates asset-task paring is invalid
    task_asset_matrix = np.array([
        [(1000, 5), (2000, 3), (1500, 4)],    # task 1: asset 1, asset 2, asset 3
        [(1200, 4), (  -1,-1), (1800, 3)],    # task 2: asset 1, asset 2, asset 3
        [(1100, 6), (2100, 4), (1600, 5)]     # task 3: asset 1, asset 2, asset 3
    ])

    # Find the minimum time period duration based on the task_asset_matrix
    min_duration = np.min(task_asset_matrix[:, :, 1][task_asset_matrix[:, :, 1] > 0])  # minimum non-zero duration

    # Sandbox for building out the scheduler
    scheduler = Scheduler(task_asset_matrix, tasks, assets, task_dependencies, weather, min_duration)
    scheduler.optimize()