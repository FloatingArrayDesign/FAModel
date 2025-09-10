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
- [] Consolidate the loops in the constraints building section
- [] Figure out how to determine which constraint is violated if the problem is infeasible
- [] Add testing
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
            Must be len(tasks) x len(assets) x 2. NOTE: The duration must be in units of scheduling periods (same as weather period length).
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
        self.period_duration = period_duration # duration of each scheduling period. Used for converting from periods to real time.

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
        
        # --- Process inputs ---

        self.num_tasks = len(self.tasks)
        self.num_assets = len(self.assets)
        self.num_periods = len(weather)  # number of scheduling periods

        # Checks for negative duration and cost in task_asset_matrix (0 cost or duration permitted)
        self.num_valid_ta_pairs = int(np.sum((self.task_asset_matrix[:,:,0] >=0) | (self.task_asset_matrix[:,:,1] >= 0))) # number of valid task-asset pairs (cost*duration >= 0)

        # --- Debug helpers ---
        # make a list of indices to help with building constraints
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
        integrality = np.ones(decision_vars.shape) # x_i are int. So set integrality to 1

        if wordy > 0:
            print("Bounds and integrality for decision variables set. Begining to build constraints...")

        # --- build the constraints ---

        # A_ub * x <= b_ub
        # A_eq * x == b_eq
        # A_lb * x >= b_lb 

        '''
        A note on constraints: There are two constraint matrices, the equality constraints (A_eq, b_eq) and the upper bound constraints (A_ub, b_ub).
        Each row in the coefficient matrices corresponds to a constraint, and each column corresponds to a decision variable. Thus the number of columns 
        is equal to the number of decision variables (P * T * A), and the number of rows is equal to the number of constraints.
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
            x_{p,t,a} for p in 0:P, t in 0:T, a in 0:A

        Constraints column (decision variable) indexing used in definitions below:
        x = x_000, ..., x_00A, x_001, ..., x_00A, ..., x_0T1, ..., x_0TA,  # task asset pairings in period 0
            x_100, ..., x_10A, x_101, ..., x_10A, ..., x_1T1, ..., x_1TA,  # task asset pairings in period 1
            ...,
            x_P00, ..., x_P0A, x_P21, ..., x_P2A, ..., x_PT1, ..., x_PTA   # task asset pairings in period p

        where:
        - p is the period index (0 to P), t is the task index (0 to T), and a is the asset index (0 to A)
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
        (x_000 + ... + x_pta) <= P
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

        (x_0jk + ... + x_pjk) = 0   # for all tasks j in range(0:t) and assets k in range(0:a) where task_asset_matrix[j, k, goal_index] <= 0
        '''
        
        mask = np.zeros(X.shape, dtype=int)

        for p in range(self.num_periods):
            mask[p,:,:] = self.task_asset_matrix[:, :, goal_index] <= 0  # Create a mask of invalid task-asset pairings where cost is negative (indicating invalid)

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
        # '''
        # This enforces task dependencies by ensuring that a task can only be assigned to a time period if all its dependencies have been completed in previous periods.
        # TODO: right now this doesnt necessarily enforce that the dependency task has been completed, just that it was assigned in a previous period. This will need to change when 
        # tasks are assigned to the number of periods their duration is.

        # The period of the dependency task must be less than the period of the current task:
        # p * x_pra  < p * x_pta    # for each task t in tasks, for all dependencies r in task_dependencies[t], for all assets a in assets, for all periods p in periods
        # '''
        
        # # Number of task dependency rows
        # A_lb_2 = np.zeros((len(self.task_dependencies), len(decision_vars)), dtype=int)
        # b_lb_2 = np.zeros(len(self.task_dependencies), dtype=int)

        # # Extract dependencies from task_dependencies dict
        # index = 0
        # for task, deps in self.task_dependencies.items():
            
        #     # if the task has dependencies, build the constraint
        #     if len(deps) > 0:
                
        #         # get task index by matching task name and index in self.tasks
        #         if task not in self.tasks:
        #             raise ValueError(f"Task '{task}' in task_dependencies not found in tasks list.")
                
        #         t = self.tasks.index(task)
    
        #         mask = np.zeros(X.shape, dtype=int)

        #         for dep in deps:
        #             # get task index by matching dependency name and index in self.tasks
        #             if dep not in self.tasks:
        #                 raise ValueError(f"Dependency task '{dep}' for task '{task}' not found in tasks list.")
        #             r = self.tasks.index(dep) # get index of dependency
                    
        #             # TODO: need to figure out how to enforce / track temporal ordering of tasks
                    
        #     A_lb_2[index, :] = mask.flatten()
        #     index += 1

        # if wordy > 2:
        #     print("A_lb_2^T:")
        #     print("            T1   T2  ") # Header for 2 tasks
        #     for i in range(A_lb_2.transpose().shape[0]):
        #         pstring = str(self.x_indices[i])
        #         for column in A_lb_2.transpose()[i]:
        #             pstring += f"{ column:5}"
        #         print(pstring)
        #     print("b_lb_2: ", b_lb_2)

        # if wordy > 0:
        #     print("Constraint 2 built.")

        # 3) assets cannot be assigned in a time period where the weather is above the maximum capacity
        # TODO: weather is disabled until this is added

        # 4) assets cannot be assigned to multiple tasks in the same time period
        '''
        this is a simplification, eventually we want to allow multiple assets per task, or parallel tasks
        Sum of all asset-period pairs must be <= 1:
        
        (x_000 + ... + x_pt0) <= 1   # for asset 0
        (x_001 + ... + x_pt1) <= 1   # for asset 1
        ...
        (x_00a + ... + x_pta) <= 1   # for asset t
        '''
        
        # num-periods * num_assets rows
        A_ub_4 = np.zeros((self.num_periods * self.num_assets, len(decision_vars)), dtype=int)
        b_ub_4 = np.ones(self.num_periods * self.num_assets, dtype=int)  # right-hand side is 1 for each asset

        index = 0
        for p in range(self.num_periods):
            for a in range(self.num_assets):
                # Create a mask for all variables for asset a
                mask = np.zeros(X.shape, dtype=int)
                mask[p, :, a] = 1             # Set all periods and tasks for asset a to 1 (so they are included in the sum)

                A_ub_4[index, :] = mask.flatten()
                index += 1

        if wordy > 2:
            print("A_ub_4^T:")
            print("            P1A1 P1A2 P2A1") # Header for 2 tasks and 2 assets example with T2A2 invalid
            for i in range(A_ub_4.transpose().shape[0]):
                pstring = str(self.x_indices[i])
                for column in A_ub_4.transpose()[i]:
                    pstring += f"{ column:5}"
                print(pstring)
            print("b_ub_4: ", b_ub_4)
 
        A_ub_list.append(A_ub_4)
        b_ub_list.append(b_ub_4)

        if wordy > 0:
            print("Constraint 4 built.")

        # 5) The total number of tasks assigned cannot be greater than the number of tasks available (NOTE: Is this necessary or is it already enforced by the fact that there t = number of tasks?)
        # TODO: enforce task limits

        # 6) The total number of assets assigned cannot be greater than the number of assets available (NOTE: Is this necessary or is it already enforced by the fact that there a = number of assets?)
        # TODO: enforce asset limits

        # 7) Ensure tasks are assigned as early as possible 
        '''
        A task cannot be assigned if it could have been assigned in an earlier period. This encourages the solver to assign tasks to the earliest possible periods.
        '''
        # TODO: implement this constraint

        # 8) All tasks must be assigned to at least one time period
        '''

        The sum of all decision variables for each task must be greater than 1, indicating all tasks were assigned at least once:

        (x_000 + ... + x_p0a) >= 1   # for task 0
        (x_010 + ... + x_p1a) >= 1   # for task 1
        ...
        (x_0t0 + ... + x_pta) >= 1   # for task t
        '''

        # num_tasks rows
        A_lb_8 = np.zeros((self.num_tasks, len(decision_vars)), dtype=int)
        b_lb_8 = np.ones(self.num_tasks, dtype=int)

        for t in range(self.num_tasks):
            # Create a mask for all variables for task t
            mask = np.zeros(X.shape, dtype=int)
            mask[:, t, :] = 1             # Set all periods and assets for task t to 1 (so they are included in the sum)
            A_lb_8[t, :] = mask.flatten()

            if wordy > 2:
                print(f"Task {t}:")
                for i in range(len(self.x_indices)):
                    print(f"  {self.x_indices[i]}: {A_lb_8[t, i]}")

        A_lb_list.append(A_lb_8)
        b_lb_list.append(b_lb_8)

        if wordy > 0:
            print("Constraint 8 built.")

        # # 9) A task must be assigned to the continuous number of time periods equal to its duration for the asset assigned to it
        # '''
        # This ensures the duration of a task-asset pair is respected. If a task has a duration of 3 periods, it must be assigned to 3 consecutive periods.

        # (x_ijk + x_(i+1)jk + ... + x_(i+d-1)jk) >= d   # for all tasks j in range(0:t) and assets k in range(0:a) where d is the duration of task j with asset k, and i is the period index. This formulation does not allow for multiple assets per task

        # '''

        # # num task-asset pairings rows
        # A_eq_9 = np.zeros((self.num_valid_ta_pairs, len(decision_vars)), dtype=int)
        # b_eq_9 = np.zeros(self.num_valid_ta_pairs, dtype=int)

        # # Loop through tasks and assets
        # pair_i = 0
        # for t in range(self.num_tasks):
        #     for a in range(self.num_assets):

        #         duration = self.task_asset_matrix[t, a, 1]  # duration of task t with asset a
        #         if duration > 0: # If valid pairing, make constraint
                    
        #             # Create a mask for all variables for task t and asset a
        #             mask = np.zeros(X.shape, dtype=int)
        #             for p in range(self.num_periods):
        #                 mask[p:p+duration, t, a] = 1
                    
        #             b_eq_9[pair_i] = self.task_asset_matrix[t, a, 1] # Duration
        #             A_eq_9[pair_i, :] = mask.flatten()
        #             pair_i += 1

        # if wordy > 0:
        #     # Print out the constraint matrix for debugging 
        #     print("A_eq_9^T:")
        #     print("            T1A1 T1A2 T2A1") # Header for 2 tasks and 2 assets example with T2A2 invalid
        #     for i in range(A_eq_9.transpose().shape[0]):
        #         pstring = str(self.x_indices[i])
        #         for column in A_eq_9.transpose()[i]:
        #             pstring += f"{ column:5}"
        #         print(pstring)
        #     print("b_eq_9: ", b_eq_9)

        # A_eq_list.append(A_eq_9)
        # b_eq_list.append(b_eq_9)

        # if wordy > 0:
        #     print("Constraint 9 built.")

        # 10) A task duration plus the first time period it is assigned to must be less than the total number of time periods available
        '''
        This ensures that a task is not assigned to a period that would cause it to exceed the total number of periods available.

        (p * x_{p,t,a} + d_{t,a} * x_{p,t,a}) <= P   # for all t in 0..T, a in 0..A, p in 0..P
        '''
        
        # num_periods rows
        A_ub_10 = np.zeros((self.num_periods, len(decision_vars)), dtype=int)
        b_ub_10 = np.ones(self.num_periods, dtype=int) * self.num_periods  

        for p in range(self.num_periods):
            # Create a mask for the period
            mask = np.zeros(X.shape, dtype=int)
            
            # Loop through pairs
            for t in range(self.num_tasks):
                for a in range(self.num_assets):
                    duration = self.task_asset_matrix[t, a, 1]  # duration of task t with asset a
                    if duration > 0:
                        mask[p, t, a] = p + duration  # Set the specific variable to i + d_jk

            A_ub_10[p, :] = mask.flatten()
            
            if wordy > 2:
                print(f"Period {p}:")
                for i in range(len(self.x_indices)):
                    print(f"  {self.x_indices[i]}: {A_ub_10[p, i]}")
                print("Upper bound limit: ", b_ub_10[p])

        A_ub_list.append(A_ub_10)
        b_ub_list.append(b_ub_10)

        if wordy > 0:
            print("Constraint 10 built.")

        # --- End Constraints ---

        if wordy > 0:
            print("All constraints built. Stacking and checking constraints...")

        # --- Assemble the SciPy Constraints ---
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

            if wordy > 1:
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
                    for t in range(self.num_tasks):
                        for a in range(self.num_assets):
                            if X_optimal[p, t, a] == 1:
                                task_name = self.tasks[t]
                                asset_name = self.assets[a]['name'] if 'name' in self.assets[a] else f"Asset {a+1}"
                                cost = self.task_asset_matrix[t, a, 0]
                                duration = self.task_asset_matrix[t, a, 1]
                                print_string += f"Asset '{asset_name}' assigned to Task '{task_name}' (Cost: {cost}, Duration: {duration})"

                    print(print_string)

        if wordy > 0:
            print("Optimization function complete.")
        

if __name__ == "__main__":

    os.system("clear") # for clearing terminal on Mac

    # A simple dummy system to test the scheduler

    # 10 weather periods = 10 time periods
    weather = [1]*5 # Three weather types for now. Example weather windows. The length of each window is equal to min_duration
    
    # Example tasks, assets, dependencies, and task_asset_matrix. Eventually try with more tasks than assets, more assets than tasks, etc.
    tasks = [
        "task1",
        "task2"
    ]
    assets = [
        {"name": "asset1", "max_weather" : 3},
        {"name": "asset2", "max_weather" : 2}
    ]

    # task dependencies
    task_dependencies = {
        "task1": ["task1"],
        "task2": []
    }

    # cost and duration tuples for each task-asset pair. -1 indicates asset-task paring is invalid
    task_asset_matrix = np.array([
        [(1000, 2), (2000, 3)],    # task 1: asset 1, asset 2
        [(1200, 5), (  -1,-1)]     # task 2: asset 1, asset 2
    ])

    # optimal assignment: task 1 with asset 1 in periods 1-2, task 2 with asset 1 in period 3

    # Find the minimum time period duration based on the task_asset_matrix
    min_duration = np.min(task_asset_matrix[:, :, 1][task_asset_matrix[:, :, 1] > 0])  # minimum non-zero duration

    # Sandbox for building out the scheduler
    scheduler = Scheduler(task_asset_matrix, tasks, assets, task_dependencies, weather, min_duration)
    scheduler.optimize()

    
    
    # # A more complex dummy system to test the scheduler (uncomment and comment out above to run)

    # # 10 weather periods = 10 time periods
    # weather = [1]*5 + [2]*1 + [3]*1 + [1]*3 # Three weather types for now. Example weather windows. The length of each window is equal to min_duration
    
    # # Example tasks, assets, dependencies, and task_asset_matrix. Eventually try with more tasks than assets, more assets than tasks, etc.
    # tasks = [
    #     "task1",
    #     "task2",
    #     "task3"
    # ]
    # assets = [
    #     {"name": "asset1", "max_weather" : 3},
    #     {"name": "asset2", "max_weather" : 2},
    #     {"name": "asset3", "max_weather" : 1},
    #     {"name": "asset4", "max_weather" : 1}
    # ]

    # # task dependencies
    # task_dependencies = {
    #     "task1": [],
    #     "task2": ["task1"],
    #     "task3": ["task1"]
    # }

    # # random cost and duration tuples for each task-asset pair. -1 indicates asset-task paring is invalid
    # task_asset_matrix = np.array([
    #     [(3000, 2), (2000, 3), (1000, 4), (4000, 5)],    # task 1: asset 1, asset 2, asset 3, asset 4
    #     [(1200, 5), (  -1,-1), (  -1,-1), (  -1,-1)],    # task 2: asset 1, asset 2, asset 3, asset 4
    #     [(2500, 3), (1500, 2), (  -1,-1), (  -1,-1)]     # task 3: asset 1, asset 2, asset 3, asset 4
    # ])

    # # optimal assignment: task 1 with asset 1 in periods 1-2, task 2 with asset 1 in period 3

    # # Find the minimum time period duration based on the task_asset_matrix
    # min_duration = np.min(task_asset_matrix[:, :, 1][task_asset_matrix[:, :, 1] > 0])  # minimum non-zero duration

    # # Sandbox for building out the scheduler
    # scheduler = Scheduler(task_asset_matrix, tasks, assets, task_dependencies, weather, min_duration)
    # scheduler.optimize()