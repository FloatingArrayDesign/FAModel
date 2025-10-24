# author: @rdavies, 9-8-2025

# Scheduler class for managing actions and tasks
# WIP, to be merged into Irma later on

'''
--- TODO List ---
- [] How to expand this to multiple assets per task?
- [] Eventually enable parallel tasks and multiple assets per task
- [] Convert input tasks and assets from dicts to Task and Asset objects
 - [] When tasks and assets are converted from lists to objects, update the type hints for task and asset list at class initialization.
- [] Add a delay cost, i.e. a cost for each time period where X = 0 <-- do we want this? might not be needed
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
    def __init__(self, task_asset_matrix : np.ndarray, tasks : list[str] = [], assets : list[dict] = [], task_dependencies = {}, dependency_types = {}, weather : list[int] = [], period_duration : float = 1, **kwargs):
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
        dependency_types : dict
            A dictionary mapping each task dependency pair to its type:
            - "finish_start" (default): dependent task starts after prerequisite finishes
            - "start_start": dependent task starts when prerequisite starts  
            - "finish_finish": dependent task finishes when prerequisite finishes
            - "start_finish": dependent task finishes when prerequisite starts
            - "offset": dependent task starts/finishes with time offset (requires offset value)
            - "same_asset": dependent task must use same asset as prerequisite
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
        self.dependency_types = dependency_types
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

        self.T = len(self.tasks)
        self.A = len(self.assets)
        self.P = len(weather)  # number of scheduling periods
        self.S = self.P        # number of start times

        # Checks for negative duration and cost in task_asset_matrix (0 cost and duration permitted)
        self.num_valid_ta_pairs = int(np.sum((self.task_asset_matrix[:,:,0] >=0) & (self.task_asset_matrix[:,:,1] >= 0))) # number of valid task-asset pairs (cost and duration >= 0)

        # --- Debug helpers ---
        # make a list of indices to help with building constraints
        self.Xta_indices = [f"Xta_[{t}][{a}]" for t in range(self.T) for a in range(self.A)]
        self.Xtp_indices = [f"Xtp_[{t}][{p}]" for t in range(self.T) for p in range(self.P)]
        self.Xap_indices = [f"Xap_[{a}][{p}]" for a in range(self.A) for p in range(self.P)]
        self.Xts_indices = [f"Xts_[{t}][{s}]" for t in range(self.T) for s in range(self.S)]
        self.X_indices = self.Xta_indices + self.Xtp_indices + self.Xap_indices + self.Xts_indices

        if wordy > 0:
            print(f"Scheduler initialized with {self.P} time periods, {self.T} tasks, {self.A} assets, and {self.S} start times.")

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

        # Decision variables:
        # Xta = task asset pairs
        # Xtp = task period pairs
        # Xap = period asset pairs
        # Xts = task start-time pairs
        num_variables = (self.T * self.A) + (self.T * self.P) + (self.A * self.P) + (self.T * self.S)  # number of decision variables

        self.Xta_start = 0  # starting index of Xta in the flattened decision variable vector
        self.Xta_end = self.Xta_start + self.T * self.A  # ending index of Xta in the flattened decision variable vector
        self.Xtp_start = self.Xta_end  # starting index of Xtp in the flattened decision variable vector
        self.Xtp_end = self.Xtp_start + self.T * self.P  # ending index of Xtp in the flattened decision variable vector
        self.Xap_start = self.Xtp_end  # starting index of Xap in the flattened decision variable vector
        self.Xap_end = self.Xap_start + self.A * self.P  # ending index of Xap in the flattened decision variable vector
        self.Xts_start = self.Xap_end  # starting index of Xts in the flattened decision variable vector
        self.Xts_end = self.Xts_start + self.T * self.S  # ending index of Xts in the flattened decision variable vector

        # Values vector: In every planning period, the value of assigning asset a to task t is the same. Constraints determine which periods are chosen.
            # Note: Intentionally using values here instead of "cost" to avoid confusion between the program 'cost' of a pairing (which could be financial cost, duration, or some other target metric for minimization) to the solver and the financial cost of a asset-task pairing.
        values = np.zeros(num_variables, dtype=int) # NOTE: enforces discrete cost and duration
        values[self.Xta_start:self.Xta_end] = self.task_asset_matrix[:, :, goal_index].flatten()  # Set the cost or duration for the task-asset pair

        # Add small penalties for later start times (Constraint 7 implementation)
        # This encourages the solver to choose earlier start times when possible
        max_task_cost = np.max(self.task_asset_matrix[:, :, goal_index])
        early_start_penalty_factor = max_task_cost * 0.001  # Very small penalty (0.1% of max cost)
        
        for t in range(self.T):
            for s in range(self.S):
                # Add small penalty proportional to start time
                # Later start times get higher penalties
                penalty = int(early_start_penalty_factor * s)
                values[self.Xts_start + t * self.S + s] = penalty

        # The rest of values (for period variables) remains zero because they do not impact cost or duration

        if wordy > 1:
            print("Values vector of length " + str(values.shape[0]) + " created")

        # lb <= x <= ub
            # Constrain decision variables to be 0 or 1
        bounds = optimize.Bounds(0, 1)  # 0 <= x_i <= 1
        integrality = np.ones(num_variables, dtype=int) # x_i are int. So set integrality to 1

        if wordy > 0:
            print("Bounds and integrality for decision variables set. Begining to build constraints...")

        # --- build the constraints ---

        # A_ub * x <= b_ub
        # A_eq * x == b_eq
        # A_lb * x >= b_lb 

        '''
        A note on constraints: There are two constraint matrices, the equality constraints (A_eq, b_eq) and the upper bound constraints (A_ub, b_ub).
        Each row in the coefficient matrices corresponds to a constraint, and each column corresponds to a decision variable. Thus the number of columns
        is equal to the number of decision variables (T*A + T*P + T*S), and the number of rows is equal to the number of constraints.
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

        The indexing of decision variables are:
            Xta = [Xta_00, ..., Xta_0A, Xta_10, ..., Xta_1A, ..., Xta_T0, ..., Xta_TA]    # task asset pairs
            Xtp = [Xtp_00, ..., Xtp_0P, Xtp_10, ..., Xtp_1P, ..., Xtp_T0, ..., Xtp_TP]    # task period pairs
            Xap = [Xap_00, ..., Xap_0P, Xap_10, ..., Xap_1P, ..., Xap_A0, ..., Xap_AP]    # asset period pairs
            Xts = [Xts_00, ..., Xts_0S, Xts_10, ..., Xts_1S, ..., Xts_T0, ..., Xts_TS]    # task start-time pairs

        The global decision variable is then:
            X = [Xta, Xtp, Xap, Xts]

            The starting indices of each section in this global variable are saved as self.Xta_start, self.Xtp_start, self.Xap_start, and self.Xts_start.
            While the values vector is only nonzero for self.Xta_start:self.Xta_end, the constraints will leverage all decision variables.

        where:
        -  t is the task index (0 to T), a is the asset index (0 to A), p is the period index (0 to P), and s is the start time index (0 to S)
        '''

        # Empty list of constraint coefficient matrices
        A_ub_list = []
        A_eq_list = []
        A_lb_list = []

        # Empty list of constraint limit vectors
        b_ub_list = []
        b_eq_list = []
        b_lb_list = []

        # 1) asset can only be assigned to a task if asset is capable of performing the task (value of pairing is non-negative)
        '''
        if task t cannot be performed by asset a, then Xta_ta = 0 

        (Xta_00 + ... + Xta_TA) = 0   # for all tasks t in range(0:T) and assets a in range(0:A) where task_asset_matrix[t, a, goal_index] <= 0
        '''

        # 1 row
        rows = []
        for t in range(self.T):
            for a in range(self.A):
                if self.task_asset_matrix[t, a, goal_index] <= 0:  # Invalid pairing
                    row = np.zeros(num_variables, dtype=int)
                    row[self.Xta_start + t * self.A + a] = 1
                    rows.append(row)
        
        if rows:  # Only create constraint if there are invalid pairings
            A_eq_1 = np.vstack(rows)
            b_eq_1 = np.zeros(A_eq_1.shape[0], dtype=int)
            
            if wordy > 1:
                print("A_eq_1^T:")
                for i in range(self.Xta_start,self.Xta_end):
                    pstring = str(self.X_indices[i])
                    for column in A_eq_1.transpose()[i]:
                        pstring += f"{ column:5}"
                    print(pstring)
                print("b_eq_1: ", b_eq_1)

            A_eq_list.append(A_eq_1)
            b_eq_list.append(b_eq_1)

        if wordy > 0:
            print("Constraint 1 built.")

        # 2) task dependencies must be respected (i.e., a task cannot start until all its dependencies have been satisfied)
        '''
        This enforces task dependencies by ensuring that a task can only be assigned to a time period if all its dependencies have been completed in previous periods.
        
        Different dependency types:
        - finish_start: Task B starts after Task A finishes  
        - start_start: Task B starts when Task A starts
        - finish_finish: Task B finishes when Task A finishes  
        - start_finish: Task B finishes when Task A starts
        - same_asset: Task B must use the same asset as Task A
        
        For finish_start dependencies (most common):
        If task t depends on task d, then task t cannot start before task d finishes.
        
        Using start times: Xts[t,s] = 1 implies Xts[d,sd] = 1 where sd + duration_d <= s
        
        Constraint: For all valid start times s for task t, if Xts[t,s] = 1, 
        then there must exist some start time sd for task d such that Xts[d,sd] = 1 
        and sd + duration_d <= s
        
        Implementation: Xts[t,s] <= sum(Xts[d,sd] for sd where sd + duration_d <= s)
        '''
        
        rows_2 = []
        vec_2 = []
        
        # Convert task names to indices for easier processing
        task_name_to_index = {task: i for i, task in enumerate(self.tasks)}
        
        for task_name, dependencies in self.task_dependencies.items():
            if task_name not in task_name_to_index:
                continue  # Skip if task not in our task list
                
            t = task_name_to_index[task_name]  # dependent task index
            
            for dep_task_name in dependencies:
                if dep_task_name not in task_name_to_index:
                    continue  # Skip if dependency not in our task list
                    
                d = task_name_to_index[dep_task_name]  # dependency task index
                
                # Get dependency type (default to finish_start)
                dep_key = f"{dep_task_name}->{task_name}"
                dep_type = self.dependency_types.get(dep_key, "finish_start")
                
                if dep_type == "finish_start":
                    # Task t cannot start until task d finishes
                    # For each possible start time s of task t
                    for s in range(self.S):
                        # Task t can start at time s only if task d has already finished
                        # Find minimum duration of task d across all possible assets
                        min_duration_d = float('inf')
                        for a_d in range(self.A):
                            duration_d = self.task_asset_matrix[d, a_d, 1]
                            if duration_d > 0:  # Valid task-asset pairing
                                min_duration_d = min(min_duration_d, duration_d)
                        
                        if min_duration_d == float('inf'):
                            continue  # No valid asset for dependency task
                        
                        # Task d must finish before time s
                        # So task d must start at latest at time (s - min_duration_d)
                        # But we need to account for the actual duration based on asset choice
                        
                        # For this constraint: if task t starts at time s, then task d must have started 
                        # and finished before time s
                        latest_start_d = s - min_duration_d
                        
                        if latest_start_d < 0:
                            # Task t cannot start at time s because task d cannot finish in time
                            row = np.zeros(num_variables, dtype=int)
                            row[self.Xts_start + t * self.S + s] = 1  # Xts[t,s] = 0 (cannot start)
                            rows_2.append(row)
                            vec_2.append(0)  # Xts[t,s] <= 0, so Xts[t,s] = 0
                        else:
                            # Task t can start at time s only if task d starts at time <= latest_start_d
                            row = np.zeros(num_variables, dtype=int)
                            row[self.Xts_start + t * self.S + s] = 1  # Xts[t,s]
                            
                            # Add all valid start times for task d
                            has_valid_dep_start = False
                            for sd in range(min(latest_start_d + 1, self.S)):  # sd from 0 to latest_start_d
                                row[self.Xts_start + d * self.S + sd] = -1  # -Xts[d,sd]
                                has_valid_dep_start = True
                            
                            if has_valid_dep_start:
                                rows_2.append(row)
                                vec_2.append(0)  # Xts[t,s] - sum(Xts[d,valid_sd]) <= 0
                
                elif dep_type == "start_start":
                    # Task t starts when task d starts (same start time)
                    for s in range(self.S):
                        row = np.zeros(num_variables, dtype=int)
                        row[self.Xts_start + t * self.S + s] = 1   # Xts[t,s]
                        row[self.Xts_start + d * self.S + s] = -1  # -Xts[d,s]
                        rows_2.append(row)
                        vec_2.append(0)  # Xts[t,s] - Xts[d,s] <= 0, so Xts[t,s] <= Xts[d,s]
                
                elif dep_type == "finish_finish":
                    # Task t finishes when task d finishes
                    # This requires both tasks to have the same end time
                    for s_t in range(self.S):
                        for a_t in range(self.A):
                            duration_t = self.task_asset_matrix[t, a_t, 1]
                            if duration_t > 0:  # Valid pairing for task t
                                end_time_t = s_t + duration_t
                                
                                # Find start times for task d that result in same end time
                                for s_d in range(self.S):
                                    for a_d in range(self.A):
                                        duration_d = self.task_asset_matrix[d, a_d, 1]
                                        if duration_d > 0:  # Valid pairing for task d
                                            end_time_d = s_d + duration_d
                                            
                                            if end_time_t == end_time_d:
                                                # If task t starts at s_t with asset a_t AND task d starts at s_d with asset a_d,
                                                # then they finish at the same time (constraint satisfied)
                                                continue
                                            else:
                                                # Prevent this combination
                                                row = np.zeros(num_variables, dtype=int)
                                                row[self.Xts_start + t * self.S + s_t] = 1     # Xts[t,s_t] 
                                                row[self.Xta_start + t * self.A + a_t] = 1     # Xta[t,a_t]
                                                row[self.Xts_start + d * self.S + s_d] = 1     # Xts[d,s_d]
                                                row[self.Xta_start + d * self.A + a_d] = 1     # Xta[d,a_d]
                                                rows_2.append(row)
                                                vec_2.append(3)  # At most 3 of these 4 can be 1 simultaneously
                
                elif dep_type == "same_asset":
                    # Task t must use the same asset as task d
                    for a in range(self.A):
                        # If both tasks can use asset a
                        if (self.task_asset_matrix[t, a, 1] > 0 and 
                            self.task_asset_matrix[d, a, 1] > 0):
                            row = np.zeros(num_variables, dtype=int)
                            row[self.Xta_start + t * self.A + a] = 1   # Xta[t,a]
                            row[self.Xta_start + d * self.A + a] = -1  # -Xta[d,a]
                            rows_2.append(row)
                            vec_2.append(0)  # Xta[t,a] - Xta[d,a] <= 0, so if t uses a, then d must use a
        
        # Build constraint matrices if we have any dependency constraints
        if rows_2:
            A_ub_2 = np.vstack(rows_2)
            b_ub_2 = np.array(vec_2, dtype=int)
            A_ub_list.append(A_ub_2)
            b_ub_list.append(b_ub_2)

        if wordy > 0:
            print("Constraint 2 built.")

        # 3) at least one asset must be assigned to each task
        '''
        Sum of all task-asset pairs must be >= 1 for each task:
        (Xta_00 + ... + Xta_0A) >= 1   # for task 0
        (Xta_10 + ... + Xta_1A) >= 1   # for task 1
        ...
        (Xta_T0 + ... + Xta_TA) >= 1   # for task T
        '''

        # num_tasks rows
        A_lb_3 = np.zeros((self.T, num_variables), dtype=int)
        b_lb_3 = np.ones(self.T, dtype=int)

        for t in range (self.T):
            # set the coefficient for each task t to one
            A_lb_3[t, (self.Xta_start + t * self.A):(self.Xta_start + t * self.A + self.A)] = 1  # Set the coefficients for the Xta variables to 1 for each task t

        if wordy > 1:
            print("A_lb_3^T:")
            print("             T1   T2") # Header for 2 tasks
            for i in range(self.Xta_start,self.Xta_end):
                pstring = str(self.X_indices[i])
                for column in A_lb_3.transpose()[i]:
                    pstring += f"{ column:5}"
                print(pstring)
            print("b_lb_3: ", b_lb_3)

        A_lb_list.append(A_lb_3)
        b_lb_list.append(b_lb_3)

        if wordy > 0:
            print("Constraint 3 built.")

        # 4) each asset can only be used by one task per time period  
        '''
        Multiple tasks can be assigned to the same asset (Xta[t1,a] = Xta[t2,a] = 1),
        but they cannot use it simultaneously in the same period.
        
        The relationship is enforced through Constraint 12:
        Xtp[t,p] + Xta[t,a] - Xap[a,p] ≤ 1 and ≥ 0
        
        This means: if Xtp[t,p] = 1 AND Xta[t,a] = 1, then Xap[a,p] = 1
        Since Xap[a,p] is binary, it can only be 1 for one reason.
        
        We ensure: Xap[a,p] ≤ 1 for each asset a, period p
        This constraint is automatically satisfied for binary variables, but we include it explicitly.
        
        The key insight: if Constraint 12 is working correctly, it should prevent conflicts
        by ensuring that if multiple tasks are assigned to the same asset and try to be
        active simultaneously, the Xap[a,p] variable relationships will prevent this.
        '''
        
        # Ensure each asset can be active in at most one context per period
        # (This is automatic for binary variables but explicit for clarity)
        A_ub_4 = np.zeros((self.A * self.P, num_variables), dtype=int)
        b_ub_4 = np.ones(self.A * self.P, dtype=int)
        '''
        row = 0
        for a in range(self.A):
            for p in range(self.P):
                A_ub_4[row, self.Xap_start + a * self.P + p] = 1  # Xap[a,p] ≤ 1
                row += 1
        '''
        # Add temporal conflict prevention for tasks that could share assets
        # For each asset that multiple tasks could use, add constraints to prevent
        # simultaneous usage by different tasks
        rows_4b = []
        bounds_4b = []
        
        for a in range(self.A):
            tasks_for_asset = [t for t in range(self.T) if self.task_asset_matrix[t, a, 1] > 0]
            
            if len(tasks_for_asset) > 1:  # Multiple tasks could use this asset
                for p in range(self.P):
                    # Create a constraint involving ALL tasks that could use this asset
                    # Σ(Xtp[t,p] for t in tasks_for_asset) + Σ(Xta[t,a] for t in tasks_for_asset) ≤ bound
                    # Logic: If tasks are assigned to asset a, at most 1 can be active in period p
                    row = np.zeros(num_variables, dtype=int)
                    
                    # Add all task-period variables for this period
                    for t in tasks_for_asset:
                        row[self.Xtp_start + t * self.P + p] = 1  # Xtp[t,p]
                    
                    # Add all task-asset variables for this asset
                    for t in tasks_for_asset:
                        row[self.Xta_start + t * self.A + a] = 1  # Xta[t,a]
                    
                    rows_4b.append(row)
                    
                    # Calculate bound: if all tasks assigned to asset a, max 1 can be active
                    num_tasks = len(tasks_for_asset)
                    max_active_when_all_assigned = 1  # Only 1 task can use asset per period
                    max_assignments = num_tasks  # All could potentially be assigned to the asset
                    bound = max_active_when_all_assigned + max_assignments
                    bounds_4b.append(bound)

        if rows_4b:
            A_ub_4b = np.vstack(rows_4b)
            b_ub_4b = np.array(bounds_4b, dtype=int)
            
            A_ub_4 = np.vstack([A_ub_4, A_ub_4b])
            b_ub_4 = np.concatenate([b_ub_4, b_ub_4b])

        if wordy > 1: 
            print("A_ub_4^T:")
            print("             P1   P2   P3   P4   P5") # Header for 5 periods
            for i in range(self.Xap_start,self.Xap_end):
                pstring = str(self.X_indices[i])
                for column in A_ub_4.transpose()[i]:
                    pstring += f"{ column:5}"
                print(pstring)
            print("b_ub_4: ", b_ub_4)
        
        A_ub_list.append(A_ub_4)
        b_ub_list.append(b_ub_4)

        if wordy > 0:
            print("Constraint 4 built.")

        # 10) A task duration plus the start-time it is assigned to must be less than the total number of time periods available
        '''
        This ensures that a task is not assigned to a period that would cause it to exceed the total number of periods available.

        (Xts * s + d_ta) <= P  # for all tasks t in range(0:T) where d is the duration of task-asset pair ta
        '''

        rows = []
        for t in range(self.T):
            for a in range(self.A):
                duration = self.task_asset_matrix[t, a, 1]  # duration of task t with asset a
                if duration > 0: # If valid pairing, make constraint
                    for s in range(self.S):
                        if s + duration > self.P:
                            row = np.zeros(num_variables, dtype=int)
                            row[self.Xts_start + t * self.S + s] = 1
                            row[self.Xta_start + t * self.A + a] = 1
                            rows.append(row)

        A_ub_10 = np.vstack(rows)
        b_ub_10 = np.ones(A_ub_10.shape[0], dtype=int)  # Each infeasible combination: Xta + Xts <= 1

        if wordy > 1:
            print("A_ub_10^T:")
            print("             T1A1 T1A2 T2A1") # Header for 3 task-asset pairs example with T2A2 invalid
            for i in range(self.Xta_start,self.Xta_end):
                pstring = str(self.X_indices[i])
                for column in A_ub_10.transpose()[i]:
                    pstring += f"{ column:5}"
                print(pstring)
            for i in range(self.Xts_start,self.Xts_end):
                pstring = str(self.X_indices[i])
                for column in A_ub_10.transpose()[i]:
                    pstring += f"{ column:5}"
                print(pstring)
            print("b_ub_10: ", b_ub_10)

        A_ub_list.append(A_ub_10)
        b_ub_list.append(b_ub_10)

        if wordy > 0:
            print("Constraint 10 built.")

        # 11) The total number of task period pairs must be greater than or equal to the number of task-start time pairs
        '''
        This ensures that the task start-time decision variable is non-zero if a task is assigned to any period.

        (Xtp_00 + ... + Xtp_TP) >= (Xts_00 + ... + Xts_TS)   # for all tasks t in range(0:T)
        '''
        """
        A_lb_11 = np.zeros((self.T, num_variables), dtype=int)
        b_lb_11 = np.ones(self.T, dtype=int) * 2
        
        for t in range(self.T):
            A_lb_11[t, (self.Xtp_start + t * self.P):(self.Xtp_start + t * self.P + self.P)] = 1
            A_lb_11[t, (self.Xts_start + t * self.S):(self.Xts_start + t * self.S + self.S)] = 1

        if wordy > 1:
            print("A_lb_11^T:")
            print("             T1   T2") # Header for 2 tasks
            for i in range(self.Xtp_start,self.Xts_end):
                pstring = str(self.X_indices[i])
                for column in A_lb_11.transpose()[i]:
                    pstring += f"{ column:5}"
                print(pstring)
            print("b_lb_11: ", b_lb_11)

        A_lb_list.append(A_lb_11)
        b_lb_list.append(b_lb_11)

        if wordy > 0:
            print("Constraint 11 built.")
        """
        # 12) The period an asset is assigned to must match the period the task in the task-asset pair is assigned to
        '''
        This ensures the chosen task and asset in a task asset pair are assigned to the same period. This means that if an asset 
        is assigned to a task, then the corresponding task-period and asset-period pairs must be equal. 

        if Xta = 1, then Xtp = Xap, else if Xta = 0, then Xtp and Xap can be anything. This requires two constriants:
        
        Xtp[t, p] - Xap[a, p] <= 1 - Xta[t, a]    --> Xtp[t, p] - Xap[a, p] + Xta[t, a] <= 1
        Xtp[t, p] - Xap[a, p] >= -(1 - Xta[t, a]) -->  Xtp[t, p] - Xap[a, p] + Xta[t, a] >= 0

        '''
        """
        A_12 = np.zeros((self.T * self.A * self.P, num_variables), dtype=int)
        b_ub_12 = np.ones(self.T * self.A * self.P, dtype=int)
        b_lb_12 = np.zeros(self.T * self.A * self.P, dtype=int)
        
        row = 0
        for t in range(self.T):
            for a in range(self.A):
                for p in range(self.P):
                    A_12[row, self.Xtp_start + t * self.P + p] = 1
                    A_12[row, self.Xap_start + a * self.P + p] = -1
                    A_12[row, self.Xta_start + t * self.A + a] = 1

                    row += 1
        """
        
        rows_ub = []
        rows_lb = []

        for t in range(self.T):
            for a in range(self.A):
                # Only create constraints for valid task-asset pairs
                if self.task_asset_matrix[t, a, 1] > 0:  # Valid pairing (duration > 0)
                    for p in range(self.P):
                        row = np.zeros(num_variables, dtype=int)
                        row[self.Xtp_start + t * self.P + p] = 1      # Xtp[t,p]
                        row[self.Xap_start + a * self.P + p] = -1     # -Xap[a,p]
                        row[self.Xta_start + t * self.A + a] = 1      # Xta[t,a]
                        
                        rows_ub.append(row.copy())  # Upper bound constraint
                        rows_lb.append(row.copy())  # Lower bound constraint
        """
        if rows_ub:
            A_ub_12 = np.vstack(rows_ub)
            b_ub_12 = np.ones(len(rows_ub), dtype=int)
            A_lb_12 = np.vstack(rows_lb)  
            b_lb_12 = np.zeros(len(rows_lb), dtype=int)
            
            A_ub_list.append(A_ub_12)
            b_ub_list.append(b_ub_12)
            A_lb_list.append(A_lb_12)
            b_lb_list.append(b_lb_12)
        """
        if wordy > 1:
            print("A_12^T:")
            for i in range(self.Xta_start,self.Xap_end):
                pstring = str(self.X_indices[i])
                for column in A_12.transpose()[i]:
                    pstring += f"{ column:5}"
                print(pstring)
            print("b_ub_12: ", b_ub_12)
            print("b_lb_12: ", b_lb_12)
        '''
        A_ub_list.append(A_12)
        b_ub_list.append(b_ub_12)
        A_lb_list.append(A_12)
        b_lb_list.append(b_lb_12)
        '''
        if wordy > 0:
            print("Constraint 12 built.")

        # 14) if a task-starttime pair is selected, the corresponding task-period pair must be selected for the period equal to the start time plus the duration of the task
        '''
        This ensures that if a task is assigned a start time, the corresponding task-period pair for the period equal to the start time plus the duration of the task is also selected.
        Xts[t, s] <= Xtp[t, s : s + d]   # for all tasks t in range(0:T) and start times s in range(0:S) where d is the duration of task t with the asset assigned to it
        '''

        # TODO: commenting out this constraint allows the optimizer to find an optimal solution

        # TODO: this is very very close. The Xtp are being assigned blocks equal to the starttime + duration. But it is causing the optimizer to fail...?
        rows_14a = []
        vec_14a = []
        rows_14b = []
        vec_14b = []
        #rows = []
        #vec = []

        # 14a) Simple start time to period mapping: Xts[t,s] <= Xtp[t,s]
        for t in range(self.T):
            for s in range(self.S):
                if s < self.P:  # Only if start time is within valid periods
                    row = np.zeros(num_variables, dtype=int)
                    row[self.Xts_start + t * self.S + s] = 1    # Xts[t,s]
                    row[self.Xtp_start + t * self.P + s] = -1   # -Xtp[t,s]
                    rows_14a.append(row)
                    vec_14a.append(0)  # Xts[t,s] - Xtp[t,s] <= 0

        '''
        for t in range(self.T):
            for a in range(self.A):
                duration = self.task_asset_matrix[t, a, 1]
                if duration > 0: # If valid pairing, make constraint
                    for s in range(min(self.S, self.P - duration + 1)):
                        row = np.zeros(num_variables, dtype=int)
                        row[self.Xta_start + t * self.A + a] = 1
                        row[self.Xts_start + t * self.S + s] = -1
                        start_idx = self.Xtp_start + t * self.P + s
                        end_idx = min(start_idx + duration, self.Xtp_start + (t + 1) * self.P)
                        row[start_idx:end_idx] = 1
                        #row[self.Xtp_start + t * self.P + s : self.Xtp_start + t * self.P + s + duration] = 1
                        rows.append(row)
                        vec.append(1)
        '''
        # 14b) Duration enforcement: if task t uses asset a and starts at s, 
        #      then it must be active for duration periods
        for t in range(self.T):
            for a in range(self.A):
                duration = self.task_asset_matrix[t, a, 1]
                if duration > 0:  # Valid task-asset pairing
                    for s in range(min(self.S, self.P - duration + 1)):  # Valid start times
                        for p in range(s, min(s + duration, self.P)):  # Each period in duration
                            row = np.zeros(num_variables, dtype=int)
                            # If Xta[t,a] = 1 AND Xts[t,s] = 1, then Xtp[t,p] = 1
                            # Constraint: Xta[t,a] + Xts[t,s] - Xtp[t,p] <= 1
                            row[self.Xta_start + t * self.A + a] = 1      # Xta[t,a]
                            row[self.Xts_start + t * self.S + s] = 1      # Xts[t,s]  
                            row[self.Xtp_start + t * self.P + p] = -1     # -Xtp[t,p]
                            rows_14b.append(row)
                            vec_14b.append(1)  # Xta[t,a] + Xts[t,s] - Xtp[t,p] <= 1

        
        #A_lb_14 = np.vstack(rows)
        #b_lb_14 = np.array(vec, dtype=int)

        if rows_14a:
            A_ub_14a = np.vstack(rows_14a)
            b_ub_14a = np.array(vec_14a, dtype=int)
            A_ub_list.append(A_ub_14a)
            b_ub_list.append(b_ub_14a)

        if rows_14b:
            A_ub_14b = np.vstack(rows_14b)
            b_ub_14b = np.array(vec_14b, dtype=int)
            A_ub_list.append(A_ub_14b)
            b_ub_list.append(b_ub_14b)

        if wordy > 1:
            print("A_lb_14^T:")
            print("           T1A1S1                   T1A2S1 ...") # Header for 3 task-asset pairs example with T2A2 invalid
            for i in range(self.Xta_start,self.Xta_end):
                pstring = str(self.X_indices[i])
                for column in A_lb_14.transpose()[i]:
                    pstring += f"{ column:5}"
                print(pstring)
            for i in range(self.Xtp_start,self.Xtp_end):
                pstring = str(self.X_indices[i])
                for column in A_lb_14.transpose()[i]:
                    pstring += f"{ column:5}"
                print(pstring)
            for i in range(self.Xts_start,self.Xts_end):
                pstring = str(self.X_indices[i])
                for column in A_lb_14.transpose()[i]:
                    pstring += f"{ column:5}"
                print(pstring)
            print("b_lb_14: ", b_ub_14)

        if wordy > 0:
            print("Constraint 14 built.")
        
        # 15) the number of task-starttime pairs must be equal to the number of tasks
        '''
        This ensures that each task is assigned a start time.

        (Xts_00 + ... + Xts_TS) = 1
        '''
        '''
        A_eq_15 = np.zeros((1, num_variables), dtype=int)
        b_eq_15 = np.array([self.T], dtype=int)

        A_eq_15[0,self.Xts_start:self.Xts_end] = 1
        '''
        A_eq_15 = np.zeros((self.T, num_variables), dtype=int)
        b_eq_15 = np.ones(self.T, dtype=int)

        for t in range(self.T):
            A_eq_15[t, (self.Xts_start + t * self.S):(self.Xts_start + t * self.S + self.S)] = 1
        
        if wordy > 1:
            print("A_eq_15^T:")
            for i in range(self.Xts_start,self.Xts_end):
                pstring = str(self.X_indices[i])
                for column in A_eq_15.transpose()[i]:
                    pstring += f"{ column:5}"
                print(pstring)
            print("b_eq_15: ", b_eq_15)
        
        A_eq_list.append(A_eq_15)
        b_eq_list.append(b_eq_15)

        if wordy > 0:
            print("Constraint 15 built.")

        # 16) Each task must be active for exactly the duration of its assigned asset
        '''
        This constraint works together with Constraint 14b to ensure proper duration handling:
        - Constraint 14b: Ensures tasks are active during their assigned duration periods (lower bound)
        - Constraint 16: Ensures tasks are active for exactly their total duration (upper bound)
        
        For each task t, the sum of Xtp periods must equal the duration of the assigned asset:
        sum(Xtp[t,p] for p in P) = sum(Xta[t,a] * duration[t,a] for a in A)
        '''
        rows_16 = []
        vec_16 = []

        for t in range(self.T):
            row = np.zeros(num_variables, dtype=int)
            # Left side: sum of all periods for task t
            for p in range(self.P):
                row[self.Xtp_start + t * self.P + p] = 1
            # Right side: subtract duration * assignment for each asset  
            for a in range(self.A):
                duration = self.task_asset_matrix[t, a, 1]
                if duration > 0:
                    row[self.Xta_start + t * self.A + a] = -duration
            
            rows_16.append(row)
            vec_16.append(0)  # sum(Xtp) - sum(duration * Xta) = 0
        
        if rows_16:
            A_eq_16 = np.vstack(rows_16)
            b_eq_16 = np.array(vec_16, dtype=int)
            A_eq_list.append(A_eq_16)
            b_eq_list.append(b_eq_16)
        
        if wordy > 0:
            print("Constraint 16 built.")

        # 17) Weather constraints: task-asset pairs cannot be assigned in periods with incompatible weather
        '''
        Assets have maximum weather conditions they can operate in (stored as 'max_weather' in asset dict).
        If the weather in period p exceeds an asset's max_weather capability, then no task can be 
        assigned to that asset in that period.
        
        For each asset a, period p, and task t:
        If weather[p] > asset[a]['max_weather'], then Xta[t,a] + Xtp[t,p] <= 1
        
        This prevents simultaneous assignment of both the task-asset pair AND the task-period pair
        when weather conditions are incompatible.
        '''
        rows_17 = []
        vec_17 = []
        
        for a in range(self.A):
            asset_max_weather = self.assets[a].get('max_weather', float('inf'))  # Default to no weather limit
            
            for p in range(self.P):
                period_weather = self.weather[p]
                
                if period_weather > asset_max_weather:
                    # Weather in period p is too severe for asset a
                    for t in range(self.T):
                        # Check if this task-asset pair is valid (positive duration and cost)
                        if (self.task_asset_matrix[t, a, 0] >= 0 and 
                            self.task_asset_matrix[t, a, 1] > 0):
                            
                            # Prevent task t from using asset a in period p due to weather
                            row = np.zeros(num_variables, dtype=int)
                            row[self.Xta_start + t * self.A + a] = 1      # Xta[t,a]
                            row[self.Xtp_start + t * self.P + p] = 1      # Xtp[t,p]
                            
                            rows_17.append(row)
                            vec_17.append(1)  # Xta[t,a] + Xtp[t,p] <= 1 (can't have both = 1)
        
        # Build constraint matrices if we have any weather constraints
        if rows_17:
            A_ub_17 = np.vstack(rows_17)
            b_ub_17 = np.array(vec_17, dtype=int)
            A_ub_list.append(A_ub_17)
            b_ub_list.append(b_ub_17)
            
            if wordy > 0:
                print(f"Constraint 17 built with {len(rows_17)} weather restrictions.")
        else:
            if wordy > 0:
                print("Constraint 17 built (no weather restrictions needed).")


        # --- End Constraints ---

        if wordy > 0:
            print("All constraints built. Stacking and checking constraints...")

        # --- Assemble the SciPy Constraints ---
        # A series of linear constraints required by the solver by stacking the constraint matrices and limits vectors
        # The number of rows in these matrices is equal to the number of constraints, so they can be vertically stacked

        # Check num columns of all constraint matrices matches number of decision variables before stacking
        for i, A in enumerate(A_ub_list):
            if A.size > 0 and A.shape[1] != num_variables:
                raise ValueError(f"Upper bound constraint matrix {i} has incorrect number of columns. Expected {num_variables}, got {A.shape[1]}.")
        for i, A in enumerate(A_eq_list):
            if A.size > 0 and A.shape[1] != num_variables:
                raise ValueError(f"Equality constraint matrix {i} has incorrect number of columns. Expected {num_variables}, got {A.shape[1]}.")
        for i, A in enumerate(A_lb_list):
            if A.size > 0 and A.shape[1] != num_variables:
                raise ValueError(f"Lower bound constraint matrix {i} has incorrect number of columns. Expected {num_variables}, got {A.shape[1]}.")

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
                for i in range(len(self.X_indices)):
                    print(f"  {self.X_indices[i]}: {int(res.x[i])}")

            if wordy > 0:
                print("Optimization successful. The following schedule was generated:")
                
            x_opt = res.x  # or whatever your result object is
            Xta = x_opt[self.Xta_start:self.Xta_end].reshape((self.T, self.A))
            Xtp = x_opt[self.Xtp_start:self.Xtp_end].reshape((self.T, self.P))
            Xap = x_opt[self.Xap_start:self.Xap_end].reshape((self.A, self.P))
            Xts = x_opt[self.Xts_start:self.Xts_end].reshape((self.T, self.S))

            for p in range(self.P):
                weather_condition = self.weather[p]
                pstring = f"Period {p} (weather {weather_condition}): "
                for t in range(self.T):
                    if Xtp[t, p] > 0: 
                        # Find assigned asset for this task
                        a_assigned = np.argmax(Xta[t, :])  # assumes only one asset per task
                        cost = self.task_asset_matrix[t, a_assigned, 0]
                        duration = self.task_asset_matrix[t, a_assigned, 1]
                        asset_name = self.assets[a_assigned].get('name', f'Asset {a_assigned}')
                        pstring += f"{asset_name} assigned to task {t} (cost: {cost}, duration: {duration}) | "
                    else:
                        pstring += " "*60 + "| "
                print(pstring)

        if wordy > 0:
            print("Optimization function complete.")
        

if __name__ == "__main__":

    os.system("clear") # for clearing terminal on Mac

    # A simple dummy system to test the scheduler with weather constraints

    # Weather periods with varying conditions (1=calm, 2=moderate, 3=severe)
    weather = [1, 1, 2, 3, 1]  # 6 time periods with different weather conditions
    
    # Example tasks, assets, dependencies, and task_asset_matrix
    tasks = [
        "task1",
        "task2"
    ]
    assets = [
        {"name": "heavy_asset", "max_weather": 3},   # Can work in all weather conditions
        {"name": "light_asset", "max_weather": 1}    # Can only work in calm weather (1)
    ]

    # task dependencies
    task_dependencies = {
        "task1": [],           # task1 has no dependencies  
        "task2": ["task1"]     # task2 depends on task1
    }
    
    # dependency types (optional - defaults to "finish_start" if not specified)
    dependency_types = {
        "task1->task2": "finish_start"  # task2 starts after task1 finishes
    }

    # cost and duration tuples for each task-asset pair. -1 indicates asset-task paring is invalid
    task_asset_matrix = np.array([
        [(2000, 2), (1000, 3)],    # task1: heavy_asset (expensive but fast), light_asset (cheap but slow)
        [(1500, 3), (-1, -1)]      # task2: both assets can do it, light_asset is cheaper
    ])

    # Expected behavior with weather constraints:
    # - light_asset can only work in periods 0,1,5 (weather=1)
    # - heavy_asset can work in any period (max_weather=3)
    # - task2 depends on task1 (finish_start dependency)

    # Find the minimum time period duration based on the task_asset_matrix
    min_duration = np.min(task_asset_matrix[:, :, 1][task_asset_matrix[:, :, 1] > 0])  # minimum non-zero duration

    # Sandbox for building out the scheduler
    scheduler = Scheduler(task_asset_matrix, tasks, assets, task_dependencies, dependency_types, weather, min_duration)
    scheduler.optimize()
    a = 2
    
    
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
    #     "task1": [],                    # task1 has no dependencies  
    #     "task2": ["task1"],            # task2 depends on task1
    #     "task3": ["task1", "task2"]    # task3 depends on both task1 and task2
    # }
    
    # # dependency types (optional - demonstrates different types)
    # dependency_types = {
    #     "task1->task2": "finish_start",  # task2 starts after task1 finishes (default)
    #     "task1->task3": "start_start",   # task3 starts when task1 starts  
    #     "task2->task3": "same_asset"     # task3 must use same asset as task2
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
    # scheduler = Scheduler(task_asset_matrix, tasks, assets, task_dependencies, dependency_types, weather, min_duration)
    # scheduler.optimize()