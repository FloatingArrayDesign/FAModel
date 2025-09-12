# Scheduler Mathematical Formulation (as implemented in scheduler.py)

This document describes the mathematical formulation of the scheduling problem solved by the `Scheduler` class, using multiple decision variables and following the numbering and naming conventions in `scheduler.py`. Incomplete constraints are marked as **TODO**.

## Sets and Indices
- $T$: Set of tasks, $t = 0, \ldots, T-1$
- $A$: Set of assets, $a = 0, \ldots, A-1$
- $P$: Set of periods, $p = 0, \ldots, P-1$
- $S$: Set of possible start periods, $s = 0, \ldots, S-1$ ($S = P$)

## Parameters
- $c_{t,a}$: Cost of assigning asset $a$ to task $t$
- $d_{t,a}$: Duration (in periods) required for asset $a$ to complete task $t$

## Decision Variables
- $X_{t,a} \in \{0,1\}$: 1 if task $t$ is assigned to asset $a$, 0 otherwise
- $X_{t,p} \in \{0,1\}$: 1 if task $t$ is active in period $p$, 0 otherwise
- $X_{a,p} \in \{0,1\}$: 1 if asset $a$ is used in period $p$, 0 otherwise 
- $X_{t,s} \in \{0,1\}$: 1 if task $t$ starts at period $s$, 0 otherwise

## Objective Function
Minimize total cost (cost is only determined by task-asset assignment):
$$
\min \sum_{t=0}^{T-1} \sum_{a=0}^{A-1} c_{t,a} X_{t,a}
$$

## Constraints

The below constraints are formulated such that they can be made into three giant matricies of upper and lower bounds and equalities. 
When added together, these are the upperbound constraint, the lower bound constraint, and the equality constraint. The solver
attempts to solve the object objective function subject to: 

subject to:   
$$
\text{1) } A_{ub} \text{ } x \text{ } \leq b_{ub} \\
\text{2) } A_{eq} \text{ } x \text{ } = b_{eq}    \\
\text{3) } A_{lb} \text{ } x \text{ } \geq b_{lb} \\
\text{4) } 0 \leq \text{ } x \text{ } \leq 1      \\
$$

### 0. Total Assignment Limit
The sum of all task-period assignments cannot exceed the number of periods:
$$
\sum_{t=0}^{T-1} \sum_{p=0}^{P-1} X_{t,p} \leq P
$$

### 1. Task-Asset Validity
Only valid task-asset pairs can be assigned:
$$
X_{t,a} = 0 \quad \forall t, a \text{ where } c_{t,a} < 0 \text{ or } d_{t,a} < 0
$$

### 2. Task Dependencies (**TODO**)
Tasks with dependencies must be scheduled after their dependencies are completed. Thus the starttime of a task must be greater than the end time of all the dependent tasks.

The general idea is that $X_{tp}(d,p) < X_{t,s}(t,s)$ where $d$ is the task that task $t$ is dependent on and where $p = s-1$

### 3. At Least One Asset Per Task
Sum of all task-asset pairs must be >= 1 for each task:
$$
\sum_{a=0}^{A-1} X_{t,a} \geq 1 \quad \forall t
$$

### 4. Asset Cannot Be Assigned to Multiple Tasks in Same Period
This means the sum of each asset-period pair in a given period must be less or equal to than 1. This prohibits multiple assets in a period. 
$$
\sum_{t=0}^{A-1} X_{a,p} \leq 1 \quad \forall p
$$


### 5. Every task must be assigned to at least one time period
Sum of all task-period pairs for each task must be >= 1:
$$
\sum_{p=0}^{P-1} X_{t,p} \geq 1 \quad \forall t
$$

### 6. The total number of assets assigned cannot be greater than the number of assets available but must be greater than the number of tasks. 
Sum of all asset-period pairs must be >= T:
$$
T \leq \sum_{a=0}^{A-1} \sum_{p=0}^{P-1} X_{a,p} \leq A
$$

### 7. Early Assignment Constraint (**TODO**)
A task cannot be assigned if it could have been assigned in an earlier period. This encourages the solver to assign tasks to the earliest possible periods.
This could be enforced with a penality multiplied by the Xts decision variable in the objective function.

### 8. All Tasks Must Be Assigned to At Least One Period
The sum of all task-period decision variables for each task must be greater than 1, indicating all tasks were assigned at least once:
$$
\sum_{p=0}^{P-1} X_{t,p} \geq 1 \quad \forall t
$$

### 9. Empty

### 10. A task duration plus the start-time it is assigned to must be less than the total number of time periods available
This ensures that a task is not assigned to a period that would cause it to exceed the total number of periods available.
$$
X_{t,s} = 0 \quad \forall t, a, s \text{ where } d_{t,a} > 0,\ s + d_{t,a} > P
$$

Note: this constraint is working, but $X_{t,p}$ is not currently being forced to match $X_{t,s}$ so it is not reflected in the final results (which check for $X_{tp} \neq 0$). If you look at the results generated you will see the $X_{t,s}$ decision variable respects this constraint, but the $X_{t,p}$ does not. Constraint 14 aims to force $X_{t,p}$ to start blocks of time assignments at $X_{t,s}$. 

### 11. The total number of task period pairs must be greater than or equal to the number of task-start time pairs
This ensures that the task start-time decision variable is non-zero if a task is assigned to any period.
$$
\sum_{p=0}^{P-1} X_{t,p} \geq \sum_{s=0}^{S-1} X_{t,s} \quad \forall t
$$

### 12. The period an asset is assigned to must match the period the task in the task-asset pair is assigned to
This ensures the chosen task and asset in a task asset pair are assigned to the same period. This means that if an asset 
is assigned to a task, then the corresponding task-period and asset-period pairs must be equal. 

if $X_{t,a} = 1$, then $X_{t,p} = X_{a,p}$, else if $X_{t,a} = 0$, then $X_{t,p}$ and $X_{a,p}$ can be anything. This requires two constriants:
$$
X_{t,p} - X_{a,p} + X_{t,a} \leq 1 \\
X_{t,p} - X_{a,p} + X_{t,a} \geq 0 \quad \forall t, a, p
$$

### 13. The total number of asset period pairs must be greater than or equal to the number of task-period pairs
This ensures that the 0 asset-period pairs solution is not selected
$$
\sum_{a=0}^{A-1} X_{a,p} \geq \sum_{t=0}^{T-1} X_{t,p} \quad \forall p
$$

### 14. If a task-starttime pair is selected, the corresponding task-period pair must be selected for the period equal to the start time plus the duration of the task (**In progress/Close**)
This ensures that if a task is assigned a start time, the corresponding task-period pair for the period equal to the start time plus the duration of the task is also selected.
`Xts[t, s] <= Xtp[t, s : s + d]`  for all tasks t in `range(0:T)` and start times s in `range(0:S)` where d is the duration of task t with the asset assigned to it

This constraint is very close. The matrix is being generated correctly, but the values in it (1, -1, 2, etc.) need to be decided so that the blocks of time are enforced only when a task asset pair is selected. Otherwise it should allow the optimization to continue. 

### 15. The number of task-starttime pairs must be equal to the number of tasks
This ensures that each task is assigned a start time.
$$
\sum_{s=0}^{S-1} X_{t,s} = 1 \quad \forall t
$$

### 16. Weather Constraints (**TODO**)
Assets cannot be assigned in periods where weather exceeds their capability.

---

**Notes:**
- $d_{t,a}$ is the duration for asset $a$ assigned to task $t$. If multiple assets are possible, $X_{t,a}$ determines which duration applies.
- This approach separates assignment, activity, and start variables for clarity and easier constraint management.
- Constraints marked **TODO** are not yet implemented in the code but are probably necessary for a truely opptimal solution.
- Constraints marked **In-progress** have code written but it is not yet working
- Constraints not marked **TODO** or **In-progress** are complete
- Constraints can be extended for parallel tasks, multiple assets per task, or other requirements as needed.
- One of the better references to understand this approach is `Irwan et al. 2017 <http://dx.doi.org/10.1016/j.cor.2015.09.010>`_
- The `scheduler.py` file also has some TODO's, which are focused on software development.