# Scheduler Mathematical Formulation (as implemented in scheduler.py)

This document describes the mathematical formulation of the scheduling problem solved by the `Scheduler` class, using multiple decision variables and following the numbering and naming conventions in `scheduler.py`.

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
- $X_{t,s} \in \{0,1\}$: 1 if task $t$ starts at period $s$, 0 otherwise

$x = [X_{t,a}  X_{t,p}  X_{t,s}] $

## Objective Function
Minimize total cost (cost is only determined by task-asset assignment):

$$
\min \sum c_{t,a} x
$$

The $c$ vector also contains 'cost' penalties for later start times $(X_{t,s})$ to prioritize tasks starting as early as they can (used to be Constraint 7)

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

### 1. Task-Asset Validity
Only valid task-asset pairs can be assigned:

$$
X_{t,a} = 0 \quad \forall t, a \text{ where } c_{t,a} < 0 \text{ or } d_{t,a} < 0
$$

### 3. At Least One Asset Per Task
Sum of all task-asset pairs must be >= 1 for each task:

$$
\sum_{a=0}^{A-1} X_{t,a} \geq 1 \quad \forall t
$$

### 15. The number of task-starttime pairs must be equal to the number of tasks
This ensures that each task is assigned exactly 1 start time.

$$
\sum_{s=0}^{S-1} X_{t,s} = 1 \quad \forall t
$$

### 10. A task cannot start in a period where its duration would exceed the maximum number of time periods
This ensures that a task is not assigned to a period that would cause it to exceed the total number of periods available.

$$
X_{t,a}[t,a] + X_{t,s} <= 1 \quad \forall t, a, s \text{ where } d_{t,a} > 0,\ s + d_{t,a} > P
$$

When a task-asset pair is assigned, then for each start time of that task, it $(X_{t,s})$ has to be zero under these conditions, where s+d>P

### 14a. A task must occupy the same period that it starts in
This ensures that the task start-time decision variable is non-zero if a task is assigned to any period.

$$
X_{t,p}[t,s] \geq X_{t,s}[t,s] \quad \forall t
$$

In every start time for each task, the corresponding period must be equal to that start time decision variable

### 14b. A task-asset assignment must be active for the duration required by that assignment

14a ensured that if a task is assigned a start time, the corresponding task-period pair for the period equal to the start time is selected.

14b ensures that if a task is assigned a start time, the number of periods equal to the duration of the task are also turned on.

$$
X_{t,a}[t,a] + X_{t,s}[t,s] - X_{t,p}[t,p] <= 1 \quad \forall t,a,s,p(s<s+d_{t,a}) 
$$

Loops through all tasks, assets, periods, and start times and assigns constraint coefficients $(A)$ to ensure the right number of periods are selected.

For example, $X_{t,a}[0,0] + X_{t,s}[0,0] - X_{t,p}[0,0] <= 1$
ensures that $X_{t,p}[0,0]$ is equal to 1 when the first task-asset pair is set to start in period 1, as well as the second period, given a task-asset duration of 2:
$X_{t,a}[0,0] + X_{t,s}[0,0] - X_{t,p}[0,1] <= 1$

This ensures that based on the task-asset pair selected and a given start time, the corresponding Xtp values are selected. However, this constraint represents a 'lower_bound', as these constraints don't have any control over the other Xtp values that it doesn't specify. Constraint 16 sets the 'upper bound'.

### 16. Each task is active for exactly its duration

This constraint ensures that the number of Xtp variables is EXACTLY equal to the task-asset duration, as 14b has some edge cases where some other Xtp values have the option to be turned on.

$$
\sum_{p=0}^{P-1} X_{t,p}[t,p] = \sum_{a=0}^{A-1} X_{t,a}[t,a]d_{t,a} 
$$

This prevents extra periods beyond the requirement from being set.

Overall, the combination of Constraints 15-14b-16 ensure tasks are active when they should be

(**TODO**) Handle cases with multiple assets per task

The problem with this constraint is that it assumes only 1 asset is assigned per task, and it uses that task-asset duration to set the constraint. But, what if asset 0 and asset 1 are both used for task 0? How is the duration set?


1. Ensure that each task only involves 1 asset in its definition. Would just require an update to Constraint 3 to an equals not a inequality. Would likely involve more dependencies to define, but it seems to be the cleanest.
   - Could set 'cases' instead of 'assets' and update Constraint 3 to equality
2. Updating the duration of the task based on the combination of assets. Would require some changes to 14b (not bad) and just an updated duration in 16.
   - Take the maximum duration (3 over 2)
   - Average them in a way (2 + 3 = 1.25)
   - Have additional duration information stored somewhere else
   - Sum them (make them sequential) like if $d_{t,a}[0,0]=2$ and $d_{t,a}[0,1]=3$, then this task duration would require 5 periods.




### 12. (NOT INCLUDED ANYMORE) An asset that a task is assigned to must also occur in the same period(s) that the task occurs in
This constraint is commented out and not used anymore because it created an adverse coupling of variables.

It was intended to ensure that a task and asset in a task-asset pair were both assigned to the same period (and updated the Xap variables accordingly)

$$
0 <= X_{t,p}[t,p] - X_{a,p}[a,p] + X_{t,a}[t,a] <= 1 
$$

However, this created a problem that when $X_{t,a}[0,0] = 1$, this constraint could set $X_{t,p}[0,p] = X_{a,p}[0,p]$, but if Asset 0 was used in another task at a different time, it may not occur in that same period.

Therefore, we have decided to eliminate the need for the Xap variable since they seem redundant after this investigation.


### 4. Asset Cannot Be Assigned to Multiple Tasks in Same Period
Ensures that if multiple tasks are assigned to the same asset, then only one can be active in any period p (not using any Xap variables)

$$
\sum_{t=0}^{T-1} X_{t,p}[t,p] + \sum_{t=0}^{T-1} X_{t,a}[t,a] <= 1 + T \quad \forall a,p 
$$

For two tasks, the right hand side will always be 3, which ensures that if the same asset is used for multiple tasks $(X_{t,a}[0,0]=1, X_{t,a}[1,0]=1)$, then only one $X_{t,p}[t,p]$ variable is turned on.



### 17. Weather Constraints
Assets cannot be assigned in periods where weather exceeds their capability.

Assets are given a 'max_weather' parameter (1, 2, or 3) and each period is initialized as a 1, 2, or 3 describing its state.

If weather[p] > asset[a]['max_weather'], then $X_{t,a}[t,a] + X_{t,p}[t,p] <= 1$

Meaning, do not turn on the task in a period whether the period's weather is greater than the maximum allowable weather capability of the asset

### 2. Task Dependencies
Tasks with dependencies must be scheduled according to their dependency rules.

We have a set of different dependency type options:
- Finish-Start: the dependent task starts after the prerequisite task finishes
- Start-Start: the dependent task starts when the prerequisite task starts
- Finish-Finish: the dependent task finishes when the prerequisite task finishes
- Same-Asset: the dependent task must use the same asset as the prerequisite task

For all valid start times s for task t, if $X_{t,s}[t,s]=1$, then there is some other start time $s_d$ for task d so that $X_{t,s}[d,s_d]=1$ and $s_d + duration <= s$

$X_{t,s}[t,s] <= \sum X_{t,s}[d,s_d]$ from $s$ to $sd+duration$

---

**Notes:**
- $d_{t,a}$ is the duration for asset $a$ assigned to task $t$. If multiple assets are possible, $X_{t,a}$ determines which duration applies.
- This approach separates assignment, activity, and start variables for clarity and easier constraint management.

- Constraints can be extended for parallel tasks, multiple assets per task, or other requirements as needed.
- One of the better references to understand this approach is `Irwan et al. 2017 <http://dx.doi.org/10.1016/j.cor.2015.09.010>`_
- The `scheduler.py` file also has some TODO's, which are focused on software development.
