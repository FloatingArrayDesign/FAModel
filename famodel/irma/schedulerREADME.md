# Scheduler Mathematical Formulation

This document describes the mathematical formulation of the scheduling problem solved by the `Scheduler` class.

## Sets and Indices
- $P$: Set of periods, $p = 1, \ldots, P$
- $T$: Set of tasks, $t = 1, \ldots, T$
- $A$: Set of assets, $a = 1, \ldots, A$
- $R$: set of task requirements/dependencies, $r =1, \ldots, R \text{ where } R < T$

## Parameters
- $v_{t,a}$: Value of assigning asset $a$ to task $t$. Can be either cost or duration depending on user input.
- $c_{t,a}$: Cost of assigning asset $a$ to task $t$
- $d_{t,a}$: Duration (in periods) required for asset $a$ to complete task $t$

## Decision Variables
- $x_{p,t,a} \in \{0,1\}$: 1 if asset $a$ is assigned to task $t$ in period $p$, 0 otherwise

## Objective Function
Minimize total cost:

$$
\min \sum_{p \in P} \sum_{t \in T} \sum_{a \in A} \left( c_{t,a} \right) x_{p,t,a}
$$

## Constraints

### 0. Total Assignment Limit
The sum of all assignments cannot exceed the number of periods:
$$
\sum_{p \in P} \sum_{t \in T} \sum_{a \in A} x_{p,t,a} \leq P
$$

### 1. Task-Asset Validity
Only valid task-asset pairs can be assigned:
$$
x_{p,t,a} = 0 \quad \forall p, t, a \text{ where } c_{t,a} < 0 \text{ or } d_{t,a} < 0
$$

### 2. Task Dependencies (**Disabled/In-progress**)
Tasks with dependencies must be scheduled after their dependencies are completed. This might need to be reworked, still figuring out the best way to enforce temporal constraints. 
$$
p * x_{p,r,a} < p * x_{p,t,a} \quad t \in T, \forall r \in R_t, p \in P, a \in A  
$$

### 3. Weather Constraints (**TODO**)
Assets cannot be assigned in periods where weather exceeds their capability.

### 4. Asset Cannot Be Assigned to Multiple Tasks in Same Period
Each asset can be assigned to at most one task in each period:
$$
\sum_{t \in T} x_{p,t,a} \leq 1 \quad \forall p \in P, a \in A
$$

### 5. Task Assignment Limit (**TODO**)
The total number of tasks assigned cannot exceed the number of tasks available.

### 6. Asset Assignment Limit (**TODO**)
The total number of assets assigned cannot exceed the number of assets available.

### 7. Early Assignment Constraint (**TODO**) 
A task cannot be assigned if it could have been assigned in an earlier period. This encourages the solver to assign tasks to the earliest possible periods.

### 8. All Tasks Must Be Assigned
Each task must be assigned at least once:
$$
\sum_{p \in P} \sum_{a \in A} x_{p,t,a} \geq 1 \quad \forall t \in T
$$

### 9. Assignment Duration (**Disabled/In-progress**)
Each task-asset pair must be assigned for exactly its required duration:
$$
\sum_{p \in P} x_{p,t,a} = d_{t,a} \quad \forall t \in T, a \in A \text{ with } d_{t,a} > 0
$$

### 10. Assignment Window
A task cannot be assigned to periods that would exceed the available time window:
$$
x_{p,t,a} = 0 \quad \forall p, t, a \text{ where } p + d_{t,a} > P
$$

---

**Notes:**
- Constraints marked **TODO** are not yet implemented in the code but are (likely) necessary for an optimal solution
- Constraints marked **Disabled/In-progress** are works in progress that if enabled cause an infeasible solution to be generated.
- Additional constraints (e.g., weather, dependencies) should be added.
- This approach isn't finalized. We may need additional decision variables if we want to have multiple objectives 
    - For example, one way to force earlier scheduling is to add a start-time decision variable that gives a penalty for later start-times
    - Did not implement this because we may not want to force earlier start-times