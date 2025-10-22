## Key Terminologies

### Actions
- **Definition**: The smallest unit of work that the system simulates.
- **Purpose**: Represents a specific action to be performed, such as transporting an anchor, installing a mooring, or deploying a WEC.
- **Examples**:
  - "Anchor installation"
  - "Mooring deployment"

### Tasks
- **Definition**: A logical group of one or more actions that are bounded by a "From-To Port" constraint.
- **Purpose**: Represents a higher-level grouping of actions that are executed together as part of a specific phase of the installation process.
- **Examples**:
  - "Anchor Installation Task"
  - "Mooring Deployment Task"

### Dependencies
- **Definition**: Logical constraints that determine the dependencies between actions.
- **Purpose**: Ensures that actions are executed in a given order based on logical and physical requirements.
- **Examples**:
  - "Anchor installation depends on anchor transport to the site."
  - "Mooring deployment depends on anchor installation."

### Action Sequencing
- **Definition**: The process of determining the sequence in which actions take place within a task.
- **Purpose**: Ensures that actions are executed in a logical and efficient order, respecting dependencies and resource constraints.

### Capabilities
- **Definition**: The specific functionality that an asset (e.g., vessel, port) can perform.
- **Purpose**: Determines which assets are suitable for specific actions.
- **Examples**:
  - A vessel with "crane" capability can install anchors.

### Metrics
- **Definition**: Quantifiable measurements of assets based on their capabilities.
- **Purpose**: Used to evaluate and compare assets for suitability and efficiency in performing actions.
- **Examples**:
  - **Speed**: The transit speed of a vessel.
  - **Capacity**: The cargo capacity of a vessel.

### Roles
- **Definition**: Functional assignments of assets in the context of actions.
- **Purpose**: Specifies how each asset contributes to the completion of an action.
- **Examples**:
  - **Carrier**: Assigned to carry specific equipment or materials.
  - **Operator**: Assigned to operate machinery or perform specialized tasks.

---