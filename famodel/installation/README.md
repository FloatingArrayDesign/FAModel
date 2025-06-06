This IO&M modeling framework was first set up by Rudy Alkarem. README by Ryan Davies

# Description

# TODO's

- [ ] Add dependencies of this code to famodel env
- [ ] Weather data module
- [ ] standardized and documented example cases
- [ ] remove fadesign dependency
- [ ] add set of default data that can be used. For example either break vesselDesc.yaml into sub vessel descriptions, or allow the code to read that if no custom vessel file is given
- [ ] deal with all the placeholder data added to make it run. Especially things like mobilization time (shouldn't this be calculated?)
- [ ] better error checking, with instructions on what users should do to resolve
- [ ] build a tree of what is called where and how pieces are connected. ex. where are actions looped through?
   - [ ]  Make a list of all possible actions that can be called

# Instructions for use

Run `python -m famodel.installation.01_unload4installation` or more generally `python -m famodel.installation.<filename>` from the FAModel directory `../../` from this directory.

NOTE: Need an environment that has a combo of FADesign and FAModel packages. (See the `famodel-env_mod4install.yaml` environment)

Created a `fadesign.patch` that contains the changes needed to FADesign to let pip pick up the conceptual module.

Need to run on the MoorPy dev branch (working as of commit 8b0bfe3f5aeb8457e499b3a53eeba79c4bb541b3).

Need numpy<2 (even though Floris gets mad about this)

# Misc. Notes

Units of time are hours unless stated otherwise. 

Naming convention per Leah: The naming convention used generally for the yamls describes the array so ILIA is individual line individual anchor (SLSA would be shared line shared anchor) and I think 1P2C represents the mooring line type (1 polyester section, 2 chain sections - so likely describing a chain-poly-chain setup) grd is I think for gridded shape (he has some options where there are clusters of turbines instead of grid).

The 01-05 steps were set up with a 3_ILIA_1P2C_grd_mini.yaml, which (I believe) is the same configuration as the Delmar_installation_staging_tool.xlsx, so it that can serve as a verification of this model. if we can't find that input configuration, then we can remake it with the instructions in the `Conceptual Desgin Tutorial.docx` file. 

# Copilot Explaination of the tool (outdated)

## Prompt: Explain the tool

The tool in the `FAModel/famodel/installation` directory is a framework for modeling the installation times and logistics of offshore wind mooring systems. Here's a detailed explanation based on the provided context:

### Purpose
The tool is designed to:
1. **Model Installation Times**: Simulate the time required for various stages of offshore wind mooring system installations, such as material transport, mobilization, and deployment.
2. **Logistics Planning**: Account for logistical constraints like vessel availability, port operations, and environmental conditions.
3. **Step-by-Step Guidance**: Provide numbered example scripts (`01_` to `05_`) that walk users through the process of using the tool.

### Key Components
1. **Example Files (`01_` to `05_`)**:
   - These files demonstrate the workflow of the tool:
     - 01_unload4installation.py: Unloads a project class from a `ConceptDesign` class.
     - 02_load4installation.py: Loads a project class and performs analysis and metrics.
     - 03_step1_materialItems.py: Creates material packages for installation.
     - 04_step2_actionItems.py: Generates action items for transport, mobilization, and installation.
     - 05_install_ex.py: Likely integrates all steps into a complete installation simulation.

2. **Core Classes**:
   - `InstallManager`: Manages the overall installation process, including scheduling events, registering vessels and ports, and running simulations.
   - `Vessel`: Represents vessels used in the installation process, with methods for mobilization, transit, and state logging.
   - `Port`: Represents ports where staging and logistics operations occur.
   - `ActionItem` and `Action`: Define tasks and dependencies for installation activities.

3. **Helper Functions**:
   - Found in install_helpers.py, these functions handle specific tasks like transporting materials, mobilizing vessels, and visualizing actions.

4. **Configuration Files**:
   - YAML files define input parameters for the simulations, such as mooring configurations and logistics setups.

5. **Spreadsheet Tool**:
   - `Delmar_installation_staging_tool.xlsm` appears to be a verification or supplementary tool for the model.

### Workflow
The tool's workflow involves:
1. **Loading Project Data**: Using the `ConceptDesign` or `Project` classes to load mooring system configurations.
2. **Creating Material Packages**: Grouping components into packages for transport and installation.
3. **Defining Actions**: Setting up tasks like mobilization, transport, and installation with dependencies.
4. **Simulating Installation**: Running the `InstallManager` to simulate the entire process, considering vessel and port operations.

### Notes from the README
- **Naming Conventions**: YAML file names describe the mooring array type (e.g., `ILIA` for individual line individual anchor, `1P2C` for chain-poly-chain setups).
- **Verification**: The `Delmar_installation_staging_tool.xlsm` can be used to verify the model's outputs.
- **Dependencies**: The tool relies on modules from `fadesign` and famodel.

### Next Steps
1. **Review Example Files**: Start with 01_unload4installation.py and follow the sequence to understand the workflow.
2. **Understand Core Classes**: Examine install_manager.py, vessel.py, and port.py to learn how the tool models logistics.
3. **Run Simulations**: Use the example YAML files and scripts to test the tool.
4. **Extend Functionality**: Based on your requirements, modify or add features to the tool.

Let me know if you need help with specific files or tasks!

## Prompt: Explain the hierarchy of function calls

The tool in the `FAModel/famodel/installation` directory follows a structured hierarchy of function calls and class interactions to model the installation process for offshore wind mooring systems. Here's the standard hierarchy and the roles of the key functions and classes:

---

### **Hierarchy of Function Calls**

1. **High-Level Workflow Management**
   - **`InstallManager`**:
     - Central class that orchestrates the entire installation process.
     - Manages vessels, ports, events, and the overall simulation timeline.
     - Key methods:
       - `registerVessel(file)`: Registers vessels for the installation process.
       - `registerPort(file)`: Registers ports for staging and logistics.
       - `scheduleEvent(time, agent, action, params)`: Schedules events in the simulation.
       - `run()`: Executes the simulation by processing scheduled events.
       - `createPkgs(project, stageMode)`: Clusters components into material packages and stages them at ports.

2. **Step-by-Step Examples**
   - **Example Scripts (`01_` to `05_`)**:
     - These scripts demonstrate the workflow:
       - 01_unload4installation.py: Unloads a project class from a `ConceptDesign` class.
       - 02_load4installation.py: Loads a project class and performs analysis.
       - 03_step1_materialItems.py: Groups components into material packages.
       - 04_step2_actionItems.py: Defines action items for transport, mobilization, and installation.
       - 05_install_ex.py: Integrates all steps into a complete simulation using `InstallManager`.

3. **Detailed Step Functions**
   - **`install_helpers.py`**:
     - Contains helper functions for specific tasks in the installation process.
     - Key functions:
       - `tranportTo_actionItem(vessel, distance2port)`: Creates transport action items.
       - `mobilizeM_actionItem(vessel, pkg)`: Creates mobilization action items.
       - `install_actionItem(vessel, pkg)`: Creates installation action items.
       - `visualizeAction(action)`: Visualizes action dependencies as a graph.

4. **Core Classes for Specific Details**
   - **`Vessel`**:
     - Represents vessels used in the installation process.
     - Handles mobilization, transit, and state logging.
   - **`Port`**:
     - Represents ports for staging and logistics operations.
     - Manages storage and staging of material packages.
   - **`Action`**:
     - Represents actions in the installation process, such as transport or deployment.
     - Tracks dependencies and durations of individual tasks.

---

### **Roles of Key Functions and Classes**

1. **High-Level Assembly of Information**
   - **`InstallManager`**:
     - Assembles all information about vessels, ports, and material packages.
     - Coordinates the simulation by scheduling and executing events.
     - Handles staging of material packages at ports.

2. **Specific Details of Installation Steps**
   - **install_helpers.py**:
     - Provides detailed logic for creating action items for transport, mobilization, and installation.
   - **`Vessel` and `Port`**:
     - Manage the state and operations of vessels and ports, respectively.
   - **`Action` and `ActionItem`**:
     - Define and manage individual tasks and their dependencies.

3. **Example Scripts**
   - Serve as a guide for using the tool step-by-step, from unloading a project to running a complete simulation.

---

### **Summary**
- **High-Level Coordination**: `InstallManager` is the central class that assembles all information and manages the workflow.
- **Detailed Steps**: Functions in install_helpers.py and classes like `Vessel`, `Port`, and `Action` handle specific details of the installation process.
- **Examples**: Scripts (`01_` to `05_`) demonstrate the workflow and how to use the tool effectively.


## Prompt: Explain the variables in the InstallManager.run function
Here is an explanation of all the variables in the `InstallManager.run` function:

```python
def run(self):
    """
    Run the installation phase simulation.
    This method processes scheduled events in the priority queue, executing actions
    and updating the current time.
    It continues until all events have been processed.
    The method also logs each event as it is processed.
    
    Parameters
    ----------
    None
        
    Returns
    -------
    None
    """
    while self.events:
        t, agent, action, params = heapq.heappop(self.events)
        self.now = t
        func = getattr(agent, action)
        completed, new_events = func(t, **params)
        self.logs.append({"time": t, "agent": agent.name, "action": action})
        for evt in new_events:
            self.scheduleEvent(*evt)
```

### Variables:

1. **`self.events`**:
   - A priority queue (list) that stores scheduled events.
   - Each event is a tuple of the form `(time, agent, action, params)`, where:
     - `time`: The time at which the event is scheduled to occur.
     - `agent`: The object (e.g., a `Vessel` or `Port`) responsible for performing the action.
     - `action`: The name of the method (as a string) to be executed by the agent.
     - `params`: A dictionary of parameters required for the action.

2. **`t`**:
   - The time at which the current event is scheduled to occur.
   - Extracted from the event tuple using `heapq.heappop(self.events)`.

3. **`agent`**:
   - The object responsible for performing the action in the current event.
   - This could be a `Vessel`, `Port`, or any other class instance registered with the `InstallManager`.

4. **`action`**:
   - The name of the method (as a string) to be executed by the `agent`.
   - For example, this could be `"mobilize"`, `"transport"`, or `"install"`.

5. **`params`**:
   - A dictionary of parameters required for the `action` method.
   - These parameters are passed to the method when it is called.

6. **`self.now`**:
   - The current simulation time.
   - Updated to the time of the current event (`t`) as the simulation progresses.

7. **`func`**:
   - A reference to the method of the `agent` that corresponds to the `action` string.
   - Obtained using `getattr(agent, action)`.

8. **`completed`**:
   - A boolean value returned by the `func` method, indicating whether the action was successfully completed.

9. **`new_events`**:
   - A list of new events generated by the `func` method.
   - Each event in this list is a tuple of the form `(time, agent, action, params)`, similar to the events in `self.events`.

10. **`self.logs`**:
    - A list that stores logs of all processed events.
    - Each log entry is a dictionary with the following keys:
      - `"time"`: The time at which the event occurred.
      - `"agent"`: The name of the agent that performed the action.
      - `"action"`: The name of the action that was performed.

11. **`evt`**:
    - A single event from the `new_events` list.
    - Each `evt` is unpacked and passed to the `self.scheduleEvent` method to schedule it for future execution.

---

### Summary of Workflow:
1. The function processes events from the `self.events` priority queue in chronological order.
2. For each event:
   - It updates the current simulation time (`self.now`).
   - It executes the specified action (`func`) on the agent with the given parameters.
   - It logs the event in `self.logs`.
   - It schedules any new events generated by the action for future execution.
3. The loop continues until all events in the queue have been processed.

# Dependencies to FAD tools
File: `01_unload4installation.py` Dependencies:
```python
from fadesign.conceptual.conceptDesign import ConceptDesign
```
File: `02_load4installation.py` Dependencies:
```python
from famodel.project import Project
from fadesign.conceptual import metrics as mtr
```
File: `03_step1materialItems` Dependencies:
```python
from famodel.project import Project
```
File: `04_step2_actionItems.py` Dependencies:
```
none
```
File: `05_install_ex.py` Dependencies:
```
none
```
File: `action.py` Dependencies:
```
none
```
File: `install_helpers.py` Dependencies:
```
none
```
File: `install_manager.py` Dependencies:
```
none
```
File: `port.py` Dependencies:
```
none
```
File: `vessel.py` Dependencies:
```
none
```
