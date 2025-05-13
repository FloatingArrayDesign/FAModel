This IO&M modeling framework was first set up by Rudy Alkarem. README by Ryan Davies

# Description


# Instructions for use


# Misc. Notes

Units of time are hours unless stated otherwise. 

Naming convention per Leah: The naming convention used generally for the yamls describes the array so ILIA is individual line individual anchor (SLSA would be shared line shared anchor) and I think 1P2C represents the mooring line type (1 polyester section, 2 chain sections - so likely describing a chain-poly-chain setup) grd is I think for gridded shape (he has some options where there are clusters of turbines instead of grid).

The 01-05 steps were set up with a 3_ILIA_1P2C_grd_mini.yaml, which (I believe) is the same configuration as the Delmar_installation_staging_tool.xlsx, so it that can serve as a verification of this model. if we can't find that input configuration, then we can remake it with the instructions in the `Conceptual Desgin Tutorial.docx` file. 

For getting things running, I changed to ______, from `FADesign/scripts/concepts/ILIA/______`. 

# Copilot Explaination of the tool
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
   - YAML files (e.g., VolturnUS-S_example.yaml, `presetLogistics.yaml`) define input parameters for the simulations, such as mooring configurations and logistics setups.

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
