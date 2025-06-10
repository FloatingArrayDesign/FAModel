This IO&M modeling framework was first set up by Rudy Alkarem. README by Ryan Davies

# Description

# TODO's

- [ ] Add dependencies of this code to famodel env
- [ ] Weather data module
- [ ] standardized and documented example cases
- [ ] remove fadesign dependency
- [X] add set of default data that can be used. For example either break vesselDesc.yaml into sub vessel descriptions, or allow the code to read that if no custom vessel file is given
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
