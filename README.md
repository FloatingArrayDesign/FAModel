# Floating Array Design Toolset

The Floating Array Design (FAD) Toolset is a collection of tools for
modeling and designing arrays of floating offshore structures. It was
originally designed for floating wind systems but has applicability
for many offshore applications.

A core part of the FAD Toolset is the Floating Array Model (FAModel),
which serves as a high-level library for efficiently
modeling a floating wind array. It combines site condition information and a 
description of the floating array design, and contains functions for evaluating
the array's behavior considering the site conditions. For example, it combines
information about site soil conditions and an array's anchor characteristics to
estimate the holding capacity of each anchor.

The library works in conjunction with the tools RAFT, MoorPy, and FLORIS to model floating
wind turbines, mooring systems, and array wakes respectively.

Layered on top of the floating array model is a set of design tools that can
be used for algorithmically adjusting or optimizing parts of the a floating
array. Specific tools existing for mooring lines, shared mooring systems, 
dynamic power cables, static power cable routing, and overall array layout.
These capabilities work with the design representation and evaluation functions
in FAModel, and they can be applied by users in various combinations to suit
different purposes. 

In addition to standalone uses of the FAD Toolset, a coupling has been made with
[Ard](https://github.com/WISDEM/Ard), a sophisticated and flexible wind farm
optimization tool. This coupling allows Ard to use certain mooring system
capabilities from FAD to perform layout optimization of floating wind farms
with Ard's more advanced layout optimization capabilities.

The FAD Toolset works with the [IEA Wind Task 49 Ontology](https://github.com/IEAWindTask49/Ontology),
which provides a standardized format for describing floating wind farm sites
and designs. 

See example use cases in our [examples](./examples/README.md) folder.

## Pre-installation Requirements
The FAD Toolset is built entirely in Python. It is recommended that users 
familiarize themselves with basic Python commands before use. 
For working with the library, it is important to understand the floating array 
model structure, which is described more [here](./famodel/README.md).


## Installation
To install the FAD Toolset itself, first clone this FAD-Toolset repository.

The dependencies required by FAD depend on how it is used. To install all
possible required dependencies, you can create a 
new python virtual environment based on the included yaml listing the required 
dependencies.

In the terminal (Anaconda Powershell Prompt), clone this repository to a 
directory of your choice, navigate into the main folder of the repository, and 
run the following command:

    conda env create -f famodel-env.yaml

This command will install all the dependencies required to run FAD.
Activate your virtual environment before using FAD with ```conda activate famodel-env```

To install the FAD Toolset package in your environment, enter the 
following in the command line from the FAD-Toolset directory.

For development use:

run ```python setup.py develop``` or ```pip install -e .``` from the command 
line in the main FAD-Toolset directory.

For non-development use:

run ```python setup.py``` or ```pip install .``` from the command line in 
the main FAD-Toolset directory.

You can test the installation by running ```pytest``` from the main FAD-Toolset directory.

<!-- FAD requires MoorPy and we currently install it separately. If you don't already have it,
you can install MoorPy with ```git clone https://github.com/NREL/MoorPy.git```
then navigate to the MoorPy folder and install with ```pip install .```.
Make sure your virtual enviroment is activated before installing MoorPy. -->


## Subpackages

The library has a core Project class for organizing information, classes for each component of an array and an evolving
collection of subpackages for specific functions. The current subpackages are:

- anchors: contains modules for anchor capacity calculations, in addition to the anchor class
- failures: contains modules for failure modeling with graph theory, and allows for enactment of a failure mode.
- seabed: contains modules for seabed bathymetry and boundary information
- design: contains various tools for performing design steps.

Please navigate into the subfolders above for additional information.

## Getting Started
The easiest way to create a FAD project is to provide the array 
information in an ontology yaml file. FAD has been designed
to work with a specific ontology yaml setup, which is described 
in detail in the [Ontology ReadMe](./famodel/ontology/README.md).

The [example driver file](./famodel/example_driver.py) creates a FAD Project 
object from a pre-set ontology file and shows the syntax and outputs of 
various capabilities. For guidance on creating your own ontology yaml file, 
it is recommended to read through the [Ontology ReadMe](./famodel/ontology/README.md), 
then either adapt one of the ontology samples or fill in the ontology template. 

The [core model readme](./famodel/README.md) describes the Project class structure, 
as well as the properties and methods of each component class. 

There are some limited helper functions to automatically fill in sections 
of a yaml from a MoorPy system or a list of platform locations. 
See [helpers](./famodel/helpers.py) for the full list of yaml writing capabilities. 


## Authors

The NREL Floating Wind Array Design team.