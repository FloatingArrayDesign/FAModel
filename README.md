# FAModel

The FAModel (or Floating Array Model) package serves as a high-level library for
modeling a floating wind array. It combines site condition information and a 
description of the floating array design, and contains functions for evaluating
the array's behavior considering the site conditions. For example, it combines
information about site soil conditions and an array's anchor characteristics to
estimate the holding capacity of each anchor.

The library works in conjunction with the tools RAFT and MoorPy to model floating
wind turbines and mooring systems, respectively.

In addition to the code, this repository defines a 
[Floating Array Ontology](https://github.com/FloatingArrayDesign/FAModel/tree/main/famodel/ontology), 
which provides a standardized description format for floating wind farms. 

An example of use of these tools to model three mooring lines over the bathymetry 
of the Humboldt lease area is shown below.

![Humboldt](famodel/seabed/images/slopeview4.PNG)


## Installation

To install FAModel itself, clone the FAModel repository and then enter the 
following in the command line from its directory.

For development use:

run ```python setup.py develop``` or ```pip install -e .``` from the command line in the main MoorPy directory.

For non-development use:

run ```python setup.py``` or ```pip install .``` from the command line in the main MoorPy directory.

The dependencies required by FAModel depend on how it is used. To install all
possible required dependencies, you can create a 
new python virtual environment based on the included yaml listing the required 
dependencies.

In the terminal (Anaconda Powershell Prompt), clone this repository to a 
directory of your choice, navigate into the main folder of the repository, and 
run the following command:

    conda env create -f famodel-env.yaml

This command will install all the dependencies required to run FAModel.


## Subpackages

The library has a core Project class for organizing information, classes for each component of an array and an evolving
collection of subpackages for specific functions. The two current subpackages are:

- anchors: contains modules for anchor capacity calculations, in addition to the anchor class
- failures: contains modules for failure modeling with graph theory, and allows for enactment of a failure mode in integrated FAModel tools such as MoorPy and RAFT.
- seabed: contains modules for seabed bathymetry and boundary information

Please navigate into the subfolders above for additional information.

## Getting Started
The easiest way to create an FAModel project is to provide the array information in a yaml file. FAModel has been designed to work with a specific ontology yaml setup, which is described in detail in the [Ontology ReadMe](./famodel/ontology/README.md).

The [example driver file](./famodel/example_driver.py) creates an FAModel project from a pre-set ontology file and shows the syntax and outputs of various capabilities. For guidance on creating your own ontology yaml file, it is recommended to read through the [Ontology ReadMe](./famodel/ontology/README.md), then either adapt one of the ontology samples or fill in the ontology template. 

There are some limited helper functions to auntomatically fill in sections of a yaml from a moorpy system or a list of platform locations. See [FAModel helpers](./famodel/helpers.py) for the full list of yaml writing capabilities. Many of these are a work in progress.


## Authors

The NREL Floating Wind Array Design team.