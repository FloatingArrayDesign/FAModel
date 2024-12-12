# -*- coding: utf-8 -*-
"""
This example is designed to walk through using the failure simulation package.
The failure probabilities and potential failure modes/effects are stored in excel 
files. The sample ones used in this example are stored in the famodel/failure
folder.

An initial failure is enacted based on the input criticality type. Then the array
watch circle is redone and a RAFT simulation is completed to determine if any failures
would result from the initial failure. If not, the next most critical failure is 
completed and the cycle continues until there are no more failure nodes possible.
"""
import os
import famodel
from famodel.failure.failureGraphs import failureGraph

# get directory and files for inputs
fam_directory = './Inputs/'
fam_yaml = 'OntologySample200m.yaml'
failure_directory = '../../famodel/failure/'
failure_dataFile = 'failureData.xlsx'
failure_dataSheet = 'bMMatch'
failure_probFile = 'failureProbabilities.xlsx'
failure_probSheet = 'Sheet3'

# switch to yaml input directory
os.chdir(fam_directory)

# initialize failure graphs
fg = failureGraph(fam_yaml)

# switch to failure directory
os.chdir(failure_directory)

# create failure graph
fg.create_failureGraph(failure_dataFile, failure_dataSheet, 
                       failure_probFile, failure_probSheet)

# run failure simulation
fg.run_failure_simulation()


