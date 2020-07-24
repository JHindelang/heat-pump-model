# heat-pump-model
Parameter estimation based thermodynamic heat pump model

## Requirements

Python 3.6

Packages and dependencies can be installed from requirements.txt using pip.

## Description of scripts

run_parameter_estimation.py (executable)

<ol>This script runs the parameter estimation for all samples specified in the input file 'input_filenames_estimation.csv'.
The sample sets specified there, are searched in the folder 'input_data'.
The parameter estimation is run 10 times, with different random seeds.
The parameters obtained such are stored to the outputf iles in folder 'output'.
The output files are named the same as the input files, with the prefix 'param'.

</ol> run_simulation.py (executable)

<ol>This script runs simulation for all sample files specified in the input file 'input_filenames_simulation.csv'.
The file containig the parameters that are used for each simulation, needs also to be specified in this input file. 
The sample files and parameter files specified there, are searched in the folder 'input_data'.
The simulation is run for all sample files and all operating points specified in the sample files.
The reuslts are stored to the output files in folder 'output'.
The output files are named the same as the sample input files, with the prefix 'sim'.

</ol> run_validation.py (executable)

<ol>This script runs the validation for all sample files specified in the input file 'input_filenames_validation.csv'.
The file containig the parameters that are used for each validation, needs also to be specified in this input file. 
The sample files and parameter files specified there, are searched in the folder 'input_data'.
The validation is run for all sample files and samples specified in the sample files.
The reuslts are stored to the output files in folder 'output'.
The output files are named the same as the sample input files, with the prefix 'validation'.
During the validation the SSE for the sample set is computed.

</ol> model_functions.py (not executable)

<ol>This file contains the model functions and is required for running the other files.
