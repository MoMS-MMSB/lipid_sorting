# LIPID SORTING IN TUBULES

## Introduction
Coming soon! Once this gets published, I can write more here :)

## Method
The pipeline to set up tubules is based on TS2CG. It will soon be added, in some form, into this directory.

Also, this directory will eventually feature the script to generate pores and restraints files, followed by an example workflow. Automating this entire process can be tricky, as longer systems will likely require use of clusters or supercomputing resources, but expect example workflows once all the pieces are in one place.

Contained here are a series of python functions used for the analysis of the resulting production runs for this system. These modules can be found in /src/modules.py.

process_trajectory.py is the example script used to run the analysis on one or several systems. Depending on the size of the system, this can take a few hours. It is recommended to remove water from the .gro and trajectory files.

# Files
example_results contain the results for the figures used in the publication. Likewise, Figures/ contains various renders of each system, as well as the tcl scripts used to generate these in vmd. 