# LIPID SORTING IN TUBULES

## Introduction
Coming soon! Once this gets published, I can write more here :)

## Method
The pipeline to set up tubules is based on TS2CG. It will soon be added, in some form, into this directory.

Also, this directory will eventually feature the script to generate pores and restraints files, followed by an example workflow. Automating this entire process can be tricky, as longer systems will likely require use of clusters or supercomputing resources, but expect example workflows once all the pieces are in one place.

Contained here are a series of python functions used for the analysis of the resulting production runs for this system. These modules can be found in /src/modules.py.

process_trajectory.py is the example script used to run the analysis on one or several systems. Depending on the size of the system, this can take a few hours. It is recommended to remove water from the .gro and trajectory files.

## Files
example_results contain the results for the figures used in the publication. Likewise, Figures/ contains various renders of each system, as well as the tcl scripts used to generate these in vmd. 

# What if I want to run my own?
The following protocol will allow you to generate *your own* tubules from scratch - all you know is your desired composition, and nothing else.

You will need: 
- Access to a terminal
- This repository; run ```git clone https://github.com/MoMS-MMSB/lipid_sorting.git``` in a terminal
- A conda environment for this repository. This gives you access to the in-house scripts generated for setup and analysis, and all their required python package dependencies. Notably, this includes Nextflow, the software used to run the workflow for initial structure generation. This is done by running ```conda env create --name lipid-sorting --file=environments.yml```
- A working GROMACS install. This entire process was developed using GROMACS/2023.1, the installation instructions for which can be found [here](https://manual.gromacs.org/documentation/2023.1/install-guide/index.html)

Once you've cloned the repository, created the environment, and entered into the directory, we can begin.

It's a good idea to create a seperate folder in which to begin. Enter that folder, and create a file called generate.str, which will contain all the important information the creation of your tubule.

Here's an example I used to create the simple POPC/POPE lipid mixtures with a radius of 10nm.

```
[Lipids List]
Domain 0
POPC    0.5  0.5   0.70
POPE    0.5  0.5   0.70
End

[Shape Data]
ShapeType Cylinder
Box 30 30 30
Thickness 2
Radius 10
End
```