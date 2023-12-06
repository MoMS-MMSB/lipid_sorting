# LIPID SORTING IN TUBULES

![](figures/Renders/POPC_POPE_r10_l10_pore/x_110_5deg_dof_notrj.gif)

## Introduction
Coming soon! Once this gets published, I can write more here :)

## What is in this Repository?
The pipeline to set up tubules is based on TS2CG. Contained in this repository
is a fork of a setup pipeline which can be used to generate tubules of desired size and 
composition.

Contained in the folder src/ is a modules file in python called pore_modules.py.
In here, you can find all the required functions to create a pore, update the GROMACS topology, and generate a .gro file for position restraints. Following creation of this file, the .itp file itself will need to be modified - for now, examples are provided; expect some further explanations in future updates. Unfortunately, modification of the .itp file for the pore restraints is the most "hands-on" step, and is highly system specific.

An example of this is found in the script **SCRIPT**

Following generation of all of these files, you will be ready to run a production run. Depending on the system size, these can take considerable amounts of time to run, especially locally. Use of dedicated computing resources such as a cluster or supercomputer are desireable, however if they are inaccessible, consider sticking to smaller systems and shorter simulations. 

Once the production run has completed, you will need to run the analysis. The resulting files will be saved as a .csv, after which you can visualise these. All relevant modules can be found in /src/analysis_modules.py.

process_trajectories.py is the example script used to run the analysis on one or several systems; here, showing the exact commands used to process the systems from the publication. It is recommended to remove water from the .gro and trajectory files.

## Files
example_results contain the results for the figures used in the publication. Likewise, Figures/ contains various renders of each system, as well as the tcl scripts used to generate these in vmd. 

## What if I want to run my own?
The following protocol will allow you to generate *your own* tubules from scratch - all you know is your desired composition, and nothing else.

If you end up using this tutorial, you should cite our book chapter (coming soon) and also TS2CG, the paper for which is found [here](https://rdcu.be/drEQr): 

You will need: 
- Access to a terminal
- This repository; run ```git clone --recurse-submodules https://github.com/MoMS-MMSB/lipid_sorting.git``` in a terminal. We call recurse-submodules as this project uses another repository we developed based on [TS2CG](https://github.com/marrink-lab/TS2CG1.1) to help with the initial tubule setup.
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


```
nextflow run ../TS2CG_Setup_Pipeline/main.nf --input ${input} --outDir ${pwd}/out
```

Based on your system and computer, the time taken can vary. The benefit of nextflow is that if any one step returns an error, it will at which steo the error occurred, and can even be resumed from the prior step by adding ```--resume``` to the command line call once the errors have been resolved.

After the workflow has run, files named "workflow_out.gro" and "workflow_out.top" will be generated.

# **TUTORIAL TO BE CONTINUED**