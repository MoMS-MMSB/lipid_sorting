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

# TUTORIAL

The following protocol will allow you to generate *your own* tubules from scratch - all you know is your desired composition, and nothing else. You can follow the exact steps in the tutorial, for which the files are included; or, you can change the composition of the tubule yourself, and try your own system. The protocol remains the same.

If you end up using this tutorial, you should cite our book chapter (coming soon) and also TS2CG, the paper for which is found [here](https://rdcu.be/drEQr): 

You will need: 
- Access to a terminal
- This repository; run ```git clone --recurse-submodules https://github.com/MoMS-MMSB/lipid_sorting.git``` in a terminal. We call recurse-submodules as this project uses another repository we developed based on [TS2CG](https://github.com/marrink-lab/TS2CG1.1) to help with the initial tubule setup.
- A conda environment for this repository. This gives you access to the in-house scripts generated for setup and analysis, and all their required python package dependencies. Notably, this includes Nextflow, the software used to run the workflow for initial structure generation. This is done by running ```conda env create --name lipid-sorting --file=environments.yml```
- A working GROMACS install. This entire process was developed using GROMACS/2023.1, the installation instructions for which can be found [here](https://manual.gromacs.org/documentation/2023.1/install-guide/index.html)

Once you've cloned the repository and created the environment, we can begin. Activate the conda environment:
 ```
 conda activate lipid-sorting
 ```

## 1. Run the pipeline
Create a new folder, and enter. For this tutorial, we will call this folder 1.initial/

Create a file called generate.str, which will contain all the important information the creation of your tubule.

Here's an example I used to create the simple POPC/POPE lipid mixtures with a radius of 10nm and length of 10nm, specified in the book chapter as the system with reduced cost. We will use this system for the sake of the tutorial, as it is probably the easiest to run on most machines, and yields results the "fastest"

```
[Lipids List]
Domain 0
POPC    0.5  0.5   0.70
POPE    0.5  0.5   0.70
End

[Shape Data]
ShapeType Cylinder
Box 30 30 10
Thickness 2
Radius 10
End
```

We can then run the nextflow pipeline, by running:


```
nextflow run ../../TS2CG-Setup-Pipeline/main.nf 
```

Based on your system and computer, the time taken can vary. The benefit of nextflow is that if any one step returns an error, it will at which steo the error occurred, and can even be resumed from the prior step by adding ```--resume``` to the command line call once the errors have been resolved.

After the workflow has run, files named "eq.gro", "index.ndx", and "topol.top" will be generated in the results/ subfolder of 1.initial/

## 2. Equilibrate (without pores)

Go up one level and make a new folder; 2.eq_nopore/

In the corresponding tutorial file, you will find an .mdp file to run the first "proper" equilibration. The example here is for 500ns, and will most definitely require a cluster of sorts. 

Run the gromacs pre-processing function gmx grompp:

```
gmx grompp -f eq-nopore.mdp -c ../1.initial/results/eq.gro -p ../1.initial/results/topol.top -n ../1.initial/results/index.ndx -o eq_nopore.tpr
```

Run the system however you can, eventually yielding "eq_nopore.gro". An example gro file has been provided in the relevant tutorial files, for those who simply wish to follow the process without running the actual simulations (a smart idea!)

## 3. Generate the pores
Let's make a new folder, 3.create_pore/

As outlined in the chapter, pores are *really* permitted via (i) their initial creation and (ii) holding them open in subsequent simulations using flat-bottomed potentials.

### **i)** Create pores
In the respective tutorial folder is a script, generate_pores.py. If you haven't already, enter into this new folder. Then run:

```
python create_pore.py -c ../2.eq_nopore/eq_nopore.gro -o create_pore_x.gro --axis x --radius 2.5 -p ../1.initial/results/topol.top -po modified_x.top --residues POPC POPE

python create_pore.py -c create_pore_x.gro -o create_pore_xy.gro --axis y --radius 2.5 -p modified_x.top -po modified_xy.top --restraints --residues POPC POPE
```

We're first generating a pore in the x- dimension, then in the y-dimension; both with radius 2.5nm. Since we're not specifying coordinates (done with the command line flag --coords), they'll automatically be set to the middle of the box - here, at 14.799, 14.799, and 05.000 for x, y, and z, respectively.

One last thing for the sake of being tidy is to use the GROMACS ```gmx editconf``` command to renumber the residues in the new gro file we've just created. Since creation of pores deletes residues and reorganises them in the gro file, naturally the gro file will be slightly unordered - this isn't necessarily a problem, but for the sake of being tidy we should run:

```
gmx editconf -f create_pore_xy.gro -resnr 1 -o init.gro
```

And with that, we have our starting structure!

### **ii)** Define flat-bottomed potentials
The definition of flat-bottomed potentials (FBPs) is the most hands-on step (and potentially, the most tedious...) of the entire process.

In order to activate an FBP in GROMACS, three steps are required:

1) Creation of a restraints.gro file (which we have just done), which defines the **geometric center of the FBP**
2) Modification of the molecule's topology file to define the **shape, size and force** of the FBP
3) Actually activating the FBP during the molecular dynamics run, by defining a flag in the .mdp file.

A more in-depth explanation on restraints files is found in [S1.](#s1-restraints-file)
In saying that, a comparison of our structure file ```init.gro```, our restraints file ```restraints.gro```, and the knowledge that we have defined the center of restraint for all pores as being the center of the box at 14.799, 14.799, and 05.000, should it rather clear what we are doing.

Modification of the .itp file needs to occur on a system-by-system basis, and occurs based on the FBP requirements for a given system. Included here is an .itp which contains only the lipids POPC and POPE, with a section added:

```
#ifdef POSRES_PL
; Flat-bottomed position restraint for each PL
[ position_restraints ]
; numatoms  functype  g   r   k
;                       (nm) (kJ mol−1nm−2)
       05      2      6  -2.5   5000
       06      2      6  -2.5   5000
       07      2      6  -2.5   5000
       08      2      6  -2.5   5000
       09      2      6  -2.5   5000
       10      2      6  -2.5   5000
       11      2      6  -2.5   5000
       12      2      6  -2.5   5000
       05      2      7  -2.5   5000
       06      2      7  -2.5   5000
       07      2      7  -2.5   5000
       08      2      7  -2.5   5000
       09      2      7  -2.5   5000
       10      2      7  -2.5   5000
       11      2      7  -2.5   5000
       12      2      7  -2.5   5000
#endif

```
Here, we are defining an FBP with an "if" statement. This if statement activates the FBPs we define - two, here; one for each dimension in which we have pores.

For a longer explanation of this section, consult [S2.](#s2-fbps-in-the-topology-file), which explains each column.

We also need to reference the correct .itp in the system topology file, as opposed to the old one. Run:
```
cp modified_xy.top system.top
```
delete the line "TK"

And replace it with `#include m3_POPC_POPE_pore.itp`
## 4. Perform production run, holding the pores open



## S1. Restraints file

## S2. FBPs in the topology file
Consider this section in the topology file for both POPC and POPE:
```
#ifdef POSRES_PL
; Flat-bottomed position restraint for each PL
[ position_restraints ]
; numatoms  functype  g   r   k
;                       (nm) (kJ mol−1nm−2)
       05      2      6  -2.5   5000
       06      2      6  -2.5   5000
       07      2      6  -2.5   5000
       08      2      6  -2.5   5000
       09      2      6  -2.5   5000
       10      2      6  -2.5   5000
       11      2      6  -2.5   5000
       12      2      6  -2.5   5000
       05      2      7  -2.5   5000
       06      2      7  -2.5   5000
       07      2      7  -2.5   5000
       08      2      7  -2.5   5000
       09      2      7  -2.5   5000
       10      2      7  -2.5   5000
       11      2      7  -2.5   5000
       12      2      7  -2.5   5000
#endif

```
Here, column one denotes the bead/atom number to which the FBP is applied. Both POPC and POPE have 12 beads, but we don't apply the FBP to the first four, which correspond to the lipid heads; this is how we allow the lipids to sort via flip-flop.

The second column defines the function type: for ```[ position_restraints ]```, function type 2 refers to a FBP, simple!

The next column refers to the **geometry** of the FBP. You will notice that we are applying two types of FBP, then; one of type 6 and one of type 7. Type 6 refers to a cylinder which has it's length in the x-dimension, and type 7 is a cylinder with length in y. For more information about possible geomtries, consult the [GROMACS manual](https://manual.gromacs.org/current/reference-manual/functions/restraints.html#flat-bottomed-position-restraints).

The center of **both** of these tubules comes from the center of geometry defined on the particles in the ```restraints.gro``` file. 

> Note: In my experience, this is where things most often go wrong: oftentimes, my systems explode because I've missed an atom or so in the restraints file, defining it's coordinates elsewhere, at which point GROMACS applies the massive FBP force penalty for the atom being out of place, destroying everything. If in doubt, make sure that every single instance of a molecule with FBPs defined has the restraint coordinates correctly specified in `restraints.gro`!

The next column defines the **radius** of the FBP. Here, we're using `-2.5` - the minus sign denotes an **inverse** FBP, keeping all particles outside of the cylinder's radius, as opposed to inside.

Lastly, we have the column denoting the force constant k(kJ/mol^-1/nm^2), applied to any 