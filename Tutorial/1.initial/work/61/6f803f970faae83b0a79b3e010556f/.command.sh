#!/bin/bash -ue
gmx grompp -c input.gro -p input.top -f eq.mdp -o eq.tpr -maxwarn 1
gmx mdrun -s eq.tpr -c eq.gro -nt 16
