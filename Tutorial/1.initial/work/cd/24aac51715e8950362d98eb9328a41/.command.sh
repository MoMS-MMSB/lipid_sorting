#!/bin/bash -ue
gmx grompp -c input.gro -p input.top -f em.mdp -o em.tpr -maxwarn 1
gmx mdrun -s em.tpr -c em.gro -nt 16 -rdd 1.6
