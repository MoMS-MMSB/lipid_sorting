#!/bin/bash -ue
gmx solvate -cs /data2/jackson/lipid_sorting/lipid_sorting_in_tubules/TS2CG-Setup-Pipeline/top/water.gro -cp input.gro -p topol.top -radius 0.21 -o solvated.gro
