#!/bin/bash -ue
gmx solvate -cs /data1/jackson/TS2CG_pipeline/top/water.gro -cp input.gro -p topol.top -radius 0.21 -o solvated.gro
