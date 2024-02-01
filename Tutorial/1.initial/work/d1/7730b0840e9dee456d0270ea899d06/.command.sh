#!/bin/bash -ue
echo -e "del 1-99
a W NA CL CA
0 & ! 1
name 1 Solvent
name 2 Solute
q" | gmx make_ndx -f input.gro
