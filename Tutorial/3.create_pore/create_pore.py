import os
from src import pore_modules
import argparse

parser = argparse.ArgumentParser(description = "Take an input gro structure, from which you can select multiple residues using -r. These residues will have restraint coordinates for x, y and z written. These are written at 0 by default, but this can be changed using the -x, -y and -z flags.")

## generation of pore
parser.add_argument("-c", "--input", help = "<.gro> input structure file", required = True)
parser.add_argument("-o", "--output", help = "<.gro> desired output .gro file (default = 'output.gro') (opt.)", default = "output.gro")
parser.add_argument("-r", "--radius", type=float, default=2)
parser.add_argument("--residues", nargs="*", type=str, help = "<str> residues to restrain (call multiple times if needed)", required = True)
parser.add_argument("-a", "--axis", help = "pore axis", default = "x")
parser.add_argument("--coords", help = "3-dimensional list of integers or floats for precising pore location. If not specified, the coordinates are set in the center of the box")
parser.add_argument("--renumber", type=bool, default=True, help="Sort the residues such that residues of the same type are grouped; set to True by default as it causes conflicts with the (current version of) topology rewriting. Gro files created in TS2CG are ordered by residue and leaflet, unfortunately we lose this specificity here.")

## generate new topology
parser.add_argument("-p", "--topology", help = "<.top> input topology file")
parser.add_argument("-po", "--top_output", help = "<.top> desired output topology file (default = 'modified.top') (opt.)", default = "modified.top")

## create restraints file
parser.add_argument("--restraints", action ='store_true', help = "Generate a corresponding restraints file.")
parser.add_argument("--restraints_out", default="restraints.gro")
args = parser.parse_args()

reader = open(args.input, "r").read()

# ASSERT STATEMENTS!
if any(residue not in reader for residue in args.residues):
    raise ValueError("One or more of your selected residues was not found in the input .gro file.")

if not args.coords:
    centerXYZ = [(vector/2) for vector in pore_modules.boxVectors(args.input)[1]]
else:
    centerXYZ = args.coords

# 1. Generate Pore
pore_modules.create_pore(args.input, args.output, args.axis, args.radius, args.residues, centerXYZ)
if args.renumber == True: 
    pore_modules.reorder_gro(args.output, pore_modules.mol_list(args.output))
    os.rename("renumbered.gro", args.output)

# 2. Generate Topology
if args.topology:
    pore_modules.rewrite_top(args.topology, pore_modules.num_unique_molecules(args.output, pore_modules.mol_list(args.output)), args.top_output)

# 3. Generate Restraints
if args.restraints:
    pore_modules.generate_restraints(args.output, args.residues, args.restraints_out, centerXYZ)