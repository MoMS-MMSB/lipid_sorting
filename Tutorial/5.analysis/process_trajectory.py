import os
import time
from src import analysis_modules
import MDAnalysis as mda
import argparse
import numpy as np 
import pandas as pd

parser = argparse.ArgumentParser(description = "Take an input gro structure, from which you can select multiple residues using -r. These residues will have restraint coordinates for x, y and z written. These are written at 0 by default, but this can be changed using the -x, -y and -z flags.")

## parse arguments
# required arguments
parser.add_argument("-f", "--trajectory", help = "<.xtc> GROMACS trajectory file", required = True)
parser.add_argument("-s", "--structure", help = "<.gro> GROMACS structure file", required=True)
parser.add_argument("-e", "--energy", help = "<.edr> GROMACS energy file", required=True)

# additional arguments
# axis
parser.add_argument("-os", "--output_sorting", help = "<.csv> leaflets sorting results in csv", default="leaflets.csv")
parser.add_argument("-or", "--output_radius", help = "<.csv> radius of system results in csv", default="radius.csv")
parser.add_argument("-of", "--output_force", help = "<.csv> force of system results in csv", default="force.csv")
args = parser.parse_args()

## generate MDAnalysis universe
u =  mda.Universe(args.structure, args.trajectory)

## calculate lipids per leaflet
print("calculating lipids per leaflet...")
analysis_modules.process_trajectory_parallel(u, output=args.output_sorting)
print("...done, saved as " + str(args.output_sorting))

## calculate radius
print("calculating radius...")
radii = np.array(analysis_modules.run_radius_parallel(u))/10 # divide by 10 to go from angstrom to nm
print("...done, saved as " + str(args.output_radius))
time = np.arange(0,len(radii))/100 # divide by 100 to get results in µs
df = pd.DataFrame({"Time (µs)":time, "Radius (nm)":radii})
df.to_csv(args.output_radius)

## calculate force
analysis_modules.calc_force_edr(args.energy, output=args.output_force)