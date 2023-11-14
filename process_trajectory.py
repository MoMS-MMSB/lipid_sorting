import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from src import modules
import MDAnalysis as mda
import matplotlib.pyplot as plt

gro = "/data1/jackson/MD/Membrane_Systems/Tubules/DOPC_POPC/r10/production/noW.gro"
trj = "/data1/jackson/MD/Membrane_Systems/Tubules/DOPC_POPC/r10/production/noW.xtc"

colours = {"DOPC" : "red", "POPC":"blue", "POPE":"orange"}

universe = mda.Universe(gro,trj)
print(universe)
modules.process_trajectory(universe, colours=False, skip=100)