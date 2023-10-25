import MDAnalysis as mda
from src import modules

u = mda.Universe("test_systems/ER8_r10.gro")
residues = modules.residue_names(u)

modules.first_and_last_atoms(u, residues)