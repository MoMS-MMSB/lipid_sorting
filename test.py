import MDAnalysis as mda
from src import modules

u = mda.Universe("test_systems/ER8_r10.gro").select_atoms("not resname W NA CL CA")
residues = modules.residue_names(u)

resname = "DOPC"
head = "NC3"
tail = "C4A"

modules.lipids_per_tubule_leaflet(u)