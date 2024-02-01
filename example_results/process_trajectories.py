import os
import time
from src import analysis_modules
import MDAnalysis as mda

root = "/data1/jackson/MD/Membrane_Systems/Tubules/"
output = "./"
if not os.path.exists(output):
    os.makedirs(output)
    print("Created output file: " + output)

DOPC_POPC_10 = mda.Universe(root + "DOPC_POPC/r10/create_pore/create_pore_out.gro", 
                 root + "DOPC_POPC/r10/40us.xtc")

DOPC_POPC_20 = mda.Universe(root + "DOPC_POPC/r20/create_pore/create_pore_out.gro", 
                 root + "DOPC_POPC/r20/40us.xtc")

DOPC_POPC_30 = mda.Universe(root + "DOPC_POPC/r30/create_pore/create_pore_out.gro", 
                 root + "DOPC_POPC/r30/40us.xtc")

POPC_POPE_10 = mda.Universe(root + "/POPC_POPE/r10/create_pore/create_pore_out.gro", 
                 root + "POPC_POPE/r10/40us.xtc")

POPC_POPE_10_l10 = mda.Universe(root + "/POPC_POPE/r10_l10/create_pore/create_pore_xy.gro", 
                 root + "POPC_POPE/r10_l10/40us.xtc")

POPC_POPE_20 = mda.Universe(root + "POPC_POPE/r20/create_pore/create_pore_out.gro", 
                 root + "POPC_POPE/r20/40us.xtc")

POPC_POPE_30 = mda.Universe(root + "POPC_POPE/r30/create_pore/create_pore_out.gro", 
                 root + "POPC_POPE/r30/40us.xtc")

# systems = [POPC_POPE_10_l10, DOPC_POPC_10, DOPC_POPC_20, DOPC_POPC_30, POPC_POPE_10, POPC_POPE_20, POPC_POPE_30]
systems = [POPC_POPE_10_l10]
system_outputs = {DOPC_POPC_10 : "DOPC_POPC_10.csv",
                   DOPC_POPC_20 : "DOPC_POPC_20.csv",
                   DOPC_POPC_30 : "DOPC_POPC_30.csv",
                   POPC_POPE_10 : "POPC_POPE_10.csv",
                   POPC_POPE_10_l10 : "POPC_POPE_10_10.csv",
                   POPC_POPE_20 : "POPC_POPE_20.csv",
                   POPC_POPE_30 : "POPC_POPE_30.csv"
                   }


for system in systems:
    print('\n\nProcessing: ' + system_outputs[system].split(".")[0] + "\n...")
    start = time.time()
    # run analysis
    analysis_modules.process_trajectory_parallel(system, output=(output + system_outputs[system]))
    end = time.time()
    print("Results saved in " + system_outputs[system])
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    print("Completed in {:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds))