from src import analysis_modules
import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
struct = "/data1/jackson/MD/Membrane_Systems/Tubules/POPC_POPE/r30/create_pore/create_pore_out.gro"
traj = "/data1/jackson/MD/Membrane_Systems/Tubules/POPC_POPE/r30/40us.xtc"

u = mda.Universe(struct, traj)

radii = np.array(analysis_modules.run_radius_parallel(u))/10
time = np.arange(0,len(radii))/100
print(np.average(radii))
df = pd.DataFrame({"Time (µs)":time, "Radius (nm)":radii})
df.to_csv("radius.PP.30.csv")
plt.plot(time, radii)
plt.xlabel('Time (µs)')
plt.ylabel('Radius (nm)')
plt.title ('Radius v. Time')
# plt.savefig("radius.10.png")
# plt.show()