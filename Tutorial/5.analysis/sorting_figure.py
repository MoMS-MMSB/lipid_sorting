import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src import analysis_modules

sorting = pd.read_csv("leaflets.csv",index_col=0)
radii = pd.read_csv("radius.csv")
force = pd.read_csv("force.csv")

colours = {"POPC": "#6a6adf",
           "POPE": "#ffb753"}

sorting_proportions = analysis_modules.df_proportions(sorting)
sorting_proportions.insert(0, "Time (µs)", np.arange(0, sorting_proportions.shape[0])/100)
print(sorting_proportions)

def inner_outer_plots(axs, df, resnames, colours):
    axs.plot(df["Time (µs)"], df[resnames[0] + " Outer"].rolling(50).mean(), color=colours[resnames[0]], linewidth=3, linestyle="-", label= resnames[0] + " Outer")
    axs.plot(df["Time (µs)"], df[resnames[0] + " Inner"].rolling(50).mean(), color=colours[resnames[0]], linewidth=3, linestyle="--", label= resnames[0] + " Inner")
    axs.plot(df["Time (µs)"], df[resnames[1] + " Outer"].rolling(50).mean(), color=colours[resnames[1]], linewidth=3, linestyle="-", label= resnames[1] + " Outer")
    axs.plot(df["Time (µs)"], df[resnames[1] + " Inner"].rolling(50).mean(), color=colours[resnames[1]], linewidth=3, linestyle="--", label= resnames[1] + " Inner")
    axs.plot(df["Time (µs)"], df[resnames[0] + " Outer"], color=colours[resnames[0]], linewidth=3, linestyle="-", label= resnames[0] + " Outer", alpha=0.15)
    axs.plot(df["Time (µs)"], df[resnames[0] + " Inner"], color=colours[resnames[0]], linewidth=3, linestyle="--", label= resnames[0] + " Inner", alpha=0.15)
    axs.plot(df["Time (µs)"], df[resnames[1] + " Outer"], color=colours[resnames[1]], linewidth=3, linestyle="-", label= resnames[1] + " Outer", alpha=0.15)
    axs.plot(df["Time (µs)"], df[resnames[1] + " Inner"], color=colours[resnames[1]], linewidth=3, linestyle="--", label= resnames[1] + " Inner", alpha=0.15)

fig, axs =  plt.subplots(3,1, figsize=(12, 16))
fig.tight_layout(pad=2.0)
fig.set_size_inches(12, 10)
plt.rc('font', size=10)          # controls default text sizes
plt.rc('axes', titlesize=14)     # fontsize of the axes title
plt.rc('axes', labelsize=10)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=10)    # fontsize of the tick labels
plt.rc('ytick', labelsize=10)    # fontsize of the tick labels

plt.subplots_adjust(hspace=0.4, wspace=0.7)

# sorting figure
inner_outer_plots(axs[0], sorting_proportions, ["POPC", "POPE"], colours)
axs[0].yaxis.set_ticks(np.arange(0.4, 0.61, 0.05))
axs[0].set_ylim((0.4, 0.6))
axs[0].set_title("Lipid sorting vs. time")
axs[0].set_ylabel("Proportion in leaflet")
axs[0].set_xlabel("time (µs)")
axs[0].legend(bbox_to_anchor=(1.1, 1.05),labels=["POPC Outer", "POPC Inner", "POPE Outer", "POPE Inner"], borderaxespad=0.)


axs[1].plot(radii["Time (µs)"], radii["Radius (nm)"].rolling(50).mean(), linewidth=3, color="#6a6adf")
axs[1].set_title("Radius vs. time")
axs[1].set_ylabel("Radius (nm)")
axs[1].set_xlabel("time (µs)")

axs[2].plot(force["Time (ns)"]/1000, force["Force (N)"].rolling(20000).mean(), color="#ffb753")
axs[2].set_title("Force vs. time")
axs[2].set_ylabel("Force (N)")
axs[2].set_xlabel("time (µs)")

plt.savefig("analysis", bbox_inches='tight')
plt.show()