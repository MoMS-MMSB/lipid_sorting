import numpy as np
import pandas as pd
from tqdm import tqdm
import MDAnalysis as mda
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
from multiprocessing import Pool
from functools import partial

def non_water_COG(universe):
    """
    computes the center of geometry of all non-water and non-ion particles in the system
    """
    return(universe.select_atoms("not resname W NA CA CL").center_of_geometry(compound="group"))


def residue_names(universe):
    """
    return a list of all residue names.
    """
    residue_names = []
    for atom in universe.atoms:
        if atom.resname not in residue_names:
            residue_names.append(atom.resname)
    
    return(residue_names)

def first_and_last_atoms(universe, list_of_residues):
    """
    take a universe and a list of residues. Create a dictionary in which each key is the residue name and each value is a list of two strings corresponding to the names of the first and last atoms of that residue.
    """
    residue_atom_dict = {}
    for residue in list_of_residues:
        residue_first_last_atoms = []
        for r in universe.residues:
            if str(r.resname) == residue:
                residue_first_last_atoms.append((r.atoms[0].name))
                residue_first_last_atoms.append((r.atoms[-1].name))
                break
    
        residue_atom_dict[residue] = residue_first_last_atoms
    return(residue_atom_dict)

def martini_lipid_tails(universe, list_of_residues):
    """
    take a universe and a list of residues. Return a dictionary where the keys are the residues and the values
    are lists of three strings; the head bead, and the last bead of each tail.
    This works only in Martini, as the tail beads end in A or B for the seperate tails.
    """
    residue_atom_dict = {}
    for residue in list_of_residues:
        residue_first_last_atoms = []
        for r in universe.residues:
            if str(r.resname) == residue:
                residue_first_last_atoms.append((r.atoms[0].name)) # head
                tail_dic = {} 
                for atom in r.atoms[:]:
                    last_char = atom.name[-1]
                    if (last_char == "A") or (last_char == "B"): #replace with is srtring
                        if last_char not in tail_dic:
                             tail_dic[last_char] = []
                        tail_dic[last_char].append(atom.name)
                for tail in tail_dic:
                    residue_first_last_atoms.append(tail_dic[tail][-1])    
                break
        # print(residue_first_last_atoms)
        residue_atom_dict[residue] = residue_first_last_atoms
    return(residue_atom_dict)
    




def distance_between_points(point1, point2):
    """Calculate the distance between two two-dimensional points"""
    assert (type(point1) == list) and (len(point1) == 2)
    assert (type(point2) == list) and (len(point2) == 2)
    return (np.sqrt(((point1[1]-point1[0])**2) + ((point2[1] - point2[0])**2)))

def lipids_per_tubule_leaflet(universe, axis="z", residue_atom_dict=False, center=False):
    """
    From an MDAnalysis universe and a dictionary containing residue names with their first and 
    last atoms, assign each residue (i.e. lipid) to either the inner or outer tubule leaflet,
    where the tubule is periodic in the dimension specified by variable 'axis'
    """

    if residue_atom_dict == False:
        residue_atom_dict = first_and_last_atoms(universe, residue_names(universe.select_atoms("not resname W NA CA CL")))
    
    if center == False:
        center = non_water_COG(universe)

    assert axis in ["x", "y", "z"], "axis is not one of 'x', 'y', or 'z'"
    assert all(isinstance(key, str) and isinstance(value, list) and len(value) == 2 and all(isinstance(item, str) for item in value) for key, value in residue_atom_dict.items()), "Invalid dictionary format"

    if axis == "x":
        tubule_center = [center[1], center[2]]
        dims = [1,2]
    elif axis == "y":
        tubule_center = [center[0], center[2]]
        dims = [0,2]
    elif axis == "z":
        tubule_center = [center[0], center[1]]
        dims = [0,1]
    
    overall_totals = []

    for resname in residue_atom_dict:
        heads = universe.select_atoms(f'resname {resname} and name {residue_atom_dict[resname][0]}')
        tails = universe.select_atoms(f'resname {resname} and name {residue_atom_dict[resname][1]}')

        outer_total = 0
        inner_total = 0
        for i in range(len(heads)):
            head_coords = [heads.positions[i][dims[0]], heads.positions[i][dims[1]]]
            tail_coords = [tails.positions[i][dims[0]], tails.positions[i][dims[1]]]
            head_dist = (distance_between_points(head_coords, tubule_center))
            tail_dist = (distance_between_points(tail_coords, tubule_center))
            if head_dist > tail_dist:
                outer_total += 1
            else:
                inner_total += 1
        overall_totals.append(outer_total)
        overall_totals.append(inner_total)
    return(overall_totals)

# def lipids_per_tubule_leaflet_parallel(frame_index, universe, axis="z"):
#     """
#     From an MDAnalysis universe and a dictionary containing residue names with their first and 
#     last atoms, assign each residue (i.e. lipid) to either the inner or outer tubule leaflet,
#     where the tubule is periodic in the dimension specified by variable 'axis'
#     """

#     universe.universe.trajectory[frame_index]

#     residue_atom_dict = first_and_last_atoms(universe, residue_names(universe.select_atoms("not resname W NA CA CL")))

#     center = non_water_COG(universe)

#     assert axis in ["x", "y", "z"], "axis is not one of 'x', 'y', or 'z'"
#     assert all(isinstance(key, str) and isinstance(value, list) and len(value) == 2 and all(isinstance(item, str) for item in value) for key, value in residue_atom_dict.items()), "Invalid dictionary format"

#     if axis == "x":
#         tubule_center = [center[1], center[2]]
#         dims = [1,2]
#     elif axis == "y":
#         tubule_center = [center[0], center[2]]
#         dims = [0,2]
#     elif axis == "z":
#         tubule_center = [center[0], center[1]]
#         dims = [0,1]
    
#     overall_totals = []

#     for resname in residue_atom_dict:
#         heads = universe.select_atoms(f'resname {resname} and name {residue_atom_dict[resname][0]}')
#         tails = universe.select_atoms(f'resname {resname} and name {residue_atom_dict[resname][1]}')

#         outer_total = 0
#         inner_total = 0
#         for i in range(len(heads)):
#             head_coords = [heads.positions[i][dims[0]], heads.positions[i][dims[1]]]
#             tail_coords = [tails.positions[i][dims[0]], tails.positions[i][dims[1]]]
#             head_dist = (distance_between_points(head_coords, tubule_center))
#             tail_dist = (distance_between_points(tail_coords, tubule_center))
#             if head_dist > tail_dist:
#                 outer_total += 1
#             else:
#                 inner_total += 1
#         overall_totals.append(outer_total)
#         overall_totals.append(inner_total)
#     return(overall_totals)

def lipids_per_tubule_leaflet_parallel(frame_index, universe, axis="z", multiple_tails = False):
    """
    From an MDAnalysis universe and a dictionary containing residue names with their first and 
    last atoms, assign each residue (i.e. lipid) to either the inner or outer tubule leaflet,
    where the tubule is periodic in the dimension specified by variable 'axis'
    """

    universe.universe.trajectory[frame_index]

    if multiple_tails:
        residue_atom_dict = martini_lipid_tails(universe, residue_names(universe.select_atoms("not resname W NA CA CL")))
    else:
        residue_atom_dict = first_and_last_atoms(universe, residue_names(universe.select_atoms("not resname W NA CA CL")))

    center = non_water_COG(universe)

    assert axis in ["x", "y", "z"], "axis is not one of 'x', 'y', or 'z'"
    # assert all(isinstance(key, str) and isinstance(value, list) and len(value) == 2 and all(isinstance(item, str) for item in value) for key, value in residue_atom_dict.items()), "Invalid dictionary format"

    if axis == "x":
        tubule_center = [center[1], center[2]]
        dims = [1,2]
    elif axis == "y":
        tubule_center = [center[0], center[2]]
        dims = [0,2]
    elif axis == "z":
        tubule_center = [center[0], center[1]]
        dims = [0,1]
    
    overall_totals = []


    # use old version for cheaper analysis...
    for resname in residue_atom_dict:
        if multiple_tails:
            tails_string = " ".join(residue_atom_dict[resname][1:])

        outer_total = 0
        inner_total = 0
        resids = universe.select_atoms(f'resname {resname}').residues
        for resid in resids:
            resid_str = str(resid).split()[-1].rstrip(">")
            
            head_com = (universe.select_atoms(f'resid {resid_str} and name {residue_atom_dict[resname][0]}').center_of_mass())
            head_coords = [head_com[dims[0]], head_com[dims[1]]]
            # print(head_coords)
            if multiple_tails:
                tail_com = (universe.select_atoms(f'resid {resid_str} and name {tails_string}').center_of_mass())
            else:
                tail_com = (universe.select_atoms(f'resid {resid_str} and name {residue_atom_dict[resname][1]}').center_of_mass())
            tail_coords = [tail_com[dims[0]], tail_com[dims[1]]]
            head_dist = (distance_between_points(head_coords, tubule_center))
            tail_dist = (distance_between_points(tail_coords, tubule_center))

            if head_dist > tail_dist:
                outer_total += 1
            else:
                inner_total += 1
        overall_totals.append(outer_total)
        overall_totals.append(inner_total)
    return(overall_totals)

# root = "/data1/jackson/MD/Membrane_Systems/Tubules/"
# DOPC_POPC_10 = mda.Universe(root + "DOPC_POPC/r10/production/noW.gro", root + "DOPC_POPC/r10/production/noW.xtc")
# resname = residue_names(DOPC_POPC_10)
# print(martini_lipid_tails(DOPC_POPC_10, resname))
# print(lipids_per_tubule_leaflet_parallel(1, DOPC_POPC_10, axis="z", multiple_tails=True))

def results_to_df(array, resnames):
    """
    Take an array and list of resnames. For each resname, make a column for both Outer and Inner
    Leaflets. The array should match this dimensionality. Return a pandas dataFrame.
    """
    headers = []
    for resname in resnames: 
        headers.append(resname + " Outer"), headers.append(resname + " Inner")

    df = pd.DataFrame(array,columns=headers)
    df.index += 1
    return(df)


def process_trajectory(trajectory, skip=1, output="dataframe.csv", axis = "z"):
    """
    Perform inner/outer tubule analysis on entire trajectory, given here
    as MDAnalysis universe. Saves as a csv named by variable 'output'
    """
    resnames = residue_names(trajectory.select_atoms("not resname W"))

    sel = trajectory.atoms.select_atoms("not resname W", updating = True)
    trajectory_output = []

    for ts in trajectory.trajectory[::skip]:
        trajectory_output.append(lipids_per_tubule_leaflet(sel, axis))
    
    df = results_to_df(trajectory_output, resnames)
    df.to_csv(output)

def process_trajectory_parallel(universe,  output="dataframe.csv", axis = "z"):
    """
    Perform inner/outer tubule analysis on entire trajectory, given here
    as MDAnalysis universe. Saves as a csv named by variable 'output'
    """
    resnames = residue_names(universe.select_atoms("not resname W"))

    run_per_frame = partial(lipids_per_tubule_leaflet_parallel, universe=universe, axis=axis)

    frame_values = np.arange(universe.trajectory.n_frames)

    with Pool(25) as worker_pool:
        result = worker_pool.map(run_per_frame, frame_values)

    df = results_to_df(result, resnames)
    df.to_csv(output)

def df_to_plot(df, resnames, rolling=1, bg=False, title=False, colours=False, x_label="Frame"):
    """
    Plot a dataframe, with lipids taking their colours from a dictionary.
    Rolling average can be applied with input variable rolling.
    """
    for resname in resnames:
        if colours:
            plt.plot(df[resname + " Outer"].rolling(rolling).mean(), color=colours[resname], linewidth=3, linestyle="-", label= resname + " Outer")
            plt.plot(df[resname + " Inner"].rolling(rolling).mean(), color=colours[resname], linewidth=3,linestyle="--", label= resname + " Inner")
            if bg:
                plt.plot(df[resname + " Outer"], color=colours[resname], linewidth=3, alpha = 0.10)
                plt.plot(df[resname + " Inner"], color=colours[resname], linewidth=3, alpha = 0.10)
        else:
            plt.plot(df[resname + " Outer"].rolling(rolling).mean(), linestyle="--", label= resname + " Outer")
            plt.plot(df[resname + " Inner"].rolling(rolling).mean(), linestyle="-.", label= resname + " Inner")
            if bg:
                plt.plot(df[resname + " Outer"], linewidth=3, alpha = 0.25)
                plt.plot(df[resname + " Inner"], linewidth=3, alpha = 0.25)
    plt.legend()
    plt.ylabel("# Lipids")
    plt.xlabel(x_label)
    if title:
        plt.title(title)

    return(plt)

def trajectory_to_plot(trajectory, axis="z", colours=False, skip=1):
    """
    Perform several of the above steps all in one, all from a single trajectory.
    Not recommended for huge trajectories, as the files are not saved.
    """
    resnames = residue_names(trajectory.select_atoms("not resname W"))

    sel = trajectory.atoms.select_atoms("not resname W", updating = True)
    trajectory_output = []

    for ts in tqdm(trajectory.trajectory[::skip]):
        trajectory_output.append(lipids_per_tubule_leaflet(sel, axis))
    
    df = results_to_df(trajectory_output, resnames)

    if colours:
        df_to_plot(df, resnames, colours)
    else:
        df_to_plot(df, resnames)

    plt.show()

def trajectory_to_plot_jupyter(trajectory, rolling=1, title=False, colours=False, skip=1):
    """
    Does the same as above, but in jupyter.
    Not recommended for huge trajectories, as the files are not saved.
    """

    resnames = residue_names(trajectory.select_atoms("not resname W"))

    sel = trajectory.atoms.select_atoms("not resname W", updating = True)
    trajectory_output = []

    with tqdm(total=len(trajectory.trajectory[::skip])) as pbar:
        for ts in trajectory.trajectory[::skip]:
            trajectory_output.append(lipids_per_tubule_leaflet(sel))
            pbar.update(1)
    
    df = results_to_df(trajectory_output, resnames)

    df_to_plot(df, resnames, rolling, title, colours)

    plt.show()


def csv_to_plot(csv, resnames, rolling=1, bg=False, title=False, colours=False, index_scaling=1, x_label="Frame", out=False):
    """
    Given a csv file and list of resnames, create dataframes and plots in accordance with
    earlier function df_to_plot()
    """
    df = pd.read_csv(csv)
    df[df.columns[0]] = df[df.columns[0]]/index_scaling
    df = df.set_index(df.columns[0])
    # print(rolling)
    df_to_plot(df, resnames, rolling, bg, title, colours, x_label)

    if out:
        plt.savefig(out)
    return(plt)

def outer_df_proportions(df):
    """
    Take a dataframe with odd columns being Outer leaflet, even columns being Inner leaflet,
    and return a dataframe with each entry as a proportion of their given leaflet.
    """
    # create inner/outer dfs
    df_outer = df.iloc[:, ::2]
    # df_inner = df.iloc[:, 1::2]
    
    # add column which is sum of all rows
    df_outer["Outer Total"] = df_outer.sum(axis=1)
    # df_inner["Inner Total"] = df_inner.sum(axis=1)

    # divide all rows by the total column
    df_outer_prop = df_outer.div(df_outer.iloc[:, -1], axis=0)
    # df_inner_prop = df_inner.div(df_inner.iloc[:, -1], axis=0)

    # remove total column
    df_outer_prop = df_outer_prop.iloc[:,:-1]
    # df_inner_prop = df_inner_prop.iloc[:,:-1]
    return(df_outer_prop)

def df_prop_to_plot(df, resnames, rolling=1, bg=False, title=False, colours=False, x_label="Frame"):
    """
    Plot a dataframe, with lipids taking their colours from a dictionary.
    Rolling average can be applied with input variable rolling.
    """
    for resname in resnames:
        if colours:
            plt.plot(df[resname + " Outer"].rolling(rolling).mean(), color=colours[resname], linewidth=3, linestyle="-", label= resname + " Outer")
            if bg:
                plt.plot(df[resname + " Outer"], color=colours[resname], linewidth=3, alpha = 0.10)
        else:
            plt.plot(df[resname + " Outer"].rolling(rolling).mean(), linestyle="--", label= resname + " Outer")
            if bg:
                plt.plot(df[resname + " Outer"], linewidth=3, alpha = 0.25)
    plt.legend()
    plt.ylabel("# Lipids")
    plt.xlabel(x_label)
    if title:
        plt.title(title)

    return(plt)

def csv_to_prop_plot(csv, resnames, rolling=1, bg=False, title=False, colours=False, index_scaling=1, x_label="Frame", out=False):
    """
    Given a csv file and list of resnames, create dataframes and plots in accordance with
    earlier function df_to_plot()
    """
    df = pd.read_csv(csv)
    df[df.columns[0]] = df[df.columns[0]]/index_scaling
    df = df.set_index(df.columns[0])
    df_prop = outer_df_proportions(df)
    # print(df_prop)
    df_prop_to_plot(df_prop, resnames, rolling, bg, title, colours, x_label)

    # if out:
    #     plt.savefig(out)
    # return(plt)