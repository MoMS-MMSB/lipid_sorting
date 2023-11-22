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

def lipids_per_tubule_leaflet_parallel(frame_index, universe, axis="z"):
    """
    From an MDAnalysis universe and a dictionary containing residue names with their first and 
    last atoms, assign each residue (i.e. lipid) to either the inner or outer tubule leaflet,
    where the tubule is periodic in the dimension specified by variable 'axis'
    """

    universe.universe.trajectory[frame_index]

    residue_atom_dict = first_and_last_atoms(universe, residue_names(universe.select_atoms("not resname W NA CA CL")))

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

def df_to_plot(df, resnames, rolling=1, title=False, colours=False, x_label="Frame"):
    """
    Plot a dataframe, with lipids taking their colours from a dictionary.
    Rolling average can be applied with input variable rolling.
    """
    for resname in resnames:
        if colours:
            plt.plot(df[resname + " Outer"].rolling(rolling).mean(), color=colours[resname], linestyle="--", label= resname + " Outer")
            plt.plot(df[resname + " Inner"].rolling(rolling).mean(), color=colours[resname], linestyle="-.", label= resname + " Inner")
        else:
            plt.plot(df[resname + " Outer"].rolling(rolling).mean(), linestyle="--", label= resname + " Outer")
            plt.plot(df[resname + " Inner"].rolling(rolling).mean(), linestyle="-.", label= resname + " Inner")
    
    plt.legend()
    plt.ylabel("# Lipids")
    plt.xlabel(x_label)
    if title:
        plt.title(title,fontsize=40)

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


def csv_to_plot(csv, resnames, rolling=1, title=False, colours=False, index_scaling=1, x_label="Frame", out=False):
    """
    Given a csv file and list of resnames, create dataframes and plots in accordance with
    earlier function df_to_plot()
    """
    df = pd.read_csv(csv)
    df[df.columns[0]] = df[df.columns[0]]/index_scaling
    df = df.set_index(df.columns[0])
    df_to_plot(df, resnames, rolling, title, colours, x_label)

    if out:
        plt.savefig(out)
    plt.show()