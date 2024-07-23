from scipy.constants import pi, Boltzmann
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
                residue_first_last_atoms.append(([r.atoms[0].name])) # head
                if residue == "CHOL": # this could be handled smarter, as "if not lipid" somehow
                    residue_first_last_atoms.append(([r.atoms[-1].name]))
                else:
                    tail_dic = {} 
                    for atom in r.atoms[:]:
                        last_char = atom.name[-1]
                        if (last_char == "A") or (last_char == "B"): #replace with is srtring
                            if last_char not in tail_dic:
                                tail_dic[last_char] = []
                            tail_dic[last_char].append(atom.name)
                    tail_first_last = []
                    for tail in tail_dic:
                        tail_first_last.append(tail_dic[tail][-1])    
                    residue_first_last_atoms.append(tail_first_last)
                break
        # print(residue_first_last_atoms)
        residue_atom_dict[residue] = residue_first_last_atoms
    return(residue_atom_dict)  

def distance_between_points(point1, point2):
    """Calculate the distance between two two-dimensional points"""
    # assert (type(point1) == list) and (len(point1) == 2)
    # assert (type(point2) == list) and (len(point2) == 2)
    return (np.sqrt(((point2[0]-point1[0])**2) + ((point2[1] - point1[1])**2)))


def lipids_per_tubule_leaflet_parallel(frame_index, universe, axis="z", cutoff=5):
    """
    From an MDAnalysis universe and a dictionary containing residue names with their first and 
    last atoms, assign each residue (i.e. lipid) to either the inner or outer tubule leaflet,
    where the tubule is periodic in the dimension specified by variable 'axis'
    """

    universe.universe.trajectory[frame_index]

    residue_atom_dict = martini_lipid_tails(universe, residue_names(universe.select_atoms("not resname W NA CA CL ION TO DO TUBE")))

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

    for resname in residue_atom_dict:
        heads = np.array(universe.select_atoms(f'resname {resname} and name {residue_atom_dict[resname][0][0]}').positions[:])
        tails = []
        for tail_name in residue_atom_dict[resname][1]:
            tails.append(universe.select_atoms(f'resname {resname} and name {tail_name}').positions[:])
        tails = np.mean(tails, axis=0)

        head_coords = [heads[:, dims[0]], heads[:, dims[1]]]
        tail_coords = [tails[:, dims[0]], tails[:, dims[1]]]

        center_array = np.array([[tubule_center[0]] * len(heads), [tubule_center[1]] * len(heads)])

        head_dists = distance_between_points(head_coords, center_array)
        tail_dists = distance_between_points(tail_coords, center_array)

        outer_total = np.sum(head_dists > tail_dists + cutoff)
        inner_total = np.sum(tail_dists > head_dists + cutoff)

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

def process_trajectory_parallel(universe,  output="dataframe.csv", ncores=4, axis = "z"):
    """
    Perform inner/outer tubule analysis on entire trajectory, given here
    as MDAnalysis universe. Saves as a csv named by variable 'output'
    """
    assert ncores > 1
    assert isinstance(ncores, int)
    
    resnames = residue_names(universe.select_atoms("not resname W ION DO TO TUBE"))
    print(resnames)

    run_per_frame = partial(lipids_per_tubule_leaflet_parallel, universe=universe, axis=axis)

    frame_values = np.arange(universe.trajectory.n_frames)

    with Pool(ncores) as worker_pool:
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

def df_proportions(df):
    """
    Take a dataframe with odd columns being Outer leaflet, even columns being Inner leaflet,
    and return a dataframe with each entry as a proportion of their given leaflet.
    """
    # create inner/outer dfs
    df_outer = df.iloc[:, ::2]
    df_inner = df.iloc[:, 1::2]
    
    # add column which is sum of all rows
    df_outer["Outer Total"] = df_outer.sum(axis=1)
    df_inner["Inner Total"] = df_inner.sum(axis=1)

    # divide all rows by the total column
    df_outer_prop = df_outer.div(df_outer.iloc[:, -1], axis=0)
    df_inner_prop = df_inner.div(df_inner.iloc[:, -1], axis=0)

    # remove total column
    df_outer_prop = df_outer_prop.iloc[:,:-1]
    df_inner_prop = df_inner_prop.iloc[:,:-1]

    df_prop = pd.merge(df_outer_prop, df_inner_prop, right_index=True, left_index=True)
    return(df_prop)

def df_prop_to_plot(df, resnames, rolling=1, bg=False, title=False, colours=False, x_label="Frame"):
    """
    Plot a dataframe, with lipids taking their colours from a dictionary.
    Rolling average can be applied with input variable rolling.
    """
    for resname in resnames:
        if colours:
            plt.plot(df[resname + " Outer"].rolling(rolling).mean(), color=colours[resname], linewidth=3, linestyle="-", label= resname + " Outer")
            plt.plot(df[resname + " Inner"].rolling(rolling).mean(), color=colours[resname], linewidth=3, linestyle="--", label= resname + " Inner")
            if bg:
                plt.plot(df[resname + " Outer"], color=colours[resname], linewidth=3, alpha = 0.10)
                plt.plot(df[resname + " Inner"], color=colours[resname], linewidth=3, alpha = 0.10)
        else:
            plt.plot(df[resname + " Outer"].rolling(rolling).mean(), linestyle="-", label= resname + " Outer")
            plt.plot(df[resname + " Inner"].rolling(rolling).mean(), linestyle="--", label= resname + " Inner")
            if bg:
                plt.plot(df[resname + " Outer"], linewidth=3, alpha = 0.25)
                plt.plot(df[resname + " Inner"], linewidth=3, alpha = 0.25)
    plt.legend()
    plt.ylabel("Proportion in Leaflet")
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
    df_prop = df_proportions(df)
    # print(df_prop)
    df_prop_to_plot(df_prop, resnames, rolling, bg, title, colours, x_label)

    if out:
        plt.savefig(out)
    return(plt)

def calc_radius(universe, axis="z"):
    center = non_water_COG(universe)

    assert axis in ["x", "y", "z"], "axis is not one of 'x', 'y', or 'z'"

    if axis == "x":
        tubule_center = [center[1], center[2]]
        dims = [1,2]
    elif axis == "y":
        tubule_center = [center[0], center[2]]
        dims = [0,2]
    elif axis == "z":
        tubule_center = [center[0], center[1]]
        dims = [0,1]

    bead_dim1 = universe.select_atoms("name PO4").positions[:, dims[0]]
    bead_dim2 = universe.select_atoms("name PO4").positions[:, dims[1]]
    beads_points = np.square(np.abs(bead_dim2 - bead_dim1))
    center_points = abs(tubule_center[1] - tubule_center[0])**2
    radii = np.sqrt(beads_points + center_points)
    return(radii)

def calc_radius_parallel(frame_index, universe, axis="z"):
    universe.universe.trajectory[frame_index]

    center = non_water_COG(universe)

    assert axis in ["x", "y", "z"], "axis is not one of 'x', 'y', or 'z'"

    if axis == "x":
        tubule_center = [center[1], center[2]]
        dims = [1,2]
    elif axis == "y":
        tubule_center = [center[0], center[2]]
        dims = [0,2]
    elif axis == "z":
        tubule_center = [center[0], center[1]]
        dims = [0,1]

    bead_dim1 = universe.select_atoms("name PO4").positions[:, dims[0]]
    bead_dim2 = universe.select_atoms("name PO4").positions[:, dims[1]]
    beads_points = np.square(np.abs(bead_dim2 - bead_dim1))
    center_points = abs(tubule_center[1] - tubule_center[0])**2
    radii = np.sqrt(beads_points + center_points)
    return(np.average(radii))

def run_radius_parallel(universe, axis="z"):
    run_per_frame = partial(calc_radius_parallel, universe=universe, axis=axis)

    frame_values = np.arange(universe.trajectory.n_frames)

    with Pool(25) as worker_pool:
        result = worker_pool.map(run_per_frame, frame_values)
    return(result)


def calc_radius_slices_parallel(frame_index, universe, axis="z"):
    universe.universe.trajectory[frame_index]

    center = non_water_COG(universe)

    assert axis in ["x", "y", "z"], "axis is not one of 'x', 'y', or 'z'"

    if axis == "x":
        axis_idx = 0
        tubule_center = [center[1], center[2]]
        dims = [1,2]
    elif axis == "y":
        axis_idx = 1
        tubule_center = [center[0], center[2]]
        dims = [0,2]
    elif axis == "z":
        axis_idx = 2
        tubule_center = [center[0], center[1]]
        dims = [0,1]
    radii_slices = []
    for i in range(10): #10 being default number of slices:
        # get max length in tubule dir
        tube_length = universe.dimensions[axis_idx]
        slice_min = tube_length*(i/10)
        slice_max = tube_length*(i/10) + tube_length/10
        bead_dim1 = universe.select_atoms(f"name PO4 and prop {axis} > {slice_min} and prop {axis} < {slice_max}").positions[:, dims[0]]
        bead_dim2 = universe.select_atoms(f"name PO4 and prop {axis} > {slice_min} and prop {axis} < {slice_max}").positions[:, dims[1]]
        beads_points = np.square(np.abs(bead_dim2 - bead_dim1))
        center_points = abs(tubule_center[1] - tubule_center[0])**2
        radii = np.sqrt(beads_points + center_points)
        radii_slices.append(np.average(radii))
    avg_slices = np.average(radii_slices)
    std_slices = np.std(radii_slices)
    radii_slices.append(avg_slices)
    radii_slices.append(std_slices)
    return radii_slices


def run_radius_slices_parallel(universe, axis="z"):
    run_per_frame = partial(calc_radius_slices_parallel, universe=universe, axis=axis)

    frame_values = np.arange(universe.trajectory.n_frames)

    with Pool(25) as worker_pool:
        result = worker_pool.map(run_per_frame, frame_values)
    return(result)

def calc_f(pressure_tube_dir, pressure_bulk, box_l_1, box_l_2):

    pressure_diff = pressure_bulk - pressure_tube_dir

    force = pressure_diff * box_l_1 * box_l_2
    return(force)

def calc_force_edr(edr, axis="z", output="energy.csv", verbose=False):
    aux = mda.auxiliary.EDR.EDRReader(edr, convert_units=False)
    all_terms = aux.get_data(["Temperature", "Box-X", "Box-Y", "Box-Z", "Pres-XX", "Pres-YY", "Pres-ZZ"])
    if axis == "x":
        box_lengths = [all_terms["Box-Y"]/1e+9, all_terms["Box-Z"]/1e+9]
        pressure_tube_direction = all_terms["Pres-XX"]
        pressure_bulk = (all_terms["Pres-YY"] + all_terms["Pres-ZZ"])/2
    elif axis == "y":
        box_lengths = [all_terms["Box-X"]/1e+9, all_terms["Box-Z"]/1e+9]
        pressure_tube_direction = all_terms["Pres-YY"]
        pressure_bulk = (all_terms["Pres-XX"] + all_terms["Pres-ZZ"])/2
    elif axis == "z":
        box_lengths = [all_terms["Box-X"]/1e+9, all_terms["Box-Y"]/1e+9]
        pressure_tube_direction = all_terms["Pres-ZZ"]
        pressure_bulk = (all_terms["Pres-XX"] + all_terms["Pres-YY"])/2
    force = calc_f(pressure_tube_direction, pressure_bulk, box_lengths[0], box_lengths[1]) *100000 # why by 100000? understand units - pascal?
    df = pd.DataFrame({"Time (ns)":all_terms["Time"]/1000,
                       "Box Dim1 (m)":box_lengths[0],
                       "Box_Dim2 (m)":box_lengths[1],
                       "Tube Pressure (bar)": pressure_tube_direction,
                       "Bulk Pressure (bar)":pressure_bulk,
                       "Force (N)":force})
    if verbose:
        df = pd.DataFrame({"Time (ns)":all_terms["Time"]/1000,
                            "Box Dim1 (m)":box_lengths[0],
                            "Box_Dim2 (m)":box_lengths[1],
                            "Tube Pressure (bar)": pressure_tube_direction,
                            "Bulk Pressure (bar)":pressure_bulk,
                            "Force (N)":force})
    else:
        df = pd.DataFrame({"Time (ns)":all_terms["Time"]/1000,
                            "Force (N)":force})        

    df.to_csv(output)

def sorting_efficiency(csv, lipids, total_lipids = 1, output=False):
    df_in = pd.read_csv(csv)
    df_in = df_in.set_index(df_in.columns[0])
    deltas = []
    for lipid in lipids:
        delta = df_in[lipid + " Outer"] - df_in[lipid + " Inner"]
        deltas.append(delta)

    deltas = np.array(deltas).T
    sorting_efficiency = (deltas[:,0] - deltas[:,1])/total_lipids
    results_array = np.column_stack((deltas, sorting_efficiency))
    headers = lipids + ["Sorting Efficiency"]
    df_out = pd.DataFrame(results_array, columns=headers)
    df_out["Sorting Efficiency"].plot()
    if output:
        df_out.to_csv(output)
    return(df_out)
