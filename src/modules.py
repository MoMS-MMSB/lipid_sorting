import numpy as np
import MDAnalysis as mda

def non_water_COG(universe):
    """
    computes the center of geometry of all non-water and non-ion particles in the system
    """
    return(universe.select_atoms("not resname W NA CA CL").center_of_geometry(compound="group"))

u = mda.Universe("./test_systems/ER8_r10.gro")

def residue_names(universe):
    """
    return a list of all residue names.
    """
    print(universe)
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
    # assert statements
    return (np.sqrt(((point1[1]-point1[0])**2) + ((point2[1] - point2[0])**2)))

def lipids_per_tubule_leaflet(universe, residue_atom_dict=False, axis="x", center=False):
    if residue_atom_dict == False:
        residue_atom_dict = first_and_last_atoms(universe, residue_names(universe))
    
    if center == False:
        center = non_water_COG(universe)

    # assert universe is mda object containing one ts?
    # assert isinstance(center, list) and len(center) == 3, "center is not a list containing three floats"
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
        residue_totals = []
        for i in range(len(heads)):
            head_coords = [heads.positions[i][dims[0]], heads.positions[i][dims[1]]]
            tail_coords = [tails.positions[i][dims[0]], tails.positions[i][dims[1]]]
            head_dist = (distance_between_points(head_coords, tubule_center))
            tail_dist = (distance_between_points(tail_coords, tubule_center))
            if head_dist > tail_dist:
                outer_total += 1
            else:
                inner_total += 1
        print(resname, outer_total, inner_total)
        residue_totals.append(outer_total)
        residue_totals.append(inner_total)
        overall_totals.append(residue_totals)
        print(overall_totals)
    # print(overall_totals)