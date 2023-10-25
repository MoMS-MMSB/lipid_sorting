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

def lipids_per_tubule_leaflet(universe, residue_atom_dict, axis="x", center=False):
    if center == False:
        center = non_water_COG(universe)
 
    # assert universe is mda object containing one ts?
    assert isinstance(center, list) and len(center) == 3 and all(isinstance(x, float) for x in center), "center is not a list containing three floats"
    assert axis in ["x", "y", "z"], "axis is not one of 'x', 'y', or 'z'"
    assert all(isinstance(key, str) and isinstance(value, list) and len(value) == 2 and all(isinstance(item, str) for item in value) for key, value in residue_atom_dict.items()), "Invalid dictionary format"

    if axis == "x":
        dim1, dim2 = center[1], center[2]
    elif axis == "y":
        dim1, dim2 = center[0], center[2]
    elif axis == "z":
        dim1, dim2 = center[0], center[1]
    
