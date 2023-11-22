### Creation of Pore:
######################################################################################
def boxVectors(gro):
    """Return dimensions of box as list of strings and list of floats"""
    with open(gro) as f:
        for line in f:
            pass
        last_line = line
    vectors_str = last_line.split()
    vectors_str = vectors_str[0:3]
    vectors_float = ([float(vector) for vector in vectors_str])
    return(vectors_str, vectors_float)

def parseGroLine(line):
    """Parse a gro line, returning variables according to fixed file format columns"""
    atomNum = line[0:5]
    atomType = line[5:10]
    resid = line[0:10]
    x = float(line[20:28])
    y = float(line[28:36])
    z = float(line[36:44])
    return(atomNum, atomType, resid, x, y, z)

def isWithinCircle(dimA, dimB, centerA, centerB, boxA, boxB, radius=2):
    """Function determining whether a given point (dimA,dimB) lies within a circle with center (centerA, centerB) with radius."""
    dx = dimA - centerA
    dy = dimB - centerB
    
    # taking care of periodic boundary conditions
    if dx > 0.5 * boxA: 
        dx -= boxA
    elif dx < 0.5 *-boxA:
        dx += boxA
    if dy > 0.5 * boxB: 
        dy -= boxB
    elif dy < 0.5 *-boxB:
        dy += boxB

    return((dx*dx + dy*dy) <= radius*radius)

def residuesInPore(input_gro, axis, pore_centerA, pore_centerB, radius, moltypes):
    """If a molecule in list moltypes is found with any atom in pore area, add to list residues_in_pore"""
    assert type(moltypes) == list
    residues_in_pore = []
    box_dims = boxVectors(input_gro)[1]
    with open(input_gro) as f:
        lines = f.readlines()[2:-1]
        for line in lines:
            atomNum, atomType, resid, x, y, z = parseGroLine(line)
            # axis definitions
            if axis == "x":
                dimA, dimB = y, z
                boxA, boxB = box_dims[1], box_dims[2]
            elif axis == "y":
                dimA, dimB == x, z
                boxA, boxB = box_dims[0], box_dims[2]
            elif axis == "z":
                dimA, dimB == x, y
                boxA, boxB = box_dims[0], box_dims[1]

            if (isWithinCircle(dimA, dimB, pore_centerA, pore_centerB, boxA, boxB, radius)) and (atomType.strip() in moltypes):
                if resid.strip() not in residues_in_pore:
                    residues_in_pore.append(resid.strip())
    return(residues_in_pore)

def deleteResiduesInList(input_gro, output_gro, residues_to_delete):
    """Delete atoms from input gro file that are within the list defined by residues_to_delete"""
    new_lines = []
    with open(input_gro) as f:
        lines = f.readlines()[:2]
        for line in lines:
            new_lines.append(line)
    with open(input_gro) as f:
        lines = f.readlines()[2:-1]
        for line in lines:
            resid = parseGroLine(line)[2]
            if resid.strip() not in residues_to_delete:
                new_lines.append(line)
    with open(input_gro) as f:
        new_lines.append(f.readlines()[-1])
    with open(output_gro, "w") as w:
        for line in new_lines:
            w.write(line)

def update_gro_particles(gro_file):
    """Take a given .gro file and edit the second line to correctly reflect the number of atoms in the system"""
    with open(gro_file) as f:
        lines = f.readlines()
        lines[1] = str(len(lines)-3) + "\n"
    
    with open(gro_file, "w") as w:
        for line in lines:
            w.write(line)

def create_pore(input_gro, output_gro, axis, radius, moltypes, centerXYZ=False):
    """
    Within a given input gro file, define a pore based on axis, pore center given as
    a list containing coordinates for x, y, z, and a list of molecules in moltypes to
    remove from the pore region.
    Coordinates are based on the axis; while all three must be provided, the two coordinates
    other than the axis specified are used for calculation of the pore.
    """
    assert (type(centerXYZ) == list) or (centerXYZ is False)
    assert type(moltypes) == list
    assert (axis == "x") or (axis == "y") or (axis =="z")
    assert (type(radius) == float) or (type(radius) == int)

    if centerXYZ == False:
        centerXYZ = [(vector/2) for vector in boxVectors(input_gro)[1]]
    
    if axis == "x":
        pore_centerA, pore_centerB = centerXYZ[1], centerXYZ[2]
    elif axis == "y":
        pore_centerA, pore_centerB = centerXYZ[0], centerXYZ[2]
    elif axis == "z":
        pore_centerA, pore_centerB = centerXYZ[0], centerXYZ[1]
    
    residues_in_pore = residuesInPore(input_gro, axis, pore_centerA, pore_centerB, radius, moltypes)
    deleteResiduesInList(input_gro, output_gro, residues_in_pore)
    update_gro_particles(output_gro)


### Rewriting Topology:
######################################################################################
def mol_list(gro):
    """Return a list of all molecules types in a given gro file."""
    mol_list = []
    with open(gro, 'r') as f:
        lines = f.readlines()[2:-11]
        for line in lines:
            if line[5:11] not in mol_list:
                mol_list.append(str(line[5:11]))
    return(mol_list)

def first_index_mol(gro, mol_list):
    """Return the index of the first molecule for each molecule type provided in mol_list"""
    assert type(mol_list) ==  list
    first_indices = []
    for mol in mol_list:
        with open(gro) as f:
            lines = f.readlines()[2:-1]
            for line in lines:
                if mol in line:
                    first_indices.append(str(line[0:11]))
                    break
    return(first_indices)

def num_unique_molecules(gro, mol_list):
    """Given a list of molecules, return the number of unique molecules for each."""
    mol_counts = []
    assert type(mol_list) ==  list
    for mol in mol_list:
        entry = []
        entry.append(mol.strip())
        last_seen = None
        count = 0
        with open(gro) as f:
            lines = f.readlines()[2:-1]
            for line in lines:
                if (mol in line) and (line[0:11] != last_seen):
                    last_seen = line[0:11]
                    count += 1
            entry.append(str(count))
        mol_counts.append(entry)
    return(mol_counts)

def rewrite_top(input_top, mol_counts, output_top):
    """Take the head #include statements from a topology, add updated molecule counts, and write."""
    writer = open(output_top, "w")
    with open(input_top, "r") as f:
        for line in f.readlines():
            if "molecules" not in line:
                writer.write(line)
            else:
                break
    writer.write("[ molecules ]\n; name    number\n")   
    for entry in mol_counts:
        writer.write("%-10s%6d" % (entry[0], int(entry[1])) + "\n")

### Generating Restraints file:
######################################################################################
def generate_restraints(input_gro, mol_types, restraints_gro="restraints.gro", restraintsXYZ=False):
    """
    Take an input gro, a list of molecules to be restrained, 
    and a given set of coordinates at which to restrain chosen molecules.
    If coordinates are not specified, they will default to half of the respective
    box coordinates, to stay in line with other programs in this module.
    """
    assert type(mol_types) == list

    if restraintsXYZ == False:
        restraintsXYZ = [(vector/2) for vector in boxVectors(input_gro)[1]]

    line_num = 0
    with open(restraints_gro, "w") as writer:
        for line in open(input_gro).readlines():
            if line_num < 2:
                writer.write(line)
            elif any(molecule in line for molecule in mol_types):
                restrained_line = [
                line[0:20],
                "{:06.3f}".format(restraintsXYZ[0]),
                "{:06.3f}".format(restraintsXYZ[1]),
                "{:06.3f}".format(restraintsXYZ[2]),
                ]
                writer.write("  ".join(restrained_line) + "\n")  
            elif line_num > 2:
                writer.write(line)
            line_num += 1