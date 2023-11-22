def boxVectors(gro):
    with open(gro) as f:
        for line in f:
            pass
        last_line = line
    vectors_str = last_line.split()
    vectors_str = vectors_str[0:3]
    vectors_float = ([float(vector) for vector in vectors_str])
    return(vectors_str, vectors_float)

def parseGroLine(line):
    atomName = line[0:5]
    atomType = line[5:10]
    resid = line[0:10]
    x = float(line[20:28])
    y = float(line[28:36])
    z = float(line[36:44])
    return(atomName, atomType, resid, x, y, z)

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
            atomName, atomType, resid, x, y, z = parseGroLine(line)
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
                    # print(line)
                    residues_in_pore.append(resid.strip())
    return(residues_in_pore)


def deleteResiduesInPore(input_gro, output_gro, residues_in_pore):
    counter = 0
    new_lines = []
    with open(input_gro) as f:
        lines = f.readlines()[:2]
        for line in lines:
            new_lines.append(line)
    with open(input_gro) as f:
        lines = f.readlines()[2:-1]
        for line in lines:
            resid = parseGroLine(line)[2]
            if resid.strip() not in residues_in_pore:
                new_lines.append(line)
                counter += 1
    with open(input_gro) as f:
        new_lines.append(f.readlines()[-1])
    with open(output_gro, "w") as w:
        for line in new_lines:
            w.write(line)

def update_gro_particles(gro_file):
    # change second line with file length -3
    print("hello")

## test zone
input_test = "/data1/jackson/MD/Membrane_Systems/Tubules/DOPC_POPC/r10/eq_nopore/eq.1.part0001.gro"
res_in_pore = residuesInPore(input_test, "x", 15, 0, 1, ["DOPC", "POPC"])
print(res_in_pore)
print("resid " + ' '.join(res_in_pore))
output_test="/data1/jackson/MD/Membrane_Systems/Tubules/DOPC_POPC/r10/pore_python_test.gro"
deleteResiduesInPore(input_test, output_test, res_in_pore)