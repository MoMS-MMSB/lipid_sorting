def boxVectors(gro):
    with open(gro) as f:
        for line in f:
            pass
        last_line = line
    vectors_str = last_line.split()
    vectors_str = vectors_str[0:3]
    vectors_float = ([float(vector) for vector in vectors_str])
    print(vectors_str, vectors_float)

def parseGroLine(line):
    atomName = line[0:5]
    atomType = line[5:10]
    resid = line[0:10]
    x = float(line[20:28])
    y = float(line[28:36])
    z = float(line[36:44])
    return(atomName, atomType, resid, x, y, z)

def isWithinCircle(dimA, dimB, centerA, centerB, radius=2):
    """Function determining whether a given point (dimA,dimB) lies within a circle with center (centerA, centerB) with radius."""
    dx = dimA - centerA
    dy = dimB - centerB
    
    # taking care of periodic boundary conditions
    if dx > 0.5 * centerA: 
        dx -= centerA
    elif dx < 0.5 *-centerA:
        dx += centerA
    if dy > 0.5 * centerB: 
        dy -= centerB
    elif dy < 0.5 *-centerB:
        dy += centerB

    return((dx*dx + dy*dy) <= radius*radius)

def residuesInPore(input_gro, axis, pore_centerA, pore_centerB, radius, moltypes):
    """If a molecule in list moltypes is found with any atom in pore area, add to list residues_in_pore"""
    residues_in_pore = []
    with open(input_gro) as f:
        lines = f.readlines()[2:-1]
        for line in lines:
            atomName, atomType, resid, x, y, z = parseGroLine(line)
            # axis definitions
            if axis == "x":
                dimA, dimB = y, z
            elif axis == "y":
                dimA, dimB == x, z
            elif axis == "z":
                dimA, dimB == x, y

            if isWithinCircle(dimA, dimB, pore_centerA, pore_centerB, radius):
                if resid not in residues_in_pore:
                    residues_in_pore.append(resid)
    print(len(residues_in_pore))

def update_gro_particles(gro_file):
    # change second line with file length -3

    
## test zone
boxVectors("/data1/jackson/MD/Membrane_Systems/Tubules/DOPC_POPC/r10/production/noW.gro")
parseGroLine("    1DOPC   NC3    1  25.485   9.803  14.523  0.0730  0.2954  0.0083")
isWithinCircle(2,2,30,30)
isWithinCircle(1,1,30,30)
residuesInPore("/data1/jackson/MD/Membrane_Systems/Tubules/DOPC_POPC/r10/production/noW.gro", "x", 15, 15, 2, "DOPC")