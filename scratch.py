# To compare to Lambda output
"""
import constants as c
import utils, p.p_bar, fft.fft, fft.n, rotation.utils
from p.pdb_parser import parseLocal
import os
from math import pi as pi
import numpy as np

import time


def run_stage(constants, indices):

    # Does the ... of ...
    # constants: scan or discrimination stage constants
    # indices: those to do
    # return: Python list of all [k_max, vox_shift]

    # Protein 'a' calculations
    # vs_a: voxel shift vector to convert from array index to voxel position
    vs_a, a_bar = p.p_bar.out(a_c, c.RHO, constants)

    if len(indices) == 0: # Then we're in the scan stage
        index_list = range(len(c.ANGLES))
        #index_list = range(200)
        vs_a, a_bar = p.p_bar.add_thickness(a_bar, vs_a, c.RHO) # Add thickness to 'a'
    else: # We're in the discrimination stage
        index_list = indices
        vs_a, a_bar = p.p_bar.add_thickness(a_bar, vs_a, c.RHO) # Add thickness to 'a'
        vs_a, a_bar = p.p_bar.add_thickness(a_bar, vs_a, c.RHO) # Add thickness to 'a'

    # Find the fft grid dimension N
    N = fft.n.find_N(a_bar, b_c, constants['eta'])

    # reshape to NxNxN
    a_grid = fft.fft.p_grid(a_bar, N)


    # Rotation and FFT

    fft_out = [] # Will hold all max fft output values and shift vectors

    for i in index_list:
        
        b_r = rotation.utils.q_rotate(b_c,c.Q_LIST[i])
        
        # Protein 'b' calculations
        vs_b, b_bar = p.p_bar.out(b_r, c.DEL_B, constants)
        b_grid = fft.fft.p_grid(b_bar, N) # Reshape to NxNxN


        k_max, v_k = fft.fft.run(a_grid, b_grid)

        vox_shift = (vs_a - vs_b) - v_k # shift to get 'b' to fit at 'a' native position

        fft_out.append([k_max, vox_shift])

    return fft_out


# ##
# Get a_c, b_c
# ##

DIR = os.path.join('..','..','..','PyMol','2hhb')
# Location of the PDB file on the local console. 
# Hard coded, but it's ok because path is also needed in AWS code.

filepath_A = os.path.join(DIR, 'A Chain - 2hhb.cif')
filepath_B = os.path.join(DIR, 'B Chain - 2hhb.cif')
# complete file path

a_c, b_c_pdb = np.array(parseLocal(filepath_A)), np.array(parseLocal(filepath_B))
# Cartesian coordinates of atoms in proteins a, b, extracted from local PDB files

# Shift geometric center of b to the origin
centroid, b_c = rotation.utils.f_centroid(b_c_pdb) # To avoid recalculation, keep b_c_pdb
#tests.utils.plot_multiple([b_c_pdb, b_c]) # f_centroid works as intended


# ##
# Scan stage
# ##

# 1: scan; 2: discrimination; 3: scan_constant; 4. discrim_constant
# Note: there is a way to get comments to show up when I hover
# constants: (eta, atom_voxels)
constants = utils.set_constants(1)

start = time.time()
stage_scan = run_stage(constants, []) # Output from scan stage
end = time.time()
print('Scan stage: ', end-start)

# Find the indices of the top 20 maximum values

# Array of only the k_max values
arr = np.array([i[0] for i in stage_scan])

# Get the indices that would sort the array in descending order
sorted_indices = np.argsort(arr)[::-1]

# Get the indices of the top 20 largest values
top_indices = sorted_indices[:20]
print('\n',f'Top 20: {top_indices}','\n\n')
for i in top_indices:
    print(stage_scan[i][0], stage_scan[i][1]) # Print k_max and voxel shift for each index

# Get the values of the top 20 largest values
#top_values = arr[top_indices]

# ##
# Discrimination stage
# ##

# Proceed with discrimination stage
constants = utils.set_constants(2)

start = time.time()
stage_discr = run_stage(constants, top_indices) # Output from discrimination stage
end = time.time()
print('Discrimination stage: ', end-start)
for i in stage_discr:
    print(i[0], i[1]) # Print k_max and voxel shift for each index


# ##
# Output values
# ##
i = np.argmax([value[0] for value in stage_discr]) # Index of max k_max in stage_discr
angles_index = top_indices[i] # Index of angles for max k_max
v_s = stage_discr[i][1] # voxel shift

e_angles = np.array(c.ANGLES[angles_index])*c.DELTA # Euler rotation angles
v_shift = np.round(
    constants['eta']*v_s - rotation.utils.e_rotate(e_angles, centroid), 
    decimals=3)
# translation vector = eta*v_s - R x centroid

output = (e_angles, v_shift)

# Rotation & translation matrix for PyMol
# Can use this rotation matrix with the PyMol line:
#   cmd.transform_selection('b_c', ['Matrix here'])
phi, theta, psi = e_angles
R = rotation.utils.R(phi, theta, psi) # Euler Matrix for the zx'z' convention
s = v_shift # Translation
matrix = [
    R[0,0],R[0,1],R[0,2],s[0],
    R[1,0],R[1,1],R[1,2],s[1],
    R[2,0],R[2,1],R[2,2],s[2],
    0,0,0,1
] # PyMol rotation/translation matrix
matrix = [float(i) for i in matrix] # Must be simple list for PyMol

# Will need these to print the PyMol script later.

print('Output values:', '\n',
    f'i: {i}', '\n',
    f'e_angles: {e_angles}', '\n',
    f'v_shift: {v_shift}', '\n',
    f'R: {R}', '\n',
    f's: {s}', '\n',
    f'matrix: {matrix}', '\n'
)
"""

# Output for comparison with Lambda - index_list=21
"""
Top 20: [ 0 15 19 13 14  7 17  6  1 18  3 11  2 10 20 12  4  9 16  5] 

Discrimination stage results - ARM
0.3961738051672549 [33  6 70]
0.18247009283637103 [-12  15  11]
0.23700138494838943 [39  1 66]
0.17969251142941725 [18  8 67]
0.21721820309070392 [61 58 46]
0.15758523084346657 [-1 40 63]
0.1690356684802913 [-14  36  35]
0.1917664877494403 [62 -7 56]
0.1536739427397961 [10 50 81]
0.1610430362684486 [32 55 19]
0.14477434517058113 [82 11 44]
0.17130308187372215 [-1 32  6]
0.16784527644874125 [64 57 64]
0.15429748142298927 [ 8 21  2]
0.16824207379259065 [32 54 18]
0.16348050566638592 [62 48 28]
0.16665488441718895 [28 46  4]
0.14040957438822402 [-5 46  6]
0.1353078942530019 [ 9 11 67]
0.1603628122504178 [66 18 76]

Discrimination stage results - Intel
[0, [0.396173805167255, [33, 6, 70]]], 
[1, [0.1536739427397961, [10, 50, 81]]], 
[2, [0.16784527644874098, [64, 57, 64]]], 
[3, [0.14477434517058116, [82, 11, 44]]], 
[4, [0.16665488441718998, [28, 46, 4]]], 
[5, [0.16036281225041765, [66, 18, 76]]], 
[6, [0.1917664877494391, [62, -7, 56]]], 
[7, [0.15758523084346657, [-1, 40, 63]]], 
[9, [0.1404095743882243, [-5, 46, 6]]],
[10, [0.1542974814229898, [8, 21, 2]]], 
[11, [0.17130308187372245, [-1, 32, 6]]], 
[12, [0.16348050566638309, [62, 48, 28]]], 
[13, [0.17969251142941814, [18, 8, 67]]], 
[14, [0.21721820309070514, [61, 58, 46]]], 
[15, [0.18247009283637086, [-12, 15, 11]]], 
[16, [0.1353078942530028, [9, 11, 67]]], 
[17, [0.16903566848029453, [-14, 36, 35]]], 
[18, [0.1610430362684486, [32, 55, 19]]], 
[19, [0.2370013849483897, [39, 1, 66]]], 
[20, [0.16824207379259107, [32, 54, 18]]]
# This particular output was pulled from Lambda. 
# Matches with Acer in number; 
# form [q, [k_max, vox_shift]] comes from Lambda.

Final output:
 0 
 [0. 0. 0.] 
 [-0.324  0.049 -0.731] 
 [[ 1. -0.  0.]
 [ 0.  1. -0.]
 [ 0.  0.  1.]] 
 [-0.324  0.049 -0.731] 
 [1.0, -0.0, 0.0, -0.324, 0.0, 1.0, -0.0, 0.049, 0.0, 0.0, 1.0, -0.731, 0.0, 0.0, 0.0, 1.0] 
"""

# Output for comparison with Lambda - index_list=200
"""
Scan stage:  15.668599128723145

 Top 20: [  0  15 172 104 100  19  32  35  26 191 153  70  43 112  80 154  41 159
 184  13] 

0.23240522096092847 [23  5 49]
0.19092891827989572 [ 6 11 44]
0.1782200274346405 [-5  8 32]
0.17812150890095704 [ 2 21 47]
0.1778259532999043 [13 48 46]
0.16807261846517074 [28  1 46]
0.1662992848588577 [48 27 20]
0.1652155809883323 [29 38 17]
0.16403335858412313 [-5  6 29]
0.161668913775701 [29 49 44]
0.16127483964096428 [ 5 30 53]
0.16078224697254398 [11 37 52]
0.1593044689672794 [-2 26 41]
0.15605335735570372 [12 22 58]
0.15516669055254623 [17  6 47]
0.1540829866820208 [48 31 20]
0.15309780134518036 [ 0 33 41]
0.14974817119991926 [28 34 10]
0.1497481711999187 [ 37 -10  32]
0.1493540970651812 [17  2 46]

Discrimination stage:  6.176717042922974
0.3961738051672549 [33  6 70]
0.18247009283637103 [-12  15  11]
0.21744494443004686 [-8 11 46]
0.14947922796194646 [ 2 30 69]
0.16024944158074406 [36 10 75]
0.23700138494838943 [39  1 66]
0.15254023604307962 [63 -8 56]
0.15701837749510822 [41 56 23]
0.20066608531865476 [-8  9 42]
0.17073622852536696 [-4  4 38]
0.18127970080481837 [53 63 40]
0.14001277704437368 [52 -1 76]
0.17011268984217184 [ 31 -22  58]
0.1733437539278106 [48  6 74]
0.1742507192851828 [25  8 68]
0.1679019617835747 [69 45 28]
0.1853043595781612 [ 56 -10  20]
0.1789556020765508 [39 49 14]
0.18201661015768586 [57 -9 19]
0.17969251142941725 [18  8 67]

Output values: 
 i: 0 
 e_angles: [0. 0. 0.] 
 v_shift: [-0.324  0.049 -0.731] 
 R: [[ 1. -0.  0.]
 [ 0.  1. -0.]
 [ 0.  0.  1.]] 
 s: [-0.324  0.049 -0.731] 
 matrix: [1.0, -0.0, 0.0, -0.324, 0.0, 1.0, -0.0, 0.049, 0.0, 0.0, 1.0, -0.731, 0.0, 0.0, 0.0, 1.0] 
"""

# Output for comparison with Lambda - index_list=2628
"""
Scan stage:  205.7004919052124

 Top 20: [ 0 1415 1517 1950 1404  661   15 2595 2010 1694 2002 1637  259  499
 1301 1913  711 1112  172  104] 


0.23240522096092847 [23  5 49]
0.21122373621883464 [ 0  7 17]
0.21063262501672722 [ 5 29 41]
0.2084652172756749 [17 11 49]
0.1969385488346306 [38 21  6]
0.19102743681358061 [47 11 14]
0.19092891827989572 [ 6 11 44]
0.18787484373568838 [13 45 35]
0.18708669546621406 [ 6 20 51]
0.18698817693252917 [36 -8 36]
0.18403262092200637 [17 30  2]
0.183342991186217 [48  6 33]
0.18314595411884677 [34 -1 41]
0.18156965757990218 [38 -4 19]
0.18038743517569059 [22 -3 44]
0.17969780543990366 [27 18  4]
0.17900817570411345 [42  1 21]
0.17881113863674575 [28 32  9]
0.1782200274346405 [-5  8 32]
0.17812150890095704 [ 2 21 47]

Discrimination stage:  6.102051258087158
0.3961738051672549 [33  6 70]
0.1790689727462204 [-1  9 24]
0.1701126898421708 [30 -7 69]
0.22305679257878933 [24 15 71]
0.16722173776554358 [35 50 13]
0.22895206740171087 [68 15 20]
0.18247009283637103 [-12  15  11]
0.15888899354468522 [ 2 40 64]
0.1744207752896907 [ 8 28 73]
0.1723801032356006 [ 1  4 16]
0.18371717020275954 [15 57 72]
0.2181251684480723 [70  8 47]
0.2098491095620523 [50 -2 59]
0.15322046006111015 [46 16 88]
0.17198330589174524 [32 -5 64]
0.17872886073721025 [40 25  4]
0.20684478681575394 [ 46 -12  51]
0.20321692538626715 [41 44 12]
0.21744494443004686 [-8 11 46]
0.14947922796194646 [ 2 30 69]

Output values: 
 i: 0 
 e_angles: [0. 0. 0.] 
 v_shift: [-0.324  0.049 -0.731] 
 R: [[ 1. -0.  0.]
 [ 0.  1. -0.]
 [ 0.  0.  1.]] 
 s: [-0.324  0.049 -0.731] 
 matrix: [1.0, -0.0, 0.0, -0.324, 0.0, 1.0, -0.0, 0.049, 0.0, 0.0, 1.0, -0.731, 0.0, 0.0, 0.0, 1.0] 
"""

# ###
# Scratch
# ###
from p.pdb_parser import parseLocal, preprocessing
import os
import numpy as np

from mmcif.io.IoAdapterCore import IoAdapterCore
from mmcif.api.PdbxContainers import DataContainer
from mmcif.api.DataCategory import DataCategory

from pathlib import Path # TESTING

def write_docked_file(
    filepath_A,
    filepath_B,
    transformed_coordinates_B,
    output_file,
):
    """
    Create a new mmCIF containing proteins A and B.

    Protein A:
        - coordinates unchanged
        - chain ID = A
        - entity ID = 1

    Protein B:
        - coordinates replaced by transformed_coordinates_B
        - chain ID = B
        - entity ID = 2

    Atom IDs are renumbered sequentially in the output.

    Parameters
    ----------
    filepath_A : str
        Filepath to protein A.

    filepath_B : str
        Filepath to protein B.

    transformed_coordinates_B : numpy.ndarray
        Array of shape (N, 3) containing the transformed
        coordinates for protein B.

    output_file : str
        Path of the output mmCIF file.
    """

    # ------------------------------------------------------------------
    # Get data containers (NOT IN CHAT SUGGESTION)
    # ------------------------------------------------------------------
    
    # Preprocessing: remove HOH, alternate sites, etc.
    # Returns preprocessed mmCIF container for proteins
    atom_site_A = preprocessing(filepath_A)
    atom_site_B = preprocessing(filepath_B)

    # ------------------------------------------------------------------
    # Verify coordinate count
    # ------------------------------------------------------------------

    transformed_coordinates_B = np.asarray(
        transformed_coordinates_B
    )

    if transformed_coordinates_B.shape != (
        atom_site_B.getRowCount(),
        3,
    ):
        raise ValueError(
            "transformed_coordinates_B must have shape "
            f"({atom_site_B.getRowCount()}, 3)"
        )

    # ------------------------------------------------------------------
    # Create output DataContainer
    # ------------------------------------------------------------------

    output_container = DataContainer("rdock_docked")

    # ------------------------------------------------------------------
    # Create _struct category
    # ------------------------------------------------------------------

    struct = DataCategory(
        "struct",
        [
            "entry_id",
            "title",
        ],
    )

    struct.append(
        [
            "rdock_docked",
            "Protein-protein rigid-body docking model generated by RDock",
        ]
    )

    output_container.append(struct)

    # ------------------------------------------------------------------
    # Create _entity category
    # ------------------------------------------------------------------

    entity = DataCategory(
        "entity",
        [
            "id",
            "type",
            "pdbx_description",
        ],
    )

    entity.append(
        [
            "1",
            "polymer",
            "Protein A",
        ]
    )

    entity.append(
        [
            "2",
            "polymer",
            "Protein B",
        ]
    )

    output_container.append(entity)

    # ------------------------------------------------------------------
    # Create _struct_asym category
    # ------------------------------------------------------------------

    struct_asym = DataCategory(
        "struct_asym",
        [
            "id",
            "entity_id",
        ],
    )

    struct_asym.append(
        [
            "A",
            "1",
        ]
    )

    struct_asym.append(
        [
            "B",
            "2",
        ]
    )

    output_container.append(struct_asym)

    # ------------------------------------------------------------------
    # Create _atom_site category
    #
    # Use the same atom_site attributes as protein A.
    # This preserves the atom-level information already present in A.
    # ------------------------------------------------------------------

    attributes = atom_site_A.getAttributeList()

    output_atom_site = DataCategory(
        "atom_site",
        attributes,
    )

    # ------------------------------------------------------------------
    # Attribute indices
    # ------------------------------------------------------------------

    atom_site_A_indices = {
        attribute: atom_site_A.getAttributeIndex(attribute)
        for attribute in attributes
    }

    atom_site_B_indices = {
        attribute: atom_site_B.getAttributeIndex(attribute)
        for attribute in attributes
        if atom_site_B.hasAttribute(attribute)
    }

    # ------------------------------------------------------------------
    # Add Protein A
    # ------------------------------------------------------------------

    output_atom_id = 1

    for row in atom_site_A.data:

        output_row = list(row)

        # New atom ID
        if "id" in atom_site_A_indices:
            output_row[atom_site_A_indices["id"]] = str(output_atom_id)

        # New chain ID
        if "label_asym_id" in atom_site_A_indices:
            output_row[
                atom_site_A_indices["label_asym_id"]
            ] = "A"

        # New entity ID
        if "label_entity_id" in atom_site_A_indices:
            output_row[
                atom_site_A_indices["label_entity_id"]
            ] = "1"

        output_atom_site.append(output_row)

        output_atom_id += 1

    # ------------------------------------------------------------------
    # Add Protein B
    # ------------------------------------------------------------------

    x_idx = atom_site_B.getAttributeIndex("Cartn_x")
    y_idx = atom_site_B.getAttributeIndex("Cartn_y")
    z_idx = atom_site_B.getAttributeIndex("Cartn_z")

    for row, (x, y, z) in zip(
        atom_site_B.data,
        transformed_coordinates_B,
    ):

        output_row = list(row)

        # New atom ID
        if "id" in atom_site_B_indices:
            output_row[
                atom_site_B_indices["id"]
            ] = str(output_atom_id)

        # New chain ID
        if "label_asym_id" in atom_site_B_indices:
            output_row[
                atom_site_B_indices["label_asym_id"]
            ] = "B"

        # New entity ID
        if "label_entity_id" in atom_site_B_indices:
            output_row[
                atom_site_B_indices["label_entity_id"]
            ] = "2"

        # Transformed coordinates
        output_row[x_idx] = f"{x:.3f}"
        output_row[y_idx] = f"{y:.3f}"
        output_row[z_idx] = f"{z:.3f}"

        output_atom_site.append(output_row)

        output_atom_id += 1

    # ------------------------------------------------------------------
    # Add atom_site to output container
    # ------------------------------------------------------------------

    output_container.append(output_atom_site)

    # ------------------------------------------------------------------
    # Write output mmCIF
    # ------------------------------------------------------------------

    io = IoAdapterCore()
    io.writeFile(
        output_file,
        [output_container],
    )


# Need to prepare:
#    filepath_A,
#    filepath_B,
#    transformed_coordinates_B,
#    output_file,

DIR = os.path.join('..','..','..','..','OneDrive','Non Sys','Projects','PyMol','2hhb')
# Location of the PDB file on the local console. 
# Hard coded, but it's ok because path is also needed in AWS code.

filepath_A = os.path.join(DIR, 'A Chain - 2hhb.cif')
filepath_B = os.path.join(DIR, 'B Chain - 2hhb.cif')
# complete file path

# transformed_coordinates_B
# Rotation comes later. Right now, I just want to test that this works
b_c_pdb = np.array(parseLocal(filepath_B))

# output_file
output_file = os.path.join('..','..','..','..','Downloads','merged_output.cif')

write_docked_file(filepath_A, filepath_B, b_c_pdb, output_file)
