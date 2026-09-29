import constants as c
import utils, p.p_bar, fft.fft, fft.n, rotation.utils
from p.pdb_parser import parseLocal, write_docked_file
import os
from math import pi as pi
import numpy as np

import time, tests.utils # FOR NOW. Remove when finished developing locally.
import tests.proteins.p_3apr as proteins # Parsed test proteins


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
        #index_list = range(100)
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

#a_c, b_c_pdb = np.array(proteins.PARSED_A), np.array(proteins.PARSED_B)
# Test: Cartesian coordinates of atoms in proteins a, b, parsed and ready to go

#a_c, b_c_pdb = tests.utils.simple_input('4') 
# Test: Simple input for testing


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


# For testing: export fft_out to txt
tests.utils.to_txt([
    float(i[0]), [int(j) for j in i[1]]
] for i in stage_scan)


# Find the indices of the top 20 maximum values

# Array of only the k_max values
#stage_scan = tests.utils.to_python('scan_out_t.txt') # For testing
arr = np.array([i[0] for i in stage_scan])

# Get the indices that would sort the array in descending order
sorted_indices = np.argsort(arr)[::-1]

# Get the indices of the top 20 largest values
top_indices = sorted_indices[:20]

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

# Plots:
#tests.utils.plot_1D(np.real([i[0] for i in stage_scan])) # All fft_out
tests.utils.plot_2D(top_indices, arr[top_indices].real) # Top only
tests.utils.plot_2D(top_indices, np.real([i[0] for i in stage_discr]))

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

# Will need 'matrix', s to print the PyMol script later.

# write the docked-proteins file
# Need to prepare:
#    filepath_A,
#    filepath_B,
#    transformed_coordinates_B,
#    output_file,

# filepath_A, filepath_B: initialized at the beginning of run.py

# transformed_coordinates_B
#transformed_B = b_c_pdb @ R.T + s # Gives overflow error when N is large. Might be version, OS or memory-dependent.
transformed_B = np.dot(b_c_pdb, R.T) + s # Equivalent to above; produces no errors.

# output_file
#output_file = os.path.join('..','..','..','..','Downloads','merged_output.cif') # Studio
output_file = os.path.join('..','..','..','Downloads','merged_output.cif') # Air

write_docked_file(filepath_A, filepath_B, b_c_pdb, output_file)
