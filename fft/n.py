from constants import R
from math import ceil, floor
import numpy as np


def find_L_a(a_bar):

    # Finds L_a, the maximum length for all 3 directions.
    # Input: a_bar. Numpy array
    # return: L_a

    # Note: it is simple now, but this may change in the future

    return np.max(a_bar.shape)


def find_L_b(r, eta):
    
    # Finds L, the maximum length for all 3 directions.
    # r: PDB coordinates for protein b. Numpy array.
    # return: L_b

    # If protein is just one atom, L_b = 2ceil(r/eta) + 1. Note the +1 accounts 
    # for the PDB point voxel. There is no +1 extra measure here; r can never 
    # exceed ceil(r/eta) 
    if len(r) == 1:
        return int(2*np.ceil(R/eta) + 1)


    # Find the max distance between any two points in the protein

    # Calculate the difference between each pair of vectors
    differences = r[:, np.newaxis, :] - r[np.newaxis, :, :]

    # Calculate the dot product of each difference vector with itself
    dot_products = np.einsum('ijk,ijk->ij', differences, differences)

    # Find the index of the maximum value in the dot_products array
    max_index = np.unravel_index(np.argmax(dot_products), dot_products.shape)

    # Extract the points corresponding to these indices
    r_i = r[max_index[0]]
    r_j = r[max_index[1]]

    # Max r_ij (the scalar)
    # KEEP BECAUSE COULD STOP HERE AND KEEP THE FORMULA VERY SIMPLE:
    #   L_b = (ceil(r_ij/eta) + 1 + 2ceil(R/eta)).astype(int)
    #   But there will potentially be extra voxels; maybe three. Will check.
    #   So maybe put these in notes and remove from here?
    #r_ij = np.sqrt(dot_products[max_index])

    
    # Rotate, translate r_i, r_j to make r_ij parallel to an axis

    # Step 1: Normalize r_ij
    r_ij = r_j - r_i
    #r_ij = r_ij/np.linalg.norm(r_ij)
    if np.all(r_ij==0): # if r_ij = [0 0 0]
        r_ij = r_ij
    else:
        r_ij = r_ij/np.linalg.norm(r_ij) # Normalize

    # Step 3: Determine the rotation axis (cross product)
    # Chose to make r_ij parallel x, but axis choice doesn't matter.
    target_axis = np.array([1,0,0]) 
    axis = np.cross(r_ij, target_axis)
    if np.all(axis==0): # if axis = [0 0 0]
        axis = axis
    else:
        axis = axis/np.linalg.norm(axis) # Normalize

    # Step 4: Calculate the rotation angle
    angle = np.arccos(np.dot(r_ij, target_axis))

    # Step 5: Construct the rotation matrix
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])
    I = np.eye(3) # 3x3 Identity matrix
    M = I + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K,K)

    # Step 6: Apply the rotation to both vectors
    r_i_prime = np.dot(M, r_i)
    r_j_prime = np.dot(M, r_j)

    # Extract the scalars that matter
    ri_p, rj_p = r_i_prime[0], r_j_prime[0] # [0] bc we align with the x-axis


    # Find the number of voxels in v_r for each point.

    v_i, v_j = np.round(ri_p/eta), np.round(rj_p/eta) # voxel positions for ri, r_j

    v_ij = np.abs(v_j - v_i) + 1 # Number of voxels occupied by r_ij

    vr_i = np.floor((np.abs(ri_p) + R)/eta) - np.abs(v_i) 
    # Number of voxels within v_r for r_i

    vr_j = np.floor((np.abs(rj_p) + R)/eta) - np.abs(v_j) 
    # Number of voxels within v_r for r_j


    # return L_b
    return int(v_ij + vr_i + vr_j)


def find_N(a_bar, b_c, eta):
    
    # returns N, number of voxels in 1D of grid
    # a_bar: Numpy array of shape (l,m,n)
    # b_c: Numpy array of shape (n,3)
    # Note: this is simple for now, but may change later

    # Find N
    L_a, L_b = find_L_a(a_bar), find_L_b(b_c, eta)
    N = L_a + L_b # Buffers are included in the L_b calculation

    return N
