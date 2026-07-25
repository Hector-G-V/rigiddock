import constants as c
from math import pi as pi
import numpy as np


"""def to_voxels(p_c, eta):

    # This function is ok, but adds inconsistency for points x/eta that land on a multiple of 0.5:
    # If the protein is not entirely in the positive or negative octant, points that land on a
    # multiple of 0.5 will stretch the protein inconsistently with a linear shift anywhere else.
    # The number of affected points is 0-2, so could leave alone, but the complexity below is big
    # compared to round(). Furthermore: the time save is negligible, but np.round is 3-4 times 
    # faster than the original function. To do all rotations, round takes 0.06s, and the original takes 0.24s

    # Performs the following algorithm on every single scalar in p_c:
    #if x > 0:
    #    vox = math.ceil(x/eta - 1/2)
    #else:
    #    vox = math.floor(x/eta + 1/2)
    
    # Convert Cartesian x values to 1D voxel positions for each element in the array
    pos_indices = p_c > 0
    neg_indices = p_c <= 0

    # Apply the ceiling operation for positive values
    pos_voxels = np.ceil(p_c[pos_indices] / eta - 1/2)
    
    # Apply the floor operation for negative values
    neg_voxels = np.floor(p_c[neg_indices] / eta + 1/2)
    
    # Initialize the output array
    voxels = np.empty_like(p_c)
    
    # Fill the output array with the computed values
    voxels[pos_indices] = pos_voxels
    voxels[neg_indices] = neg_voxels
    
    return voxels.astype(int)"""


def to_voxels(p_c, eta):
    # For every point in p_c, this function chooses the voxel whose xyz 
    # coordinates are closest to the point (pc_x, pc_y,pc_z)/eta.
    # Return: p, array of voxel coordinates.
    
    return np.round(p_c/eta).astype(int)


def atom_vox_variable(p_c, constants):

    # Finds p, all the voxels of the protein p_c.
    # Voxels in atom vary with position.
    # p_c: PDB coordinates of the protein. Numpy array.
    # constants: eta, v_candidates, v_points

    # Initialize the variables
    d = constants # Shorten dict variable name for clarity in math below
    eta = d['eta']
    
    # PDB points converted to voxel position
    ctr_vox = to_voxels(p_c, eta) # These are the center voxels of their atom voxels

    # Candidate voxel points shifted to each atom's voxel center position
    v_p = ctr_vox[:, np.newaxis, :] + d['v_candidates'] # v_p: protein voxel points

    # Distance b/n protien coord and center of voxel candidate: p_c - v_c
    # v_c = v_p*eta is now a Cartesian coordinate; no longer a voxel position.
    # Note: 'p_c - v_c' is 'coord-adjusted' in previous function iterations.
    pc_vc = (p_c[:,np.newaxis,:] - eta*v_p).reshape(-1,3)

    # Test which voxels are inside the protein
    inside_test = np.sum(pc_vc ** 2, axis=1)<=c.R*c.R # (p_c - v_c)^2 <= r^2

    # Select candidates that are inside the protein
    v_inside = v_p.reshape(-1,3)[inside_test]

    # return
    #p = np.unique(v_inside, axis=0) # np.unique orders the output array
    p = v_inside # Finding unique takes more time than without!

    return p


def atom_vox_constant(p_c, constants):

    # Finds p, all the voxels of the protein p_c. 
    # Voxels in atom do not depend on atom position.
    # p_c: PDB coordinates of the protein. Numpy array.
    # constants: eta, atom_voxels

    # Initialize the variables
    d = constants # Shorten dict variable name for clarity in math below
    eta = d['eta']

    # Convert p_c to voxel position
    ctr_vox = to_voxels(p_c, eta) # These are the center voxels of their atom voxels    

    # Add atom_voxels to each point in out_p
    p = (ctr_vox[:, np.newaxis, :] + d['atom_voxels']).reshape(-1, 3)
    # Note: deleting duplicates via np.unique costs time; it doesn't save time.
    
    return p


def find_p(p_c, constants):

    # Finds p, all the voxels of the protein p_c
    # p_c: PDB coordinates of the protein. Numpy array
    # constants: eta, atom_voxels

    # Initialize the variables
    if 'v_candidates' in constants: # If constants==(eta, v_candidates)
        p = atom_vox_variable(p_c, constants)
    elif 'atom_voxels' in constants:  # If constants==(eta, const_atom_voxels)
        p = atom_vox_constant(p_c, constants)

    return p


"""def find_p_bar(p, inside_value):

    # This function uses the zero values to populate p_bar. But after studying the breakdown
    # of real proteins, I found that most of the grid is empty. Therefore, I changed the 
    # funciton to use the non-zero values.

    # Notes to transfer: The algorithm that uses the multi-dim array to find p_bar; 
    # if I use an (n,3) array, I will run into the same problem that slowed down 
    # computation: I will search a list ~N^2 times!

    # Finds p_bar
    # p: protein voxel positions. Numpy array of shape (n,3)
    # inside_value: c.RHO for protein 'a', c.DEL_B for 'b.' Float.
    # return p_bar, Numpy array of shape (l,m,n)

    # Find the max and min values of p along xyz
    p_max, p_min = np.max(p, axis=0), np.min(p, axis=0) 
    # p_min is the shift vector v_s, needed to convert from p_bar array index
    # to voxel position

    # Shape of the array that will hold p_bar
    shape = p_max - p_min + 1 # +1 adds the origin voxel, which is excluded

    p_bar = np.zeros(shape) # Initialize p_bar array

    # Convert p voxel positions to a tuple of arrays for advanced indexing
    indices_tuple = tuple((p-p_min).T) # p-p_min shifts voxels to the 1st Octant

    # Populate p_bar
    p_bar[indices_tuple] = inside_value # Populate p_bar with the inside value (rho or delta)


    # To find the outside voxels, must the pad p_bar array with zeros
    # Algorithm: if p_bar value is 0, then all non-diagonal neighbors are surface

    # Define the padding width for each axis
    pad_width = ((1, 1), (1, 1), (1, 1))  # Padding 1 element on both sides for each axis

    # Pad the array with 0
    p_padded = np.pad(p_bar, pad_width, mode='constant', constant_values=0)

    out_p = np.array(np.where(p_padded == 0)).T # Indices for outside p (the zeros)

    n_offsets = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
    # This is 'non-diagonal neighbors' in my language

    # Add n_offsets to each point in out_p
    n = out_p[:, np.newaxis, :] + n_offsets

    # Flatten the indices for advanced indexing
    n_flat = n.reshape(-1, 3)

    # Unpack the indices for advanced indexing
    indices = (n_flat[:, 0], n_flat[:, 1], n_flat[:, 2])

    # Get the dimensions of p_padded
    max_dims = np.array(p_padded.shape) - 1

    # Create a boolean mask to filter out out-of-bound indices
    valid_mask = (
        (indices[0] >= 0) & (indices[0] <= max_dims[0]) &
        (indices[1] >= 0) & (indices[1] <= max_dims[1]) &
        (indices[2] >= 0) & (indices[2] <= max_dims[2])
    )

    # Filter indices using the valid_mask
    filtered_indices = (indices[0][valid_mask], indices[1][valid_mask], indices[2][valid_mask])

    # Check where p_padded is not 0 and set those positions to 1
    p_padded[filtered_indices] = np.where(p_padded[filtered_indices] != 0, 1, p_padded[filtered_indices])

    # return
    p_bar = p_padded[1:-1, 1:-1, 1:-1] # Remove padding

    # p_min is the shift vector to needed convert from array index to voxel position
    return p_min, p_bar 
"""

def find_p_bar(p, inside_value):

    # Finds p_bar
    # p: protein voxel positions. Numpy array of shape (n,3)
    # inside_value: c.RHO for protein 'a', c.DEL_B for 'b.' Float.
    # return p_bar, Numpy array of shape (l,m,n)

    # Find the max and min values of p along xyz
    p_max, p_min = np.max(p, axis=0), np.min(p, axis=0) 
    # p_min is the shift vector v_s, needed to convert from p_bar array index
    # to voxel position

    # Shape of the array that will hold p_bar
    shape = p_max - p_min + 1 # +1 adds the origin voxel, which is excluded

    p_bar = np.zeros(shape) # Initialize p_bar array

    # Convert p voxel positions to a tuple of arrays for advanced indexing
    indices_tuple = tuple((p-p_min).T) # p-p_min shifts voxels to the 1st Octant

    # Populate p_bar
    p_bar[indices_tuple] = 1 # Populate p_bar with 1


    # FROM CHATGPT
    
    # Create shifted versions of p_bar
    shift_xp = np.pad(p_bar[1:, :, :], ((0, 1), (0, 0), (0, 0)), mode='constant')
    shift_xm = np.pad(p_bar[:-1, :, :], ((1, 0), (0, 0), (0, 0)), mode='constant')
    
    shift_yp = np.pad(p_bar[:, 1:, :], ((0, 0), (0, 1), (0, 0)), mode='constant')
    shift_ym = np.pad(p_bar[:, :-1, :], ((0, 0), (1, 0), (0, 0)), mode='constant')
    
    shift_zp = np.pad(p_bar[:, :, 1:], ((0, 0), (0, 0), (0, 1)), mode='constant')
    shift_zm = np.pad(p_bar[:, :, :-1], ((0, 0), (0, 0), (1, 0)), mode='constant')

    # Check the condition for updating p_bar
    condition = (
        (p_bar == 1) &
        (shift_xp == 1) & (shift_xm == 1) &
        (shift_yp == 1) & (shift_ym == 1) &
        (shift_zp == 1) & (shift_zm == 1)
    )

    # Update p_bar where the condition is met
    p_bar[condition] = inside_value

    # p_min is the shift vector to needed convert from array index to voxel position
    return p_min, p_bar 


def add_thickness(p_bar, v_s, inside_value):

    # Adds thickness to the protein. 
    # p_bar: protein to be given thickness. Numpy array of shape (l,m,n)
    # v_s: shift vector to be updated. Numpy array.
    # inside_value: c.RHO for protein 'a', or c.DEL_B for protein 'b'
    # return: p_bar with added thickness and its updated shift vector v_s'

    # If the protein is 'b', then do not add thickness
    if inside_value == c.DEL_B:
        return v_s, p_bar

    # If not b, then must be a. Add thickness.
    # Algorithm:
    #   For every value 1: all zero-value, non-diagonal neighbors are set to 1.

    # To add thickness, must the pad p_bar array with zeros
    # Define the padding width for each axis
    pad_width = ((1, 1), (1, 1), (1, 1))  # Padding 1 element on both sides for each axis

    # Pad the array with 0
    p_padded = np.pad(p_bar, pad_width, mode='constant', constant_values=0)

    val_one = np.array(np.where(p_padded == 1)).T # Indices for voxel value=1

    n_offsets = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
    # This is 'non-diagonal neighbors' in my language

    # Add n_offsets to each point in val_one
    n = val_one[:, np.newaxis, :] + n_offsets

    # Flatten the indices for advanced indexing
    n_flat = n.reshape(-1, 3)

    # Unpack the indices for advanced indexing
    indices = (n_flat[:, 0], n_flat[:, 1], n_flat[:, 2])

    # Get the dimensions of p_padded
    max_dims = np.array(p_padded.shape) - 1

    # Create a boolean mask to filter out out-of-bound indices
    valid_mask = (
        (indices[0] >= 0) & (indices[0] <= max_dims[0]) &
        (indices[1] >= 0) & (indices[1] <= max_dims[1]) &
        (indices[2] >= 0) & (indices[2] <= max_dims[2])
    )

    # Filter indices using the valid_mask
    filtered_indices = (indices[0][valid_mask], indices[1][valid_mask], indices[2][valid_mask])

    # Check where p_padded is 0 and set those positions to 1
    p_padded[filtered_indices] = np.where(p_padded[filtered_indices] == 0, 1, p_padded[filtered_indices])

    # return
    # Unlike the p_bar calculation, the padding now stays.
    # v_s is now shifted by [-1,-1,-1] 
    return (v_s - [1,1,1], p_padded)


def out(p_c, inside_value, constants):
    
    p = find_p(p_c, constants)
    
    v_s, p_bar = find_p_bar(p, inside_value) 
    # v_s: shift vector to convert from voxel position to array index

    return v_s, p_bar
