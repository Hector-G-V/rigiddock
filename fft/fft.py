import numpy as np

def p_grid(p_bar, N):

    # Puts p_bar on the grid by padding with zeros to shape NxNxN
    # p_bar: a_bar or b_bar. Numpy array of shape (l,m,n)
    # N: Number of voxels along 1D of grid.
    # return: array of shape NxNxN

    pad = N-np.array(p_bar.shape)

    # Define the padding width for each axis
    pad_width = ((0, pad[0]), (0, pad[1]), (0, pad[2]))  # Padding for shape NxNxN

    # Pad the array with 0 for shape NxNxN
    p_padded = np.pad(p_bar, pad_width, mode='constant', constant_values=0)

    return p_padded


def kk_algorithm(a_bar, b_bar):

    # Input: a, b values, reshaped for np.fft input. Numpy arrays.
    # returns: [max value, max value index]

    # Note: pyfftw can potientially improve the computation time by 2-5 fold.
    # The algorithm below is fine, given the scan step. But pyfftw will be the 
    # next step toward optimization. See ChatGPT suggestion.
    #   Check: does this install on aws???
    #   Windows: Try: build a new python 3.11 environment, to install distutils

    # "ortho" does not make a difference in computational time
    A = np.conjugate(
        np.fft.ifftn(a_bar, norm="ortho")
    )
    B = np.fft.ifftn(b_bar, norm="ortho")

    C = A*B

    c = np.fft.fftn(C, norm="ortho").real # Im(c) = 0, always
    # c = np.round(c,10) # Obsolete: Essential for floating-point error in complex numbers

    # The line below seems to pull only the first index of the max value. 
    # I believe it's ok, given that the proteins will be large,
    # and there will likely be only one maximum. 
    max_index = np.unravel_index(np.argmax(c), c.shape)
    
    return [c[max_index], np.array(max_index)]
    
    
def eval_shift(v_shift, N):
    
    # This function evaluates the shift that is output by kk_algorithm, which 
    # only takes steps in the (-) direction. The fft algorithm loops back onto
    # the grid, but vectors do not.

    # Algorithm: for each of the xyz coordinates, choose the smallest step
    # direction between the (-) step and the (+) step.

    # v_shift: voxel shift vector given by kk_algorithm. Nupmy array of ints
    # N: number of voxels along one direction of the grid cube
    # return: optimal voxel shift vector

    u = np.array([], dtype = 'int64') # Will hold the new vector

    for i in range(3):
        v = v_shift[i]

        if abs(v) <= abs(v-N):
            u = np.append(u, int(v))
        else:
            u = np.append(u, int(v-N))

    return u


def run(a_grid, b_grid):

    # kk_max: the max value returned by kk_algorithm
    # v_k: the shift vector returned by kk_algorithm
    kk_max, v_k = kk_algorithm(a_grid, b_grid)

    # Evaluate v_k for the smallest step in xyz directions
    v_k = eval_shift(v_k, len(a_grid))

    # (a_c_vox-b_c_vox)-v_k: voxel shift vector for the original position of b
    return kk_max, v_k
