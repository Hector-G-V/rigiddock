# If any functions here get too big or complicated,
# then make the function a separate module.

import constants as c
import math
from math import pi as pi
from numpy import sin as Sin
from numpy import cos as Cos


# Set stage constants
def set_constants(stage):
    # Sets stage constants
    # stage:
    #   1: scan; 2: discrimination; 3: scan_constant; 4. discrim_constant
    # return constants: (eta, atom_voxels)

    if stage == 1: # scan stage
        return {
            'eta': c.ETA_S, 
            'v_candidates': c.V_CANDIDATES_S}
    elif stage == 2: # discrimination stage
        return {
            'eta': c.ETA_D, 
            'v_candidates': c.V_CANDIDATES_D}
    elif stage == 3: # scan stage, constant atom voxels
        return {
            'eta': c.ETA_S, 'atom_voxels': c.CONST_ATOM_VOXELS_S}
    elif stage == 4: # discrimination stage, constant atom voxels
        return {'eta': c.ETA_D, 'atom_voxels': c.CONST_ATOM_VOXELS_D}


# Voxel Functions
def atom_voxels(r, eta):

    # Produces the voxel candidates of an atom given r and eta. 
    # Produces a txt file, and returns arrays for testing or manipulation.
    # The output is to be pasted into constants.py

    v_candidates = [] # voxel candidates

    # Find v_r, the 'voxel radius' along an axis. 
    # Will be used to produce the cube of voxel candidates.
    v_r = math.ceil(r/eta)

    # Populate v_candidates
    for i in range(-v_r, v_r+1):
        for j in range(-v_r, v_r+1):
            for k in range(-v_r, v_r+1):
                v_candidates.append([i,j,k])

    # Write the output to txt files
   
    candidates_file = "v_candidates.txt"
    with open(candidates_file, 'w') as f:
        line = '['
        for sublist in v_candidates:
            sublist_str = str(sublist)
            if len(line) + len(sublist_str) <= 80:
                if line[-1] == '[':
                    line += sublist_str
                else:
                    line += ', ' + sublist_str
            else:
                f.write(line + '\n')
                line = '    ' + sublist_str
        if line != '[':
            f.write(line + '\n]')
        else:
            f.write('[]\n')

    # Optional: return for testing or manipulation
    return (v_candidates)
    

# Rotation Functions
def q(phi, theta, psi):
    # Part of the function series to find all the quaternions for step size Delta

    # Euler to quaternion conversion.
    # Note : this function is for the zx'z' convention.
    # return: quaternion q derived from Euler angles

    return [
        Cos((phi + psi)/2) * Cos(theta/2),
        Cos((phi - psi)/2) * Sin(theta/2),
        Sin((phi - psi)/2) * Sin(theta/2),
        Sin((phi + psi)/2) * Cos(theta/2)
    ]


def euler_to_q(delta):

    # Part of the function series to find all the quaternions for step size Delta    

    # Finds all quaternions from every Euler angle combination.
    # Finds unique quaternions only — those that produce a unique rotation.
    # Note: this function is for the zx'z' convention.
    # return: list of all quaternions 

    tolerance = 1e-10 # Tolerance value for floating points in the while loop

    qlist = [] # List of quaternions from all Euler angle combinations
    
    angles = [] # List of Euler angles that pertain to the quaternions. 
    # Keeping only i,j,k to save storage space.


    # First find n_2p, the number of deltas that fit inside 2pi
    n_0, r = divmod(2*pi, c.DELTA) # Quotient, remainder of 2pi/delta

    if r == 0:
       n_2p = int(n_0 - 1)
    else:
        n_2p = int(n_0)

    # Repeat for n_p, the number of deltas that fit inside pi
    n_0, r = divmod(pi, c.DELTA) # Quotient, remainder of pi/delta
    
    if r == 0:
       n_p = int(n_0 - 1)
    else:
        n_p = int(n_0)

    # theta = 0, psi = 0    
    for i in range(n_2p):
        phi = i*delta
        qlist.append(q(phi,0,0))
        angles.append([i,0,0])

    # theta = pi, psi = 0
    for i in range(n_2p):
        phi = i*delta
        qlist.append(q(phi,pi,0))
        angles.append([i,pi/delta,0])

    # delta <= theta < 2pi
    for i in range(n_2p):
        for j in range(1,n_p):
            for k in range(n_2p):
                phi,theta,psi = i*delta, j*delta, k*delta
                qlist.append(q(phi,theta,psi))
                angles.append([i,j,k])

    # Now export to a txt file

    # First replace numbers close to zero with 0.0 in each sublist
    qlist = [[0.0 if abs(num) < tolerance else num for num in sublist] for sublist in qlist]

    # Remove negative sign from -0.0 values
    #qlist = [0.0 if num == -0.0 else num for num in qlist]

    # Specify the file path
    file_path = "qlist.txt"

    # Open the file in write mode
    with open(file_path, "w") as file:
        # Write each item in the list to the file
        for item in qlist:
            file.write(str(item) + ',' + '\n')

    print("List exported to", file_path)

    # Define the filename for the output text file
    filename = 'angles.txt'

    # Open the file in write mode
    with open(filename, 'w') as file:
        # Iterate over angles list, formatting and writing to file
        for i in range(0, len(angles), 8):
            # Join 8 sublists into a single line
            line = ','.join(map(str, angles[i:i+8]))
            file.write(line + ',' + '\n')

    print("List exported to", filename)
