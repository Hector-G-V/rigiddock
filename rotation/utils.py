import numpy as np
from numpy import sin as Sin
from numpy import cos as Cos

def f_centroid(p_c):
    
    # Function shifts the protein's geometric center (centroid) to the origin.
    # For protien rotation.

    # Do I need to do this function? Yes:
    #   There are symmetry benefits.
    #   There is predictability: rotation doesn't shift the centroid point.

    # p_c: protein coordinates pulled directly from the PDB file. Python List.
    # Return: [centriod, shifted protein coordinates]
    # Return the centroid AND the shifted p_c array; 
    # will need the new array immediately and the centroid later.

    v = np.array(p_c) # Convert to array to do the task
    
    centroid = np.mean(v, axis=0)
    p_c_shifted = np.subtract(v, centroid)

    return [centroid, p_c_shifted]


def q_rotate(b_c, q):
    
    # Rotates all 3-vector p in b_c using quaternion q
    # q: rotation quaternion
    # Output: all rotated p 3-vector from quaternion multiplication qpq^-1
    
    # The function below does the following quaternion multiplication on each point in b_c:
    #2*(p[1]*(q[1]*q[2] - q[0]*q[3]) + p[2]*(q[0]*q[2] + q[1]*q[3])) + p[0]*(q[0]**2 + q[1]**2 - q[2]**2 - q[3]**2), 
    #2*(p[0]*(q[1]*q[2] + q[0]*q[3]) + p[2]*(q[2]*q[3] - q[0]*q[1])) + p[1]*(q[0]**2 - q[1]**2 + q[2]**2 - q[3]**2), 
    #2*(p[0]*(q[1]*q[3] - q[0]*q[2]) + p[1]*(q[0]*q[1] + q[2]*q[3])) + p[2]*(q[0]**2 - q[1]**2 - q[2]**2 + q[3]**2)
    
    q0, q1, q2, q3 = q
    p0, p1, p2 = b_c.T  # p: b_c point to be rotated. p.T to separate the components
    q0q0, q1q1, q2q2, q3q3 = q0**2, q1**2, q2**2, q3**2
    q0q1, q0q2, q0q3 = q0*q1, q0*q2, q0*q3
    q1q2, q1q3, q2q3 = q1*q2, q1*q3, q2*q3

    r0 = 2 * (p1 * (q1q2 - q0q3) + p2 * (q0q2 + q1q3)) + p0 * (q0q0 + q1q1 - q2q2 - q3q3)
    r1 = 2 * (p0 * (q1q2 + q0q3) + p2 * (q2q3 - q0q1)) + p1 * (q0q0 - q1q1 + q2q2 - q3q3)
    r2 = 2 * (p0 * (q1q3 - q0q2) + p1 * (q0q1 + q2q3)) + p2 * (q0q0 - q1q1 - q2q2 + q3q3)

    # Return qpq1(b_c, q)
    #return np.stack((r0, r1, r2), axis=-1)
    return np.round(np.stack((r0, r1, r2), axis=-1), decimals=3)


def R(phi, theta, psi):

    # Euler Matrix for the zx'z' convention
    # phi, theta, psi: radians

    return np.array([
        [
            Cos(phi)*Cos(psi) - Cos(theta)*Sin(phi)*Sin(psi), 
            -Cos(theta)*Cos(psi)*Sin(phi) - Cos(phi)*Sin(psi),
            Sin(theta)*Sin(phi)
        ], 
        [
            Cos(psi)*Sin(phi) + Cos(theta)*Cos(phi)*Sin(psi), 
            Cos(theta)*Cos(phi)*Cos(psi) - Sin(phi)*Sin(psi),
            -Cos(phi)*Sin(theta)
        ],
        [
            Sin(theta)*Sin(psi), 
            Cos(psi)*Sin(theta), 
            Cos(theta)
        ]
    ])


def e_rotate(angles, p):
    
    # Rotates point p using euler angles
    # angles: Python List or numpy array
    # p: numpy array
    
    phi, theta, psi = angles

    # Euler Matrix for the zx'z' convention

    matrix = R(phi, theta, psi)

    if np.array(p).shape==(3,):
        return np.round(np.matmul(matrix, p),decimals=3)
    else:
        return [
            np.round(np.matmul(matrix,p_i), decimals=3) for p_i in p
        ]
