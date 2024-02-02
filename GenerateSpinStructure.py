from Constants import sep, z_axis_separation
import numpy as np

class skyrmion:

    zaxis_vector = np.array([0, 0, 1])
    xaxis_vector = np.array([1, 0, 0])
    wavevector_1 = 1
    wavevector_2 = 1
    wavevector_3 = 1

    def __init__(self, dim_x: int, dim_y: int, dim_z: int, wavevector_wavelength: float):
        self.latticesize_xaxis = np.linspace(-dim_x, dim_x, (2*dim_x)+1)
        self.latticesize_yaxis = np.linspace(-dim_y, dim_y, (2*dim_y)+1)
        self.latticesize_zaxis = np.linspace(-dim_z, dim_z, (2*dim_z)+1)
        self.wavevector_wavelength = wavevector_wavelength


"""
After specifying a lattice size, an atom's spin orientation at a given coordinate is determined by a fourier transform 
of a helical spin structure. 
"""

zaxis_vector = np.array([0, 0, 1])
xaxis_vector = np.array([1, 0, 0])

latticesize_xaxis = np.linspace(-5, 5, 11)
latticesize_yaxis = np.linspace(-5, 5, 11)
latticesize_zaxis = np.linspace(-5, 5, 11)


def Spin_Helix(a, b, k, pos):
    """Use of a fourier transform to map spin orientations"""
    A = a * np.exp(1j * np.dot(k, pos))
    B = b * np.exp(-1j * np.dot(k, pos))
    s_d = [B]
    return s_d


def Rot_Mat(phi, theta, psi, r):
    """Defining the rotation matrix using Euler angles"""

    R = np.zeros(shape=(3, 3))
    R[0, 0] = (np.cos(psi)*np.cos(phi)) - (np.cos(theta)*np.sin(phi)*np.sin(psi))
    R[0, 1] = (np.cos(psi) * np.sin(phi)) + (np.cos(theta) * np.cos(phi) * np.sin(psi))
    R[0, 2] = np.sin(psi)*np.sin(theta)
    R[1, 0] = -(np.sin(psi) * np.cos(phi)) - (np.cos(theta) * np.sin(phi) * np.cos(psi))
    R[1, 1] = -(np.sin(psi) * np.sin(phi)) + (np.cos(theta) * np.cos(phi) * np.cos(psi))
    R[1, 2] = np.cos(psi) * np.sin(theta)
    R[2, 0] = np.sin(theta)*np.sin(phi)
    R[2, 1] = -np.sin(theta)*np.cos(phi)
    R[2, 2] = np.cos(theta)
    result = np.dot(R, r)
    return result


"""
Populating a results list, with each element representing a lattice point and it's associated spin
"""

k1 = np.array([np.pi, np.pi, np.pi]) # ---- > modulation vectors
a = np.array([1, 0, 0])
b = np.array([0, 0, 1])

A = []
x_length, y_length, z_length = latticesize_xaxis, latticesize_yaxis, latticesize_zaxis

"""
Run over entire structure to generate spin surface
"""

for x, y, z in np.ndindex(len(x_length), len(y_length), len(z_length)):
    position_in_lattice = [x_length[x], y_length[y], z_length[z]]
    position_in_lattice_scaled = [x_length[x], y_length[y]]
    position_in_lattice_scaled = [i * sep for i in position_in_lattice_scaled[:2]]
    position_in_lattice_scaled.append(z_length[z] * z_axis_separation)

    print('position of point in lattice:', position_in_lattice_scaled)

# for x in range(0, len(x_length)):
#     for y in range(0, len(y_length)):
#         for z in range(0,len(z_length)):
#             dir = [x_length[k], y_length[j], z_length[z]]
#             pos = [x[k], y[j]]
#             pos = [i * sep for i in pos[:2]]


#             scale = 3.18e-10
#             pos.append(z[p]*scale)
#             print(pos)
#
#             direction = Spin_Helix(a, b, k1, dir)
#
#             spintest = sep*(direction[0][2].real)
#             A.append([pos[0], pos[1], pos[2], 0, 0, spintest])

#Sphere_s_e = create_sphere(pos[0],pos[1],pos[2],rvalue)



