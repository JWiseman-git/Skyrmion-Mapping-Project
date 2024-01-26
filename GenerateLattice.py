from ArrowGenerator import Arrow3D
from Constants import spinmag, mu_0, mu_B, psi, sep, sep_star
import numpy as np
import math
import cmath

"""
After specifing a lattice size, an atom's spin orientation at a given coordinate is determined by a Fourier transform 
of a helical spin structure. Following constructing the resulting dipole (mag) field at a given position is calculable. 
The determination of a local field is used for determining muon precession. 
"""

def uv(vec):
    """normalising wavector"""
    return vec/np.linalg.norm(vec)

def fourier(a,b,k,pos):
    """Use of a fourier series to map spin orientations"""
    A = (a) * np.exp(1j * np.dot(k, pos))
    B = (b) * np.exp(-1j * np.dot(k, pos))
    s_d = [B]
    return s_d

def r_mod(r0,r1):
    """Function to determine |r| between the field generating dipole and the point being considered"""
    vector = []
    zip_object = zip(r0, r1)
    for r0_i, r1_i in zip_object:
        vector.append(r1_i - r0_i)
    modulus = np.abs(np.sqrt(((vector[0]) ** 2) + ((vector[1]) ** 2) + ((vector[2]) ** 2)))
    return modulus

def mu_s(sx,sy,sz):
    """Function to produce m list from the direction of the spin vector S"""
    moment = [sx, sy, sz]
    moment = [i * spinmag for i in moment]
    return moment

def dipoletensor(r0,r1,abs_r):
    """Function to produce the Dipole Tensor"""
    r = []
    zip_object = zip(r0, r1)
    for r0_i, r1_i in zip_object:
        r.append(r1_i - r0_i)

    prefactor = 1e-7/(abs_r**3)
    matrix = np.zeros(shape = (3,3))
    #Populating the matrix
    matrix[0, 0] = -1 + (3*((r[0])**2))/(abs_r**2)
    matrix[0, 1] = (3*r[0]*r[1])/(abs_r**2)
    matrix[0, 2] = (3*r[0]*r[2])/(abs_r**2)
    matrix[1, 0] = (3*r[1]*r[0])/(abs_r**2)
    matrix[1, 1] = -1 + (3*((r[1])**2))/(abs_r**2)
    matrix[1, 2] = (3*r[1]*r[2])/(abs_r**2)
    matrix[2, 0] = (3*r[2]*r[0])/(abs_r**2)
    matrix[2, 1] = (3*r[2]*r[1])/(abs_r**2)
    matrix[2, 2] = -1 + (3*((r[2])**2))/(abs_r**2)
    tensor = matrix*prefactor
    return tensor

def B_field(dtensor,moment):
    """Calculating the resultant magnetic field vector"""
    B = np.dot(dtensor, moment).reshape(1,3)
    B = B.tolist()
    return B

def analytical_sol(moment,r0,r1,abs_r):
    """Analytical solution of the dipole field equation"""

    r = []
    zip_object = zip(r0, r1)
    for r0_i, r1_i in zip_object:
        r.append(r1_i - r0_i)

    r_unit = np.array(r)/abs_r
    m_dot_R = np.dot(moment, r_unit)
    moment = np.array([moment])

    prefactor = 1e-7/(abs_r**3)
    a = (3 * m_dot_R * r_unit)
    b = moment
    B = prefactor*(a-b)
    return B


def create_sphere(cx,cy,cz, r, resolution=360):
    '''
    Generating a sphere with center (cx, cy, cz) and radius r
    '''
    phi = np.linspace(0, 2*np.pi, 2*resolution)
    theta = np.linspace(0, np.pi, resolution)

    theta, phi = np.meshgrid(theta, phi)

    r_xy = r*np.sin(theta)
    x = cx + np.cos(phi) * r_xy
    y = cy + np.sin(phi) * r_xy
    z = cz + r * np.cos(theta)

    return np.stack([x,y,z])

'''
Populating a results list, with each element representing a lattice point and it's associated spin
'''

k1 = np.array([np.pi,np.pi,np.pi])
a = np.array([1, 0, 0])
b = np.array([0, 0, 1])

A = []
x = np.linspace(-5,5,11)
y = np.linspace(-5,5,11)
z = np.linspace(-5,5,11)

#for l in range(0, len(z)):
for k in range(0, len(x)):
    for j in range(0, len(y)):
        for p in range(0,len(z)):
            dir = [x[k], y[j], z[p]]
            pos = [x[k], y[j]]
            pos = [i * 4.7e-10 for i in pos[:2]]
            scale = 3.18e-10
            pos.append(z[p]*scale)
            print(pos)

            direction = fourier(a, b, k1, dir)

            spintest = sep*(direction[0][2].real)
            A.append([pos[0], pos[1], pos[2], 0, 0, spintest])

#Sphere_s_e = create_sphere(pos[0],pos[1],pos[2],rvalue)

#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax.set_xlim3d(-max(x)*sep, max(x)*sep) # Change to meet the new conditions
#ax.set_ylim3d(-max(y)*sep, max(y)*sep)
#ax.set_zlim3d(-max(z)*sep, max(z)*sep)
#ax.plot3D(pos[0], pos[1], pos[2], marker='x', markersize=10, color='red')
#for i in range(0, len(A)):
#    ax.plot3D(A[i][0], A[i][1], A[i][2], marker='o', markersize=1, color='black')
#    ax.arrow3D(A[i][0], A[i][1], A[i][2], A[i][3], A[i][4], A[i][5], mutation_scale=10,
#               arrowstyle="-|>", linestyle='solid', color ='black')
#
#ax.plot_surface(Sphere_s_e[0], Sphere_s_e[1], Sphere_s_e[2], rstride=1, cstride=1, color='r', alpha=0.3, linewidth=0)
#totalB_ang = cmath.phase(totalB[2])
#dir_mag_x = totalB_mag*math.cos(totalB_ang)
#dir_mag_y = totalB_mag*math.sin(totalB_ang)
#dir_mag_z = totalB_mag*math.sin(totalB_ang)
#ax.arrow3D(pos[0], pos[1], pos[2], dir_mag_x, 0, dir_mag_z, mutation_scale=10,
           #arrowstyle="-|>", linestyle='solid', color ='r')
#Plotting the analytical results
#ax.xaxis.pane.fill = False
#ax.xaxis.pane.set_edgecolor('white')
#ax.yaxis.pane.fill = False
#ax.yaxis.pane.set_edgecolor('white')
#ax.zaxis.pane.fill = False
#ax.zaxis.pane.set_edgecolor('white')
#ax.grid(False)
#ax.w_zaxis.line.set_lw(0.)
#ax.set_zticks([])
#ax.set_xlabel("’$\AA$’")
#plt.tight_layout()
#plt.show()
#cor = np.array(direction)
        #tol = 1e-14
        #cor.real[abs(cor.real) < tol] = 0.0
        #cor.imag[abs(cor.imag) < tol] = 0.0
        #direction = cor.tolist()

