from Constants import sep, spinmag
from GenerateSpinStructure import r_mod

"""
After specifying a position within the lattice the local field is calculated using numerical and analytical methods.
Difference between analytical and numerical solutions calculated. A convergence test is also performed between the two
"""

position_cartesian = [2.2, 0, 0]
position_scaled = [i * 1e-10 for i in position_cartesian]
rs = []
moments = []
Bf = []
ana_Bf = []
Bf = []
rvalue = 10*sep

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

for i in range(0, len(A)):
    rs.append(r_mod(position_scaled, A[i][:3]))
    s_f = spinmag/sep
    moments.append([A[i][3]*s_f, A[i][4]*s_f, A[i][5]*s_f])

I_vals = []
for i in range(0, len(A)):
    if rs[i] <= rvalue:
        I_vals.append(i)
    else:
        i = i

for i in range(0, len(I_vals)):
    j = I_vals[i]
    t = dipoletensor(pos, A[j][:3], rs[j])
    m_a = np.array(moments[j]).reshape((3, 1))
    [v] = B_field(t, m_a)
    Bf.append(v)

totalB = [sum(y) for y in zip(*Bf)]

totalB_mag = (np.sqrt((totalB[0])*np.conj(totalB[0]) + (totalB[1])*np.conj(totalB[1])
                      + (totalB[2])*np.conj(totalB[2])))
totalB_mag = totalB_mag.real

for i in range(0,len(I_vals)):
    j = I_vals[i]
    m_a = np.array(moments[j])

    loc = np.array(A[j][:3])

    [r] = analytical_sol(m_a,pos, loc, rs[j])

    ana_Bf.append(r)

totalB_ana = [sum(y) for y in zip(*ana_Bf)]
totalB_ana_mag = (np.sqrt((totalB_ana[0])*np.conj(totalB_ana[0]) + (totalB_ana[1])*np.conj(totalB_ana[1])
                      + (totalB_ana[2])*np.conj(totalB_ana[2])))
totalB_ana_mag = totalB_ana_mag.real

print('Numerical B vector',totalB)
print('Analytical B vector',totalB_ana)
print('Analytical solution:', totalB_ana_mag)
print("Numerical Solution:", totalB_mag)
print('Difference:', totalB_mag - totalB_ana_mag)
