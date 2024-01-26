from Constants import sep, spinmag
from GenerateLattice import r_mod

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
