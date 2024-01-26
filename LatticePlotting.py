import matplotlib.pyplot as plt

"""
Plotting of a given lattice and 'x' - the point at which the local magnetic field strength is being measured.
Arrows indicate the orientation of individual atomic spins. 
"""

Sphere_s_e = create_sphere(pos[0],pos[1],pos[2],rvalue)

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