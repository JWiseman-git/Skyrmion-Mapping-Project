import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import matplotlib as mpl
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
import datetime as datetime
import math
csfont = {'fontname':'Times New Roman'}
'''
Script to plot the underlying helix structure of a crystal lattice
'''

todaydate = '{}'.format(datetime.datetime.now().date())
date_time_obj = datetime.datetime.strptime(todaydate, '%Y-%m-%d').strftime('%m_%d_%y')
date = date_time_obj
csfont = {'fontname':'Times New Roman'}
spinmag = 2.6*9.274010e-24
mu_0 = 4e-7*np.pi
mu_B = 9.274010e-24

class Arrow3D(FancyArrowPatch):

    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1+dx, y1+dy, z1+dz)
        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)

    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs)

def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    """Adding 3d arrow to Axes3D class."""
    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)
setattr(Axes3D,'arrow3D',_arrow3D)


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


def fourier(a, k, pos):
    """Use of a fourier series to map spin orientations"""
    a_xy = np.array([a[0], a[1], 0])
    a_z = np.array([0, 0, a[2]])
    fr = a_z*np.cos(np.dot(k, pos)+np.pi) + a_xy*np.sin(np.dot(k, pos)+np.pi)
    return fr


def ferro(a, k, pos):


    """Use of a fourier series to map spin orientations"""
    fe = (1/2)*((a * np.exp(1j *(np.dot(k, pos)))) + (np.conj(a) * np.exp(-1j * (np.dot(k, pos)))))
    return fe


spinmag = 0.4*9.274010e-24
sep = 4.558e-10
q_p = 93.5
q_mag = 2*np.pi/q_p # unit cells suggested in Lan Trans.
x = np.linspace(-100, 0, 101)
y = np.linspace(0, 100, 101)
z = [0]

size = '{}x{}x{}'.format(2*max(x), 2*max(y), 2*max(z))
q = [q_mag, 0, 0]
q_init = np.array(q)
a = np.array([0, 1, 1])


A = []
rot1 = -np.pi/3

for k in range(0, len(x)):
    for j in range(0, len(y)):
        for p in range(0, len(z)):

            dir = np.array([x[k], y[j], z[p]])
            pos = [x[k], y[j], z[p]]
            trans_a = Rot_Mat(rot1+np.pi , 0, 0, a)
            trans_q = Rot_Mat(rot1, 0, 0, q_init)
            direction = fourier(trans_a, trans_q, dir)
            xspin = direction[0].real
            yspin = direction[1].real
            zspin = direction[2].real
            A.append([pos[0], pos[1], pos[2], xspin, yspin, zspin])

#Q2

B = []
rot2 = np.pi

for k in range(0, len(x)):
    for j in range(0, len(y)):
        for p in range(0, len(z)):

            dir = np.array([x[k], y[j], z[p]])
            pos = [x[k], y[j], z[p]]

            trans_b = Rot_Mat(rot2+np.pi, 0, 0, a)
            trans_q2 = Rot_Mat(rot2, 0, 0, q_init)
            direction = fourier(trans_b, trans_q2, dir)
            xspin = direction[0].real
            yspin = direction[1].real
            zspin = direction[2].real
            B.append([pos[0], pos[1], pos[2], xspin, yspin, zspin])

#Q1

C = []
rot3 = 2*np.pi/6
for k in range(0, len(x)):
    for j in range(0, len(y)):
        for p in range(0, len(z)):

            dir = np.array([x[k], y[j], z[p]])
            pos = [x[k], y[j], z[p]]
            trans_c = Rot_Mat(rot3+np.pi, 0, 0, a)
            trans_q3 = Rot_Mat(rot3, 0, 0, q_init)
            direction = fourier(trans_c, trans_q3, dir)
            xspin = direction[0].real
            yspin = direction[1].real
            zspin = direction[2].real
            C.append([pos[0], pos[1], pos[2], xspin, yspin, zspin])

D = []
rot3 = 0
qferro = np.array([0, 0, 0])
b = np.array([0, 0, 1])
for k in range(0, len(x)):
    for j in range(0, len(y)):
        for p in range(0, len(z)):

            dir = np.array([x[k], y[j], z[p]])
            pos = [x[k], y[j], z[p]]
            direction = fourier(b, qferro, dir)
            xspin = direction[0].real
            yspin = direction[1].real
            zspin = direction[2].real
            D.append([pos[0], pos[1], pos[2], xspin, yspin, 1])
R = []
for i in range(0, len(A)):
    Ar = np.array(A[i][3:])
    Br = np.array(B[i][3:])
    Cr = np.array(C[i][3:])
    Dr = np.array(D[i][3:])
    Total = Dr + Ar + Br + Cr
    print(Total)
    mod = np.sqrt((Total[0].real ** 2) + (Total[1].real ** 2) + (Total[2].real ** 2))
    Total = Total.tolist()
    R.append([A[i][0], A[i][1], A[i][2], Total[0]
              /mod, Total[1]/mod, Total[2]/mod])


fig = plt.figure(figsize=(12, 12))
ax = plt.axes(projection='3d')
ax.set_xlim3d(min(x), max(x))
ax.set_ylim3d(min(y), max(y))
ax.set_zlim3d(0, 1)

#ax.arrow3D(-82, 0, 0,  0, 165, 0, mutation_scale=4,
 #                 arrowstyle="-|>", linestyle='solid', color='red')
for i in range(0, len(R)):
    #if R[i][5] >= 0.7:
    #    ax.arrow3D(R[i][0], R[i][1], R[i][2],  R[i][3],  R[i][4], R[i][5], mutation_scale=4,
     #             arrowstyle="-|>", linestyle='solid', color='red')
    #if  R[i][5] < 0.7:
     #   ax.arrow3D(R[i][0], R[i][1], R[i][2], 4 * R[i][3], 4 * R[i][4], R[i][5], mutation_scale=4,
                            #arrowstyle="-|>", linestyle='solid', color='black')
    if R[i][5] >= 0.8:
        ax.arrow3D(R[i][0], R[i][1], R[i][2], R[i][3], R[i][4], R[i][5], mutation_scale=4,
                  arrowstyle="-|>", linestyle='solid', color='midnightblue')
    elif .8 > R[i][5] >= 0.5:
        ax.arrow3D(R[i][0], R[i][1], R[i][2], R[i][3], R[i][4], R[i][5], mutation_scale=8,
                   arrowstyle="-|>", linestyle='solid', color='mediumblue')
    elif 0.5 > R[i][5] >= 0.3:
        ax.arrow3D(R[i][0], R[i][1], R[i][2], R[i][3], R[i][4], R[i][5], mutation_scale=8,
                   arrowstyle="-|>", linestyle='solid', color='cornflowerblue')
    elif 0.3 > R[i][5] >= 0.1:
        ax.arrow3D(R[i][0], R[i][1], R[i][2],R[i][3], R[i][4], R[i][5], mutation_scale=8,
                   arrowstyle="-|>", linestyle='solid', color='mediumspringgreen')
    elif 0.1 > R[i][5] >= -0.2:
        ax.arrow3D(R[i][0], R[i][1], R[i][2], R[i][3], R[i][4], R[i][5], mutation_scale=8,
                   arrowstyle="-|>", linestyle='solid', color='yellow')
    elif -0.2 > R[i][5] >= -0.7:
        ax.arrow3D(R[i][0], R[i][1], R[i][2], R[i][3], R[i][4], R[i][5], mutation_scale=8,
                  arrowstyle="-|>", linestyle='solid', color='darkorange')
    elif -0.7 > R[i][5] >= -1.5:
        ax.arrow3D(R[i][0], R[i][1], R[i][2], R[i][3], R[i][4], R[i][5], mutation_scale=8,
                   arrowstyle="-|>", linestyle='solid', color='red')

ax.xaxis.pane.fill = False
ax.xaxis.pane.set_edgecolor('white')
ax.yaxis.pane.fill = False
ax.yaxis.pane.set_edgecolor('white')
ax.zaxis.pane.fill = False
ax.zaxis.pane.set_edgecolor('white')
ax.grid(False)
ax.w_zaxis.line.set_lw(0.)
ax.w_xaxis.line.set_lw(1)
ax.w_yaxis.line.set_lw(1)
ax.set_zticks([])
ax.set_xlim([min(x), max(x)])
ax.set_ylim([0, max(y)])
ax.view_init(elev=90, azim=0)
ax.text2D(0.17, 0.8, r"(a)", transform=ax.transAxes,**csfont, fontsize = 14)
ax.set_xlabel(r'$\it{x}/\it{a}$',**csfont, fontsize = 14)
ax.set_ylabel(r'$\it{y}/\it{a}$',**csfont, fontsize = 14)

name = "Square Lattice trial {}".format(q_p, size)
#plt.savefig('Final_Graphs/%s.png' % name, dpi = 500)
plt.tight_layout()
plt.show()






