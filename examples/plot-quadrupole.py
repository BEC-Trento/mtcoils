#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 03-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>

import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.patches import Circle
# from mpl_toolkits.mplot3d import Axes3D
# from plotpatch3d import pathpatch_2d_to_3d, pathpatch_translate

from coils import Coil, CoilSystem

from scipy.optimize import curve_fit
from scipy.constants import h, atomic_mass, pi
mass = 23*atomic_mass

from ruamel.yaml import YAML
yaml = YAML(typ='safe')

plt.rcParams['image.aspect'] = 1
# plt.rcParams['figure.figsize'] = (6,6)

data = 'magnetic-trap--bec1.yaml'

with open(data) as f:
    D = yaml.load(f)

print(D)

# D['pinch']['position'] = [2.415e-2, 0, 0]

# coil_system = CoilSystem([D[k] for k in ['quadrupole_down', 'quadrupole_up']])
coil_system = CoilSystem([D[k] for k in ['quadrupole_down', 'quadrupole_up', 'pinch']])
# coil_system = CoilSystem([D[k] for k in ['compensation_1', 'compensation_2', 'pinch']])
# coil_system = CoilSystem(D.values())

zz = np.linspace(-50,50, 21)*1e-6
xx = 10 * zz
zz = xx
X,Y,Z = np.meshgrid(xx, zz, zz, indexing='ij')

# zz = np.linspace(-10, 10, 21)*1e-2
# X,Y,Z = np.meshgrid(zz, zz, zz, indexing='ij')
ix0 = len(zz)//2

Bquad = coil_system.calc_B(X,Y,Z, 100)


# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
#
# ax.quiver(X,Y,Z,*Bd)
# Qd.plot_coil(ax)
#
#
fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111, projection='3d')

S = (slice(None, None, 2),)*3

ax.quiver(X[S], Y[S], Z[S], *Bquad[(slice(None),)+S], length=1e-5)
coil_system.plot_coils(ax, color=['C0', 'C0', 'r', 'r', 'g'], )

X1 = X[:, :, 0]
Y1 = Y[:, :, 0]


B_xy = Bquad[:,:,:,ix0]  # slice at z = 0
abs_B_xy = np.linalg.norm(B_xy, axis=0)

fig, ax = plt.subplots()
C = ax.contourf(X1, Y1, abs_B_xy, )
ax.quiver(X1, Y1, B_xy[0, ...], B_xy[1, ...], cmap='viridis')
plt.colorbar(C, ax=ax)
ax.set_aspect(1)

print(zz[ix0])

#
abs_B = np.linalg.norm(Bquad, axis=0)

Vx = abs_B[:, ix0, ix0]
Vy = abs_B[ix0, :, ix0]
Vz = abs_B[ix0, ix0, :]


# def Vpot2D(XY, cx, cy, mx, my, b0, ravel=True):
#     x, y = XY
#     v = cx*(x - mx)**2 + cy*(y - my)**2 + b0
#     if ravel:
#         v = v.ravel()
#     return v
#
#
# p0 = (1e7, 1e7, 0, 0, 0)
# p, cov = curve_fit(Vpot2D, (X1.ravel(), Y1.ravel()), abs_B_xy.ravel(), p0)
# p = {k: v for k, v in zip(('cx', 'cy', 'mx', 'my', 'b0'), p)}
# print(p)
# fitted = Vpot2D((X1, Y1), ravel=False, **p)
# fitx = fitted[:, ix0]
# fity = fitted[ix0, :]


def Vpot3D(XYZ, cx, cy, cz, mx, my, mz, b0, ravel=True):
    x, y, z = XYZ
    v = cx*(x - mx)**2 + cy*(y - my)**2 + cz*(z - mz)**2 + b0
    if ravel:
        v = v.ravel()
    return v


p0 = (1e7, 1e7, 1e7, 0, 0, 0, 0)
p, cov = curve_fit(Vpot3D,
    (X.ravel(), Y.ravel(), Z.ravel()), abs_B.ravel(), p0)
p = {k: v for k, v in zip(('cx', 'cy', 'cz', 'mx', 'my', 'mz', 'b0'), p)}
print(p)
fitted = Vpot3D((X, Y, Z), ravel=False, **p)
fitx = fitted[:, ix0, ix0]
fity = fitted[ix0, :, ix0]
fitz = fitted[ix0, ix0, :]

fig, ax = plt.subplots()
plotting = {'x': (Vx, fitx),
            'y': (Vy, fity),
            'z': (Vz, fitz),
            }
for label, (V, fitt) in plotting.items():
    omega2 = 2*h*0.7e6*p[f'c{label}'] / mass
    print(f"Trap freq {label}: {np.sqrt(omega2)/2/pi:.2f} Hz")
    l, = ax.plot(zz*1e2, V, 'o', label=f'axis {label}')
    ax.plot(zz*1e2, fitt, '-', color=l.get_color())
# ax.set_ylim(top=100)
ax.grid()
ax.legend()
plt.show()
