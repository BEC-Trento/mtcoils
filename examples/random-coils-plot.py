#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 03-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>
"""
FIXME: plotting random coils the 3D quiver shows that the coil patch is not at the right place.
It must be an issue with the patch, not with the field itself, since the contourf on xy looks ok.
"""
import sys
sys.path.insert(0, '..')

import numpy as np
import matplotlib.pyplot as plt

from numpy.random import random

from mtcoils.coils import CoilSystem


plt.rcParams['image.aspect'] = 1


def random_coil_dict(Rmax=0.2):
    d0 = {
        # a small thin coil
        'Rint': random()*Rmax,
        'Nr': 2,
        'Nz': 2,
        'lr': 1e-2,
        'lz': 1e-2,
        'position': (random(3) - 0.5)*Rmax,
        'normal': (random(3) - 0.5)
    }
    return d0.copy()


Ncoils = 1
D = [random_coil_dict() for _ in range(Ncoils)]

# D = [{
#     # a small thin coil
#     'Rint': 5e-2,
#     'Nr': 3,
#     'Nz': 3,
#     'lr': 1e-2,
#     'lz': 1e-2,
#     'position': (0.07,0,0),
#     'normal': (0.43102544,  0.44422806, -0.00386388)
# }]
print(D)
coil_system = CoilSystem(D)

zz = np.linspace(-10, 10, 21)*1e-2
X,Y,Z = np.meshgrid(zz, zz, zz, indexing='ij')
ix0 = len(zz)//2

Bquad = coil_system.calc_B(X, Y, Z, 100)  # 100 A

# 3D plot
fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111, projection='3d')

S = (slice(None, None, 2),)*3

ax.quiver(X[S], Y[S], Z[S], *Bquad[(slice(None),)+S], length=2e-4)
coil_system.plot_coils(ax, color=[f'C{j}' for j in range(len(D))], )

# 2D contour and field on xy plane
X1 = X[:, :, 0]
Y1 = Y[:, :, 0]


B_xy = Bquad[:,:,:,ix0]  # slice at z = 0
abs_B_xy = np.linalg.norm(B_xy, axis=0)

fig, ax = plt.subplots()
C = ax.contourf(X1*1e2, Y1*1e2, abs_B_xy, )
ax.quiver(X1*1e2, Y1*1e2, B_xy[0, ...], B_xy[1, ...], cmap='viridis')
plt.colorbar(C, ax=ax)
ax.set(xlabel='x [cm]', ylabel='y [cm]')

ax.set_aspect(1)

plt.show()
