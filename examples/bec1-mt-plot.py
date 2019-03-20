#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 03-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>
import sys
sys.path.append('..')

import numpy as np
import matplotlib.pyplot as plt

from mtcoils.coils import CoilSystem
from mtcoils.setup import setup_bec1 as D


plt.rcParams['image.aspect'] = 1


print(D)
coil_system = CoilSystem(D.values())

zz = np.linspace(-10, 10, 21)*1e-2
X,Y,Z = np.meshgrid(zz, zz, zz, indexing='ij')
ix0 = len(zz)//2

Bquad = coil_system.calc_B(X, Y, Z, 100)  # 100 A

# 3D plot
fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111, projection='3d')

S = (slice(None, None, 2),)*3

ax.quiver(X[S], Y[S], Z[S], *Bquad[(slice(None),)+S], length=2e-5)
coil_system.plot_coils(ax, color=['C0', 'C0', 'r', 'r', 'g'], )

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
