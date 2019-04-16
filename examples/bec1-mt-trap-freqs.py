#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Create: 03-2019 - Carmelo Mordini <carmelo> <carmelo.mordini@unitn.it>
# import sys
# sys.path.insert(0, '..')

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.constants import pi, h, atomic_mass

from mtcoils.coils import CoilSystem
from mtcoils.setup import setup_bec1 as D

import pprint

mass = 23*atomic_mass  # ^23 Na

plt.rcParams['image.aspect'] = 1

pprint.pprint(D)
coil_system = CoilSystem(D.values())

# small region around 0
zz = np.linspace(-50,50, 21)*1e-6
xx = 10 * zz
X,Y,Z = np.meshgrid(xx, zz, zz, indexing='ij')

# zz = np.linspace(-10, 10, 21)*1e-2
# X,Y,Z = np.meshgrid(zz, zz, zz, indexing='ij')
ix0 = len(zz)//2

Bquad = coil_system.calc_B(X,Y,Z, 45)  # choose one current

X1 = X[:, :, 0]
Y1 = Y[:, :, 0]

B_xy = Bquad[:,:,:,ix0]  # slice at z = 0, where Bz = 0
abs_B_xy = np.linalg.norm(B_xy, axis=0)

fig, ax = plt.subplots()
C = ax.contourf(X1*1e6, Y1*1e6, abs_B_xy, )
ax.quiver(X1*1e6, Y1*1e6, B_xy[0, ...], B_xy[1, ...], cmap='viridis')
plt.colorbar(C, ax=ax)
ax.set(xlabel='x [um]', ylabel='y [um]')
# ax.set_aspect(1)

print(zz[ix0])


abs_B = np.linalg.norm(Bquad, axis=0)

# Shift by gravity to include sag
gravity = 9.806
B_gravity = mass * gravity * Z /(h * 0.7e6)

abs_B += B_gravity

# slices along x, y, z axes
Vx = abs_B[:, ix0, ix0]
Vy = abs_B[ix0, :, ix0]
Vz = abs_B[ix0, ix0, :]


# parabolic fit to B field modulus
def Vpot3D(XYZ, cx, cy, cz, mx, my, mz, b0, ravel=True):
    x, y, z = XYZ
    v = cx*(x - mx)**2 + cy*(y - my)**2 + cz*(z - mz)**2 + b0
    if ravel:
        v = v.ravel()
    return v


p0 = (1e7, 1e7, 1e7, 0, 0, 0, 0)
p, cov = curve_fit(Vpot3D,
    list(map(np.ravel, (X, Y, Z))), abs_B.ravel(), p0)
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

z_sag = - mass * gravity / (2*h*0.7e6*p['cz'])
z_sag_fit = p['mz']
print(f"grav. sag: {z_sag*1e6:.3f} um (fitted: {z_sag_fit*1e6:.3f} um)")
ax.axvline(z_sag*1e2, ls='--', color=l.get_color())
# ax.set_ylim(top=100)
ax.grid()
ax.legend()
plt.show()
