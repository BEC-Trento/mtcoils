---
title: "`mtcoils` math"
author: Carmelo Mordini
---
# mtcoils
Static B field calculator for arbitrary configurations of coils

## The math

Each coil is defined by

- the (inner) radius $R_{in}$
- a number of radial and axial turns $N_r$ and $N_z$
- the radial and axial wire thickness $l_r$, $l_z$

And its field is obtained by summing $N_r \times N_z$ coils of radius $R_{in} + l_r(i + 1/2)$ $(i = 0 \ldots N_r-1)$ posititoned at $z_j = l_z(j - (N_z-1)/2)$ $(j = 0 \ldots N_z-1)$. The composite coil is placed at the origin, the plane $z=0$ corresponding to the middle of the block. For the radial and azimuthal components of the field we use Eq. 24-25 in [1].

Arbitrary configurations are specified by a position vector $r_0$ and a unit vector $\hat n$ normal to the plane of the coil. The field at a point $r$ will lay on the plane specified by $x = r - r_0$ and $\hat n$ and will be given by

$$\left\{
    \begin{aligned}
    x &= z\, \hat n + \rho\, \hat u \\
    \vec B &= B_z\,\hat n + B_r\,\hat u
    \end{aligned}
\right.$$
