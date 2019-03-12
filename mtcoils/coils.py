import numpy as np
from functools import partial
from scipy.special import ellipk, ellipe
from scipy.constants import pi, mu_0

from matplotlib.patches import Wedge

from .plotpatch3d import pathpatch_2d_to_3d, cycle_kwargs


class Coil():
    '''
    ref: https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20010038494.pdf
    '''

    def __init__(self, Rint, Nr, Nz, lr, lz, position=(0, 0, 0), normal=(0, 0, 1)):
        """A single coil
         Args:
             Rint: internal radius.
             Nr: number of radial turns
             Nz: number of axial turns
             lr: radial wire thickness
             lz: axial wire thickness
             position (tuple(float)): a 3-tuple or array-like, the coordinates of
                the coil center
             normal (tuple(float)): a 3-tuple or array-like, components of the
                vector normal to the coil. This defines the current sign according to
                the right hand's rule.

        """
        self.position = position
        self.normal = normal

        rs = Rint + lr * (np.arange(Nr) + 0.5)
        ds = lz * (np.arange(-Nz / 2, Nz / 2) + 0.5)
        self.rs, self.ds = np.meshgrid(rs, ds, indexing='ij')
        self.coil_shape = self.rs.shape
        # print(self.coil_shape, (Nr, Nz))
        # assert self.coil_shape == (Nr, Nz)

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, position):
        assert len(position) == 3
        self._position = np.asarray(position)

    @property
    def normal(self):
        return self._normal

    @normal.setter
    def normal(self, normal):
        assert len(normal) == 3
        self._normal = np.asarray(normal)
        self._normal = self._normal / np.linalg.norm(self._normal)

    def _get_coords(self, rho, z):
        # this will raise a ValueError if rho and z have wrong shape
        rho, z = np.broadcast_arrays(rho, z)
        rho, z, rr, dd = map(partial(np.broadcast_to, shape=self.coil_shape + rho.shape),
                             [rho, z,
                              self.rs.reshape(self.rs.shape + (1,) * rho.ndim),
                              self.ds.reshape(self.rs.shape + (1,) * rho.ndim),
                              ])

        return rho, z, rr, dd

    def _simple_Bz2(self, rho, z, I, R):
        alpha2 = R**2 + rho**2 + z**2 - 2 * R * rho
        beta2 = R**2 + rho**2 + z**2 + 2 * R * rho
        m = 1 - alpha2 / beta2
        Ek = ellipe(m)
        Kk = ellipk(m)
        Bz = 1e4 * mu_0 * I / pi / \
            (2 * alpha2 * np.sqrt(beta2)) * \
            ((R**2 - rho**2 - z**2) * Ek + alpha2 * Kk)
        return Bz

    def _simple_Br2(self, rho, z, I, R):
        alpha2 = R**2 + rho**2 + z**2 - 2 * R * rho
        beta2 = R**2 + rho**2 + z**2 + 2 * R * rho
        gamma = np.where(rho != 0, z / rho, 0)
        m = 1 - alpha2 / beta2
        Ek = ellipe(m)
        Kk = ellipk(m)
        Br = 1e4 * mu_0 * I / pi / \
            (2 * alpha2 * np.sqrt(beta2)) * gamma * \
            ((R**2 + rho**2 + z**2) * Ek - alpha2 * Kk)
        return Br

    def calc_Bfield1(self, r, I,):
        # Keep it simple
        assert len(r) == 3
        r = np.asarray(r)
        X = r - self.position
        z = np.dot(X, self.normal)
        hat_u = (X - z * self.normal)
        rho = np.linalg.norm(hat_u)
        if rho > 0:
            hat_u /= rho

        rho, z, rr, dd = self._get_coords(rho, z)
        Br = self._simple_Br2(rho, z - dd, I, rr).sum(axis=(0, 1))
        Bz = self._simple_Bz2(rho, z - dd, I, rr).sum(axis=(0, 1))
        return Br * hat_u + Bz * self.normal

    def calc_B(self, x, y, z, I,):
        x, y, z = np.broadcast_arrays(x, y, z)
        addshape = (-1,) + (1,) * x.ndim
        r = np.concatenate([q[np.newaxis, ...] for q in (x, y, z)], axis=0)
        X = r - self.position.reshape(addshape)
        # print('X,r', X.shape, r.shape)
        z = np.sum(X * self.normal.reshape(addshape), axis=0)
        hat_u = (X - z[None, ...] * self.normal.reshape(addshape))
        rho = np.linalg.norm(hat_u, axis=0)
        # print(rho.shape)
        # if rho > 0:
        hat_u = np.where(rho > 0, hat_u / rho, hat_u)

        rho, z, rr, dd = self._get_coords(rho, z)
        Br = self._simple_Br2(rho, z - dd, I, rr).sum(axis=(0, 1))
        Bz = self._simple_Bz2(rho, z - dd, I, rr).sum(axis=(0, 1))
        return Br * hat_u + Bz * self.normal.reshape(addshape)

    def plot_coil(self, ax, **patch_kw):
        rr = self.rs.max()
        width = self.rs.ptp()
        p = Wedge((0,0), rr, 0, 360, width=width, **patch_kw)
        ax.add_patch(p)
        pathpatch_2d_to_3d(p, normal=self.normal, translate=self.position)
        # pathpatch_translate(p, self.position)


class CoilSystem():
    def __init__(self, coils_list):
        self.coils = [Coil(**item) for item in coils_list]

    def calc_B(self, x, y, z, I):
        Bfields = [coil.calc_B(x, y, z, I) for coil in self.coils]
        return np.sum(Bfields, axis=0)

    def plot_coils(self, ax, **patch_kw):
        cc = cycle_kwargs(**patch_kw)
        if len(cc) < len(self.coils):
            print("WARNING: missing keyword arguments in patch_kw list. Some coil will not be plotted")
        for d, coil in zip(cc, self.coils):
            coil.plot_coil(ax, **d)
