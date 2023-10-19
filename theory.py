import numpy as np
import scipy.optimize as opt
from uncertainties import ufloat


def vdp_function(rs, ra, rb):
    return np.exp(-np.pi*ra/rs) + np.exp(-np.pi*rb/rs) - 1


def vdp_function_prime(rs, ra, rb):
    return np.pi/rs**2*(ra*np.exp(-np.pi*ra/rs) + rb*np.exp(-np.pi*rb/rs))


def sheet_resistance_vdp(ra, rb):
    r0 = np.pi * (ra + rb) / (2 * np.log(2))

    rs_v = opt.newton(vdp_function, r0.n, fprime=vdp_function_prime, args=(ra.n, rb.n))
    rs = ufloat(rs_v, r0.s)
    return rs


def hall_coefficient(db, dr):
    return dr/db


def density_mobility_majority(rs, rh, thickness):
    pass
