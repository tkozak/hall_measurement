import numpy as np
import scipy.optimize as opt
import uncertainties
from uncertainties import umath


def vdp_function(rs, ra, rb):
    return np.exp(-np.pi*ra/rs) + np.exp(-np.pi*rb/rs) - 1


def vdp_function_prime(rs, ra, rb):
    return np.pi/rs**2*(ra*np.exp(-np.pi*ra/rs) + rb*np.exp(-np.pi*rb/rs))


# def sheet_resistance_vdp(ra, rb):
#     r0 = np.pi * (ra + rb) / (2 * np.log(2))
#     rs = opt.newton(vdp_function, r0, fprime=vdp_function_prime, args=(ra, rb))
#     return rs

def sheet_resistance_vdp(ra, rb, delta=5e-4):
    zi1 = 2*np.log(2)/np.pi/(ra+rb)
    while True:
        if isinstance(ra, uncertainties.UFloat):
            yi = umath.exp(-np.pi * zi1 * ra) + umath.exp(-np.pi * zi1 * rb)
            dzi = - (1 - yi) / np.pi / (ra * umath.exp(-np.pi * zi1 * ra) + rb * umath.exp(-np.pi * zi1 * rb))
        else:
            yi = np.exp(-np.pi * zi1 * ra) + np.exp(-np.pi * zi1 * rb)
            dzi = - (1 - yi) / np.pi / (ra * np.exp(-np.pi * zi1 * ra) + rb * np.exp(-np.pi * zi1 * rb))
        zi = zi1 + dzi
        if (zi-zi1)/zi < delta:
            break
        zi1 = zi
    return 1/zi


def hall_coefficient(dr, db):
    return dr/db


def density_mobility_majority(rs, rh, thickness):
    pass
