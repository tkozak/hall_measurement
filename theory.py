import numpy as np
import scipy.optimize as opt


def vdp_function(rs, ra, rb):
    return np.exp(-np.pi*ra/rs) + np.exp(-np.pi*rb/rs) - 1


def vdp_function_prime(rs, ra, rb):
    return np.pi/rs**2*(ra*np.exp(-np.pi*ra/rs) + rb*np.exp(-np.pi*rb/rs))


def sheet_resistance_vdp(ra, rb, ra_s=0, rb_s=0):
    r0 = np.pi * (ra + rb) / (2 * np.log(2))
    rs_s0 = np.pi/(2 * np.log(2)) * np.sqrt(ra_s**2 + rb_s**2)
    rs = opt.newton(vdp_function, r0, fprime=vdp_function_prime, args=(ra, rb))
    return rs, rs_s0

