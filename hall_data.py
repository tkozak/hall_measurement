import numpy as np
import theory
from uncertainties import ufloat

ec = 1.602e-19

class VdpData:
    def __init__(self, current, voltage):
        self.i = current
        self.u = voltage

        if self.i.shape[0] == 2:   # symmetric current measurement
            self.rc = np.squeeze(np.diff(self.u, axis=0) / np.diff(self.i, axis=0))  # resistance without dc offset
            self.offset = self.u[1, :] - self.rc * self.i[1, :]  # save dc offset
        else:
            raise ValueError("Current array must have 2 rows")

        # rc_s = 0.02 * rc  # estimate error as 2% of the calculated resistance value
        # rc = unumpy.uarray(rc, 0.02*rc)

        ra_n = 0.5 * (self.rc[0] + self.rc[1])  # average two symmetric cases
        # estimate error as the difference between one value and the mean
        self.ra = ufloat(ra_n, np.abs(self.rc[0] - ra_n))
        # possible we can add some systematic error

        rb_n = 0.5 * (self.rc[2] + self.rc[3])  # average two symmetric cases
        # estimate error as the difference between one value and the mean
        self.rb = ufloat(rb_n, np.abs(self.rc[2] - rb_n))
        self.rs = theory.sheet_resistance_vdp(self.ra, self.rb)


class HallData:
    def __init__(self, current, voltage, field):
        self.i = current
        self.u = voltage
        self.b = field

        if self.i.shape[0] == 2:  # symmetric current measurement
            self.rc = np.squeeze(np.diff(self.u, axis=0) / np.diff(self.i, axis=0))  # resistance without dc offset
            self.offset = self.u[1, :] - self.rc * self.i[1, :]  # save dc offset
        else:
            raise ValueError("Current array must have 2 rows")

        dr_v = np.array([self.rc[0] - self.rc[1], self.rc[2] - self.rc[3]])
        dr_n = 0.5 * (dr_v[0] + dr_v[1])  # average both resistance differences
        dr = ufloat(dr_n, np.abs(dr_v[0] - dr_n))  # estimate error as difference between dr
        # perhaps too crude, as this is seldom equal, but the average seems to be quite reproducible

        db = ufloat((self.b[0] - self.b[1]) * 1e-4, 0.01 * np.sqrt(2))  # estimate error as 100 G for one measurement
        self.rh = theory.hall_coefficient(db, dr)


class DataPoint:
    def __init__(self, current_set_point, temp_set_point):
        self.current = current_set_point
        self.temp = temp_set_point
        self.thickness = None
        self.vdp = None
        self.hall = None
        self.rho = None
        self.r_hall = None
        self.n = None
        self.mu = None

    def __str__(self):
        if isinstance(self.temp, str):
            return f'I = {self.current:.2e} A, T = {self.temp}: {self.rho} \u2126m'
        else:
            return f'I = {self.current:.2e} A, T = {self.temp:.1f} K: {self.rho} \u2126m'

    def __repr__(self):
        if isinstance(self.temp, str):
            return f'I{self.current:.2e} T{self.temp}: R{self.rho.n:.3e}'
        else:
            return f'I{self.current:.2e} T{self.temp:.1f}: R{self.rho.n:.3e}'


    def set_thickness(self, thickness):
        self.thickness = thickness
        self.recalculate()

    def set_vdp(self, vdp: VdpData):
        self.vdp = vdp
        self.recalculate()

    def set_hall(self, hall: HallData):
        self.hall = hall
        self.recalculate()

    def recalculate(self):
        if self.thickness is not None:
            if self.vdp is not None:
                self.rho = self.vdp.rs * self.thickness
            if self.hall is not None:
                self.r_hall = self.hall.rh * self.thickness
                self.n = 1./(ec * self.r_hall)
        if self.vdp is not None and self.hall is not None:
            self.mu = self.hall.rh/self.vdp.rs
