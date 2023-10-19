import numpy as np
import theory
from uncertainties import ufloat


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
    def __init__(self, current_setpoint, temp_setpoint):
        self.thickness = None
        self.current = current_setpoint
        self.temp = temp_setpoint
        self.vdp = None
        self.hall = None

