import numpy as np
import theory
from uncertainties import ufloat


def length_unit_factor(length_units):
    if length_units == 'm':
        return 1.0
    elif length_units == 'cm':
        return 100.


def mean(vl):
    m = ufloat(0., 0.)
    for v in vl:
        m += v
    m /= len(vl)
    return m


def mean2(vl):
    s1 = 0
    s2 = 0
    for v in vl:
        s1 += v.n/v.s**2
        s2 += 1/v.s**2

    return ufloat(s1/s2, np.sqrt(1/s2))


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
        self.rh = theory.sheet_hall_coefficient(db, dr)


class DataPoint:
    def __init__(self, current_set_point, temp_set_point):
        self.current = current_set_point
        self.temp = temp_set_point
        self.vdp = None
        self.hall = None
        self.rs = None
        self.rh = None
        self.rho = None
        self.r_hall = None
        self.n = None
        self.mu = None

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        if isinstance(self.temp, str):
            temp_str = f'RT'
        else:
            temp_str = f'T{self.temp:.1f}'
        if self.current is not None:
            current_str = f'I{self.current:.2e}'
        else:
            current_str = ''
        s = ''
        if self.rs is not None:
            s += 'V'
        if self.rh is not None:
            s += 'H'
        return temp_str+current_str+':'+s

    def set_data(self, vdp: VdpData, hall: HallData):
        self.vdp = vdp
        self.rs = vdp.rs
        self.hall = hall
        self.rh = hall.rh

    def set_results(self, rs, rh):
        self.rs = rs
        self.rh = rh

    def recalculate(self, thickness):
        if self.rs is not None:
            self.rho = theory.resistivity(self.rs, thickness)
        if self.rh is not None:
            self.r_hall = theory.hall_coefficient(self.rh, thickness)
            self.n = theory.density(self.rh, thickness)
        if self.rs is not None and self.rh is not None:
            self.mu = theory.mobility(self.rs, self.rh)


class DataList(list):
    def __init__(self, name, thickness):
        super().__init__()
        self.name = name
        self.thickness = thickness

    def collect_by_variable(self, var):
        if var == 'temp':
            x_list = [x.temp for x in self]
            x_set = set([x.temp for x in self])
        else:
            ValueError('Unknown (not implemented) variable '+var)
            return
        new_list = DataList(self.name+'-'+var, self.thickness)
        for xj in x_set:
            dp = DataPoint(None, xj)
            rs1 = []
            rh1 = []
            for i, xi in enumerate(x_list):
                if xi == xj:
                    rs1.append(self[i].rs.n)
                    rh1.append(self[i].rh.n)
            rs1 = np.array(rs1)
            rh1 = np.array(rh1)
            if var == 'temp':
                dp.set_results(ufloat(np.mean(rs1), np.std(rs1)),
                               ufloat(np.mean(rh1), np.std(rh1)))
            new_list.append(dp)
        return new_list

    def set_thickness(self, thickness):
        self.thickness = thickness
        self.recalculate()

    def recalculate(self):
        for x in self:
            x.recalculate(self.thickness)

    def to_file(self, filename, length_units='cm'):
        fl = length_unit_factor(length_units)

        # print header

        # print data
        for x in self:
            pass




