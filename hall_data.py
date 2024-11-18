import csv
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
        s1 += v.n / v.s ** 2
        s2 += 1 / v.s ** 2

    return ufloat(s1 / s2, np.sqrt(1 / s2))


class VdpData:
    def __init__(self, current, voltage):
        self.i = current
        self.u = voltage

        if self.i.shape[0] == 2:  # symmetric current measurement
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

        self.dr_v = np.array([self.rc[0] - self.rc[1], self.rc[2] - self.rc[3]])
        dr_n = 0.5 * (self.dr_v[0] + self.dr_v[1])  # average both resistance differences
        self.dr = ufloat(dr_n, np.abs(self.dr_v[0] - dr_n))  # estimate error as difference between dr
        # perhaps too crude, as this is seldom equal, but the average seems to be quite reproducible

        db = ufloat((self.b[0] - self.b[1]) * 1e-4, 0.01 * np.sqrt(2))  # estimate error as 100 G for one measurement
        self.rh = theory.sheet_hall_coefficient(db, self.dr)


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
        return temp_str + current_str + ':' + s

    def set_data(self, vdp: VdpData, hall: HallData):
        self.vdp = vdp
        if self.vdp is not None:
            self.rs = vdp.rs
        self.hall = hall
        if self.hall is not None:
            self.rh = hall.rh

    def set_results(self, rs, rh):
        self.rs = rs
        self.rh = rh

    def recalculate(self, thickness):
        if self.rs is not None:
            self.rho = theory.resistivity(self.rs, thickness)
        if self.rh is not None:
            self.r_hall = theory.hall_coefficient(self.rh, thickness)
            if self.rh.n == 0:
                self.n = ufloat(np.nan, 0.)
            else:
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
            x_tol = 0.1
            x_list = []
            for x in self:
                if type(x.temp) == str:
                    x_list.append(x.temp)
                else:
                    x_list.append(round(x.temp, 0))
            x_set = set(x_list)
        else:
            ValueError('Unknown (not implemented) variable ' + var)
            return
        new_list = DataList(self.name + '-' + var, self.thickness)
        for xj in x_set:
            dp = DataPoint(None, xj)
            rs1 = []
            rh1 = []
            for i, xi in enumerate(x_list):
                if type(xj) == str:
                    if xi == xj:
                        if self[i].rs is not None:
                            rs1.append(self[i].rs.n)
                        if self[i].rh is not None:
                            rh1.append(self[i].rh.n)
                else:
                    if abs(xi - xj) < x_tol:
                        if self[i].rs is not None:
                            rs1.append(self[i].rs.n)
                        if self[i].rh is not None:
                            rh1.append(self[i].rh.n)
            rs1 = np.array(rs1)
            rh1 = np.array(rh1)
            if rh1.size > 0:
                rh = ufloat(np.mean(rh1), np.std(rh1))
            else:
                rh = None
            if rs1.size > 0:
                rs = ufloat(np.mean(rs1), np.std(rs1))
            else:
                rs = None
            if var == 'temp':
                if rh is not None:
                    if rh.s > rh.n:
                        rh = ufloat(np.mean(rh1), 0.8*np.mean(rh1))
                dp.set_results(rs, rh)
            new_list.append(dp)
        # sort
        if var == 'temp':
            new_list.sort(key=lambda z: z.temp)
        return new_list


    def set_thickness(self, thickness):
        self.thickness = thickness
        self.recalculate()

    def recalculate(self):
        for x in self:
            x.recalculate(self.thickness)

    def table_csv(self, filename, length_units='cm'):
        fl = length_unit_factor(length_units)

        header1 = ['Temperature', 'Current']
        header2 = ['K', 'A']
        if self[0].rs is not None:
            header1.extend(['Sheet resistance', 'SR error', 'Resistivity', 'R error'])
            header2.extend(['\u2126', '\u2126', '\u2126' + length_units, '\u2126' + length_units])

        if self[0].rh is not None:
            header1.extend(['Reduced Hall coefficient', 'RHC error','Hall coefficient', 'HC error','Density', 'D error',
                   'Mobility', 'M error'])
            header2.extend([length_units + '\u00b2/C', length_units + '\u00b2/C', length_units + '\u00b3/C',
                   length_units + '\u00b3/C',
                   length_units + '\u207b\u00b3', length_units + '\u207b\u00b3', length_units + '\u00b2/Vs',
                   length_units + '\u00b2/Vs'])


        with open(filename, 'w', encoding='utf8', newline='') as f:
            w = csv.writer(f, delimiter=',')
            w.writerow(['Name', self.name])
            w.writerow(['Thickness', f'{self.thickness * 1e9:.0f} nm'])
            w.writerow(header1)
            w.writerow(header2)

            # print data
            for x in self:
                if isinstance(x.temp, str):
                    temp_str = 'RT'
                else:
                    temp_str = f'{x.temp:.1f}'
                if x.current is not None:
                    current_str = f'{x.current:.3e}'
                else:
                    current_str = '-'

                data = [temp_str, current_str]
                if x.rs is not None:
                    data.extend([f'{x.rs.n:.4e}', f'{x.rs.s:.1e}', f'{x.rho.n * fl:.4e}', f'{x.rho.s * fl:.1e}'])
                if x.rh is not None:
                    data.extend([f'{x.rh.n * fl ** 2:.4e}', f'{x.rh.s * fl ** 2:.1e}',
                        f'{x.r_hall.n * fl ** 3:.4e}', f'{x.r_hall.s * fl ** 3:.1e}',
                        f'{x.n.n / fl ** 3:.4e}', f'{x.n.s / fl ** 3:.1e}',
                        f'{x.mu.n * fl ** 2:.4e}', f'{x.mu.s * fl ** 2:.1e}'])
                    header1
                w.writerow(data)

    def report_txt(self, filename, length_units='cm'):
        with open(filename, 'w', encoding='utf8') as f:
            f.write(f'Name: {self.name}\n')
            f.write(f'Thickness: {self.thickness * 1e9:.0f} nm\n\n')
            for x in self:
                if isinstance(x.temp, str):
                    temp_str = 'RT'
                else:
                    temp_str = f'{x.temp:.1f}'
                if x.current is not None:
                    current_str = f'{x.current:.3e}'
                else:
                    current_str = '-'
                f.write(f'T = {temp_str}, I = {current_str} A\n')
                f.write('-------------------------\n')
                if x.vdp is not None:
                    f.write('    VDP:     R_12/43       R_34/21       R_23/14       R_41/32\n')
                    f.write('      ')
                    for rc in x.vdp.rc:
                        f.write(f'{rc:14.4e}')
                    f.write('\n')
                    ra_str = f'{x.vdp.ra}'
                    rb_str = f'{x.vdp.rb}'
                    f.write('                               R_A                         R_B\n')
                    f.write(' ' * (34 - len(ra_str)) + ra_str + ' ' * (28 - len(rb_str)) + rb_str + '\n')
                    f.write('                               R_S\n')
                    rs_str = f'{x.vdp.rs}'
                    f.write(' ' * (34 - len(rs_str)) + rs_str + '\n\n')

                if x.hall is not None:
                    f.write('    Hall: R_13/24(+)    R_13/24(-)    R_24/31(+)    R_24/31(-) \n')
                    f.write('      ')
                    for rc in x.hall.rc:
                        f.write(f'{rc:14.4e}')
                    f.write('\n')
                    f.write('                          dR_13/24                    dR_24/31\n')
                    f.write('      ')
                    for dr in x.hall.dr_v:
                        f.write(f'{dr:28.4e}')
                    f.write('\n')
                    f.write('                                B+                          B-\n')
                    f.write('      ')
                    for b in x.hall.b:
                        f.write(f'{b * 1e-4:28.4f}')
                    f.write('\n')
                    f.write('                               r_H\n')
                    rh_str = f'{x.hall.rh}'
                    f.write(' ' * (34 - len(rh_str)) + rh_str + '\n\n')
