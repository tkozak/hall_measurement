import re
import numpy as np
from uncertainties import unumpy, umath

import theory


def sig_round(x, p): return float(f'{x:.{p - 1}e}')


def read_current_voltage_columns(lines):
    values = np.genfromtxt(lines, delimiter=',', usecols=(2, 3, 4, 5))
    current = values[[0, 2], :]
    voltage = values[[1, 3], :]
    current = current[:, [0, 2, 1, 3]]
    voltage = voltage[:, [0, 2, 1, 3]]
    return current, voltage


def load_csv_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    data = {'title': lines[1].split(',')[0], 'setpoint': [], 'vdp': [], 'hall': []}
    i = 2
    while i < len(lines):
        if re.match('Thickness', lines[i]):
            t_str = lines[i].split(',')[1]
            data['thickness'] = float(t_str[:-2]) * 1e-6
        if re.match('Point #', lines[i]) and re.match('Hall', lines[i - 1]):
            current, voltage = read_current_voltage_columns(lines[i + 1: i + 5])
            field_str_list = lines[i+6].split(',')
            field = np.array([float(field_str_list[2]), float(field_str_list[4])])
            data['hall'].append({'current': current,
                                 'voltage': voltage,
                                 'field': field})
            i += 12
        if re.search('Current 1', lines[i]) and not re.match('Hall', lines[i - 2]):
            # correct +- and ++ in line 2
            lines[i + 1] = re.sub('\+\+', '+', lines[i + 1])
            lines[i + 1] = re.sub('\+-', '-', lines[i + 1])
            current, voltage = read_current_voltage_columns(lines[i: i + 4])
            temp_str = lines[i].split(',')[6]
            try:
                temp_setpoint = float(temp_str)
            except ValueError:
                temp_setpoint = temp_str.strip()

            current_setpoint = sig_round(np.mean(-0.5*np.diff(current, axis=0)), 2)
            data['setpoint'].append({'current': current_setpoint,
                                     'temp': temp_setpoint})
            data['vdp'].append({'current': current,
                                'voltage': voltage})
            i += 8
        if re.match('AVERAGE', lines[i]):
            break
        i += 1
    return data


def process_data(data):
    data['vdp_dc_offset'] = []
    data['hall_dc_offset'] = []
    data['rs'] = []
    data['rh'] = []
    for k in range(0, len(data['setpoint'])):

        i = data['vdp'][k]['current']
        u = data['vdp'][k]['voltage']
        # r = u/i
        rc = np.squeeze(np.diff(u, axis=0)/np.diff(i, axis=0))   # resistance without dc offset
        data['vdp_dc_offset'].append(u[1, :] - rc * i[1, :]) # save dc offset
        # data['reversal_error'] = u[1, :] + u[2, :]
        rc = unumpy.uarray(rc, 0.02*rc)   # estimate error as 2% of the calculated resistance value
        ra = 0.5 * (rc[0] + rc[1])   # ra
        rb = 0.5 * (rc[2] + rc[3])
        rs = theory.sheet_resistance_vdp(ra, rb)
        data['rs'].append(rs)
        # rho = rs*data['thickness']*100  # convert to Ohm*cm
        # print(f'{ra:.2e}, {rb:.2e}, {rs:.2e}, {rho:.2e}')

        i = data['hall'][k]['current']
        u = data['hall'][k]['voltage']
        b = data['hall'][k]['field']
        rc = np.squeeze(np.diff(u, axis=0) / np.diff(i, axis=0))  # resistance without dc offset
        data['hall_dc_offset'].append(u[1, :] - rc * i[1, :])  # save dc offset

        rc = unumpy.uarray(rc, 0.0002*np.abs(rc))   # estimate error as 0.02% of the calculated resistance value (???)
        b = unumpy.uarray(b, 100)   # estimate error as 100 Gauss (??)
        dra = rc[0] - rc[1]
        drb = rc[2] - rc[3]
        dr = 0.5*(dra + drb)
        db = (b[0] - b[1]) * 1e-4   # change in magnetic field (T)
        rh = theory.hall_coefficient(dr, db)
        data['rh'].append(rh)

    return data
