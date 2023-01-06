import re
import numpy as np

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
    for k in range(0, len(data['setpoint'])):
        i = data['vdp'][k]['current']
        u = data['vdp'][k]['voltage']
        r = u/i
        rc = np.squeeze(np.diff(u, axis=0)/np.diff(i, axis=0))   # resistance without dc offset
        data['dc_offset'] = u[1, :] - rc * i[1, :]  # save dc offset
        rc_s = 0.02*rc   # estimate error as 2% of the calculated resistance value
        ra = 0.5 * (rc[0] + rc[1])   # ra
        ra_s = 0.5 * np.sqrt(rc_s[0]**2 + rc_s[1]**2)
        rb = 0.5 * (rc[2] + rc[3])
        rb_s = 0.5 * np.sqrt(rc_s[2]**2 + rc_s[3]**2)
        rs, rs_s = theory.sheet_resistance_vdp(ra, rb, ra_s, rb_s)
        rho = rs*data['thickness']*100
        print(f'{ra:.2e}, {rb:.2e}, {rho:.2e}')

    return data
