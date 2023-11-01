import re
import numpy as np
import hall_data


def sig_round(x, p): return float(f'{x: .{p - 1}e}')


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
    title = lines[1].split(',')[0]
    sp = []
    vd = []
    hd = []
    i = 2
    th = 0.
    while i < len(lines):
        if re.match('Thickness', lines[i]):
            t_str = lines[i].split(',')[1]
            try:
                th = float(t_str[:-2]) * 1e-6  # thickness
            except ValueError:  # possible empty line with no value - do nothing
                pass
        if re.match('Point #', lines[i]) and re.match('Hall', lines[i - 1]):
            current, voltage = read_current_voltage_columns(lines[i + 1: i + 5])
            field_str_list = lines[i+6].split(',')
            field = np.array([float(field_str_list[2]), float(field_str_list[4])])
            hd.append({'current': current, 'voltage': voltage, 'field': field})
            i += 12
        if re.search('Current 1', lines[i]) and not re.match('Hall', lines[i - 2]):
            # correct +- and ++ in line 2
            lines[i + 1] = re.sub(r'\+\+', '+', lines[i + 1])
            lines[i + 1] = re.sub(r'\+-', '-', lines[i + 1])
            current, voltage = read_current_voltage_columns(lines[i: i + 4])
            temp_str = lines[i].split(',')[6]
            try:
                temp_sp = float(temp_str)
            except ValueError:
                temp_sp = temp_str.strip()

            current_sp = sig_round(np.mean(-0.5*np.diff(current, axis=0)), 2)
            sp.append({'current': current_sp, 'temp': temp_sp})
            vd.append({'current': current, 'voltage': voltage})
            i += 8
        if re.match('AVERAGE', lines[i]):
            break
        i += 1

    # process data
    data_list = hall_data.DataList(title, th)
    for k, set_point in enumerate(sp):
        data_point = hall_data.DataPoint(set_point['current'], set_point['temp'])
        data_point.set_data(hall_data.VdpData(vd[k]['current'], vd[k]['voltage']),
                            hall_data.HallData(hd[k]['current'], hd[k]['voltage'], hd[k]['field']))
        data_list.append(data_point)

    return data_list
