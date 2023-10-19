import re
import numpy as np
import hall_data


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
    title = lines[1].split(',')[0]
    sp = []
    vd = []
    hd = []
    i = 2
    while i < len(lines):
        if re.match('Thickness', lines[i]):
            t_str = lines[i].split(',')[1]
            th = float(t_str[:-2]) * 1e-6  # thickness
        if re.match('Point #', lines[i]) and re.match('Hall', lines[i - 1]):
            current, voltage = read_current_voltage_columns(lines[i + 1: i + 5])
            field_str_list = lines[i+6].split(',')
            field = np.array([float(field_str_list[2]), float(field_str_list[4])])
            hd.append({'current': current,
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
            sp.append({'current': current_setpoint, 'temp': temp_setpoint})
            vd.append({'current': current, 'voltage': voltage})
            i += 8
        if re.match('AVERAGE', lines[i]):
            break
        i += 1

    # process data
    for k in range(0, len(sp)):
        data = hall_data.VdpData(vd[k]['current'], vd[k]['voltage'])
        print(data.rs)

        data = hall_data.HallData(hd[k]['current'], hd[k]['voltage'], hd[k]['field'])
        print(data.rh)

        #TODO put all together

    return data

