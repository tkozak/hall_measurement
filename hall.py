import mmr_data
import argparse
import glob


def main():

    parser = argparse.ArgumentParser(description='Process MMR hall measurements')
    parser.add_argument('output_name', type=str, help='core name of output files')
    parser.add_argument('files', type=str, nargs='+', help='path to input files')
    parser.add_argument('-u', '--unit', type=str, default='cm', help='length units for output')
    parser.add_argument('-t', '--temperature', help='collect data points by temperature',
                        action="store_true")
    parser.add_argument('-i', '--current', help='collect data points by current',
                        action="store_true")

    args = parser.parse_args()

    files = []
    for p in args.files:
        files.extend(glob.glob(p))
    print('List of loaded files:')
    for i, f in enumerate(files):
        print(f'  {f}')

    all_data = None
    for f in files:
        data = mmr_data.load_csv_file(f)
        if all_data is None:
            all_data = data
        else:
            all_data.extend(data)

    all_data.report_txt(args.output_name + '_report.txt')

    data = all_data
    if args.temperature:
        data = data.collect_by_variable('temp')
    if args.current:
        data = data.collect_by_variable('current')

    data.recalculate()
    data.table_csv(args.output_name + '.csv', length_units=args.unit)


if __name__ == "__main__":
    main()