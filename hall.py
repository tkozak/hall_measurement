import mmr_data


def main():
    data = mmr_data.load_csv_file('testing/CuO-03-11 October 2022-1506.csv')
    print(data)
    data.report_txt('testing/CuO_report.txt')

    x = data.collect_by_variable('temp')
    x.recalculate()
    x.table_csv('testing/CuO.csv', length_units='cm')



if __name__ == "__main__":
    main()