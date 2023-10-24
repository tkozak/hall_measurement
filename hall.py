import mmr_data


def main():
    data = mmr_data.load_csv_file('testing/CuO-03-11 October 2022-1506.csv')
    print(data)

    x = data.collect_by_variable('temp')
    x.recalculate()
    print(x[0].rho)
    print(x[0].mu)


if __name__ == "__main__":
    main()