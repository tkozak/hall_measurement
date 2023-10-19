import mmr_data


def main():
    data = mmr_data.load_csv_file('testing/CuO-03-11 October 2022-1506.csv')
    print(data)
    print(data[0])
    print(data[0].vdp.rs)

    print(data[0].rho)
    print(data[0].r_hall)
    print(data[0].n)
    print(data[0].mu)


if __name__ == "__main__":
    main()