import mmr_data




def main():
    data = mmr_data.load_csv_file('testing/CuO-03-11 October 2022-1506.csv')
    data = mmr_data.process_data(data)

    print(data['rs'])
    print(data['rh'])

    # data = mmr_data.load_csv_file('testing/HfSiBCN_Y3-09 February 2017-666.csv')


if __name__ == "__main__":
    main()