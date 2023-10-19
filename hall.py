import mmr_data


def print_report(filename, data):
    


def main():
    data = mmr_data.load_csv_file('testing/CuO-03-11 October 2022-1506.csv')
    data = mmr_data.process_data(data)

    print(data['rs'])
    print(data['rh'])




if __name__ == "__main__":
    main()