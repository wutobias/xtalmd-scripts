import os


def cif_temperature():
    temp_data = dict()
    for f in os.listdir('../supercellbuilding/CIF/'):
        cif_file = open('../supercellbuilding/CIF/' + f, 'r')
        temp_data[str(f.split(".")[0])] = 298.15
        for line in cif_file.readlines():
            if line.startswith("_cell_measurement_temperature"):
                temp_data[str(f.split(".")[0])] = float(line.split()[1].split("(")[0])
    return temp_data

