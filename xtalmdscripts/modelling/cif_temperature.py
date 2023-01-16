import os
import pandas as pd

# This script will extract the temperature for MD simulation from database (CIF files)
def cif_temperature():
    temp_data = dict()
    for f in os.listdir('../supercellbuilding/CIF/'):
        cif_file = open('../supercellbuilding/CIF/' + f, 'r')
        temp_data[str(f.split(".")[0])] = "None"
        for line in cif_file.readlines():
            if line.startswith("_cell_measurement_temperature"):
                temp_data[str(f.split(".")[0])] = float(line.split()[1].split("(")[0])
    return temp_data


# This function will extract the temperature for database temperature distribution
def cif_temperature():
    file_name = []
    compound = []
    temperature = []
    year = []
    i = 0
    for f in os.listdir('../supercellbuilding/CIF/'):
        cif_file = open('../supercellbuilding/CIF/' + f, 'r')
        file_name.append(str(f.split(".")[0]))
        temperature.append(None)
        compound.append(None)
        year.append(None)
        for line in cif_file.readlines():
            if line.startswith("_chemical_formula_sum"):
                compound[i] = line.split("'")[1]
            if line.startswith("_cell_measurement_temperature"):
                temperature[i] = float(line.split()[1].split("(")[0])
            if line.startswith("_journal_year"):
                year[i] = line.split(" ")[-1].strip("\n")
        i += 1
    s1 = pd.Series(file_name)
    s2 = pd.Series(compound)
    s3 = pd.Series(temperature)
    s4 = pd.Series(year)
    data = pd.concat([s1, s2, s3, s4], names=['File name', 'Compound', 'Temperature(K)', 'Year'], axis=1)
    return data

data = cif_temperature()
df = pd.DataFrame(data)
print(df.describe())
f = open('./temperature.csv', 'w')  # Open file as append mode
df.to_csv(f, index=False, header=['File name', 'Compound', 'Temperature(K)', 'Year'])
f.close()
