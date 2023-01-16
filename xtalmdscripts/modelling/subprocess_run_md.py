import os
from subprocess import run
from cif_temperature import cif_temperature
temp_data = cif_temperature()

# This script will run the subprocess to run run_md.py to run MD simulation with used-defined parameter
count = 0
for f in os.listdir('../build_system/MM/'):
    if f.endswith("parsley_1.0.0.xml") or f.endswith("parsley_1.1.0.xml") or f.endswith("parsley_1.2.0.xml") or f.endswith("parsley_1.3.0.xml") or f.endswith("sage_2.0.0.xml"):
        print(f)
        run(['/home/yu-tang/anaconda3/envs/xtalmd/bin/python',
             './run_md.py', 'xml',
             '-i', '../build_system/MM/' + f,
             '-p', './MIN/' + f.strip(".xml") + '_min.pdb',
             '--temperature=' + str(temp_data[f.split("_")[0]]),
             '--nanoseconds', '2',
             '--replicates', '1',
             '--prefix', f.strip(".xml") + '_md'
             ])
        count += 1

print("Run " + str(count) + " md parsley/sage files success")
