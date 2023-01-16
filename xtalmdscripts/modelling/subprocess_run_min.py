import os
from subprocess import run, Popen
import csv
import time


# This script will run the subprocess Popen to control run_min.py to run minimization in a batch
# Write the minimization result to CSV files
with open('./minimization_results.csv', 'w', newline='') as outcsv:
    writer = csv.DictWriter(outcsv, fieldnames=['COD ID', 'RMSD', 'RMSD (Remove Periodicity)', 'RMSD20', 'RMSD20 MAX',
                                                'RMSD20 MIN', 'RMSD20 Mean', 'method', 'Original Energy',
                                                'Minimized Energy (OpenMM)', 'Minimized Energy (Cell Minimization)',
                                                'Original Box Vectors', 'Minimized Box Vectors'])
    writer.writeheader()


# Run minimization in a batch parallelly with memory control
count = 0
for f in os.listdir('../build_system/MM/'):
    if f.endswith("parsley_1.0.0.xml") or f.endswith("parsley_1.1.0.xml") or f.endswith(
            "parsley_1.2.0.xml") or f.endswith("parsley_1.3.0.xml") or f.endswith("sage_2.0.0.xml"):
        if not os.path.exists('./MIN/' + f.strip(".xml") + '_min.pdb'):
            Popen(['/home/yu-tang/anaconda3/envs/xtalmd/bin/python',
                   './run_min.py', 'xml',
                   '-i', '../build_system/MM/' + f,
                   '-p', '../build_system/MM/' + f.strip(".xml") + '.pdb',
                   '-pre', '/MIN/' + f.strip(".xml") + '_min',
                   '--method', 'L-BFGS-B'
                   ])
            count += 1
            # Memory control
            if ((count + 1) % 10) == 0:
                time.sleep(300)
print("Run " + str(count) + " minimization files success")

