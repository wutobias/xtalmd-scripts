import os
from subprocess import run, Popen
import csv
count = 0

with open('./minimization_results.csv', 'w', newline='') as outcsv:
    writer = csv.DictWriter(outcsv, fieldnames = ['COD ID', 'RMSD','RMSD (Remove Periodicity)','RMSD20','RMSD20 MAX','RMSD20 MIN','RMSD20 Mean','method','Original Energy',
        'Minimized Energy (OpenMM)', 'Minimized Energy (Cell Minimization)','Original Box Vectors','Minimized Box Vectors'])
    writer.writeheader()


for f in os.listdir('../build_system/MM/'):
    if f.endswith("parsley_1.0.0.xml") or f.endswith("parsley_1.1.0.xml") or f.endswith("parsley_1.2.0.xml") or f.endswith("parsley_1.3.0.xml") or f.endswith("sage_2.0.0.xml"):
        print(f)

        run(['/home/yu-tang/anaconda3/envs/xtalmd/bin/python',
             './run_min.py', 'xml',
             '-i', '../build_system/MM/' + f,
             '-p', '../build_system/MM/' + f.strip(".xml") + '.pdb',
             '-pre', '/MIN/' + f.strip(".xml") + '_min',
             '--method', 'L-BFGS-B'
             ])

        count += 1



print("Run " + str(count) + " minimization parsley/sage files success")
