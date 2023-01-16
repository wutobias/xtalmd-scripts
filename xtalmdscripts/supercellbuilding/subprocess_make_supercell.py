#from subprocess import check_output
import os
from subprocess import run

# This script will run the subprocess to control make_supercell.py and build supercell
# System with specifying size of supercell from previous downloaded CIF files
count = 0
for f in os.listdir('CIF'):
    print(f)
    run(['/home/yu-tang/anaconda3/envs/xtalmd/bin/python', './make_supercell.py', '-i', './CIF/' + f, '-pre',
         './supercell/' + f.split(".")[0] + "_supercell", '--a_min_max', '2', '--b_min_max',  '2', '--c_min_max', '2'])
    count += 1
print("Generate " + str(count) + " supercell files success")

