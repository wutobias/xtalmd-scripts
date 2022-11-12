#from subprocess import check_output
import os
from subprocess import run

count = 0
for f in os.listdir('CIF'):
    print(f)
    run(['/home/yu-tang/anaconda3/envs/xtalmd/bin/python', './make_supercell.py', '-i', './CIF/' + f, '-pre',
         './supercell/' + f.split(".")[0] + "_supercell", '--a_min_max', '2', '--b_min_max',  '2', '--c_min_max', '2'])
    count += 1
print("Generate " + str(count) + " supercell files success")

