import os
from subprocess import run

# This script will run the subprocess to control build_system.py and build Molecular Mechanics Atomic-Level Systems
# System with specifying size of supercell and openforcefield version from previous downloaded CIF files
# The size of supercell is 0.9*2.2 ~ 2nm ~ 20A
count = 0
for f in os.listdir('../supercellbuilding/CIF/'):
    print(f)
    run(['/home/yu-tang/anaconda3/envs/xtalmd/bin/python', './build_system.py',
         '-i', '../supercellbuilding/CIF/' + f,
         '-pre', './MM/' + f.split(".")[0] + "_parsley_1.0.0",
         '-ff', 'parsley',
         '--version', '1.0.0',
         '--nbcutoff', '0.9',
         '--axislengthfactor', '2.2'])
    run(['/home/yu-tang/anaconda3/envs/xtalmd/bin/python', './build_system.py',
         '-i', '../supercellbuilding/CIF/' + f,
         '-pre', './MM/' + f.split(".")[0] + "_parsley_1.1.0",
         '-ff', 'parsley',
         '--version', '1.1.0',
         '--nbcutoff', '0.9',
         '--axislengthfactor', '2.2'])
    run(['/home/yu-tang/anaconda3/envs/xtalmd/bin/python', './build_system.py',
         '-i', '../supercellbuilding/CIF/' + f,
         '-pre', './MM/' + f.split(".")[0] + "_parsley_1.2.0",
         '-ff', 'parsley',
         '--version', '1.2.0',
         '--nbcutoff', '0.9',
         '--axislengthfactor', '2.2'])
    run(['/home/yu-tang/anaconda3/envs/xtalmd/bin/python', './build_system.py',
         '-i', '../supercellbuilding/CIF/' + f,
         '-pre', './MM/' + f.split(".")[0] + "_parsley_1.3.0",
         '-ff', 'parsley',
         '--version', '1.3.0',
         '--nbcutoff', '0.9',
         '--axislengthfactor', '2.2'])
    run(['/home/yu-tang/anaconda3/envs/xtalmd/bin/python', './build_system.py',
         '-i', '../supercellbuilding/CIF/' + f,
         '-pre', './MM/' + f.split(".")[0] + "_sage_2.0.0",
         '-ff', 'sage',
         '--version', '2.0.0',
         '--nbcutoff', '0.9',
         '--axislengthfactor', '2.2'])
    count += 1
print("Generate " + str(count) + " molecular mechanics parsley 1.0.0 / parsley 1.1.0 / parsley 1.2.0 / parsley 1.3.0 / sage parsley 2.0.0 files success")
