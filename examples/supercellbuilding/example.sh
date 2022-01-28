# Default 3x3x3 lattice
make_supercell -i 1101307.cif -o 1101307.pdb
# 5x5x5 lattice
make_supercell -i 1101307.cif -o 1101307.pdb --a_min_max -2 2 --b_min_max -2 2 --c_min_max -2 2

