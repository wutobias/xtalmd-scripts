# Default 3x3x3 lattice
make_supercell -i ../data/adipamide/1101307.cif \
               -pre 1101307_3x3x3
make_supercell -i ../data/oxamide/1473522.cif   \
               -pre 1473522_3x3x3
# 5x5x5 lattice
make_supercell -i ../data/adipamide/1101307.cif \
               -pre 1101307_5x5x5 \
               --a_min_max -2 2 \
               --b_min_max -2 2 \
               --c_min_max -2 2
make_supercell -i ../data/oxamide/1473522.cif \
               -pre 1473522_5x5x5 \
               --a_min_max -2 2 \
               --b_min_max -2 2 \
               --c_min_max -2 2
