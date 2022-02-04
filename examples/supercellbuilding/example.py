#!/usr/bin/env python3

from xtalmdscripts.supercellbuilding import make_supercell

### Set a,b,c limits for 3x3x3
a_min_max = [0,2]
b_min_max = [0,2]
c_min_max = [0,2]

### Load the cif file as gemmi structure
strc = make_supercell.parse_cif("../data/adipamide/1101307.cif")
replicated_mol_list = make_supercell.generate_replicated_mol_list(
    strc,
    a_min_max,
    b_min_max,
    c_min_max,
    )
pdb_str = make_supercell.get_pdb_str(
    replicated_mol_list,
    strc,
    a_min_max,
    b_min_max,
    c_min_max
    )
with open("./strc.pdb", "w") as fopen:
    fopen.write(pdb_str)
