## Supercell building

The script `make_supercell.py` can generate supercells for usage in PBC molecular dynamics simulations. The input structure is a cif file (for instance from the CSD database) and the output is a PDB file.

### Usage

To build a 3x3x3 supercell lattice type:

```make_supercell.py -i 1101307.cif -o 1101307.pdb -a -1 1 -b -1 1 -c -1 1```

The `-a`, `-b` and `-c` arguments tell the program where it should start (first integer) and stop (second integer) with replicating the lattice along the a, b or c direction, respectively.

Type `make_supercell.py --help` for more help.

