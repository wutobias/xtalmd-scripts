## Supercell building

The script `make_supercell.py` can generate supercells for usage in PBC molecular dynamics simulations. The input structure is a cif file (for instance from the CSD database) and the output is a PDB file.

### Usage

To build a 3x3x3 supercell lattice type:

```make_supercell.py -i 1101307.cif -pre adipamide -a -1 1 -b -1 1 -c -1 1```

The `-a`, `-b` and `-c` arguments tell the program where it should start (first integer) and stop (second integer) with replicating the lattice along the a, b or c direction, respectively.

Type `make_supercell.py --help` for more help.

### Output

#### PDB file

The output `{prefix}.pdb` contains the atomic coordinates and box coordinates of the supercell. The residue numbers are increasing over the whole lattice. The residue name is consistent for identical molecules across the whole lattice. Residues are named as `M{X}`, where `X` goes from 0 to whatever number of unique molecules are on the lattice.

#### CSV file

Useful for post processing. This csv file `{prefix}.csv` contains one row for each molecule on the lattice. The first column is a unique int identifier for each molecule on the lattice, second column is a unique identifier for each molecule in each unit cell. Third to fifth columns are the fractional coordinates of the unit cell origin of the unit cell the molecule is in.

#### json file

Useful for post processing. This json file `{prefix}.json` can be used to reconstruct the rdkit object used for generating the supercell pdb file. Something like the following can be used to get the rdkit object.

```
from rdkit import Chem
with open("adipamide.json", "r") as fopen:
	mol=Chem.JSONToMols(fopen.read())
```
