## Table Formatting

The script `make_tables.py` reads data from md trajs, cif and others and
incorporates it into a single excel sheet.

## Usage

Run the script using

`make_tables.py -i input.json -o data.xlsx`

All the input data must be specified in the file `input.json`. The
formatting should be pretty self explanatory. Meaning of the input
variables:

- FORCEFIELD: This is a list containing data computed with all the different force fields.

	- XTAL_ENERGY: File containing output from  `gmx energy`  for the MD of the xtal lattice. File must be in xvg format. Will only read potential energy from this file.

	- XTAL_CELL: File containing output from  `gmx energy`  for the MD of the xtal lattice. File must be in xvg format. Will only read cell data from this file.

	- GAS_ENERGY: File containing output from  `gmx energy`  for the MD of the gas phase molecule. File must be in xvg format. Will only read potential energy from this file.

	- NAME: Name of the force field.

	- UNITCELLS_A: Number of unit cells in direction a.

	- UNITCELLS_B: Number of unit cells in direction b.

	- UNITCELLS_C: Number of unit cells in direction c.

	- MOLS_PER_CELL: Number of molecules per unit cell.

- EXPT: This is a list containing data from experiments

	- CIF: cif file

	- SUBLIMATION_ENTHALPY: The sublimation enthalpy

	- SUBLIMATION_ENTHALPY_ERROR: The sublimation enthalpy error

