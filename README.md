# xtalmd scripts

This is a package for building atomic-level systems supercell lattices, energy minimizing supercell lattices, running MD simulations of supercell lattices and analysis of the simulations.

## Installation

### Clone the repository

`git clone https://github.com/wutobias/xtalmd-scripts`

### Setup conda environment

Setup a conda environment to get the dependencies right with `xtalmd.yml`:

```
conda env create -f xtalmd.yml
```

Activate it:
```
conda activate xtalmd
```

If you want to run energy minimizations or lattice simulations, also install the package `ray`:
```
pip install ray['default']
```

### Install the package

`python setup.py install`


## Dependencies

Almost all dependencies should be covered in the `xtalmd.yml` file. All other dependencies are listed below.

### xyz2mol
We make heavy usage of `xyz2mol`, which is developed in the Jensen group and can be found [here](https://github.com/jensengroup/xyz2mol).

### Charmm Force Field
In order to be able to use the cgenff CHARMM force field, the programs [`cgenff`](https://silcsbio.com/) and [`charmm`](https://academiccharmm.org/) must be in `PATH`. We recommend using `cgenff` version 2.5.1 and `charmm` version 46b1, both can be obtained free of charge for academic users.
The force field files for CgenFF should be obtained [here](https://www.charmm.org/archive/charmm/resources/charmm-force-fields/#charmm).

### Oplsaa Force Field
The program `BOSS` must be installed and available in `PATH` in order to use this force field. Make sure that the environment variable `BOSSdir` is set and points to the directory of the BOSS installation. The `BOSS` program is available for free [here](http://zarbi.chem.yale.edu/software.html). After that install the program `ligpargen` from [here](https://github.com/Isra3l/ligpargen).

### OpenEye Toolkits
Some of the optional system building routines can only be used together with OpenEye toolkit, which requires a license.

## Examples

A bunch of examples for the different things this package can do, have a look at the [examples](examples/README.md).

## Authors

Tobias Hüfner (UC San Diego)