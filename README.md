## xtalmd scripts

This is a collection of scripts for the xtal MD project in the Gilson Lab.

## Notes

### xyz2mol
This script makes heavy usage of `xyz2mol`, which is developed in the Jensen group and can be found [here](https://github.com/jensengroup/xyz2mol).

### Charmm Force Field
In order to be able to use the cgenff CHARMM force field, the programs `cgenff` and `charmm` must be in `PATH`.

### Oplsaa Force Field
The program `BOSS` must be installed and available in `PATH` in order to use this force field. Make sure that the environment variable `BOSSdir` is set and points to the directory of the BOSS installation. The `BOSS` program is available for free [here](http://zarbi.chem.yale.edu/software.html). After that install the program `ligpargen` from [here](https://github.com/Isra3l/ligpargen).

## Installation

### Clone the repository

`git clone https://github.com/wutobias/xtalmd-scripts`

### Setup conda environment

Setup a conda environment to get the dependencies right with `xtalmd.yml`. Just run

```
conda env create -f xtalmd.yml
```

If you want to run xtal minimization or xtal md on a cluster, better also install ray using pip
```
pip install ray['default']
```

### Install the package

`python setup.py install`

