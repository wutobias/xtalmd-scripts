## xtalmd scripts

This is a collection of scripts for the xtal MD project in the Gilson Lab.

## Notes

### xyz2mol
This script makes heavy usage of `xyz2mol`, which is developed in the Jensen group and can be found [here](https://github.com/jensengroup/xyz2mol).

### Charmm Force Field
In order to be able to use the cgenff CHARMM force field, the programs `cgenff` and `charmm` must be in `PATH`.

### Oplsaa Force Field
In order to be able to use the OplsAA force field, there must be conda environment called `ligpargen` that has the program `LigParGen` installed. Furthermore, the program `BOSS` must be installed and available in `PATH`.
Instructions on how to install `LigParGen` and get a `BOSS` license can be found [here](http://zarbi.chem.yale.edu/ligpargen/). A yml file to setup the `ligpargen` conda environment is in the top level folder of this GitHub repo.

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

