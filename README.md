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

## Instruction 
1. run `python xtalmdscripts/supercellbuilding/COD_import.py`

Download list of desired structures from COD using list of COD IDs in `COD_ID_List.txt` and `COD_import.py` script. This downloads CIF files from COD and converts to PDB format using pybel. These are output in the `xtalmdscripts/supercellbuilding/CIF` and `xtalmdscripts/supercellbuilding` directories.

2. run `python xtalmdscripts/supercellbuilding/make_supercell.py`

Create supercell from CIF file. This is done in the `xtalmdscripts/supercellbuilding/make_supercell.py` script and output in the `xtalmdscripts/supercellbuilding/supercell` directory. This step can be skipped because step3 would build the supercell and molecular mechanics atomic-level systems at the same time.

3. run `python xtalmdscripts/build_system/subprocess_build_system.py`

Build supercell and molecular mechanics atomic-level systems at the same time from CIF file via subprocess to control `xtalmdscripts/build_system/build_system.py`. This cell is created to be large enought to satisfy that periodic boundary conditions are greater than 2.0 nm. This is done in the `xtalmdscripts/build_system/subprocess_build_system.py` script and output in the `xtalmdscripts/build_system/MM` directory. 

4. run `python xtalmdscripts/modelling/subprocess_run_min.py`

Paramaterize using openFF a in openMM via subprocess to control `xtalmdscripts/modelling/subprocess_run_min.py`. Perform energy minimization by 'L-BFGS-B' or 'trust-constr' (with period boundary constraint) and record the initial and final energy values as well as the RMSD20 between initial and final state. The energy v.s. iterations for each minimization would also be plotted after minimization. This is done in the `xtalmdscripts/modelling/subprocess_run_min.py` script and output in the `xtalmdscripts/modelling/MIN` and `xtalmdscripts/modelling/minimization_results.csv`. 

5. run `python xtalmdscripts/modelling/subprocess_run_md.py`

Perform MD simulation by previous energy minimized structure for 2 nanoseconds and calculate centroids, RMSDF and B-factor after MD simulation. This is done in the `xtalmdscripts/modelling/subprocess_run_md.py` script and output in the `xtalmdscripts/modelling/run_0`. 

## Examples for fine-tuning each parameter

A bunch of examples for the different things this package can do, have a look at the [examples](examples/README.md).


## Current Issues 

`Independent study report.pdf` have already listed the current progresses and issues that need to be solved

1. Figure out different bonds can contribute to different energy for different forcefields and
why different forcefield exists so different minimized energy with a relatively similar
RMSD20.

2. Further B-factor data analysis to estimate the performance of OpenFF is by observing how
atoms move in MD simulation. This can help us estimate the performance of OpenFF in the
crystal structure prediction.

## Authors

Tobias Hüfner (UC San Diego)

