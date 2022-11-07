# Examples

This is to illustrate some use cases of this package. The examples provided here are based on crystal structures of [adipamide](https:/en.wikipedia.org/wiki/Adipamide) and [oxamide](https:/en.wikipedia.org/wiki/Oxamide).

## Building supercells (command line)

Here we build a default 3x3x3 supercell lattice:

```
make_supercell -i data/adipamide.cif -pre make_supercell/adipamide_3x3x3
make_supercell -i data/oxamide.cif   -pre make_supercell/oxamide_3x3x3
```

Next, let's build a 5x5x5 lattice:

```
make_supercell -i data/adipamide.cif -pre make_supercell/adipamide_5x5x5 --a_min_max -2 2 --b_min_max -2 2 --c_min_max -2 2
make_supercell -i data/oxamide.cif   -pre make_supercell/oxamide_5x5x5 --a_min_max -2 2 --b_min_max -2 2 --c_min_max -2 2
```

The PDB files with atomic coordinates for each of the lattices can be found in the folder `make_supercell`  that we have built with the above commands. In addition, there are `.json` and `.csv` files for each lattice. The `.json` files can be parsed to rdkit to obtain an rdkit molecule object for the whole lattice. For instance:

```python
from rdkit import Chem
with open("make_supercell/oxamide_5x5x5.json", "r") as fopen:
    rdmol = Chem.JSONToMols(fopen.read())[0]
```

The `.csv` file contains human readable information on the unit cell positions of each molecule:

```
#Mol_idx,mol_in_unitcell,unitcell_in_supercell_a,unitcell_in_supercell_b,unitcell_in_supercell_c
0,0,-1,-1,-1
1,1,-1,-1,-1
2,0,-1,-1,0
3,1,-1,-1,0
4,0,-1,-1,1
...
```

The `Mol_idx` column iterates through all molecules in the lattice and assigns them a unique integer identifier. This identifier corresponds to the position of the molecules in the list of molecules we get from `Chem.GetMolFrags(rdmol, asMols=True)` using the above rdkit molecule parsed from `.json` file. The column `mol_in_unitcell` assigns each molecule in each unit cell a unique identifier. For instance, compare `make_supercell/adipamide_3x3x3.csv`  and `make_supercell/oxamide_3x3x3.csv` to see that Adipamide has two molecules per unitcell, whereas Oxamide has only one. The last three columns in the `.csv` file (`unitcell_in_supercell_a`, `unitcell_in_supercell_b`, `unitcell_in_supercell_c`) correspond to the fractional positions of the unitcell that a given molecule is in.

By default, `make_supercell` uses the `xyz2mol` for building molecule topologies. However, sometimes `xyz2mol` cannot build a lattice and we will have to use the OpenEye toolkit for this (note that OpenEye toolkit requires a valid license located at `$OE_LICENSE`). In order to do that, simply re-run the above commands with `--use_openeye`. For instance:

```
make_supercell -i data/adipamide.cif -pre make_supercell/adipamide_3x3x3_oe --use_openeye
make_supercell -i data/oxamide.cif -pre make_supercell/oxamide_3x3x3_oe --use_openeye
```

We can also use the OpenEye toolkit to add hydrogen atoms. But check carefully, the protonation states are not assigned according to the crystal environment but to an aq. solution environment.


```
make_supercell -i data/adipamide.cif -pre make_supercell/adipamide_3x3x3_ah --use_openeye --addhs
make_supercell -i data/oxamide.cif -pre make_supercell/oxamide_3x3x3_ah --use_openeye --addhs
```

If you want to adjust protonation states according to the crystal environment, use the `--n_protonation_attempts` option to tell the program how many attempts it should take on finding protonation states (defaults to `0`). Be careful, this is currently not very efficient, not well tested and results might be wrong. Also note that it will take some time to run.

```
make_supercell -i data/adipamide.cif -pre make_supercell/adipamide_3x3x3_p1000 --n_protonation_attempts 1000
make_supercell -i data/oxamide.cif -pre make_supercell/oxamide_3x3x3_p1000 --n_protonation_attempts 1000
```

Sometimes CIF files contain an unusal space group that can't be interpreted by the underlying `gemmi` library which builds the initial unitcell. In these cases the resulting lattice will be incorrect. However, usually the symmetry operations that are explicitely provided in the CIF file are correct. In that case, activate `--use_symmetry_operations` to build the lattice directly from the symmetry operations instead of the spacegroup (which also would be a mapping to the set of symmetry operations) in the CIF file directly. This is usually the safest way to build a supercell lattice.

```
make_supercell -i data/adipamide.cif -pre make_supercell/adipamide_3x3x3_op --use_symmetry_operations
```

Some crystals need have bound water molecules that are resolved in the crystal structure. With the option `--addwater`, one can add a fixed number of water molecules to the unit cell. Default is `--addwater 0` (do not add any water molecules).

## Building supercells (python)

A short example for building a 3x3x3 supercell in python:

```python
from xtalmdscripts.supercellbuilding import make_supercell

### Set a,b,c limits for 3x3x3
a_min_max = [0,2]
b_min_max = [0,2]
c_min_max = [0,2]

### Load the cif file as gemmi structure
strc = make_supercell.parse_cif("data/adipamide.cif")
### `replicated_mol_list` is a list of rdkit molecule objects.
### One for each molecule in the lattice.
replicated_mol_list, _, _ = make_supercell.generate_replicated_mol_list(
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
with open("./make_supercell/adipamide_3x3x3_py.pdb", "w") as fopen:
    fopen.write(pdb_str)
```

## Building Molecular Mechanics Atomic-Level Systems

We can build atomic-level systems using different force fields: GAFF 1.X, GAFF 2.X, OpenFF 1.X, OpenFF 2.X, CgenFF, OPLS. For that we will use the script `build_system`, which inherits all the options from the previous `make_supercell`. The script `build_system` can build the lattice structure as a pdb file and generate the accompanying force field as an OpenMM xml file.

Note that for the CgenFF force field to work, you must have installed the Charmm/CgenFF force field files stored under `toppar_c36_jul21` as well as the programs `CgenFF` and `CHARMM`. For the opls-AA forcefield to work, you must have installed the programs `LigParGen` and `BOSS`. See [here](../README.md) for instructions on how to get these.

```
build_system -i data/oxamide.cif -pre build_system/oxamide_gaff1   -ff gaff1
build_system -i data/oxamide.cif -pre build_system/oxamide_parsley -ff parsley

build_system -i data/adipamide.cif -pre build_system/adipamide_gaff1   -ff gaff1
build_system -i data/adipamide.cif -pre build_system/adipamide_parsley -ff parsley
```

In constrast to `make_supercell`, we do not want to specify a fixed number of replicates for each unit cell dimension. Instead we define the nonbonded cutoff `--nbcutoff` and the `--axislengthfactor`, which tells the program how long in multiples of the nonbonded cutoff, each crystallographic axis should be.

Note that by default, we save force field `.xml` files not only for the supercell lattice, but also for each unique monomer molecule in the lattice (`XXX_monomer0.xml`). In addition, we have a `.json` file that contains the supercell lattice as big rdkit molecule. It can be parsed to python again using:

```python
from rdki import Chem
with open("build_system/oxamide_cgenff.json", "r") as fopen:
    rdmol = Chem.JSONToMols(fopen.read())[0]
```

## Running Lattice Minimization

After the building the lattice, we should first minimize the lattice. This can be done using the script `run_xtal_min`. In principal, it has two run modes: `run_xtal_min xml` for normal command line based input or `run_xtal_min yaml` for high-throughput runs. We will start with the command line runs:

First, let's optimize the joint coordinates of atomic positions and lattice vectors using 10 steps of L-BFGS-B. Note that in reality you want to run this longer than 10 steps.

```
run_xtal_min xml -i build_system/oxamide_gaff1.xml -p build_system/oxamide_gaff1.pdb -pre xtal_min/oxamide_gaff1_min --method L-BFGS-B --steps 10
```

Now look at the minimized structure `xtal_min/oxamide_gaff1_min.pdb` and compare it with the initial structure from `build_system/oxamide/oxamide_gaff1.pdb`. What has changed?

As an alternative to the above minimization method, we can optimize the atom positions and box vectors seperately and alternate between the two.

```
run_xtal_min xml -i build_system/oxamide_gaff1.xml -p build_system/oxamide_gaff1.pdb -pre xtal_min/oxamide_gaff1_min_alt --method L-BFGS-B --steps 10 --alternating
```

Normally, the minimization operates on the 6 degrees of freedom in unit cell matrix. However, one can also use axes lengths and angles for the minization variables. In that case, add the `--use_lengths_and_angles` option to the command above.

In general, the option `--method` can be used to define the main optimization algorithm in place. Valid choices are Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA, SLSQP, trust-const, dogleg, trust-ncg, trust-exact, trust-krylov. See the [scipy documentation](https:/docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html) for more information.

Now let's do the same thing for a larger number of lattice structures. The commands are all stored in the [`input.yaml`](xtal_min/input.yaml) Go to the folder `xtal_min` first, and then run the following command. But be careful, this will start a lot of processes. If you want to limit the number of available CPUs, first start the ray server manually `ray start --head --port=6831 --redis-password=5241590000000000 --num-cpus=16 --num-gpus=4` and add the line `"ray_host" : "localhost:6381"` to the `input.yaml` file.

```
run_xtal_min yaml -i input.yaml
```

## Running Lattice Simulations

To run MD simulations of the supercell lattice, we use `run_xtal_md`. In accordance with the minimization script `run_xtal_min`, we can run this in two modes, either `run_xtal_md xml` from the command line or `run_xtal_md yaml` as high-throughput.

The command line arguments are very self-explaining. Just run `run_xtal_md xml -h` to see what is available. For instance, to run a 2 nanosecond lattice simulation of oxamide starting from the minimized structure:

```
run_xtal_md xml -i build_system/oxamide_gaff1.xml -p xtal_min/oxamide_gaff1_min.pdb --temperature=298.15 --nanoseconds 2 --replicates 2 --prefix oxamide_gaff1
```

For the high-throughput work, have a look at the [`input.yaml`](xtal_min/input.yaml). This works sort of similarly as for the previous `run_xtal_min yaml` case.

## Data Analysis

The data analysis interface works only through `.yaml` files. This is necessary, since we need to combine lots of different information, which is easier in a structured file format than on the command line. Have a look into the file `analysis/input.yaml` to see what information goes where. The variables should be self-explaining. The final analysis is summarized and saved in an Excel `.xlsx` table. Run the analaysis using:

```
make_tables -i analysis/input.yaml -o table.xlsx
```