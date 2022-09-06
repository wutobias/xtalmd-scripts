# Examples

This is to illustrate some use cases of this package. The examples provided here are based on crystal structures of [adipamide](https://en.wikipedia.org/wiki/Adipamide) and [oxamide](https://en.wikipedia.org/wiki/Oxamide).

## Building supercells (command line)

Here we build a default 3x3x3 supercell lattice:

```
make_supercell -i data/adipamide/1101307.cif -pre buildsystem/1101307_3x3x3
make_supercell -i data/oxamide/1473522.cif -pre buildsystem/1473522_3x3x3
```

Next, let's build a 5x5x5 lattice:

```
make_supercell -i data/adipamide/1101307.cif -pre supercellbuilding/1101307_5x5x5 --a_min_max -2 2 --b_min_max -2 2 --c_min_max -2 2
make_supercell -i data/oxamide/1473522.cif -pre supercellbuilding/1473522_5x5x5 --a_min_max -2 2 --b_min_max -2 2 --c_min_max -2 2
```

By default, `make_supercell` uses the `xyz2mol` for building molecule topologies. However, sometimes that doesn't work and we will have to use the OpenEye toolkit for this. In order to do that, simply re-run the above commands with `--use_openeye`. For instance:

```
make_supercell -i data/adipamide/1101307.cif -pre supercellbuilding/1101307_3x3x3_oe --use_openeye
make_supercell -i data/oxamide/1473522.cif -pre supercellbuilding/1473522_3x3x3_oe --use_openeye
```

We can also use OpenEye toolkits to add hydrogen atoms. But check carefully, the protonation states are not assigned according to the crystal environment but to an aq. solution environment.

```
make_supercell -i data/adipamide/1101307.cif -pre supercellbuilding/1101307_3x3x3_ah --use_openeye --addhs
make_supercell -i data/oxamide/1473522.cif -pre supercellbuilding/1473522_3x3x3_ah --use_openeye --addhs
```

If you want to adjust protonation states according to the crystal environment, use the `--n_protonation_attempts` option to tell the program how many attempts it should take on finding protonation states (defaults to `0`). Be careful, this is currently not very efficient and might take some time to run.

```
make_supercell -i data/adipamide/1101307.cif -pre supercellbuilding/1101307_3x3x3_p5 --n_protonation_attempts 5
make_supercell -i data/oxamide/1473522.cif -pre supercellbuilding/1473522_3x3x3_p5 --n_protonation_attempts 5
```

Sometimes CIF files contain an unusal space group. In these cases the resulting lattice will be incorrect. However, usually the symmetry operations that are explicitely given in the CIF file are correct. In that case, activate `--use_symmetry_operations` to build the lattice from the symmetry operations in the CIF file directly. This is usually the safest way to build a supercell lattice.

```
make_supercell -i data/adipamide/1101307.cif -pre supercellbuilding/1101307_3x3x3_op --use_symmetry_operations
make_supercell -i data/oxamide/1473522.cif -pre supercellbuilding/1473522_3x3x3_op --use_symmetry_operations
```

Some crystals need have bound water molecules that are resolved in the crystal structure. With the option `--addwater`, one can add a fixed number of water molecules to the unit cell. Default is `--addwater 0` (do not add any water molecules).

## Building supercells (python)

A short example in python can be found [here](supercellbuilding/example.py)

## Building Molecular Mechanics Atomic-Level Systems

We can build atomic-level systems using different force fields: GAFF 1.X, GAFF 2.X, OpenFF 1.X, OpenFF 2.X, CgenFF, OPLS. For that we will use the script `build_system`, which inherits all the options from the previous `make_supercell`. The script `build_system` can build the lattice structure as a pdb file and generate the accompanying force field as an OpenMM xml file.

Note that for the CgenFF force field to work, you must have the Charmm/CgenFF force field files installed.

```
build_system -i data/oxamide/1473522.cif -pre build_system/oxamide/1473522_gaff1   -ff gaff1
build_system -i data/oxamide/1473522.cif -pre build_system/oxamide/1473522_gaff2   -ff gaff2
build_system -i data/oxamide/1473522.cif -pre build_system/oxamide/1473522_parsley -ff parsley
build_system -i data/oxamide/1473522.cif -pre build_system/oxamide/1473522_sage    -ff sage
build_system -i data/oxamide/1473522.cif -pre build_system/oxamide/1473522_cgenff  -ff cgenff --toppar toppar_c36_jul21
build_system -i data/oxamide/1473522.cif -pre build_system/1473522_oplsaa  -ff oplsaa

build_system -i data/adipamide/1101307.cif -pre build_system/adipamide/1101307_gaff1   -ff gaff1
build_system -i data/adipamide/1101307.cif -pre build_system/adipamide/1101307_gaff2   -ff gaff2
build_system -i data/adipamide/1101307.cif -pre build_system/adipamide/1101307_parsley -ff parsley
build_system -i data/adipamide/1101307.cif -pre build_system/adipamide/1101307_sage    -ff sage
build_system -i data/adipamide/1101307.cif -pre build_system/adipamide/1101307_cgenff  -ff cgenff --toppar toppar_c36_jul21
build_system -i data/adipamide/1101307.cif -pre build_system/adipamide/1101307_oplsaa  -ff oplsaa
```

In constrast to `make_supercell`, we do not want to specify a fixed number of replicates for each unit cell dimension. Instead we define the nonbonded cutoff `--nbcutoff` and the `--axislengthfactor`, which tells the program how long in multiples of the nonbonded cutoff, each crystallographic axis should be.

Note that by default, we save force field `.xml` files not only for the supercell lattice, but also for each unique monomer molecule in the lattice (`XXX_monomer0.xml`). In addition, we have a `.json` file that contains the supercell lattice as big rdkit molecule. It can be parsed to python again using:

```
import Chem
with open(1473522_cgenff.json, "r") as fopen:
    rdmol = Chem.JSONToMols(fopen.read())[0]
```

## Running Lattice Minimization

After the building the lattice, we should first minimize the lattice. This can be done using the script `run_xtal_min`. In principal, it has two run modes: `run_xtal_min xml` for normal command line based input or `run_xtal_min yaml` for high-throughput runs. We will start with the command line runs:

First, let's optimize the joint coordinates of atomic positions and lattice vectors using 100 steps of BFGS.

```
run_xtal_min xml -i build_system/oxamide/1473522_oplsaa.xml -p build_system/oxamide/1473522_oplsaa.pdb -pre xtal_min/oxamide_oplsaa_min --method BFGS --steps 100
```

As an alternative, we can optimize the atom positions and box vectors seperately and alternate between the two.

```
-i build_system/oxamide/1473522_oplsaa.xml -p build_system/oxamide/1473522_oplsaa.pdb -pre xtal_min/oxamide_oplsaa_min_alt --method BFGS --steps 100 --alternating
```

In general, the option `--method` can be used to define the main optimization algorithm in place. Valid choices are Nelder-Mead, Powell, CG, BFGS, Newton-CG, L-BFGS-B, TNC, COBYLA, SLSQP, trust-const, dogleg, trust-ncg, trust-exact, trust-krylov. See the [scipy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html) for more information.

Now let's do the same thing for a larger number of lattice structures. The commands are all stored in the [`input.yaml`](xtal_min/input.yaml) Go to the folder `xtal_min` first, and then run the following command. But be careful, this will start a lot of processes. If you want to limit the number of available CPUs, first start the ray server manually `ray start --head --port=6831 --redis-password=5241590000000000 --num-cpus=16 --num-gpus=4` and add the line `"ray_host" : "localhost:6381"` to the `input.yaml` file.

```
run_xtal_min yaml -i input.yaml
```

## Running Lattice Simulations

To run MD simulations of the supercell lattice, we use `run_xtal_md`. In accordance with the minimization script `run_xtal_min`, we can run this in two modes, either `run_xtal_md xml` from the command line or `run_xtal_md yaml` as high-throughput.

The command line arguments are very self-explaining. Just run `run_xtal_md xml -h` to see what is available. 

For the high-throughput work, have a look at the [`input.yaml`](xtal_min/input.yaml). This works similarly as for the previous `run_xtal_min yaml` case.

## Data Analysis

The data analysis interface works only through `.yaml` files. This is necessary, since we need to combine lots of different information, which is easier in a structured file format than on the command line. Have a look into the file `analysis/input.yaml` to see what information goes where. The variables should be self-explaining. The final analysis is summarized and saved in an Excel `.xlsx` table. Run the analaysis using:

```
make_tables -i input.yaml -o table.xlsx
```