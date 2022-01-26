#!/use/bin/env python3

import gemmi
from .supercellbuilding import make_supercell
import openmm
from openmm import unit
from openmm import app
from rdkit import Chem
import numpy as np

def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Python script for building openmm xml files for xtal md simulations with different force fields."
        )

    parser.add_argument(
        '--input', 
        "-i", 
        type=str, 
        help="Input cif file", 
        required=True
        )

    parser.add_argument(
        '--output', 
        "-o", 
        type=str, 
        help="output xml file", 
        required=True
        )

    parser.add_argument(
        '--forcefield', 
        "-ff", 
        type=str, 
        help="force field name",
        choices=[
        	"Gaff1", 
        	"Gaff2", 
        	"Parsley",
        	"Sage",
        	"Cgenff",
        	"Oplsaa"
        	],
        required=True
        )

    parser.add_argument(
        '--nbcutoff', 
        "-nb", 
        type=float, 
        help="nonbonded cutoff distance in nm", 
        default="0.8",
        required=False
        )

    return parser.parse_args()


def build_system_gaff(
    replicated_mol_list,
    pdb_path,
    version="2.11"):

    from openff.toolkit.topology import Molecule
    from openmmforcefields.generators import GAFFTemplateGenerator
    from simtk.openmm.app import ForceField
    from simtk.openmm.app import PDBFile

    forcefield  = ForceField()
    offmol_list = [Molecule.from_rdkit(rdmol) for rdmol in replicated_mol_list]
    N_mol       = len(replicated_mol_list)
    gaff        = GAFFTemplateGenerator(
        molecules=offmol_list, 
        forcefield=f'gaff-{version}'
        )
    forcefield.registerTemplateGenerator(gaff.generator)
    pdbfile  = PDBFile(pdb_path)
    topology = pdbfile.getTopology()
    positions = pdbfile.getPositions()
    boxvectors = topology.getPeriodicBoxVectors()
    system   = forcefield.createSystem(
        topology = topology,
        nonbondedMethod=app.PME, 
        constraints=app.HBonds,
        removeCMMotion=True,        
    )

    return system


def build_system_off(
    replicated_mol_list,
    pdb_path,
    version="1.3.1"):

    from openff.toolkit.topology import Molecule
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator
    from simtk.openmm.app import ForceField
    from simtk.openmm.app import PDBFile

    forcefield  = ForceField()
    offmol_list = [Molecule.from_rdkit(rdmol) for rdmol in replicated_mol_list]
    N_mol       = len(replicated_mol_list)
    gaff        = SMIRNOFFTemplateGenerator(
        molecules=offmol_list, 
        forcefield=f"openff-{version}"
        )
    forcefield.registerTemplateGenerator(gaff.generator)
    pdbfile  = PDBFile(pdb_path)
    topology = pdbfile.getTopology()
    positions = pdbfile.getPositions()
    boxvectors = topology.getPeriodicBoxVectors()
    system   = forcefield.createSystem(
        topology = topology,
        nonbondedMethod=app.PME, 
        constraints=app.HBonds,
        removeCMMotion=True,        
    )

    return system


def build_system_cgenff(
    replicated_mol_list,
    pdb_path,
    version="all36_cgenff"):

    from .data import cgenff
    import sys
    import parmed as pmd

    toppar_dir = f"{cgenff.__path__}/toppar_c36_jul21"

    ### topology
    rtf_path = f"{toppar_dir}/top_{version}.rtf"
    ### parameters
    prm_path = f"{toppar_dir}/par_{version}.prm"
    ### mol2 path
    mol2_path = f"O{pdb_path.replace(".pdb", ".mol2")}"
    ### stream file path
    str_path   = f"O{pdb_path.replace(".pdb", ".str")}"
    ### psf file path
    psf_path   = f"O{pdb_path.replace(".pdb", ".psf")}"
    ### Number of molecules
    N_mol      = len(replicated_mol_list)

    ### Convert to mol2 file using obabel
    sys.call([
        "obabel",
        "-ipdb",
        pdb_path,
        "-omol2",
        mol2_path,
        ])

    ### Obtain cgenff stream file using cgenff
    sys.call([
        "cgenff",
        "-p",
        rtf_path,
        prm_path,
        "-b",
        "-f",
        str_path,
        mol2_path
        ])

    charmm_inp = f"""
! read topology and parameter files

! protein topology and parameter
open read card unit 10 name {toppar_dir}/top_all36_prot.rtf
read  rtf card unit 10

open read card unit 20 name {toppar_dir}/par_all36m_prot.prm
read para card unit 20 flex

! nucleic acids
open read card unit 10 name {toppar_dir}/top_all36_na.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir}/par_all36_na.prm
read para card unit 20 append flex

! carbohydrates
open read card unit 10 name {toppar_dir}/top_all36_carb.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir}/par_all36_carb.prm
read para card unit 20 append flex

! lipids
open read card unit 10 name {toppar_dir}/top_all36_lipid.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir}/par_all36_lipid.prm
read para card unit 20 append flex

! CGENFF
open read card unit 10 name {toppar_dir}/top_all36_cgenff.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir}/par_all36_cgenff.prm
read para card unit 20 append flex

stream {str_path}

! Generate sequence
! read sequence mol {N_mol}
read sequence mol 1
generate mol first none last none setup

open write unit 10 card name {psf_path}
write psf  unit 10 card

stop
    """

    with open("make_psf.inp", "w") as fopen:
        fopen.write(charmm_inp)

    ### Obtain psf file using charmm
    sys.call([
        "charmm",
        "-i",
        charmm_inp
        ])

    psf_pmd  = pmd.load_file(psf_path)
    psf_pmd *= N_mol
    pmd.Structure.write_psf(
        psf_pmd, 
        psf_path)

    pdbfile


def main():

    args = parse_arguments()

    ### We want to build super cell that is longer than 2.5 times the 
    ### nonbonded cutoff in each direction
    min_length_a = args.nbcutoff * unit.nanometer * 2.5
    min_length_b = args.nbcutoff * unit.nanometer * 2.5
    min_length_c = args.nbcutoff * unit.nanometer * 2.5

    doc  = gemmi.cif.read(args.input)[0]
    strc = gemmi.make_small_structure_from_block(doc)

    uc_length_a = strc.cell.a * unit.angstrom
    uc_length_b = strc.cell.b * unit.angstrom
    uc_length_c = strc.cell.c * unit.angstrom

    a_min_max = [0, np.ceil(min_length_a / uc_length_a)]
    b_min_max = [0, np.ceil(min_length_b / uc_length_b)]
    c_min_max = [0, np.ceil(min_length_c / uc_length_c)]

    a_len = np.max(a_min_max) - np.min(a_min_max) + 1.
    b_len = np.max(b_min_max) - np.min(b_min_max) + 1.
    c_len = np.max(c_min_max) - np.min(c_min_max) + 1.

    strc_write               = gemmi.Structure()
    strc_write.spacegroup_hm = strc.spacegroup_hm
    strc_write.cell          = gemmi.UnitCell(
        strc.cell.a * a_len,
        strc.cell.b * b_len,
        strc.cell.c * c_len,
        strc.cell.alpha,
        strc.cell.beta,
        strc.cell.gamma
    )

    ### Build the supercell as a set of rdkit molecules
    replicated_mol_list = make_supercell.generate_replicated_mol_list(
        strc,
        a_min_max,
        b_min_max,
        c_min_max,
        )

    with open("./strc.pdb", "w") as fopen:
        fopen.write(make_supercell.get_pdb_block(replicated_mol_list, strc_write))
    with open("./monomer.pdb", "w") as fopen:
        fopen.write(Chem.MolToPDBBlock(replicated_mol_list[0]))

    if args.forcefield.lower() == "gaff1":

        system = build_system_gaff(
            "./strc.pdb",
            replicated_mol_list,
            version="1.81")

    elif args.forcefield.lower() == "gaff2":

        system = build_system_gaff(
            "./strc.pdb",
            replicated_mol_list,
            version="2.11")

    elif args.forcefield.lower() == "parsley":

        system = build_system_off(
            "./strc.pdb",
            replicated_mol_list,
            version="1.3.1")

    elif args.forcefield.lower() == "sage":

        system = build_system_off(
            "./strc.pdb",
            replicated_mol_list,
            version="2.0.0")

    elif args.forcefield.lower() == "cgenff":



    ### Set nonbonded cutoff
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nbforce = forces['NonbondedForce']
    nbforce.setCutoffDistance(args.nbcutoff * unit.nanometer)

    with open(args.output, "w") as fopen:
        fopen.write(openmm.XmlSerializer.serialize(system))

def entry_point():

    main()

if __name__ == "__main__":

    entry_point()