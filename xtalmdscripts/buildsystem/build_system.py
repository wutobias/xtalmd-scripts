#!/use/bin/env python3

import gemmi
from xtalmdscripts.supercellbuilding import make_supercell
import openmm
from openmm import unit
from openmm import app
from rdkit import Chem
import numpy as np

from pkg_resources import resource_filename

oplsaa_xml_builder_path = resource_filename("xtalmdscripts.data", "oplsaa/build_xml.sh")


def parse_arguments():

    """
    Parse command line arguments.
    """

    import argparse

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
        '--prefix', 
        "-pre", 
        type=str, 
        help="Output prefix xml and pdb files.", 
        default="xtal_min",
        required=False
        )

    parser.add_argument(
        '--forcefield', 
        "-ff", 
        type=str, 
        help="force field name",
        choices=[
            "gaff1",
            "gaff2",
            "parsley",
            "sage",
            "cgenff",
            "oplsaa"
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

    parser.add_argument(
        '--toppar', 
        "-tp", 
        type=str, 
        help="path to topology and parameter files for charm", 
        required=False
        )

    return parser.parse_args()


def OPLS_LJ(system):

    """
    Helper function to get the OPLS combination rules.
    See http://zarbi.chem.yale.edu/ligpargen/openMM_tutorial.html
    """

    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    lorentz = openmm.CustomNonbondedForce(
        '4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)')
    lorentz.setNonbondedMethod(lorentz.CutoffPeriodic)
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    lorentz.setUseLongRangeCorrection(True)
    LJset = {}
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        LJset[index] = (sigma, epsilon)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(
            index, charge, sigma, epsilon * 0
        )
    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        # ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
        # FORCE
        lorentz.addExclusion(p1, p2)
        if eps._value != 0.0:
            #print p1,p2,sig,eps
            sig14 = np.sqrt(LJset[p1][0] * LJset[p2][0])
            eps14 = np.sqrt(LJset[p1][1] * LJset[p2][1]) * 0.5
            ### eps is already scaled
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
    system.addForce(lorentz)
    return system


def build_system_gaff(
    replicated_mol_list,
    pdb_path,
    version="2.11"):

    """
    Build openmm system for gaff force field.
    """

    from rdkit import Chem
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
        removeCMMotion=False,
    )

    return system


def build_system_off(
    replicated_mol_list,
    pdb_path,
    version="1.3.1"):

    """
    Build openmm system for openff force field.
    """

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
        removeCMMotion=False,
    )

    return system


def build_system_cgenff(
    replicated_mol_list,
    pdb_path,
    toppar_dir_path):

    """
    Build openmm system for cgenff force field.
    """

    import subprocess
    import parmed as pmd
    from rdkit import Chem

    from simtk.openmm.app import ForceField
    from simtk.openmm.app import PDBFile, CharmmPsfFile, CharmmParameterSet

    ### cgenff topology
    rtf_path = f"{toppar_dir_path}/top_all36_cgenff.rtf"
    ### cgenff parameters
    prm_path = f"{toppar_dir_path}/par_all36_cgenff.prm"

    unique_mapping, rdmol_list_unique  = make_supercell.get_unique_mapping(replicated_mol_list)
    unique_mol_idxs = set(unique_mapping.values())
    unique_mol_idxs = sorted(unique_mol_idxs)

    ### number replicatations necessary to build
    ### the the replicated mol list from the unique list
    N_unique_replicates = int(len(replicated_mol_list)/len(unique_mol_idxs))
    assert (len(replicated_mol_list) % len(unique_mol_idxs)) == 0

    ### This is the parmed.Structure instance that will
    ### store all topology information
    psf_pmd = pmd.Structure()
    ### Stores the path to all charmm str, topology and parameter files.
    params_path_list = [
        f'{toppar_dir_path}/top_all36_prot.rtf',
        f'{toppar_dir_path}/par_all36m_prot.prm',
        f'{toppar_dir_path}/top_all36_na.rtf',
        f'{toppar_dir_path}/par_all36_na.prm',
        f'{toppar_dir_path}/top_all36_carb.rtf',
        f'{toppar_dir_path}/par_all36_carb.prm',
        f'{toppar_dir_path}/top_all36_lipid.rtf',
        f'{toppar_dir_path}/par_all36_lipid.prm',
        f'{toppar_dir_path}/top_all36_cgenff.rtf',
        f'{toppar_dir_path}/par_all36_cgenff.prm',
    ]
    for mol_idx in unique_mol_idxs:

        mol = rdmol_list_unique[mol_idx]

        mi = mol.GetAtomWithIdx(0).GetMonomerInfo()
        resname = mi.GetResidueName()

        pdb_path_monomer = f"./mol_{mol_idx}.pdb"

        with open(pdb_path_monomer, "w") as fopen:
            fopen.write(
                Chem.MolToPDBBlock(mol)
                )

        ### mol2 path
        mol2_path_monomer = pdb_path_monomer.replace(".pdb", ".mol2")
        ### stream file path
        str_path_monomer  = pdb_path_monomer.replace(".pdb", ".str")
        params_path_list.append(str_path_monomer)
        ### psf file path
        psf_path_monomer  = pdb_path_monomer.replace(".pdb", ".psf")

        ### Convert to mol2 file using obabel
        subprocess.run([
            "obabel",
            "-ipdb",
            f"./mol_{mol_idx}.pdb",
            "--title",
            resname,
            "-omol2",
            f"-O{mol2_path_monomer}",
            ])

        ### Obtain cgenff stream file using cgenff
        subprocess.run([
            "cgenff",
            "-p",
            rtf_path,
            prm_path,
            "-b",
            "-f",
            str_path_monomer,
            mol2_path_monomer
            ])

        charmm_inp = f"""
! read topology and parameter files

! protein topology and parameter
open read card unit 10 name {toppar_dir_path}/top_all36_prot.rtf
read  rtf card unit 10

open read card unit 20 name {toppar_dir_path}/par_all36m_prot.prm
read para card unit 20 flex

! nucleic acids
open read card unit 10 name {toppar_dir_path}/top_all36_na.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir_path}/par_all36_na.prm
read para card unit 20 append flex

! carbohydrates
open read card unit 10 name {toppar_dir_path}/top_all36_carb.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir_path}/par_all36_carb.prm
read para card unit 20 append flex

! lipids
open read card unit 10 name {toppar_dir_path}/top_all36_lipid.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir_path}/par_all36_lipid.prm
read para card unit 20 append flex

! CGENFF
open read card unit 10 name {toppar_dir_path}/top_all36_cgenff.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir_path}/par_all36_cgenff.prm
read para card unit 20 append flex

stream {str_path_monomer}

! Generate sequence
read sequence {resname} 1
generate {resname} first none last none setup

open write unit 10 card name {psf_path_monomer}
write psf  unit 10 card

stop
"""
        with open("make_psf.inp", "w") as fopen:
            fopen.write(charmm_inp)

        ### Obtain psf file using charmm
        subprocess.run([
            "charmm",
            "-i",
            "make_psf.inp"
            ])

        psf_pmd += pmd.load_file(psf_path_monomer)

    ### Replicate the unique molecules
    ### in order to match the full `replicated_mol_list`
    psf_pmd *= N_unique_replicates
    psf_pmd.write_psf("./strc.psf")

    ### Build the omm system
    psffile  = CharmmPsfFile("./strc.psf")
    topology = psffile.topology
    pdbfile  = PDBFile(pdb_path)
    boxvectors = pdbfile.topology.getPeriodicBoxVectors()
    positions  = pdbfile.positions

    params  = CharmmParameterSet(*params_path_list)
    psffile.loadParameters(params)

    a = boxvectors[0]
    b = boxvectors[1]
    c = boxvectors[2]

    a_len = np.linalg.norm(a)
    b_len = np.linalg.norm(b)
    c_len = np.linalg.norm(c)

    a_normalized = a/a_len
    b_normalized = b/b_len
    c_normalized = c/c_len

    alpha = np.arccos(
        np.dot(
            c_normalized,
            b_normalized
            )
        ) * unit.radian

    beta = np.arccos(
        np.dot(
            a_normalized,
            c_normalized
            )
        ) * unit.radian

    gamma = np.arccos(
        np.dot(
            a_normalized,
            b_normalized
            )
        ) * unit.radian

    psffile.setBox(
        a_len,
        b_len,
        c_len,
        alpha,
        beta,
        gamma,
    )

    system = psffile.createSystem(
        params,
        nonbondedMethod=app.PME, 
        constraints=app.HBonds,
        removeCMMotion=False,
    )

    return system


def build_system_oplsaa(
    replicated_mol_list,
    pdb_path,
    version="CM1A"):

    """
    Build openmm system for opls aa force field.
    """

    import subprocess
    from simtk.openmm.app import ForceField
    from simtk.openmm.app import PDBFile

    set_lbcc = 0

    if version == "CM1A-LBCC":
        set_lbcc = 1
    elif version == "CM1A":
        set_lbcc = 0
    else:
        raise ValueError(
            f"version={version} not understood.")

    unique_mapping, rdmol_list_unique = make_supercell.get_unique_mapping(replicated_mol_list)
    unique_mol_idxs = set(unique_mapping.values())
    unique_mol_idxs = sorted(unique_mol_idxs)

    xml_file_list   = list()
    for mol_idx in unique_mol_idxs:

        mol = rdmol_list_unique[mol_idx]

        pdb_path_monomer = f"mol_{mol_idx}.pdb"
        with open(pdb_path_monomer, "w") as fopen:
            fopen.write(
                Chem.MolToPDBBlock(mol)
                )

        mi = mol.GetAtomWithIdx(0).GetMonomerInfo()
        resname = mi.GetResidueName().rstrip().lstrip()

        subprocess.run([
            oplsaa_xml_builder_path,
            pdb_path_monomer,
            resname,
            str(Chem.GetFormalCharge(mol)),
            str(set_lbcc)
            ])

        xml_file_list.append(
            f"{resname}.xml"
            )

    pdbfile    = PDBFile(pdb_path)
    topology   = pdbfile.getTopology()
    boxvectors = topology.getPeriodicBoxVectors()
    positions  = pdbfile.positions
    forcefield = app.ForceField(*xml_file_list)

    system = forcefield.createSystem(
        topology=topology,
        nonbondedMethod=app.PME, 
        constraints=app.HBonds,
        removeCMMotion=False,
    )

    system = OPLS_LJ(system)

    return system


def main():

    """
    Build openmm system for given ff and save as xml.
    """

    args = parse_arguments()

    ### We want to build super cell that is longer than 2.5 times the 
    ### nonbonded cutoff in each direction
    min_length_a = args.nbcutoff * unit.nanometer * 2.5
    min_length_b = args.nbcutoff * unit.nanometer * 2.5
    min_length_c = args.nbcutoff * unit.nanometer * 2.5

    strc = make_supercell.parse_cif(args.input)

    uc_length_a = strc.cell.a * unit.angstrom
    uc_length_b = strc.cell.b * unit.angstrom
    uc_length_c = strc.cell.c * unit.angstrom

    a_min_max = [0, np.ceil(min_length_a / uc_length_a)]
    b_min_max = [0, np.ceil(min_length_b / uc_length_b)]
    c_min_max = [0, np.ceil(min_length_c / uc_length_c)]

    a_len = np.max(a_min_max) - np.min(a_min_max) + 1.
    b_len = np.max(b_min_max) - np.min(b_min_max) + 1.
    c_len = np.max(c_min_max) - np.min(c_min_max) + 1.

    print(f"Building super cell with with {a_len} x {b_len} x {c_len}")

    ### Build the supercell as a list of rdkit molecules
    ### ================================================
    replicated_mol_list = make_supercell.generate_replicated_mol_list(
        strc,
        a_min_max,
        b_min_max,
        c_min_max,
        )

    ### Write pdb file
    ### ==============
    pdb_str = make_supercell.get_pdb_str(
        replicated_mol_list, 
        strc,
        a_min_max,
        b_min_max,
        c_min_max
        )

    prefix = args.prefix
    if prefix.endswith(".xml"):
        prefix = prefix.replace(".xml", "")

    with open(f"./{prefix}.pdb", "w") as fopen:
        fopen.write(pdb_str)

    if args.forcefield.lower() == "gaff1":

        system = build_system_gaff(
            replicated_mol_list,
            f"./{prefix}.pdb",
            version="1.81"
            )

    elif args.forcefield.lower() == "gaff2":

        system = build_system_gaff(
            replicated_mol_list,
            f"./{prefix}.pdb",
            version="2.11"
            )

    elif args.forcefield.lower() == "parsley":

        system = build_system_off(
            replicated_mol_list,
            f"./{prefix}.pdb",
            version="1.3.1"
            )

    elif args.forcefield.lower() == "sage":

        system = build_system_off(
            replicated_mol_list,
            f"./{prefix}.pdb",
            version="2.0.0"
            )

    elif args.forcefield.lower() == "cgenff":

        system = build_system_cgenff(
            replicated_mol_list,
            f"./{prefix}.pdb",
            args.toppar
            )

    elif args.forcefield.lower() == "oplsaa":

        system = build_system_oplsaa(
            replicated_mol_list,
            f"./{prefix}.pdb",
            version="CM1A-LBCC"
            )

        ### Set nonbonded cutoff special for oplsaa
        forces = {system.getForce(index).__class__.__name__: system.getForce(
            index) for index in range(system.getNumForces())}
        nbforce = forces['CustomNonbondedForce']
        nbforce.setCutoffDistance(args.nbcutoff * unit.nanometer)

    else:
        raise ValueError(
            f"force field {args.forcefield} not understood")

    ### Set nonbonded cutoff
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nbforce = forces['NonbondedForce']
    nbforce.setCutoffDistance(args.nbcutoff * unit.nanometer)

    with open(f"./{prefix}.xml", "w") as fopen:
        fopen.write(openmm.XmlSerializer.serialize(system))

def entry_point():

    main()

if __name__ == "__main__":

    entry_point()