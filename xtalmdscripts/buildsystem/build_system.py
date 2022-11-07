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
tip3p_xml_path = resource_filename("xtalmdscripts.data", "tip3p/tip3p.xml")
tip3p_opls_xml_path = resource_filename("xtalmdscripts.data", "tip3p/tip3p_opls.xml")


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
        '--version', 
        "-ve", 
        type=str, 
        help="Specifier for force field version. For instance, 1.81 for GAFF1, 1.3.1 for Parsley or CM1A-LBCC for oplsaa.",
        default="",
        required=False,
        )

    parser.add_argument(
        '--nbcutoff', 
        "-nb", 
        type=float, 
        help="nonbonded cutoff distance in nm", 
        default=0.9,
        required=False
        )

    parser.add_argument(
        '--axislengthfactor', 
        "-ax", 
        type=float, 
        help="Each axis is add least this value times the nbcutoff in length.", 
        default=5.0,
        required=False
        )

    parser.add_argument(
        '--toppar', 
        "-tp", 
        type=str, 
        help="path to topology and parameter files for charm", 
        required=False
        )

    parser.add_argument(
        '--addhs', 
        "-ah", 
        action='store_true',
        help="Remove any existing hydrogen and add protonate molecule. Requires OpenEye Toolkits.", 
        required=False,
        default=False,
        )

    parser.add_argument(
        '--addwater', 
        "-aw", 
        type=int, 
        help="Number of water molecules to add", 
        required=False,
        default=0,
        )

    parser.add_argument(
        '--use_symmetry_operations', 
        "-op", 
        action='store_true',
        help="Use symmetry operations in cif file instead of space group.", 
        required=False,
        default=False,
        )
    
    parser.add_argument(
        '--n_protonation_attempts', 
        "-np", 
        type=int, 
        help="Number of attempts to compute protonatation states in  unit cell.", 
        required=False,
        default=0,
        )

    parser.add_argument(
        '--use_openeye', 
        "-oe", 
        action='store_true',
        help="Use openeye-toolkit for topology building. Otherwise use xyz2mol.", 
        required=False,
        default=False,
        )

    parser.add_argument(
        '--use_tip3p', 
        "-tp3", 
        action='store_true',
        help="If there is water in the system, make it Tip3p.", 
        required=False,
        default=False,
        )

    parser.add_argument(
        '--charge_method', 
        "-cm", 
        type=str, 
        help="Charge method used for partial charge assignment. Choice `default` uses standard method for forcefield.", 
        choices=[
            "default",
            "am1bcc", 
            "am1elf10", 
            "am1-mulliken", 
            "gasteiger"
            ],
        required=False,
        default="default"
        )

    return parser.parse_args()


def OPLS_LJ(system, CutoffPeriodic=True):

    """
    Helper function to get the OPLS combination rules.
    See http://zarbi.chem.yale.edu/ligpargen/openMM_tutorial.html

    Scaling factors for 1-4 interactions are 0.5 for both LJ and charge.
    These are assumed to be already correctly scaled when the system
    is parsed.
    """

    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    lorentz = openmm.CustomNonbondedForce(
        '4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)')
    if CutoffPeriodic:
        lorentz.setNonbondedMethod(lorentz.CutoffPeriodic)
    else:
        lorentz.setNonbondedMethod(lorentz.NoCutoff)
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
            sig14 = np.sqrt(LJset[p1][0] * LJset[p2][0])
            eps14 = np.sqrt(LJset[p1][1] * LJset[p2][1]) * 0.5
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
    system.addForce(lorentz)
    return system


def build_system_gaff(
    replicated_mol_list,
    pdb_path,
    version="2.11",
    use_tip3p=False):

    """
    Build openmm system for gaff force field.
    """

    from rdkit import Chem
    from openff.toolkit.topology import Molecule
    from openmmforcefields.generators import GAFFTemplateGenerator
    from simtk.openmm.app import ForceField
    from simtk.openmm.app import PDBFile

    params = Chem.SmilesParserParams()
    params.removeHs = False
    wat_rdmol = Chem.MolFromSmiles("[H]O[H]", params)
    offmol_list = list()
    has_water = False
    wat_list  = list()
    for rdmol_idx, rdmol in enumerate(replicated_mol_list):
        if use_tip3p:
            if rdmol.HasSubstructMatch(wat_rdmol):
                has_water = True
                wat_list.append(rdmol_idx)
            else:
                offmol_list.append(
                    Molecule.from_rdkit(rdmol)
                    )
        else:
            offmol_list.append(
                Molecule.from_rdkit(rdmol)
                )

    if use_tip3p and has_water:
        forcefield  = ForceField(tip3p_xml_path)
    else:
        forcefield  = ForceField()
    gaff = GAFFTemplateGenerator(
        molecules=offmol_list, 
        forcefield=f'gaff-{version}'
        )
    forcefield.registerTemplateGenerator(gaff.generator)
    pdbfile  = PDBFile(pdb_path)
    topology = pdbfile.getTopology()
    positions = pdbfile.getPositions()
    boxvectors = topology.getPeriodicBoxVectors()
    if use_tip3p and has_water:
        for res in topology.residues():
            if res.index in wat_list:
                res.name = "HOH"
                found_one_H = False
                for atm in res.atoms():
                    if atm.element.atomic_number == 8:
                        atm.name = "O"
                    elif atm.element.atomic_number == 1 and found_one_H:
                        atm.name = "H2"
                    else:
                        atm.name = "H1"
                        found_one_H = True

    if boxvectors == None:
        nonbondedMethod = app.NoCutoff
    else:
        nonbondedMethod = app.PME
    system = forcefield.createSystem(
        topology = topology,
        nonbondedMethod=nonbondedMethod,
        constraints=None,
        removeCMMotion=False,
        rigidWater=use_tip3p and has_water,
    )

    return system


def build_system_off(
    replicated_mol_list,
    pdb_path,
    version="1.3.1",
    use_tip3p=False):

    """
    Build openmm system for openff force field.
    """

    from openff.toolkit.topology import Molecule
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator
    from simtk.openmm.app import ForceField
    from simtk.openmm.app import PDBFile

    params = Chem.SmilesParserParams()
    params.removeHs = False
    wat_rdmol = Chem.MolFromSmiles("[H]O[H]", params)
    offmol_list = list()
    has_water = False
    wat_list  = list()
    for rdmol_idx, rdmol in enumerate(replicated_mol_list):
        if use_tip3p:
            if rdmol.HasSubstructMatch(wat_rdmol):
                has_water = True
                wat_list.append(rdmol_idx)
            else:
                offmol_list.append(
                    Molecule.from_rdkit(rdmol)
                    )
        else:
            offmol_list.append(
                Molecule.from_rdkit(rdmol)
                )

    if use_tip3p and has_water:
        forcefield  = ForceField(tip3p_xml_path)
    else:
        forcefield  = ForceField()
    openff      = SMIRNOFFTemplateGenerator(
        molecules=offmol_list, 
        forcefield=f"openff-{version}"
        )
    forcefield.registerTemplateGenerator(openff.generator)
    pdbfile  = PDBFile(pdb_path)
    topology = pdbfile.getTopology()
    positions = pdbfile.getPositions()
    boxvectors = topology.getPeriodicBoxVectors()
    if use_tip3p and has_water:
        for res in topology.residues():
            if res.index in wat_list:
                res.name = "HOH"
                found_one_H = False
                for atm in res.atoms():
                    if atm.element.atomic_number == 8:
                        atm.name = "O"
                    elif atm.element.atomic_number == 1 and found_one_H:
                        atm.name = "H2"
                    else:
                        atm.name = "H1"
                        found_one_H = True
    if boxvectors == None:
        nonbondedMethod = app.NoCutoff
    else:
        nonbondedMethod = app.PME
    system = forcefield.createSystem(
        topology = topology,
        nonbondedMethod=nonbondedMethod,
        constraints=None,
        removeCMMotion=False,
        rigidWater=use_tip3p and has_water,
    )

    return system


def build_system_cgenff(
    replicated_mol_list,
    pdb_path,
    toppar_dir_path,
    use_tip3p=False):

    """
    Build openmm system for cgenff force field.
    """

    import os
    import subprocess
    import glob
    import parmed as pmd
    from rdkit import Chem
    from .utils import get_params_path_list, get_charmm_inp_header

    from simtk.openmm.app import ForceField
    from simtk.openmm.app import PDBFile, CharmmPsfFile, CharmmParameterSet

    ### cgenff topology
    rtf_path_cgenff = f"{toppar_dir_path}/top_all36_cgenff.rtf"
    ### cgenff parameters
    prm_path_cgenff = f"{toppar_dir_path}/par_all36_cgenff.prm"

    unique_mapping, rdmol_list_unique  = make_supercell.get_unique_mapping(
        replicated_mol_list,
        stereochemistry=False
        )

    ### Stores the names to all charmm str, topology and parameter files.
    params_path_list  = get_params_path_list(toppar_dir_path)
    charmm_inp_header = get_charmm_inp_header(toppar_dir_path)

    to_remove_list = list()
    str_list = list()
    pmd_list = list()
    failed   = False
    already_processed = list()

    for mol_idx in unique_mapping:
        unique_idx = unique_mapping[mol_idx]
        if unique_idx in already_processed:
            continue
        else:
            already_processed.append(
                unique_idx
                )
        rdmol = rdmol_list_unique[unique_idx]

        mi = rdmol.GetAtomWithIdx(0).GetMonomerInfo()
        resname = mi.GetResidueName()

        pdb_path_monomer = f"mol_{mol_idx}.pdb"

        with open(pdb_path_monomer, "w") as fopen:
            fopen.write(
                Chem.MolToPDBBlock(rdmol)
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
            f"mol_{mol_idx}.pdb",
            "--title",
            resname,
            "-omol2",
            f"-O{mol2_path_monomer}",
            ])

        ### Obtain cgenff stream file using cgenff
        subprocess.run([
            "cgenff",
            "-p",
            rtf_path_cgenff,
            prm_path_cgenff,
            "-b",
            "-f",
            str_path_monomer,
            mol2_path_monomer
            ])

        charmm_inp = f"""
{charmm_inp_header}

stream {str_path_monomer}

! Generate sequence
bomlev -1
read sequence {resname} 1
generate {resname} first none last none setup
bomlev 0

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
            "make_psf.inp",
            "-o",
            "/dev/null"
            ])

        to_remove_list.extend([
            str_path_monomer,
            pdb_path_monomer,
            mol2_path_monomer,
            psf_path_monomer,
            "make_psf.inp",
            "strc.psf"
            ])

        if os.path.exists(psf_path_monomer):
            try:
                pmd_list.append(
                    pmd.load_file(psf_path_monomer)
                    )
                with open(str_path_monomer, "r") as fopen:
                    str_list.append(fopen.read())
            except:
                failed = True
                break
        else:
            failed = True
            break

    if failed:
        ### Clean up
        for filename in to_remove_list:
            if os.path.exists(filename):
                os.remove(filename)
        raise RuntimeError(
            "Could not build structure with cgenff."
            )
        
    ### This is the parmed.Structure instance that will
    ### store all topology information
    psf_pmd = pmd.Structure()
    for replicated_mol_idx in range(len(replicated_mol_list)):
        psf_pmd += pmd_list[unique_mapping[replicated_mol_idx]]
    psf_pmd.write_psf("strc.psf")

    ### Build the omm system
    psffile  = CharmmPsfFile("strc.psf")
    topology = psffile.topology
    pdbfile  = PDBFile(pdb_path)
    boxvectors = pdbfile.topology.getPeriodicBoxVectors()
    positions  = pdbfile.positions

    params  = CharmmParameterSet(*params_path_list)
    psffile.loadParameters(params)

    if boxvectors == None:
        nonbondedMethod = app.NoCutoff
    else:
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

        nonbondedMethod = app.PME
    system = psffile.createSystem(
        params,
        nonbondedMethod=nonbondedMethod, 
        constraints=None,
        removeCMMotion=False,
        rigidWater=use_tip3p,
    )

    ### Clean up
    for filename in to_remove_list:
        if os.path.exists(filename):
            os.remove(filename)

    return system, str_list


def build_system_oplsaa(
    replicated_mol_list,
    pdb_path,
    version="CM1A",
    use_tip3p=False):

    """
    Build openmm system for opls aa force field.
    """

    import os
    import subprocess
    from simtk.openmm.app import ForceField
    from simtk.openmm.app import PDBFile
    from simtk.openmm.app import modeller
    from rdkit import Chem

    set_lbcc = 0

    if version == "CM1A-LBCC":
        set_lbcc = 1
    elif version == "CM1A":
        set_lbcc = 0
    else:
        raise ValueError(
            f"version={version} not understood.")

    unique_mapping, rdmol_list_unique = make_supercell.get_unique_mapping(
        replicated_mol_list,
        stereochemistry=False
        )
    already_processed = list()

    xml_file_list = list()
    to_remove_list = list()
    failed = False

    params = Chem.SmilesParserParams()
    params.removeHs = False
    wat_rdmol = Chem.MolFromSmiles("[H]O[H]", params)
    has_water = False
    wat_list  = list()
    for mol_idx in unique_mapping:
        unique_idx = unique_mapping[mol_idx]
        if unique_idx in already_processed:
            continue
        else:
            already_processed.append(
                unique_idx
                )
        rdmol = rdmol_list_unique[unique_idx]

        if use_tip3p:
            if rdmol.HasSubstructMatch(wat_rdmol):
                has_water = True
                for mol_idx_query in unique_mapping.keys():
                    if unique_mapping[mol_idx_query] == unique_idx:
                        if not mol_idx_query in wat_list:
                            wat_list.append(mol_idx_query)
                continue

        mi = rdmol.GetAtomWithIdx(0).GetMonomerInfo()
        resname = mi.GetResidueName().rstrip().lstrip()

        pdb_path_monomer = f"mol_{unique_idx}.pdb"
        mol_path_monomer = f"mol_{unique_idx}.mol"
        with open(pdb_path_monomer, "w") as fopen:
            fopen.write(
                Chem.MolToPDBBlock(rdmol)
                )
        with Chem.SDWriter(mol_path_monomer) as sdwriter:
            sdwriter.write(rdmol)

        subprocess.run([
            oplsaa_xml_builder_path,
            pdb_path_monomer,
            resname,
            str(int(Chem.GetFormalCharge(rdmol))),
            str(set_lbcc)
            ])

        xml_file_list.append(
            f"{resname}.openmm.xml"
            )

        try:
            pdbfile_renamed  = PDBFile(f"{resname}.openmm.pdb")
            pdbfile_original = PDBFile(f"mol_{unique_idx}.pdb")

            assert pdbfile_renamed.topology.getNumAtoms() == pdbfile_original.topology.getNumAtoms()

            N_atoms = pdbfile_renamed.topology.getNumAtoms()
        except:
            failed = True
            to_remove_list.extend([
                f"{resname}.openmm.pdb",
                f"{resname}.openmm.xml",
                pdb_path_monomer,
                mol_path_monomer,
            ])
            break

        to_remove_list.extend([
            f"{resname}.openmm.pdb",
            f"{resname}.openmm.xml",
            pdb_path_monomer,
            mol_path_monomer,
        ])

        import xml.etree.ElementTree as ET
        
        tree = ET.parse(f"{resname}.openmm.xml")
        root = tree.getroot()
        ### We will have double atom types and names
        ### if we have multiple molecules. So force them
        ### to be different between molecules (even when
        ### parameters are identical). Note that charges 
        ### might be slightly different even between 
        ### enantiomeric molecules.
        class_dict = dict()
        for r1 in root:
            if r1.tag == "AtomTypes":
                N_atoms = 0
                for r2 in r1:
                    _class = r2.attrib["class"]
                    class_dict[_class] = N_atoms
                    N_atoms += 1
            for r2 in r1:
                if r2.tag == "Improper" and r1.tag == "PeriodicTorsionForce":
                    ### Gather atom indices of atoms in improper dihedral
                    allatomidx = [
                        class_dict[r2.attrib["class1"]],
                        class_dict[r2.attrib["class2"]],
                        class_dict[r2.attrib["class3"]],
                        class_dict[r2.attrib["class4"]]
                    ]
                    ### First atom must be central atom. Re-order.
                    ### First find central atom.
                    bonds_count = [0,0,0,0]
                    for bond in pdbfile_original.topology.bonds():
                        id1 = bond.atom1.index 
                        id2 = bond.atom2.index
                        if id1 in allatomidx and id2 in allatomidx:
                            i1 = allatomidx.index(id1)
                            i2 = allatomidx.index(id2)
                            bonds_count[i1] += 1
                            bonds_count[i2] += 1
                    centralatom_idx = None
                    for atom, count in enumerate(bonds_count):
                        if count == 3:
                            centralatom_idx = atom
                    if centralatom_idx == None:
                        raise ValueError(
                            "Could not determine Improper central atom."
                            )
                    ### Find old and new class values and
                    ### keys.
                    for _class in class_dict:
                        if class_dict[_class] == allatomidx[centralatom_idx]:
                            _class1_old = r2.attrib["class1"]
                            _class1_new = _class
                            _class_idx  = f"class{centralatom_idx+1:d}"
                    ### Write them.
                    r2.attrib["class1"]   = _class1_new
                    r2.attrib[_class_idx] = _class1_old

                for attrib in r2.attrib:
                    if any([attrib.startswith(c) for c in ["class", "name", "type"]]):
                        r2.attrib[attrib]  += f"_{unique_idx}"
                for r3 in r2:
                    for attrib in r3.attrib:
                        if any([attrib.startswith(c) for c in ["class", "name", "type"]]):
                            r3.attrib[attrib]  += f"_{unique_idx}"

        with open(f"{resname}.openmm.xml", "wb") as fopen:
            fopen.write(ET.tostring(root))

    if failed:
        ### Clean up
        for filename in to_remove_list:
            if os.path.exists(filename):
                os.remove(filename)
        raise RuntimeError(
            "Could not build structure with oplsaa."
            )
    pdbfile    = PDBFile(pdb_path)
    topology   = pdbfile.getTopology()
    boxvectors = topology.getPeriodicBoxVectors()
    positions  = pdbfile.positions
    if use_tip3p and has_water:
        xml_file_list.append(tip3p_opls_xml_path)
    forcefield = app.ForceField(*xml_file_list)

    if use_tip3p and has_water:
        for res in topology.residues():
            if res.index in wat_list:
                res.name = "HOH"
                found_one_H = False
                for atm in res.atoms():
                    if atm.element.atomic_number == 8:
                        atm.name = "O"
                    elif atm.element.atomic_number == 1 and found_one_H:
                        atm.name = "H2"
                    else:
                        atm.name = "H1"
                        found_one_H = True

    if boxvectors == None:
        nonbondedMethod  = app.NoCutoff
        lorentz_periodic = False
    else:
        nonbondedMethod  = app.PME
        lorentz_periodic = True
    system = forcefield.createSystem(
        topology = topology,
        nonbondedMethod=nonbondedMethod,
        constraints=None,
        removeCMMotion=False,
        rigidWater=use_tip3p and has_water,
    )

    system = OPLS_LJ(system, lorentz_periodic)

    ### Clean up
    for filename in to_remove_list:
        if os.path.exists(filename):
            os.remove(filename)

    return system


def main():

    """
    Build openmm system for given ff and save as xml.
    """

    args = parse_arguments()

    strc, atom_crds_ortho, atom_num = make_supercell.parse_cif(
        args.input, 
        args.use_symmetry_operations
        )

    a_length = strc.cell.a
    b_length = strc.cell.b
    c_length = strc.cell.c

    alpha    = strc.cell.alpha * np.pi / 180.
    beta     = strc.cell.beta  * np.pi / 180.
    gamma    = strc.cell.gamma * np.pi / 180.

    ### First transform box to reduced form, which is used by openmm.
    ### This code is from openmm.app.internal.unitcell.computePeriodicBoxVectors
    ### See also https://github.com/openmm/openmm/blob/master/wrappers/python/openmm/app/internal/unitcell.py
    a = [a_length, 0, 0]
    b = [b_length*np.cos(gamma), b_length*np.sin(gamma), 0]
    cx = c_length*np.cos(beta)
    cy = c_length*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)
    cz = np.sqrt(c_length*c_length-cx*cx-cy*cy)
    c = [cx, cy, cz]

    # If any elements are very close to 0, set them to exactly 0.

    for i in range(3):
        if abs(a[i]) < 1e-6:
            a[i] = 0.0
        if abs(b[i]) < 1e-6:
            b[i] = 0.0
        if abs(c[i]) < 1e-6:
            c[i] = 0.0
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)

    # Make sure they're in the reduced form required by OpenMM.
    c = c - b*np.round(c[1]/b[1])
    c = c - a*np.round(c[0]/a[0])
    b = b - a*np.round(b[0]/a[0])

    uc_length_a = a[0] * unit.angstrom
    uc_length_b = b[1] * unit.angstrom
    uc_length_c = c[2] * unit.angstrom

    ### We want to build super cell that is longer than args.axislengthfactor
    ### times the nonbonded cutoff in each direction.
    min_length_a = args.nbcutoff * unit.nanometer * args.axislengthfactor
    min_length_b = args.nbcutoff * unit.nanometer * args.axislengthfactor
    min_length_c = args.nbcutoff * unit.nanometer * args.axislengthfactor

    a_min_max = [0, np.ceil(min_length_a / uc_length_a)]
    b_min_max = [0, np.ceil(min_length_b / uc_length_b)]
    c_min_max = [0, np.ceil(min_length_c / uc_length_c)]

    a_len = np.max(a_min_max) - np.min(a_min_max) + 1.
    b_len = np.max(b_min_max) - np.min(b_min_max) + 1.
    c_len = np.max(c_min_max) - np.min(c_min_max) + 1.

    print(f"Building super cell with with {a_len} x {b_len} x {c_len}")

    ### Build the supercell as a list of rdkit molecules
    ### ================================================
    replicated_mol_list, mol_identifies, unitcell_in_supercell_fracs = make_supercell.generate_replicated_mol_list(
        cell=strc.cell,
        atom_crds_ortho=atom_crds_ortho,
        atom_num=atom_num,
        a_min_max=a_min_max,
        b_min_max=b_min_max,
        c_min_max=c_min_max,
        addhs=args.addhs,
        addwater=args.addwater,
        N_iterations_protonation=args.n_protonation_attempts,
        use_openeye=args.use_openeye
        )

    ### Write pdb file
    ### ==============
    strc_write               = gemmi.Structure()
    strc_write.spacegroup_hm = "P1"
    strc_write.cell          = strc.cell
    pdb_str = make_supercell.get_pdb_str(
        replicated_mol_list, 
        strc_write,
        a_min_max,
        b_min_max,
        c_min_max
        )

    prefix = args.prefix
    with open(f"{prefix}.pdb", "w") as fopen:
        fopen.write(pdb_str)
    with open(f"{prefix}.csv", "w") as fopen:
        info_str = make_supercell.get_supercell_info_str(
            mol_identifies, 
            unitcell_in_supercell_fracs
            )
        fopen.write(info_str)
    with open(f"{prefix}.json", "w") as fopen:
        json_str = make_supercell.get_replicated_mol_list_json(replicated_mol_list)
        fopen.write(json_str)

    unique_mapping, rdmol_list_unique = make_supercell.get_unique_mapping(
        replicated_mol_list,
        stereochemistry=False
        )
    for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
        pdb_str = Chem.MolToPDBBlock(rdmol)
        with open(f"{prefix}_monomer{rdmol_idx}.pdb", "w") as fopen:
            fopen.write(pdb_str)

    monomer_sys_list = list()
    if args.forcefield.lower() == "gaff1":

        if args.version == "":
            args.version="1.81"
        else:
            if not args.version in ["1.4", "1.7", "1.8", "1.81"]:
                raise ValueError(
                    "Version not understood."
                    )

        system = build_system_gaff(
            replicated_mol_list,
            f"{prefix}.pdb",
            version=args.version,
            use_tip3p=args.use_tip3p,
            )
        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            monomer_sys_list.append(
                build_system_gaff(
                    [rdmol],
                    f"{prefix}_monomer{rdmol_idx}.pdb",
                    version=args.version,
                    use_tip3p=args.use_tip3p,
                    )
                )

    elif args.forcefield.lower() == "gaff2":

        if args.version == "":
            args.version="2.1"
        else:
            if not args.version in ["2.1", "2.11"]:
                raise ValueError(
                    "Version not understood."
                    )

        system = build_system_gaff(
            replicated_mol_list,
            f"{prefix}.pdb",
            version=args.version,
            use_tip3p=args.use_tip3p,
            )
        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            monomer_sys_list.append(
                build_system_gaff(
                    [rdmol],
                    f"{prefix}_monomer{rdmol_idx}.pdb",
                    version=args.version,
                    use_tip3p=args.use_tip3p,
                    )
                )

    elif args.forcefield.lower() == "parsley":

        if args.version == "":
            args.version="1.3.1"
        else:
            if not args.version in ["1.0.0", "1.0.1", "1.1.0", "1.1.1", "1.2.0", "1.3.0", "1.3.1"]:
                raise ValueError(
                    "Version not understood."
                    )

        system = build_system_off(
            replicated_mol_list,
            f"{prefix}.pdb",
            version=args.version,
            use_tip3p=args.use_tip3p,
            )
        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            monomer_sys_list.append(
                build_system_off(
                    [rdmol],
                    f"{prefix}_monomer{rdmol_idx}.pdb",
                    version=args.version,
                    use_tip3p=args.use_tip3p,
                    )
                )

    elif args.forcefield.lower() == "sage":

        if args.version == "":
            args.version="2.0.0"
        else:
            if not args.version in ["2.0.0"]:
                raise ValueError(
                    "Version not understood."
                    )

        system = build_system_off(
            replicated_mol_list,
            f"{prefix}.pdb",
            version=args.version,
            use_tip3p=args.use_tip3p,
            )
        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            monomer_sys_list.append(
                build_system_off(
                    [rdmol],
                    f"{prefix}_monomer{rdmol_idx}.pdb",
                    version=args.version,
                    use_tip3p=args.use_tip3p,
                    )
                )

    elif args.forcefield.lower() == "cgenff":

        if args.toppar == None:
            raise ValueError(f"When building cgenff ff, must provide --toppar")

        system, str_list = build_system_cgenff(
            replicated_mol_list,
            f"{prefix}.pdb",
            args.toppar,
            use_tip3p=args.use_tip3p,
            )
        for idx, str_ in enumerate(str_list):
            with open(f"{prefix}_monomer{idx}.str", "w") as fopen:
                fopen.write(str_)

        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            monomer_system, _ = build_system_cgenff(
                [rdmol],
                f"{prefix}_monomer{rdmol_idx}.pdb",
                args.toppar,
                use_tip3p=args.use_tip3p,
                )
            monomer_sys_list.append(monomer_system)


    elif args.forcefield.lower() == "oplsaa":

        if args.version == "":
            args.version="CM1A"
        else:
            if not args.version in ["CM1A", "CM1A-LBCC"]:
                raise ValueError(
                    "Version not understood."
                    )

        system = build_system_oplsaa(
            replicated_mol_list,
            f"{prefix}.pdb",
            version=args.version,
            use_tip3p=args.use_tip3p,
            )
        ### Set nonbonded cutoff special for oplsaa
        forces = {system.getForce(index).__class__.__name__: system.getForce(
            index) for index in range(system.getNumForces())}
        nbforce = forces['CustomNonbondedForce']
        nbforce.setCutoffDistance(args.nbcutoff * unit.nanometer)
        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            monomer_sys_list.append(
                build_system_oplsaa(
                    [rdmol],
                    f"{prefix}_monomer{rdmol_idx}.pdb",
                    version=args.version,
                    use_tip3p=args.use_tip3p,
                    )
                )

    else:
        raise ValueError(
            f"force field {args.forcefield} not understood")

    ### Set nonbonded cutoff
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nbforce = forces['NonbondedForce']
    nbforce.setCutoffDistance(args.nbcutoff * unit.nanometer)

    ### Change charge method
    if args.charge_method.lower() != "default":
        from openff.toolkit.topology import Molecule
        partial_charges_list = list()
        for rdmol in rdmol_list_unique:
            offmol = Molecule.from_rdkit(rdmol)
            offmol.assign_partial_charges(
                partial_charge_method=args.charge_method.lower(),
                )
            partial_charges_list.append(
                offmol.partial_charges
                )

        atom_count = 0
        old_q_list = list()
        new_q_list = list()
        for mol_idx in unique_mapping:
            unique_idx = unique_mapping[mol_idx]
            partial_charges = partial_charges_list[unique_idx]
            N_atoms = len(partial_charges)
            for i in range(N_atoms):
                q_old, sig, eps = nbforce.getParticleParameters(atom_count)
                q_new = partial_charges[i]
                old_q_list.append(q_old)
                new_q_list.append(q_new)                
                nbforce.setParticleParameters(
                    atom_count, 
                    q_new, 
                    sig, 
                    eps
                    )
                atom_count += 1

        for i in range(nbforce.getNumExceptions()):
            p1, p2, q_prod, sig, eps = nbforce.getExceptionParameters(i)
            q1 = old_q_list[p1]
            q2 = old_q_list[p2]
            f  = q_prod/(q1*q2)

            q1_new = new_q_list[p1]
            q2_new = new_q_list[p2]

            q_prod = f * q1_new * q2_new

            nbforce.setExceptionParameters(i, p1, p2, q_prod, sig, eps)

        for sys_idx, system in enumerate(monomer_sys_list):
            forces = {system.getForce(index).__class__.__name__: system.getForce(
                index) for index in range(system.getNumForces())}
            nbforce = forces['NonbondedForce']

            old_q_list = list()
            new_q_list = list()
            unique_idx = unique_mapping[sys_idx]
            partial_charges = partial_charges_list[unique_idx]
            N_atoms = len(partial_charges)
            for i in range(N_atoms):
                q_old, sig, eps = nbforce.getParticleParameters(i)
                q_new = partial_charges[i]
                old_q_list.append(q_old)
                new_q_list.append(q_new)                
                nbforce.setParticleParameters(
                    i, 
                    q_new, 
                    sig, 
                    eps
                    )

            for i in range(nbforce.getNumExceptions()):
                p1, p2, q_prod, sig, eps = nbforce.getExceptionParameters(i)
                q1 = old_q_list[p1]
                q2 = old_q_list[p2]
                f  = q_prod/(q1*q2)

                q1_new = new_q_list[p1]
                q2_new = new_q_list[p2]

                q_prod = f * q1_new * q2_new

                nbforce.setExceptionParameters(i, p1, p2, q_prod, sig, eps)

    with open(f"{prefix}.xml", "w") as fopen:
        fopen.write(openmm.XmlSerializer.serialize(system))

    for sys_idx, system in enumerate(monomer_sys_list):
        with open(f"{prefix}_monomer{sys_idx}.xml", "w") as fopen:
            fopen.write(openmm.XmlSerializer.serialize(system))


def entry_point():

    main()

if __name__ == "__main__":

    entry_point()