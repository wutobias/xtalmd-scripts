#!/use/bin/env python3

from pkg_resources import resource_filename

oplsaa_xml_builder_path = resource_filename("xtalmdscripts.data", "oplsaa/build_xml.sh")

AMBER_PROTEIN_FF  = ["ff14SB", "ff19SB", "fb15", "ff14SBonlysc"]
CHARMM_PROTEIN_FF = ["charmm36", "charmm36m"]
OPLS_PROTEIN_FF   = ["OPLS"]
WATER_FF          = ["TIP3P", "TIP3PFB", "OPC", "NONE"]

ALL_PROTEIN_FF = AMBER_PROTEIN_FF[:] + CHARMM_PROTEIN_FF[:] + OPLS_PROTEIN_FF[:]


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
        default="xtal",
        required=False
        )

    parser.add_argument(
        '--forcefield', 
        "-ff", 
        type=str, 
        help="Force field name. If `offxml` is selected, a offxml file must be provided via --offxml flag.",
        choices=[
            "gaff1",
            "gaff2",
            "parsley",
            "sage",
            "cgenff",
            "oplsaa",
            "offxml",
            "amber",
            "charmm",
            "opls"
            ],
        required=True
        )

    parser.add_argument(
        "--water_model",
        "-wm",
        type=str,
        help="Water model to use in simulation with protein forcefield. Will only affect simulations with protein forcefields.",
        choices=WATER_FF,
        required=False,
        default="NONE",
        )

    parser.add_argument(
        "--rigid_water",
        "-rwa",
        help="Use rigid water molecules.", 
        action='store_true',
        default=False,
        required=False
        )

    parser.add_argument(
        "--offxml",
        "-of",
        type=str,
        required=False,
        )

    parser.add_argument(
        "--use_openfftk",
        "-offtk",
        help="Use OpenFF Toolkit to generate system. Useful for non-standard offxml files.", 
        action='store_true',
        default=False,
        required=False
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
        '--nolongrange', 
        "-nolr", 
        help="Deactivate any analytical treatment of longrange treatment of LJ interactions.",
        action='store_true',
        default=False,
        required=False
        )

    parser.add_argument(
        '--longrange', 
        "-lr", 
        type=str, 
        help="Long range method to use in xtal simulation.", 
        choices=[
            "Ewald",
            "PME",
            "LJPME"
            ],
        default="PME",
        required=False
        )

    parser.add_argument(
        '--switching_function', 
        "-sf", 
        help="Use switching function to let LJ potential smoothly go to zero at cutoff.", 
        action='store_true',
        default=False,
        required=False
        )

    parser.add_argument(
        '--switching_cutoff', 
        "-sc", 
        type=float, 
        help="Switching function cutoff in nm.", 
        default=0.8,
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
        '--removewater', 
        "-rw", 
        action='store_true',
        help="Remove water molecules. This will be executed prior to water placement.", 
        required=False,
        default=False,
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


def rewrite_oplsaa_xml(inpath, outpath, topology, suffix):

    import xml.etree.ElementTree as ET
    
    tree = ET.parse(inpath)
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
                for bond in topology.bonds():
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
                    r2.attrib[attrib]  += f"_{suffix}"
            for r3 in r2:
                for attrib in r3.attrib:
                    if any([attrib.startswith(c) for c in ["class", "name", "type"]]):
                        r3.attrib[attrib]  += f"_{suffix}"

    with open(outpath, "wb") as fopen:
        fopen.write(
            ET.tostring(root)
            )


def remove_water_bonded_forces(system, topology):

    import openmm

    for f_idx, force in enumerate(system.getForces()):
        if isinstance(force, openmm.HarmonicBondForce):
            bond_force = force
            bond_force_idx = f_idx
        elif isinstance(force, openmm.HarmonicAngleForce):
            angle_force = force
            angle_force_idx = f_idx
    atom_list   = list(topology.atoms())
    new_bond_force = openmm.HarmonicBondForce()
    new_angle_force = openmm.HarmonicAngleForce()

    new_bond_force.setForceGroup(bond_force.getForceGroup())
    new_angle_force.setForceGroup(angle_force.getForceGroup())
    for idx in range(bond_force.getNumBonds()):
        parms = bond_force.getBondParameters(idx)
        p1,p2 = parms[:2]
        is_water = False
        for _pidx in [p1,p2]:
            if atom_list[_pidx].residue.name in ["HOH", "WAT", "TIP", "TIP3"]:
                is_water = True
                break
        if not is_water:
            new_bond_force.addBond(*parms)

    for idx in range(angle_force.getNumAngles()):
        parms = angle_force.getAngleParameters(idx)
        p1,p2,p3 = parms[:3]
        is_water = False
        for _pidx in [p1,p2,p3]:
            if atom_list[_pidx].residue.name in ["HOH", "WAT", "TIP", "TIP3"]:
                is_water = True
                break
        if not is_water:
            new_angle_force.addAngle(*parms)

    for force_idx in sorted([bond_force_idx, angle_force_idx], reverse=True):
        system.removeForce(force_idx)
    system.addForce(new_bond_force)
    system.addForce(new_angle_force)

    return system


def label_water(rdmol=None, charmm=False):

    from rdkit import Chem

    has_water = False
    wat_rdmol = Chem.MolFromSmarts(
        "[#1+0X1:1][#8H2+0:2][#1+0X1:3]"
        )
    if rdmol.HasSubstructMatch(wat_rdmol):
        has_water = True
        matches = rdmol.GetSubstructMatches(wat_rdmol)
        for match in matches:
            for wat_idx, rdmol_idx in enumerate(match):
                mi = Chem.AtomPDBResidueInfo()
                mi.SetIsHeteroAtom(True)
                if charmm:
                    mi.SetResidueName("TIP".ljust(3))
                else:
                    mi.SetResidueName("HOH".ljust(3))
                tag = wat_rdmol.GetAtomWithIdx(wat_idx).GetAtomMapNum()
                atom = rdmol.GetAtomWithIdx(rdmol_idx)
                if tag == 1:
                    mi.SetName("H1".ljust(4))
                if tag == 3:
                    mi.SetName("H2".ljust(4))
                if tag == 2:
                    if charmm:
                        mi.SetName("OH2".ljust(4))
                    else:
                        mi.SetName("O".ljust(4))
                atom.SetMonomerInfo(mi)
    return rdmol, has_water


def OPLS_LJ(system, CutoffPeriodic=True):

    """
    Helper function to get the OPLS combination rules.
    See http://zarbi.chem.yale.edu/ligpargen/openMM_tutorial.html

    Scaling factors for 1-4 interactions are 0.5 for both LJ and charge.
    These are assumed to be already correctly scaled when the system
    is parsed.
    """

    import openmm
    import numpy as np

    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    geometric_mixing = openmm.CustomNonbondedForce(
        '4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)')
    if CutoffPeriodic:
        geometric_mixing.setNonbondedMethod(geometric_mixing.CutoffPeriodic)
    else:
        geometric_mixing.setNonbondedMethod(geometric_mixing.NoCutoff)
    geometric_mixing.addPerParticleParameter('sigma')
    geometric_mixing.addPerParticleParameter('epsilon')
    geometric_mixing.setCutoffDistance(nonbonded_force.getCutoffDistance())
    geometric_mixing.setUseLongRangeCorrection(True)
    LJset = {}
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        LJset[index] = (sigma, epsilon)
        geometric_mixing.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(
            index, charge, sigma, epsilon * 0
        )
    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        # ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
        # FORCE
        geometric_mixing.addExclusion(p1, p2)
        if eps._value != 0.0:
            sig14 = np.sqrt(LJset[p1][0] * LJset[p2][0])
            eps14 = np.sqrt(LJset[p1][1] * LJset[p2][1]) * 0.5
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
    system.addForce(geometric_mixing)
    return system


def get_cell_parameters(boxvectors):

    from openmm import unit

    if boxvectors == None:
        a_len = 1. * unit.nanometer
        b_len = 1. * unit.nanometer
        c_len = 1. * unit.nanometer
        alpha = 90. * unit.degree
        beta  = 90. * unit.degree
        gamma = 90. * unit.degree
    else:
        from openmm.app.internal.unitcell import computeLengthsAndAngles
        cell_parameters = computeLengthsAndAngles(boxvectors)
        a_len = cell_parameters[0] * unit.nanometer
        b_len = cell_parameters[1] * unit.nanometer
        c_len = cell_parameters[2] * unit.nanometer
        alpha = cell_parameters[3] * unit.degree
        beta  = cell_parameters[4] * unit.degree
        gamma = cell_parameters[5] * unit.degree

    a_len = a_len.value_in_unit(unit.angstrom)
    b_len = b_len.value_in_unit(unit.angstrom)
    c_len = c_len.value_in_unit(unit.angstrom)
    alpha = alpha.value_in_unit(unit.degree)
    beta  = beta.value_in_unit(unit.degree)
    gamma = gamma.value_in_unit(unit.degree)

    return a_len, b_len, c_len, alpha, beta, gamma


def get_cryst1_header(boxvectors):

    from openmm import app
    from openmm import unit

    if boxvectors == None:
        nonbondedMethod = app.NoCutoff
        a_len = 1. * unit.nanometer
        b_len = 1. * unit.nanometer
        c_len = 1. * unit.nanometer
        alpha = 90. * unit.degree
        beta  = 90. * unit.degree
        gamma = 90. * unit.degree
    else:
        nonbondedMethod = app.NoCutoff
        from openmm.app.internal.unitcell import computeLengthsAndAngles
        cell_parameters = computeLengthsAndAngles(boxvectors)
        a_len = cell_parameters[0] * unit.nanometer
        b_len = cell_parameters[1] * unit.nanometer
        c_len = cell_parameters[2] * unit.nanometer
        alpha = cell_parameters[3] * unit.radian
        beta  = cell_parameters[4] * unit.radian
        gamma = cell_parameters[5] * unit.radian

    a_len = a_len.value_in_unit(unit.angstrom)
    b_len = b_len.value_in_unit(unit.angstrom)
    c_len = c_len.value_in_unit(unit.angstrom)
    alpha = alpha.value_in_unit(unit.degree)
    beta  = beta.value_in_unit(unit.degree)
    gamma = gamma.value_in_unit(unit.degree)
    cryst1_header = f"CRYST1{a_len:9.3f}{b_len:9.3f}{c_len:9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P1"

    return cryst1_header, nonbondedMethod


def get_mapping_dict(replicated_mol_list, replicated_mol_list_new):

    import numpy as np
    from rdkit import Chem
    import copy

    min_trans_len_cutoff = 0.5

    mapping_dict = dict()
    replicated_mol_list_new_crds_dict = dict()
    for mol_original_idx, mol_original in enumerate(replicated_mol_list):
        min_trans_len = 999999999999.
        pos = mol_original.GetConformer(0).GetPositions()
        _crds_o = list()
        for atom in mol_original.GetAtoms():
            if atom.GetAtomicNum() > 1:
                _crds_o.append(
                    pos[atom.GetIdx()].tolist()
                    )
        _crds_o = np.array(_crds_o)
        if _crds_o.ndim > 1:
            com_o = np.mean(_crds_o, axis=0)
        else:
            com_o = _crds_o
        for mol_new_idx, mol_new in enumerate(replicated_mol_list_new):
            if mol_new_idx in replicated_mol_list_new_crds_dict:
                com_n = replicated_mol_list_new_crds_dict[mol_new_idx]
            else:
                pos = mol_new.GetConformer(0).GetPositions()
                _crds_n= list()
                for atom in mol_new.GetAtoms():
                    if atom.GetAtomicNum() > 1:
                        _crds_n.append(
                            pos[atom.GetIdx()].tolist()
                            )
                _crds_n = np.array(_crds_n)
                if _crds_o.ndim > 1:
                    com_n = np.mean(_crds_n, axis=0)
                else:
                    com_n = _crds_n
                replicated_mol_list_new_crds_dict[mol_new_idx] = com_n
            trans_len = np.linalg.norm(com_n-com_o)
            if trans_len < min_trans_len:
                min_trans_len = trans_len
                mapping_dict[mol_original_idx] = mol_new_idx
        if min_trans_len > min_trans_len_cutoff:
            e_str  = f"Best match for {Chem.MolToSmiles(mol_original)} is "
            e_str += f"{Chem.MolToSmiles(replicated_mol_list_new[mapping_dict[mol_original_idx]])} "
            e_str += f"with displacment {min_trans_len}"

            raise ValueError(e_str)

    mapped_idxs = set(list(mapping_dict.values()))
    assert len(mapping_dict) == len(replicated_mol_list)
    assert len(mapping_dict) == len(mapped_idxs)

    return mapping_dict


def topology_to_rdmol(topology, positions=None):

    from rdkit import Chem
    from rdkit import Geometry
    from openmm import unit
    
    rdmol_e  = Chem.EditableMol(Chem.Mol())
    map_dict = dict()
    nonetype_dict = dict()

    for atom in topology.atoms():
        if atom.element == None:
            idx = rdmol_e.AddAtom(Chem.Atom(0))
            map_dict[atom.index] = idx
            for _atom in atom.residue.atoms():
                if _atom.index == atom.index:
                    continue
                nonetype_dict[atom.index] = _atom.index
                break
        else:
            idx = rdmol_e.AddAtom(
                Chem.Atom(
                    atom.element.atomic_number
                    )
                )
            map_dict[atom.index] = idx
        
    bonded_list = list()
    for bond in topology.bonds():
        rdmol_e.AddBond(
            map_dict[bond.atom1.index],
            map_dict[bond.atom2.index],
            Chem.BondType.UNSPECIFIED
        )
        _b = sorted([bond.atom1.index, bond.atom2.index])
        bonded_list.append(_b)
    for idx in nonetype_dict:
        _b = sorted([idx, nonetype_dict[idx]])
        if _b not in bonded_list:
            rdmol_e.AddBond(
                map_dict[idx],
                map_dict[nonetype_dict[idx]],
                Chem.BondType.UNSPECIFIED
            )

    rdmol = rdmol_e.GetMol()
    if not isinstance(positions, type(None)):
        conformer = Chem.Conformer()
        for atom in topology.atoms():
            xyz = positions[atom.index]
            if unit.is_quantity(xyz):
                xyz = xyz.value_in_unit(unit.angstrom)
            conformer.SetAtomPosition(
                map_dict[atom.index], 
                Geometry.Point3D(*xyz)
            )
        rdmol.AddConformer(conformer)

    rdmol_list = Chem.GetMolFrags(rdmol, asMols=True)
        
    return rdmol_list        


def build_system_amber(
    replicated_mol_list,
    prefix="xtal",
    boxvectors=None,
    version="ff14SB",
    water_model="tip3p",
    rigid_water=False,
    ):

    import subprocess
    from openmm import unit
    from openmm import app
    from openmm.app import AmberInpcrdFile, AmberPrmtopFile, PDBFile
    from xtalmdscripts.supercellbuilding import make_supercell

    N_mols = len(replicated_mol_list)
    for rdmol_idx in range(N_mols):
        rdmol, _ = label_water(
            replicated_mol_list[rdmol_idx]
            )
        replicated_mol_list[rdmol_idx] = rdmol

    cryst1_header, nonbondedmethod = get_cryst1_header(boxvectors)
    if boxvectors == None:
        tleap_box_str = ""
    else:
        tleap_box_str = "set mol box { 10. 10. 10. }"

    pdb_block = make_supercell.get_pdb_block(
        replicated_mol_list, 
        header=cryst1_header
        )
    with open(f"{prefix}.pdb", "w") as fopen:
        fopen.write(pdb_block)

    with open("cpptraj.in", "w") as fopen:
        fopen.write(f"""
parm {prefix}.pdb
loadcrd {prefix}.pdb name CRD
prepareforleap crdset CRD name Final pdbout {prefix}-cpptraj.pdb terbymol noh
go
""")

    with open("cpptraj.out", "w") as fopen:
        subprocess.call([
            "cpptraj",
            "-i",
            "cpptraj.in",
            ],
            stdout=fopen
            )

    if water_model.lower() == "tip3p":
        water_leaprc = "leaprc.water.tip3p"
    elif water_model.lower() == "opc":
        water_leaprc = "leaprc.water.opc"
    elif water_model.lower() == "tip3pfb":
        water_leaprc = "leaprc.water.fb3"
    else:
        raise ValueError(
            f"Water model {water_model} not known for this forcefield"
            )

    if version.lower() == "ff14sb":
        protein_leaprc = "leaprc.protein.ff14SB"
    elif version.lower() == "ff19sb":
        protein_leaprc = "leaprc.protein.ff19SB"
    elif version.lower() == "fb15":
        protein_leaprc = "leaprc.protein.fb15"
    elif version.lower() == "ff14sbonlysc":
        protein_leaprc = "leaprc.protein.ff14SBonlysc"
    else:
        raise ValueError(
            f"Protein forcefield {version} not known."
            )

    with open("tleap.in", "w") as fopen:
        fopen.write(f"""
set default nocenter on
source {protein_leaprc}
source {water_leaprc}
mol = loadpdb {prefix}-cpptraj.pdb
{tleap_box_str}
savepdb mol {prefix}-tleap.pdb
saveamberparm mol {prefix}.prmtop {prefix}.inpcrd
quit
""")
    with open("tleap.out", "w") as fopen:
        subprocess.call([
            "tleap",
            "-s",
            "-f",
            "tleap.in",
            ],
            stdout=fopen
            )

    prmtop = AmberPrmtopFile(f'{prefix}.prmtop')
    inpcrd = AmberInpcrdFile(f'{prefix}.inpcrd')

    system = prmtop.createSystem(
        nonbondedMethod=nonbondedmethod,
        constraints=None,
        removeCMMotion=True,
        rigidWater=rigid_water,
    )

    tleap_pdb = ""
    with open(f"{prefix}-tleap.pdb", "r") as fopen:
        for line in fopen:
            if line.startswith("CRYST1"):
                continue
            else:
                tleap_pdb += line
    with open(f"{prefix}.pdb", "w") as fopen:
        fopen.write(cryst1_header + "\n" + tleap_pdb)

    if boxvectors != None:
        prmtop.topology.setPeriodicBoxVectors(
            boxvectors
            )
        system.setDefaultPeriodicBoxVectors(
            *boxvectors
            )

    replicated_mol_list_new = topology_to_rdmol(
        prmtop.topology, 
        inpcrd.getPositions(asNumpy=True)
        )
    assert len(replicated_mol_list) == len(replicated_mol_list_new)

    mapping_dict = get_mapping_dict(
        replicated_mol_list, 
        replicated_mol_list_new
        )

    if rigid_water:
        system = remove_water_bonded_forces(
            system, 
            prmtop.topology
            )

    return system, mapping_dict, replicated_mol_list_new


def build_system_gaff(
    replicated_mol_list,
    prefix="xtal",
    boxvectors=None,
    version="2.11",
    water_model="tip3p",
    rigid_water=False):

    """
    Build openmm system for gaff force field.
    """

    use_water_ff = True
    if water_model.lower() == "none":
        use_water_ff = False

    from rdkit import Chem
    from openff.toolkit.topology import Molecule
    from openmmforcefields.generators import GAFFTemplateGenerator
    from openmm.app import ForceField
    from openmm.app import PDBFile
    from openmm import app
    from xtalmdscripts.supercellbuilding import make_supercell

    offmol_list = list()
    for rdmol_idx, rdmol in enumerate(replicated_mol_list):
        rdmol, has_water = label_water(
            rdmol
            )
        if has_water:
            ### If not use water ff,
            ### parameterize as small molecule
            if not use_water_ff:
                offmol_list.append(
                    Molecule.from_rdkit(rdmol)
                    )
        else:
            offmol_list.append(
                Molecule.from_rdkit(rdmol)
                )

    if use_water_ff:
        if water_model.lower() == "tip3p":
            xml_file_path = "amber/tip3p_standard.xml"
        elif water_model.lower() == "opc":
            xml_file_path = "amber/opc_standard.xml"
        elif water_model.lower() == "tip3pfb":
            xml_file_path = "amber/tip3pfb_standard.xml"
        else:
            raise ValueError(
                f"Water model {water_model.lower()} not known."
                )
        forcefield  = ForceField(xml_file_path)
    else:
        forcefield  = ForceField()

    gaff = GAFFTemplateGenerator(
        molecules=offmol_list, 
        forcefield=f'gaff-{version}'
        )
    forcefield.registerTemplateGenerator(gaff.generator)

    cryst1_header, nonbondedmethod = get_cryst1_header(boxvectors)
    pdb_block = make_supercell.get_pdb_block(
        replicated_mol_list, 
        header=cryst1_header
        )
    with open(f"{prefix}.pdb", "w") as fopen:
        fopen.write(pdb_block)

    pdbfile  = PDBFile(f"{prefix}.pdb")
    topology = pdbfile.getTopology()

    modeller = app.Modeller(
        topology, 
        pdbfile.getPositions()
        )
    modeller.addExtraParticles(forcefield)
    topology  = modeller.getTopology()
    positions = modeller.getPositions()
    with open(f"{prefix}.pdb", "w") as fopen:
        PDBFile.writeFile(
            topology,
            positions,
            fopen
            )

    replicated_mol_list_new = topology_to_rdmol(topology, positions)
    assert len(replicated_mol_list) == len(replicated_mol_list_new)

    mapping_dict = get_mapping_dict(
        replicated_mol_list, 
        replicated_mol_list_new
        )

    system = forcefield.createSystem(
        topology = topology,
        nonbondedMethod=nonbondedmethod,
        constraints=None,
        removeCMMotion=True,
        rigidWater=rigid_water,
    )

    if rigid_water:
        system = remove_water_bonded_forces(
            system, 
            topology
            )

    return system, mapping_dict, replicated_mol_list_new


def build_system_off(
    replicated_mol_list,
    prefix="xtal",
    boxvectors=None,
    version="1.3.1",
    offxml=None,
    water_model="tip3p",
    rigid_water=False,
    use_openfftk=False):

    """
    Build openmm system for openff force field.
    """

    use_water_ff = True
    if water_model.lower() == "none":
        use_water_ff = False

    from openff.toolkit.topology import Molecule
    from openmm.app import PDBFile
    from openmm import app
    from xtalmdscripts.supercellbuilding import make_supercell

    offmol_list = list()
    for rdmol_idx, rdmol in enumerate(replicated_mol_list):
        rdmol, has_water = label_water(
            rdmol
            )
        if has_water:
            ### If not use water ff,
            ### parameterize as small molecule
            if not use_water_ff:
                offmol_list.append(
                    Molecule.from_rdkit(rdmol)
                    )
        else:
            offmol_list.append(
                Molecule.from_rdkit(rdmol)
                )

    cryst1_header, nonbondedmethod = get_cryst1_header(boxvectors)

    pdb_block = make_supercell.get_pdb_block(
        replicated_mol_list, 
        header=cryst1_header
        )
    with open(f"{prefix}.pdb", "w") as fopen:
        fopen.write(pdb_block)

    if offxml != None:
        forcefield_name = offxml
    else:
        forcefield_name = f"openff-{version}"

    if use_openfftk:
        if not rigid_water:
            raise ValueError(
                "When using openff toolkit, rigid_water must be set to `False`"
                )

        from openff.toolkit.typing.engines.smirnoff import ForceField
        from openff.toolkit.topology import Topology

        forcefield_name_list = [forcefield_name]

        pdbfile  = PDBFile(f"{prefix}.pdb")
        topology = pdbfile.getTopology()

        if use_water_ff:
            if water_model.lower() == "tip3p":
                xml_file_path = "tip3p.offxml"
            elif water_model.lower() == "opc":
                xml_file_path = "opc.offxml"
            else:
                raise ValueError(
                    f"Water model {water_model.lower()} not known."
                    )
            forcefield_name_list.append(xml_file_path)

        smi_list = [Chem.MolToSmiles(mol) for mol in replicated_mol_list]

        unique_mol_list = list()
        unique_smi_list = list()
        for mol, smi in zip(replicated_mol_list, smi_list):
            if smi in unique_smi_list:
                continue
            else:
                unique_smi_list.append(smi)
                unique_mol_list.append(mol)    

        off_topology = Topology.from_openmm(
            topology,
            unique_molecules=[Molecule.from_rdkit(mol) for mol in unique_mol_list]
        )

        if boxvectors == None:
            off_topology.box_vectors = [2., 2., 2.] * unit.nanometer
        else:
            off_topology.box_vectors = boxvectors
        forcefield = ForceField(
            forcefield_name_list, 
            load_plugins=True,
            )
        interchange = forcefield.create_interchange(
            off_topology
            )

        system = interchange.to_openmm(
            combine_nonbonded_forces=False, 
            add_constrained_forces=True,
            )
        topology = interchange.to_openmm_topology()
        positions = topology.getPositions()
    else:
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator
        from openmm.app import ForceField

        if use_water_ff:
            if water_model.lower() == "tip3p":
                xml_file_path = "amber/tip3p_standard.xml"
            elif water_model.lower() == "opc":
                xml_file_path = "amber/opc_standard.xml"
            else:
                raise ValueError(
                    f"Water model {water_model.lower()} not known."
                    )
            forcefield  = ForceField(xml_file_path)
        else:
            forcefield  = ForceField()
        openff = SMIRNOFFTemplateGenerator(
            molecules=offmol_list, 
            forcefield=forcefield_name,
            )
        forcefield.registerTemplateGenerator(openff.generator)

        pdbfile  = PDBFile(f"{prefix}.pdb")
        topology = pdbfile.getTopology()

        modeller = app.Modeller(
            topology, 
            pdbfile.getPositions()
            )
        modeller.addExtraParticles(forcefield)
        topology  = modeller.getTopology()
        positions = modeller.getPositions()

        system = forcefield.createSystem(
            topology = topology,
            nonbondedMethod=nonbondedmethod,
            constraints=None,
            removeCMMotion=True,
            rigidWater=rigid_water
        )

    if rigid_water:
        system = remove_water_bonded_forces(
            system, 
            topology
            )

    with open(f"{prefix}.pdb", "w") as fopen:
        PDBFile.writeFile(
            topology,
            positions,
            fopen
            )

    replicated_mol_list_new = topology_to_rdmol(topology, positions)
    assert len(replicated_mol_list) == len(replicated_mol_list_new)

    mapping_dict = get_mapping_dict(
        replicated_mol_list, 
        replicated_mol_list_new
        )

    return system, mapping_dict, replicated_mol_list_new


def build_system_charmm(
    replicated_mol_list,
    prefix="xtal",
    boxvectors=None,
    version="charmm36",
    rigid_water=False,
    toppar_dir_path="./toppar",):

    """
    Build charmm system. Note: We strictly assume that residues
    are already correctly labled.
    """

    cleanup = True

    version = version.lower()

    valid_versions = CHARMM_PROTEIN_FF[:]
    valid_versions.append("cgenff")
    valid_versions.append("opls")

    if not version in valid_versions:
        raise ValueError(
            f"Version with name {version} not supported."
            )

    import os
    import subprocess
    import glob
    import parmed as pmd
    from rdkit import Chem
    from xtalmdscripts.supercellbuilding import make_supercell
    from .utils import get_params_path_list_charmm36, get_charmm_inp_header_charmm36
    from .utils import get_params_path_list_charmm36m, get_charmm_inp_header_charmm36m
    from .utils import get_params_path_list_opls, get_charmm_inp_header_opls

    import openmm
    from openmm.app import ForceField
    from openmm.app import PDBFile, CharmmPsfFile, CharmmParameterSet
    from openmm import app

    cgenff = False
    opls   = False
    ### Stores the names to all charmm str, topology and parameter files.
    if version == "charmm36":
        params_path_list  = get_params_path_list_charmm36(toppar_dir_path)
        charmm_inp_header = get_charmm_inp_header_charmm36(toppar_dir_path)
    elif version == "charmm36m":
        params_path_list  = get_params_path_list_charmm36m(toppar_dir_path)
        charmm_inp_header = get_charmm_inp_header_charmm36m(toppar_dir_path)
    elif version == "cgenff":
        params_path_list  = get_params_path_list_charmm36m(toppar_dir_path)
        charmm_inp_header = get_charmm_inp_header_charmm36m(toppar_dir_path)
        cgenff = True
    elif version == "opls":
        params_path_list  = get_params_path_list_opls(toppar_dir_path)
        charmm_inp_header = get_charmm_inp_header_opls(toppar_dir_path)
        opls = True
    else:
        raise ValueError(
            f"Version with name {version} not supported."
            )

    ### cgenff topology
    rtf_path_cgenff = f"{toppar_dir_path}/top_all36_cgenff.rtf"
    ### cgenff parameters
    prm_path_cgenff = f"{toppar_dir_path}/par_all36_cgenff.prm"

    unique_mapping, rdmol_list_unique = make_supercell.get_unique_mapping(
        replicated_mol_list,
        stereochemistry=False
        )

    to_remove_list = list()
    pmd_list = list()
    failed   = False
    already_processed_list = list()
    for mol_idx in unique_mapping:
        unique_idx = unique_mapping[mol_idx]
        already_processed = False
        if unique_idx in already_processed_list:
            already_processed = True
        else:
            already_processed_list.append(
                unique_idx
                )
        basename_monomer = f"m{unique_idx}"

        rdmol = replicated_mol_list[mol_idx]

        if not cgenff:
            ### We have to check to for presence of ACE or NME first.
            ### We assume that everything is correctly labled. If we
            ### find ACE or NME, then we merge it with its adjacent
            ### residue and apply a patch in charmm.
            ### Also, we want to figure out if N or C terminal ending
            ### is protonated. This will also require different patches.
            has_ACE = False
            has_NME = False
            has_CPROT = False
            has_NPROT2 = False
            has_NPROT3 = False
            is_cyclic = False
            N_term_resname = ""
            C_term_resname = ""
            sequence_dict = dict()
            for atom in rdmol.GetAtoms():
                mi = atom.GetMonomerInfo()
                resid = mi.GetResidueNumber()
                resname = mi.GetResidueName().rstrip().lstrip()
                atomname = mi.GetName().rstrip().lstrip()
                if resname == "ACE":
                    has_ACE = True
                elif resname == "NME":
                    has_NME = True

                if atomname == "HO1":
                    has_CPROT = True
                elif atomname == "HN2":
                    has_NPROT2 = True
                elif atomname == "HN3":
                    has_NPROT3 = True
                sequence_dict[resid] = resname

            if (not has_NPROT2) and (not has_NPROT3) and (not has_CPROT):
                if (not has_ACE) and (not has_NME):
                    is_cyclic = True

            N_term_resid = 1
            C_term_resid = len(sequence_dict)
            if has_ACE or has_NME:
                for atom in rdmol.GetAtoms():
                    mi = atom.GetMonomerInfo()
                    resid = mi.GetResidueNumber()
                    resname = mi.GetResidueName()
                    if resname == "ACE":
                        resname = sequence_dict[resid+1]
                        N_term_resid = 1
                    elif resname == "NME":
                        resname = sequence_dict[resid-1]
                        resid   = resid - 1
                        C_term_resid = C_term_resid - 1
                        if has_ACE:
                            resid = resid - 1
                            C_term_resid = C_term_resid - 1
                    else:
                        if has_ACE:
                            resid = resid - 1
                    mi.SetResidueNumber(resid)
                    mi.SetResidueName(resname)
                    atom.SetMonomerInfo(mi)

            N_term_resname = sequence_dict[N_term_resid]
            C_term_resname = sequence_dict[C_term_resid]

        if cgenff:
            _, has_water = label_water(
                rdmol,
                charmm=True
                )
        else:
            rdmol, has_water = label_water(
                rdmol,
                charmm=True
                )
        with open(f"{basename_monomer}.pdb", "w") as fopen:
            fopen.write(
                Chem.MolToPDBBlock(rdmol)
                )

        if cgenff and not has_water and not already_processed:

            mi = rdmol.GetAtomWithIdx(0).GetMonomerInfo()
            resname = mi.GetResidueName()

            params_path_list.append(f"{basename_monomer}.str")

            ### Convert to mol2 file using obabel
            subprocess.run([
                "obabel",
                "-ipdb",
                f"{basename_monomer}.pdb",
                "--title",
                resname,
                "-omol2",
                f"-O{basename_monomer}.mol2",
                ])

            ### Obtain cgenff stream file using cgenff
            subprocess.run([
                "cgenff",
                "-p",
                rtf_path_cgenff,
                prm_path_cgenff,
                "-b",
                "-f",
                f"{basename_monomer}.str",
                f"{basename_monomer}.mol2"
                ])

            charmm_inp_header += f"\nstream {basename_monomer}.str"

        if has_water:
            chain_name = "TIP3"
            patch_name = ["none", "none"]
            trajout_str = f"""
trajout {basename_monomer}-cpptraj.pdb
trajout {basename_monomer}-cpptraj.cor segmask :HOH,WAT,TIP3 TIP3
trajout {basename_monomer}-{mol_idx:d}-cpptraj.pdb
trajout {basename_monomer}-{mol_idx:d}-cpptraj.cor segmask :HOH,WAT,TIP3 TIP3
"""
        else:
            if cgenff or is_cyclic:
                patch_name = ["none", "none"]
            elif opls:
                patch_name = ["nter", "cter"]
                ### N terminal
                if has_ACE:
                    if N_term_resname == "PRO":
                        patch_name[0] = "acp"
                    else:
                        patch_name[0] = "ace"
                elif has_NPROT3:
                    if N_term_resname == "GLY":
                        patch_name[0] = "glyp"
                    else:
                        patch_name[0] = "nter"
                elif has_NPROT2:    
                    if N_term_resname == "PRO":
                        patch_name[0] = "prop"
                    else:
                        ### It seems like opls does not support
                        ### neutral n terminus?
                        patch_name[0] = "nter"
                else:
                    patch_name[0] = "nter"

                ### C terminal
                if has_NME:
                    patch_name[1] = "ct3"
                elif has_CPROT:
                    if C_term_resname == "GLY":
                        patch_name[1] = "ctpg"
                    else:
                        patch_name[1] = "ctp"
                else:
                    patch_name[1] = "cter"
            else:
                patch_name = ["nter", "cter"]
                
                ### N terminal
                if has_ACE:
                    ### Dipeptide
                    if len(sequence_dict) < 4:
                        if N_term_resname == "PRO":
                            patch_name[0] = "acpd"
                        else:
                            patch_name[0] = "aced"
                    else:
                        if N_term_resname == "PRO":
                            patch_name[0] = "acp"
                        else:
                            patch_name[0] = "ace"
                elif has_NPROT3:
                    if N_term_resname == "GLY":
                        patch_name[0] = "glyp"
                    else:
                        patch_name[0] = "nter"
                elif has_NPROT2:
                    if N_term_resname == "GLY":
                        patch_name[0] = "ngne"
                    elif N_term_resname == "PRO":
                        patch_name[0] = "prop"
                    else:
                        patch_name[0] = "nneu"
                else:
                    patch_name[0] = "nter"

                ### C terminal
                if has_NME:
                    patch_name[1] = "ct3"
                elif has_CPROT:
                    if C_term_resname == "ASP":
                        patch_name[1] = "cneu"
                    else:
                        patch_name[1] = "pcte"
                else:
                    patch_name[1] = "cter"
            chain_name = "PROA"
            trajout_str = f"""
trajout {basename_monomer}-cpptraj.pdb
trajout {basename_monomer}-cpptraj.cor segmask !:HOH,WAT,TIP3 PROA
trajout {basename_monomer}-{mol_idx:d}-cpptraj.pdb
trajout {basename_monomer}-{mol_idx:d}-cpptraj.cor segmask !:HOH,WAT,TIP3 PROA
"""
        ### Convert to regular pdb naming to charmm pdb naming
        ### ${MMTSB}/perl/convpdb.pl
        if "MMTSB" in os.environ:
            convpdb = f"{os.environ['MMTSB']}/perl/convpdb.pl"
        else:
            convpdb = "convpdb.pl"
        
        with open(f"{basename_monomer}-mmtsb.pdb", "w") as fopen:
            subprocess.run([
                convpdb,
                "-out",
                "charmm22",
                f"{basename_monomer}.pdb",
                ],
                stdout=fopen
                )

        with open("cpptraj.in", "w") as fopen:
            fopen.write(f"""
parm {basename_monomer}-mmtsb.pdb
trajin {basename_monomer}-mmtsb.pdb
strip @H=
{trajout_str}
go
""")
        with open("cpptraj.out", "w") as fopen:
            subprocess.call([
                "cpptraj",
                "-i",
                "cpptraj.in"
                ],
                stdout=fopen
                )
        if has_water:
            with open(f"{basename_monomer}-cpptraj.pdb", "r") as fopen:
                instr = fopen.read().replace("TIP ", "TIP3")
            with open(f"{basename_monomer}-cpptraj.pdb", "w") as fopen:
                fopen.write(instr)
            with open(f"{basename_monomer}-cpptraj.cor", "r") as fopen:
                instr = fopen.read().replace("TIP ", "TIP3")
            with open(f"{basename_monomer}-cpptraj.cor", "w") as fopen:
                fopen.write(instr)
            with open(f"{basename_monomer}-{mol_idx:d}-cpptraj.pdb", "r") as fopen:
                instr = fopen.read().replace("TIP ", "TIP3")
            with open(f"{basename_monomer}-{mol_idx:d}-cpptraj.pdb", "w") as fopen:
                fopen.write(instr)
            with open(f"{basename_monomer}-{mol_idx:d}-cpptraj.cor", "r") as fopen:
                instr = fopen.read().replace("TIP ", "TIP3")
            with open(f"{basename_monomer}-{mol_idx:d}-cpptraj.cor", "w") as fopen:
                fopen.write(instr)

        to_remove_list.extend([
             "cpptraj.in",
             "cpptraj.out",
             "charmm.out",
            f"{basename_monomer}.pdb",
            f"{basename_monomer}.mol2",
            f"{basename_monomer}.str",
            f"{basename_monomer}-mmtsb.pdb",
            f"{basename_monomer}-cpptraj.pdb",
            f"{basename_monomer}-cpptraj.cor",
            f"{basename_monomer}-{mol_idx:d}-cpptraj.pdb",
            f"{basename_monomer}-{mol_idx:d}-cpptraj.cor",
            f"{basename_monomer}-charmm.psf",
            f"{basename_monomer}-charmm.crd",
            f"{basename_monomer}-charmm.pdb",
            f"make_psf-{unique_idx}.inp",
            ])

        if already_processed:
            continue
        
        if cgenff:
            resid_str = "1"
        else:
            resid_str = ""

        charmm_inp = f"""
bomlev -2

{charmm_inp_header}

open read card unit 10 name {basename_monomer}-cpptraj.cor
read sequence {chain_name} {resid_str} coor card unit 10 resid
generate {chain_name} setup warn first {patch_name[0]} last {patch_name[1]}

open read unit 10 card name {basename_monomer}-cpptraj.cor
read coor unit 10 card

! The order seems to be important here
ic params
ic generate
ic build

open write unit 10 card name {basename_monomer}-charmm.psf
write psf  unit 10 card

open write card unit 10 name {basename_monomer}-charmm.crd
write coor unit 10 card

open write card unit 10 name {basename_monomer}-charmm.pdb
write coor pdb  unit 10 official

stop
"""
        with open(f"make_psf-{unique_idx}.inp", "w") as fopen:
            fopen.write(charmm_inp)

        ### Obtain psf file using charmm
        subprocess.run([
            "charmm",
            "-i",
            f"make_psf-{unique_idx}.inp",
            "-o",
            "charmm.out"
            ])

        if os.path.exists(f"{basename_monomer}-charmm.psf"):
            try:
                p_psf = pmd.load_file(
                    f"{basename_monomer}-charmm.psf"
                    )
                pmd_list.append(p_psf)
            except:
                failed = True
                break
        else:
            failed = True
            break

    if failed:
        raise RuntimeError(
            f"Could not build structure with {version}."
            )

    ### This is the parmed.Structure instance that will
    ### store all topology information
    psf_pmd = pmd.Structure()
    crd_pmd = pmd.Structure()
    for mol_idx in range(len(replicated_mol_list)):
        unique_idx = unique_mapping[mol_idx]
        basename_monomer = f"m{unique_idx}"
        crd_pmd += pmd.load_file(
            f"{basename_monomer}-{mol_idx:d}-cpptraj.pdb"
            )
        psf_pmd += pmd_list[unique_idx]
    psf_pmd.write_psf(f"strc.psf")
    crd_pmd.save(f"strc.crd", overwrite=True)
    to_remove_list.append("strc.psf")
    to_remove_list.append("strc.crd")
    to_remove_list.append("strc-final.crd")
    to_remove_list.append("strc-final.pdb")
    to_remove_list.append("strc-final.psf")
    to_remove_list.append("make_psf-final.inp")

    charmm_inp = f"""
bomlev -2

{charmm_inp_header}

open read unit 10 card name strc.psf
read psf  unit 10 card

open read unit 10 card name strc.crd
read coor unit 10 card

! The order seems to be important here
! Different than for single molecule case
ic generate
ic params
ic build
hbuild

ic build

open write unit 10 card name strc-final.psf
write psf  unit 10 card

open write card unit 10 name strc-final.crd
write coor unit 10 card

open write card unit 10 name strc-final.pdb
write coor pdb  unit 10 official

stop
"""
    with open("make_psf-final.inp", "w") as fopen:
        fopen.write(charmm_inp)

    ### Obtain psf file using charmm
    subprocess.run([
        "charmm",
        "-i",
        "make_psf-final.inp",
        "-o",
        "charmm.out"
        ])

    psf_pmd.coordinates = pmd.load_file(
        "./strc-final.pdb"
        ).coordinates

    cryst1_header, nonbondedmethod = get_cryst1_header(boxvectors)
    if opls:
        if boxvectors == None:
            geometric_mixing_periodic = False
        else:
            geometric_mixing_periodic = True

    ### Build the omm system
    psffile = CharmmPsfFile("strc-final.psf")
    params  = CharmmParameterSet(*params_path_list)
    psffile.loadParameters(params)

    topology = psffile.topology
    positions = psf_pmd.coordinates[:]

    system = psffile.createSystem(
        params,
        nonbondedMethod=nonbondedmethod, 
        constraints=None,
        removeCMMotion=True,
        rigidWater=rigid_water,
    )

    if opls:
        system = OPLS_LJ(system, geometric_mixing_periodic)

    ### Clean up
    if cleanup:
        for filename in to_remove_list:
            if os.path.exists(filename):
                os.remove(filename)

    if rigid_water:
        system = remove_water_bonded_forces(
            system, psffile.topology
            )

    if boxvectors != None:
        psffile.topology.setPeriodicBoxVectors(
            boxvectors
            )
        system.setDefaultPeriodicBoxVectors(
            *boxvectors
            )

    with open(f"{prefix}.pdb", "w") as fopen:
        PDBFile.writeFile(
            topology,
            positions,
            fopen
            )

    replicated_mol_list_new = topology_to_rdmol(topology, positions)
    assert len(replicated_mol_list) == len(replicated_mol_list_new)

    mapping_dict = get_mapping_dict(
        replicated_mol_list, 
        replicated_mol_list_new
        )

    return system, mapping_dict, replicated_mol_list_new


def build_system_oplsaa(
    replicated_mol_list,
    prefix="xtal",
    boxvectors=None,
    version="CM1A",
    water_model="tip3p",
    rigid_water=False,):

    """
    Build openmm system for opls aa force field.
    """

    use_water_ff = True
    if water_model.lower() == "none":
        use_water_ff = False

    import os
    import subprocess
    from openmm.app import ForceField
    from openmm.app import PDBFile
    from openmm import app
    from rdkit import Chem
    from xtalmdscripts.supercellbuilding import make_supercell

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
    already_processed_list = list()
    xml_file_list = list()
    to_remove_list = list()
    failed = False
    for mol_idx in unique_mapping:
        unique_idx = unique_mapping[mol_idx]
        if unique_idx in already_processed_list:
            continue
        else:
            already_processed_list.append(
                unique_idx
                )
        rdmol = rdmol_list_unique[unique_idx]

        if use_water_ff:
            rdmol, has_water = label_water(
                rdmol,
                charmm=True
                )
        else:
            rdmol, has_water = label_water(
                rdmol,
                )

        basename_monomer = f"mol_{unique_idx}"

        mi = rdmol.GetAtomWithIdx(0).GetMonomerInfo()
        resname = mi.GetResidueName().rstrip().lstrip()

        with open(f"{basename_monomer}.pdb", "w") as fopen:
            fopen.write(
                Chem.MolToPDBBlock(rdmol)
                )
        with Chem.SDWriter(f"{basename_monomer}.mol") as sdwriter:
            sdwriter.write(rdmol)

        if use_water_ff and has_water:
            continue

        subprocess.run([
            oplsaa_xml_builder_path,
            f"{basename_monomer}.pdb",
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
                f"{basename_monomer}.pdb",
                f"{basename_monomer}.mol",
            ])
            break

        to_remove_list.extend([
            f"{resname}.openmm.pdb",
            f"{resname}.openmm.xml",
            f"{basename_monomer}.pdb",
            f"{basename_monomer}.mol",
        ])

        rewrite_oplsaa_xml(
            f"{resname}.openmm.xml", 
            f"{resname}.openmm.xml",
            pdbfile_original.topology,
            unique_idx
            )

    if failed:
        ### Clean up
        for filename in to_remove_list:
            if os.path.exists(filename):
                os.remove(filename)
        raise RuntimeError(
            "Could not build structure with oplsaa."
            )

    cryst1_header, nonbondedmethod = get_cryst1_header(boxvectors)
    if boxvectors == None:
        geometric_mixing_periodic = False
    else:
        geometric_mixing_periodic = True

    pdb_block = make_supercell.get_pdb_block(
        replicated_mol_list, 
        header=cryst1_header
        )
    with open(f"{prefix}.pdb", "w") as fopen:
        fopen.write(pdb_block)

    if use_water_ff:
        if water_model.lower() == "tip3p":
            with open("./tip3p_oplsaa.xml", "w") as fopen:
                fopen.write("""<ForceField>
 <AtomTypes>
  <Type name="tip3p-O" class="tip3p-O" element="O" mass="15.99943"/>
  <Type name="tip3p-H" class="tip3p-H" element="H" mass="1.007947"/>
 </AtomTypes>
 <Residues>
  <Residue name="HOH">
   <Atom name="O" type="tip3p-O" charge="-0.834"/>
   <Atom name="H1" type="tip3p-H" charge="0.417"/>
   <Atom name="H2" type="tip3p-H" charge="0.417"/>
   <Bond atomName1="O" atomName2="H1"/>
   <Bond atomName1="O" atomName2="H2"/>
  </Residue>
 </Residues>
 <HarmonicBondForce>
  <Bond type1="tip3p-O" type2="tip3p-H" length="0.09572" k="462750.4"/>
 </HarmonicBondForce>
 <HarmonicAngleForce>
  <Angle type1="tip3p-H" type2="tip3p-O" type3="tip3p-H" angle="1.82421813418" k="836.8"/>
 </HarmonicAngleForce>
 <NonbondedForce coulomb14scale="0.5" lj14scale="0.5">
  <UseAttributeFromResidue name="charge"/>
  <Atom type="tip3p-O" sigma="0.31507524065751241" epsilon="0.635968"/>
  <Atom type="tip3p-H" sigma="1" epsilon="0"/>
 </NonbondedForce>
</ForceField>
""")
            xml_file_list.append(
                "./tip3p_oplsaa.xml"
                )
        else:
            raise ValueError(
                f"Water model {water_model.lower()} not known."
                )
    forcefield = app.ForceField(*xml_file_list)

    pdbfile  = PDBFile(f"{prefix}.pdb")
    topology = pdbfile.getTopology()
    modeller = app.Modeller(
        topology, 
        pdbfile.getPositions()
        )
    modeller.addExtraParticles(forcefield)
    topology  = modeller.getTopology()
    positions = modeller.getPositions()

    system = forcefield.createSystem(
        topology = topology,
        nonbondedMethod=nonbondedmethod,
        constraints=None,
        removeCMMotion=True,
        rigidWater=rigid_water,
    )

    if rigid_water and has_water:
        system = remove_water_bonded_forces(
            system, 
            topology
            )

    system = OPLS_LJ(system, geometric_mixing_periodic)

    ### Clean up
    for filename in to_remove_list:
        if os.path.exists(filename):
            os.remove(filename)

    with open(f"{prefix}.pdb", "w") as fopen:
        PDBFile.writeFile(
            topology,
            positions,
            fopen
            )

    replicated_mol_list_new = topology_to_rdmol(topology, positions)
    assert len(replicated_mol_list) == len(replicated_mol_list_new)

    mapping_dict = get_mapping_dict(
        replicated_mol_list, 
        replicated_mol_list_new
        )

    return system, mapping_dict, replicated_mol_list_new


def main():

    """
    Build openmm system for given ff and save as xml.
    """

    import gemmi
    from xtalmdscripts.supercellbuilding import make_supercell
    import openmm
    from openmm import unit
    from openmm import app
    from rdkit import Chem
    import numpy as np

    args = parse_arguments()

    AMBER_PROTEIN_FF  = ["ff14SB", "ff19SB", "fb15", "ff14SBonlysc"]
    CHARMM_PROTEIN_FF = ["charmm36", "charmm36m"]
    OPLS_PROTEIN_FF   = ["OPLS"]

    ALL_PROTEIN_FF = AMBER_PROTEIN_FF[:] + CHARMM_PROTEIN_FF[:] + OPLS_PROTEIN_FF[:]

    ALL_PROTEIN_FF    = [ff.lower() for ff in ALL_PROTEIN_FF]
    AMBER_PROTEIN_FF  = [ff.lower() for ff in AMBER_PROTEIN_FF]
    CHARMM_PROTEIN_FF = [ff.lower() for ff in CHARMM_PROTEIN_FF]
    OPLS_PROTEIN_FF   = [ff.lower() for ff in OPLS_PROTEIN_FF]

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

    print(f"Building super cell with {a_len} x {b_len} x {c_len}")

    ### Build the supercell as a list of rdkit molecules
    ### ================================================
    is_protein = False
    if args.forcefield.lower() in ["amber", "opls", "charmm"]:
        is_protein = True

    replicated_mol_list, mol_identifies, unitcell_in_supercell_fracs = make_supercell.generate_replicated_mol_list(
        cell=strc.cell,
        atom_crds_ortho=atom_crds_ortho,
        atom_num=atom_num,
        a_min_max=a_min_max,
        b_min_max=b_min_max,
        c_min_max=c_min_max,
        addhs=args.addhs,
        addwater=args.addwater,
        removewater=args.removewater,
        N_iterations_protonation=args.n_protonation_attempts,
        use_openeye=args.use_openeye,
        label_residues=is_protein
        )

    unique_mapping, rdmol_list_unique = make_supercell.get_unique_mapping(
        replicated_mol_list,
        stereochemistry=False
        )

    a_len = np.max(a_min_max) - np.min(a_min_max) + 1.
    b_len = np.max(b_min_max) - np.min(b_min_max) + 1.
    c_len = np.max(c_min_max) - np.min(c_min_max) + 1.

    from openmm.app.internal.unitcell import computePeriodicBoxVectors
    _boxvectors = computePeriodicBoxVectors(
        strc.cell.a * a_len * unit.angstrom,
        strc.cell.b * b_len * unit.angstrom,
        strc.cell.c * c_len * unit.angstrom,
        strc.cell.alpha * unit.degree,
        strc.cell.beta  * unit.degree,
        strc.cell.gamma * unit.degree,
        )
    boxvectors = np.zeros((3,3), dtype=float) * unit.nanometer
    for i in range(3):
        boxvectors[i,0] = _boxvectors[i][0]
        boxvectors[i,1] = _boxvectors[i][1]
        boxvectors[i,2] = _boxvectors[i][2]

    monomer_sys_list = list()
    if args.forcefield.lower() == "gaff1":

        if args.version == "":
            args.version="1.81"
        else:
            if not args.version in ["1.4", "1.7", "1.8", "1.81"]:
                raise ValueError(
                    "Version not understood."
                    )

        system, mapping_dict, replicated_mol_list_new = build_system_gaff(
            replicated_mol_list,
            prefix=args.prefix,
            boxvectors=boxvectors,
            version=args.version,
            water_model=args.water_model,
            rigid_water=args.rigid_water
            )
        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            system_monomer, _, _ = build_system_gaff(
                [rdmol],
                prefix=f"{args.prefix}_monomer{rdmol_idx}",
                boxvectors=None,
                version=args.version,
                water_model=args.water_model,
                rigid_water=args.rigid_water
                )
            monomer_sys_list.append(system_monomer)

    elif args.forcefield.lower() == "gaff2":

        if args.version == "":
            args.version="2.1"
        else:
            if not args.version in ["2.1", "2.11"]:
                raise ValueError(
                    "Version not understood."
                    )

        system, mapping_dict, replicated_mol_list_new = build_system_gaff(
            replicated_mol_list,
            prefix=args.prefix,
            boxvectors=boxvectors,
            version=args.version,
            water_model=args.water_model,
            rigid_water=args.rigid_water
            )
        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            system_monomer, _, _ = build_system_gaff(
                [rdmol],
                prefix=f"{args.prefix}_monomer{rdmol_idx}",
                boxvectors=None,
                version=args.version,
                water_model=args.water_model,
                rigid_water=args.rigid_water
                )
            monomer_sys_list.append(system_monomer)

    elif args.forcefield.lower() == "parsley":

        if args.version == "":
            args.version="1.3.1"
        else:
            if not args.version in ["1.0.0", "1.0.1", "1.1.0", "1.1.1", "1.2.0", "1.3.0", "1.3.1"]:
                raise ValueError(
                    "Version not understood."
                    )

        system, mapping_dict, replicated_mol_list_new = build_system_off(
            replicated_mol_list,
            prefix=f"{args.prefix}",
            boxvectors=boxvectors,
            version=args.version,
            water_model=args.water_model,
            rigid_water=args.rigid_water,
            use_openfftk=args.use_openfftk
            )
        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            system_monomer, _, _ = build_system_off(
                [rdmol],
                prefix=f"{args.prefix}_monomer{rdmol_idx}",
                boxvectors=None,
                version=args.version,
                water_model=args.water_model,
                rigid_water=args.rigid_water,
                use_openfftk=args.use_openfftk
                )
            monomer_sys_list.append(system_monomer)

    elif args.forcefield.lower() == "sage":

        if args.version == "":
            args.version="2.0.0"
        else:
            if not args.version in ["2.0.0"]:
                raise ValueError(
                    "Version not understood."
                    )

        system, mapping_dict, replicated_mol_list_new = build_system_off(
            replicated_mol_list,
            prefix=f"{args.prefix}",
            boxvectors=boxvectors,
            version=args.version,
            water_model=args.water_model,
            rigid_water=args.rigid_water,
            use_openfftk=args.use_openfftk
            )
        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            system_monomer, _, _ = build_system_off(
                [rdmol],
                prefix=f"{args.prefix}_monomer{rdmol_idx}",
                boxvectors=None,
                version=args.version,
                water_model=args.water_model,
                rigid_water=args.rigid_water,
                use_openfftk=args.use_openfftk
                )
            monomer_sys_list.append(system_monomer)

    elif args.forcefield.lower() == "offxml":

        if not args.offxml:
            raise ValueError(
                "Must also provide path to offxml file."
                )

        system, mapping_dict, replicated_mol_list_new = build_system_off(
            replicated_mol_list,
            prefix=f"{args.prefix}",
            boxvectors=boxvectors,
            offxml=args.offxml,
            water_model=args.water_model,
            rigid_water=args.rigid_water,
            use_openfftk=args.use_openfftk
            )
        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            monomer_system, _, _ = build_system_off(
                [rdmol],
                prefix=f"{args.prefix}_monomer{rdmol_idx}",
                boxvectors=None,
                offxml=args.offxml,
                water_model=args.water_model,
                rigid_water=args.rigid_water,
                use_openfftk=args.use_openfftk
                )
            monomer_sys_list.append(monomer_system)

    elif args.forcefield.lower() == "cgenff":

        if args.toppar == None:
            raise ValueError(f"When building cgenff ff, must provide --toppar")

        if args.water_model.lower() not in ["tip3p"]:
            raise ValueError(
                f"Water model {water_model} not available for cgenff."
                )

        system, mapping_dict, replicated_mol_list_new = build_system_charmm(
            replicated_mol_list,
            prefix=f"{args.prefix}",
            boxvectors=boxvectors,
            version="cgenff",
            rigid_water=args.rigid_water,
            toppar_dir_path=args.toppar,
            )

        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            monomer_system, _, _ = build_system_charmm(
                [rdmol],
                prefix=f"{args.prefix}_monomer{rdmol_idx}",
                boxvectors=boxvectors,
                version="cgenff",
                rigid_water=args.rigid_water,
                toppar_dir_path=args.toppar,
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

        system, mapping_dict, replicated_mol_list_new = build_system_oplsaa(
            replicated_mol_list,
            prefix=args.prefix,
            boxvectors=boxvectors,
            version=args.version,
            water_model=args.water_model,
            rigid_water=args.rigid_water
            )
        ### Set nonbonded cutoff special for oplsaa
        forces = {system.getForce(index).__class__.__name__: system.getForce(
            index) for index in range(system.getNumForces())}
        nbforce = forces['CustomNonbondedForce']
        nbforce.setCutoffDistance(args.nbcutoff * unit.nanometer)
        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            monomer_system, _, _ = build_system_oplsaa(
                [rdmol],
                prefix=f"{args.prefix}_monomer{rdmol_idx}",
                boxvectors=None,
                version=args.version,
                water_model=args.water_model,
                rigid_water=args.rigid_water
                )
            monomer_sys_list.append(monomer_system)

    elif args.forcefield.lower() == "amber":

        if args.water_model.lower() not in ["tip3p", "opc", "tip3pfb"]:
            raise ValueError(
                f"Water model {args.water_model} unknown for this forcefield."
                )

        _AMBER_PROTEIN_FF = [ff.lower() for ff in AMBER_PROTEIN_FF]
        if args.version.lower() not in _AMBER_PROTEIN_FF:
            raise ValueError(
                f"Forcefield version {args.version.lower()} unknown."
                )

        system, mapping_dict, replicated_mol_list_new = build_system_amber(
            replicated_mol_list,
            prefix = f"{args.prefix}",
            boxvectors=boxvectors,
            version=args.version,
            water_model=args.water_model,
            rigid_water=args.rigid_water,
            )
        
        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            monomer_system, _, _ = build_system_amber(
                [rdmol],
                prefix=f"{args.prefix}_monomer{rdmol_idx}",
                boxvectors=None,
                version=args.version,
                water_model=args.water_model,
                rigid_water=args.rigid_water,
                )
            monomer_sys_list.append(
                monomer_system
                )

    elif args.forcefield.lower() == "charmm":

        if args.water_model.lower() not in ["tip3p"]:
            raise ValueError(
                f"Water model {args.water_model} unknown for this forcefield."
                )

        _CHARMM_PROTEIN_FF = [ff.lower() for ff in CHARMM_PROTEIN_FF]
        if args.version.lower() not in _CHARMM_PROTEIN_FF:
            raise ValueError(
                f"Forcefield version {args.version.lower()} unknown."
                )

        system, mapping_dict, replicated_mol_list_new = build_system_charmm(
            replicated_mol_list,
            prefix = f"{args.prefix}",
            boxvectors=boxvectors,
            version=args.version,
            rigid_water=args.rigid_water,
            toppar_dir_path=args.toppar,
            )

        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            monomer_system, _, _ = build_system_charmm(
                [rdmol],
                prefix=f"{args.prefix}_monomer{rdmol_idx}",
                boxvectors=None,
                version=args.version,
                rigid_water=args.rigid_water,
                toppar_dir_path=args.toppar,
                )
            monomer_sys_list.append(monomer_system)

    elif args.forcefield.lower() == "opls":

        if args.water_model.lower() not in ["tip3p"]:
            raise ValueError(
                f"Water model {args.water_model} unknown for this forcefield."
                )

        _OPLS_PROTEIN_FF = [ff.lower() for ff in OPLS_PROTEIN_FF]
        if args.version.lower() not in _OPLS_PROTEIN_FF:
            raise ValueError(
                f"Forcefield version {args.version.lower()} unknown."
                )

        system, mapping_dict, replicated_mol_list_new = build_system_charmm(
            replicated_mol_list,
            prefix = f"{args.prefix}",
            boxvectors=boxvectors,
            version="opls",
            rigid_water=args.rigid_water,
            toppar_dir_path=args.toppar,
            )

        for rdmol_idx, rdmol in enumerate(rdmol_list_unique):
            monomer_system, _, _ = build_system_charmm(
                [rdmol],
                prefix=f"{args.prefix}_monomer{rdmol_idx}",
                boxvectors=None,
                version="opls",
                rigid_water=args.rigid_water,
                toppar_dir_path=args.toppar,
                )
            monomer_sys_list.append(monomer_system)

    else:
        raise ValueError(
            f"force field {args.forcefield} not understood")

    ### Write supercell info and rdkit json here after all the
    ### processing is done. This is necessary since some build system
    ### routines will alter the atom and/or residue order and we
    ### therefore must update the supercell info and rdki object.
    unique_mapping, rdmol_list_unique = make_supercell.get_unique_mapping(
        replicated_mol_list,
        stereochemistry=False
        )

    import copy
    mol_identifies_new = copy.deepcopy(mol_identifies)
    unitcell_in_supercell_fracs_new = copy.deepcopy(unitcell_in_supercell_fracs)
    unique_mapping_new = copy.deepcopy(unique_mapping)
    rdmol_list_unique_new = copy.deepcopy(rdmol_list_unique)
    for old_idx, new_idx in mapping_dict.items():
        mol_identifies_new[new_idx] = mol_identifies[old_idx]
        unitcell_in_supercell_fracs_new[new_idx] = unitcell_in_supercell_fracs[old_idx]
        unique_mapping_new[new_idx] = unique_mapping[old_idx]
        rdmol_list_unique_new[unique_mapping_new[new_idx]] = rdmol_list_unique[unique_mapping[old_idx]]
    mol_identifies = mol_identifies_new
    unitcell_in_supercell_fracs = unitcell_in_supercell_fracs_new
    unique_mapping = unique_mapping_new
    rdmol_list_unique = rdmol_list_unique_new

    with open(f"{args.prefix}.csv", "w") as fopen:
        info_str = make_supercell.get_supercell_info_str(
            mol_identifies, 
            unitcell_in_supercell_fracs
            )
        fopen.write(info_str)
    with open(f"{args.prefix}.json", "w") as fopen:
        json_str = make_supercell.get_replicated_mol_list_json(
            replicated_mol_list_new
            )
        fopen.write(json_str)

    ### Set nonbonded cutoff
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nbforce = forces['NonbondedForce']
    
    ### Sets cutoff distance
    nbforce.setCutoffDistance(
        args.nbcutoff * unit.nanometer
        )

    ### Sets nonbonded method

    opls_ff = False
    if args.forcefield.lower() in ["oplsaa", "opls"]:
        opls_ff = True

    from openmm import app
    nbmethod_dict = {
        "Ewald".lower()             : 3, #app.Ewald
        "PME".lower()               : 4, #app.PME
        "LJPME".lower()             : 5, #app.LJPME
    }
    nbforce.setNonbondedMethod(
        nbmethod_dict[args.longrange.lower()]
        )
    if args.nolongrange:
        nbforce.setNonbondedMethod(
            nbmethod_dict["PME".lower()]
            )
        nbforce.setUseDispersionCorrection(False)
        if opls_ff:
            geometric_mixing = forces["CustomNonbondedForce"]
            geometric_mixing.setUseLongRangeCorrection(False)
    else:
        if args.longrange.lower() == "LJPME".lower():
            nbforce.setUseDispersionCorrection(False)
            ### If we are using oplsaa, we must substract out 
            ### the LB mixing rules interactions.
            if opls_ff:
                geometric_mixing = forces["CustomNonbondedForce"]
                for i in range(geometric_mixing.getNumParticles()):
                    sig, eps = geometric_mixing.getParticleParameters(i)
                    q, _, _  = nbforce.getParticleParameters(i)
                    nbforce.setParticleParameters(
                        i, 
                        q,
                        sig * unit.nanometer, 
                        eps * unit.kilojoule_per_mole
                    )

                geometric_mixing.setNonbondedMethod(2)
                geometric_mixing.setUseLongRangeCorrection(False)
                geometric_mixing.setEnergyFunction("""
Ugeo - Ulb;

Ugeo= 4*epsilon*((sigmageo/r)^12-(sigmageo/r)^6);
Ulb = 4*epsilon*((sigmaari/r)^12-(sigmaari/r)^6);

sigmageo=sqrt(sigma1*sigma2);
sigmaari=(sigma1+sigma2)*0.5;
epsilon=sqrt(epsilon1*epsilon2);
""")

        else:
            nbforce.setUseDispersionCorrection(True)
            if opls_ff:
                geometric_mixing = forces["CustomNonbondedForce"]
                geometric_mixing.setUseLongRangeCorrection(True)

    ### Sets Switching function
    if args.switching_function:
        nbforce.setUseSwitchingFunction(True)
        nbforce.setSwitchingDistance(
            args.switching_cutoff * unit.nanometer
            )
        if opls_ff:
            geometric_mixing = forces["CustomNonbondedForce"]
            geometric_mixing.setUseSwitchingFunction(True)
            geometric_mixing.setSwitchingDistance(
                args.switching_cutoff * unit.nanometer
                )
    else:
        nbforce.setUseSwitchingFunction(False)
        if opls_ff:
            geometric_mixing = forces["CustomNonbondedForce"]
            geometric_mixing.setUseSwitchingFunction(False)

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

        for mol_idx in unique_mapping:
            unique_idx  = unique_mapping[mol_idx]
            mono_system = monomer_sys_list[unique_idx]

            forces = {mono_system.getForce(index).__class__.__name__: mono_system.getForce(
                index) for index in range(mono_system.getNumForces())}
            nbforce = forces['NonbondedForce']

            old_q_list = list()
            new_q_list = list()
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

    with open(f"{args.prefix}.xml", "w") as fopen:
        fopen.write(openmm.XmlSerializer.serialize(system))

    for sys_idx, mono_system in enumerate(monomer_sys_list):
        with open(f"{args.prefix}_monomer{sys_idx}.xml", "w") as fopen:
            fopen.write(openmm.XmlSerializer.serialize(mono_system))


def entry_point():

    main()

if __name__ == "__main__":

    entry_point()