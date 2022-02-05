#!/usr/bin/env python3

import gemmi
import copy
import numpy as np
from scipy.spatial import distance
import string
from rdkit import Chem
from rdkit.Geometry import Point3D

from . import xyz2mol

import argparse

def parse_arguments():

    """
    Parse command line arguments.
    """

    parser = argparse.ArgumentParser(
        description="Python script for building supercell that can be used with PBC in MD simulation."
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
        help="Output pdb file", 
        required=True
        )

    parser.add_argument(
        '--a_min_max', 
        "-a", 
        type=int, 
        help="Minimum and maximum unit cell replicates along direction `a`", 
        required=False,
        default=[-1,1],\
        nargs=2
        )

    parser.add_argument(
        '--b_min_max', 
        "-b", 
        type=int, 
        help="Minimum and maximum unit cell replicates along direction `b`", 
        required=False,
        default=[-1,1],\
        nargs=2
        )

    parser.add_argument(
        '--c_min_max', 
        "-c", 
        type=int, 
        help="Minimum and maximum unit cell replicates along direction `c`", 
        required=False,
        default=[-1,1],\
        nargs=2
        )

    return parser.parse_args()


def make_P1(strc):

    """
    Generate the P1 cell. Return tuple with atomic coordinates (in Ang) and
    atomic numbers of all atoms in P1 cell.
    """

    import networkx as nx

    strc.setup_cell_images()
    atom_crds_ortho = list()
    atom_num  = list()

    for site in strc.get_all_unit_cell_sites():
        pos = strc.cell.orthogonalize(site.fract)
        atom_crds_ortho.append([pos.x, pos.y, pos.z])
        atom_num.append(site.element.atomic_number)
    N_atoms = len(atom_num)

    found_bond = True
    ### Terminate if we haven't found any bonds.
    while found_bond:
        ### Update the frac coordinates
        atom_crds_frac = list()
        for atm_idx in range(N_atoms):
            frac = strc.cell.fractionalize(
                gemmi.Position(
                    *atom_crds_ortho[atm_idx]
                    )
                )
            atom_crds_frac.append(frac.tolist())

        ### Get the "pre-molecules" (adjacency matrix with not much chemistry)
        acmatrix, mol = xyz2mol.xyz2AC(
            atom_num,
            atom_crds_ortho,
            0,
            )
        ### Twice the total number of bonds in the pre-molecules
        n_bond2_new  = np.sum(acmatrix)
        n_bond2_best = n_bond2_new

        ### Find the disconnected graphs from the adjacency matrix
        G           = nx.convert_matrix.from_numpy_matrix(acmatrix)
        G_node_list = list(nx.connected_components(G))
        ### Translate molecules to neighboring unit cells in + direction
        ### and check if we can form new bonds. If yes, update `atom_crds_ortho`
        found_bond = False
        for g in G_node_list:
            for a in [0,-1]:
                for b in [0,-1]:
                    for c in [0,-1]:
                        atom_crds_ortho_cp = copy.deepcopy(atom_crds_ortho)
                        for atm_idx in g:
                            frac = copy.deepcopy(atom_crds_frac[atm_idx])
                            frac[0] += a
                            frac[1] += b
                            frac[2] += c
                            ortho = strc.cell.orthogonalize(
                                gemmi.Fractional(*frac)
                            ).tolist()
                            atom_crds_ortho_cp[atm_idx] = ortho
                        acmatrix, mol = xyz2mol.xyz2AC(
                            atom_num,
                            atom_crds_ortho_cp,
                            0,
                            )
                        n_bond2_new = np.sum(acmatrix)
                        if n_bond2_new > n_bond2_best:
                            atom_crds_ortho = copy.deepcopy(atom_crds_ortho_cp)
                            n_bond2_best = n_bond2_new
                            found_bond = True

    mol_list = Chem.GetMolFrags(
        xyz2mol.xyz2mol(
            atom_num,
            atom_crds_ortho,
            charge=0,
            allow_charged_fragments=True
            )[0], 
        asMols=True
    )
    N_mol           = len(mol_list)
    atom_crds_ortho = list()
    atom_num        = list()
    for mol_idx in range(N_mol):
        mol  = mol_list[mol_idx]
        conf_pos = mol.GetConformer().GetPositions()
        frac_crds = list()
        for pos in conf_pos:
            frac_crds.append(
                strc.cell.fractionalize(
                    gemmi.Position(
                        *pos
                    )
                ).tolist()
            )
        frac_crds = np.array(frac_crds)
        for i in [-1,0,1]:
            for j in [-1,0,1]:
                for k in [-1,0,1]:
                    frac_crds_cp = np.copy(frac_crds)
                    
                    frac_crds_query = np.copy(frac_crds_cp[0])
                    frac_crds_query[0] += i
                    frac_crds_query[1] += j
                    frac_crds_query[2] += k
                    
                    overlap = False
                    
                    if np.all((0. < frac_crds_query) * (1. > frac_crds_query)):
                        
                        for a in [0,-1]:
                            for b in [0,-1]:
                                for c in [0,-1]:
                                    frac_crds_query[0] += a
                                    frac_crds_query[1] += b
                                    frac_crds_query[2] += c

                                    ortho = strc.cell.orthogonalize(
                                        gemmi.Fractional(*frac_crds_query)
                                    ).tolist()

                                    if len(atom_crds_ortho) > 0:
                                        dists = distance.cdist([ortho], atom_crds_ortho)
                                        valids = np.where(dists < 0.01)[0]
                                        if valids.size > 0:
                                            overlap = True

                                    frac_crds_query[0] -= a
                                    frac_crds_query[1] -= b
                                    frac_crds_query[2] -= c

                        frac_crds_cp[:,0] += i
                        frac_crds_cp[:,1] += j
                        frac_crds_cp[:,2] += k

                        if not overlap:
                            ortho_crds = list()
                            for atm_idx, frac in enumerate(frac_crds_cp):
                                ortho = strc.cell.orthogonalize(
                                    gemmi.Fractional(*frac)
                                ).tolist()

                                atom_crds_ortho.append(ortho)
                                rd_atom = mol.GetAtomWithIdx(atm_idx)
                                atom_num.append(rd_atom.GetAtomicNum())

    atom_crds_ortho = np.array(atom_crds_ortho)
    atom_num        = np.array(atom_num)

    return atom_crds_ortho, atom_num


def make_supercell(
    strc,
    atom_crds_ortho, 
    atom_num,
    a_min, a_max,
    b_min, b_max,
    c_min, c_max):

    """
    Generate supercell based specified parameters. Assuming that `atom_crds_ortho`
    and `atom_num` are for P1 cell. See method `make_P1`. Returns list of rdkit mol
    objects with all molecules in supercell.
    """

    a_replicate = np.arange(a_min,a_max+1, dtype=int)
    b_replicate = np.arange(b_min,b_max+1, dtype=int)
    c_replicate = np.arange(c_min,c_max+1, dtype=int)

    mol_list = Chem.GetMolFrags(
        xyz2mol.xyz2mol(
            atom_num.tolist(), 
            atom_crds_ortho,
            charge=0,
            allow_charged_fragments=True)[0], 
        asMols=True
    )

    N_mol               = len(mol_list)
    replicated_mol_list = list()
    all_atom_positions  = list()
    unitcell_identifier = list()
    mol_identifies      = list()
    unit_cell_count     = 0
    for a in a_replicate:
        for b in b_replicate:
            for c in c_replicate:
                for mol_idx in range(N_mol):
                    mol      = copy.deepcopy(mol_list[mol_idx])
                    conf     = mol.GetConformer(0)
                    conf_pos = conf.GetPositions()
                    N_atoms  = mol.GetNumAtoms()
                    
                    for atom_idx in range(N_atoms):
                        pos      = conf_pos[atom_idx]                    
                        frac_pos = strc.cell.fractionalize(
                            gemmi.Position(*pos)
                        )
                        frac_pos.x += a
                        frac_pos.y += b
                        frac_pos.z += c
                        
                        conf_pos_abc = strc.cell.orthogonalize(frac_pos).tolist()
                            
                        conf.SetAtomPosition(
                            atom_idx,
                            Point3D(*conf_pos_abc)
                        )

                    replicated_mol_list.append(mol)
                    conf = mol.GetConformer(0)
                    all_atom_positions.extend(
                        conf.GetPositions().tolist()
                    )
                    mol_identifies.append(mol_idx)
                    unitcell_identifier.append(unit_cell_count)
                unit_cell_count += 1

    N_mol = len(replicated_mol_list)
    for mol_idx in range(N_mol):
        mol = replicated_mol_list[mol_idx]
        mi  = Chem.AtomPDBResidueInfo()
        mi.SetIsHeteroAtom(True)
        mi.SetResidueName(f'M{mol_identifies[mol_idx]}'.ljust(3))
        mi.SetResidueNumber(mol_idx + 1)
        mi.SetOccupancy(1.0)
        mi.SetTempFactor(0.0)
        #mi.SetChainId(f"{string.ascii_uppercase[mol_identifies[mol_idx]]}")
        
        atom_counts_dict = dict()
        for atom in mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            atomic_ele = atom.GetSymbol()
            if not atomic_num in atom_counts_dict:
                atom_counts_dict[atomic_num] = 1
            else:
                atom_counts_dict[atomic_num] += 1
            mi.SetName(
                f"{atomic_ele[0]}{atom_counts_dict[atomic_num]}".ljust(4)
                )
            atom.SetMonomerInfo(mi)

        replicated_mol_list[mol_idx] = mol

    return replicated_mol_list


def get_unique_mapping(mol_list):

    """
    Get unique mapping dict and list of unique rdkit mol objects in list of rdkit mol objects.
    """

    N_mol = len(mol_list)

    smiles_list = [Chem.MolToSmiles(mol, isomericSmiles=True) for mol in mol_list]
    smiles_list_unique = set(smiles_list)
    smiles_list_unique = list(smiles_list_unique)

    rdmol_list_unique  = list()
    unique_mapping     = dict()
    for smiles_unique_idx, smiles_unique in enumerate(smiles_list_unique):
        found_unique = False
        for mol_idx in range(N_mol):
            mol    = mol_list[mol_idx]
            smiles = smiles_list[mol_idx]
            if smiles == smiles_unique:
                if not found_unique:
                    rdmol_list_unique.append(mol)
                    found_unique = True
                else:
                    unique_mapping[mol_idx] = smiles_unique_idx

    return unique_mapping, rdmol_list_unique


def equalize_rdmols(mol_list):

    """
    Get list of rdkit mol objects in which all chemically idential mol objects
    have identical topology and pdb monomer info. Only difference are coordinates.
    """

    import copy

    unique_mapping, rdmol_list_unique = get_unique_mapping(mol_list)

    for mol_idx in unique_mapping:
        ### This is the molecule that holds the correct coordinates
        ### and pdb monomer info.
        mol_1 = mol_list[mol_idx]
        ### This is the molecule that holds the correct names, ordering, etc...
        mol_2 = copy.deepcopy(rdmol_list_unique[unique_mapping[mol_idx]])
        match = mol_1.GetSubstructMatch(mol_2)
        conf_1     = mol_1.GetConformer(0)
        conf_pos_1 = conf_1.GetPositions()
        conf_2     = mol_2.GetConformer(0)
        conf_pos_2 = conf_2.GetPositions()
        for mol2_idx, mol1_idx in enumerate(match):

            pos = conf_pos_1[mol1_idx]
            conf_2.SetAtomPosition(
                mol2_idx,
                Point3D(*pos)
            )

        N_atoms = mol_1.GetNumAtoms()
        for atm_idx in range(N_atoms):
            atom   = mol_2.GetAtomWithIdx(atm_idx)
            ### Note, we cannot `copy.copy(mi_original)`
            ### or `copy.copy(mol_1.GetAtomWithIdx(atm_idx))`
            mi = mol_1.GetAtomWithIdx(atm_idx).GetMonomerInfo()
            mi.SetResidueName(f'M{unique_mapping[mol_idx]}'.ljust(3))
            atom.SetMonomerInfo(mi)

        mol_list[mol_idx] = mol_2

    return mol_list


def generate_replicated_mol_list(
    strc, 
    a_min_max,
    b_min_max,
    c_min_max):

    """
    Generate rdkit mol object list for molecules in supercell. supercell is generated
    according input parameters.
    """

    atom_crds_ortho, atom_num = make_P1(strc)
    replicated_mol_list = make_supercell(
        strc,
        atom_crds_ortho, 
        atom_num,
        a_min_max[0], a_min_max[1],
        b_min_max[0], b_min_max[1],
        c_min_max[0], c_min_max[1],
        )

    replicated_mol_list = equalize_rdmols(replicated_mol_list)

    return replicated_mol_list    


def get_pdb_block(
    replicated_mol_list, 
    strc_write):

    """
    Get pdb block as str. strc_write is gemmi structure object and must reflect
    the dimensions of the supercell.
    """

    ### Combine all rdmols in a single big rdmol
    N_mol   = len(replicated_mol_list)
    mol_new = Chem.Mol()
    for mol_idx in range(N_mol):
        mol = copy.deepcopy(replicated_mol_list[mol_idx])
        mol_new = Chem.CombineMols(mol_new, mol)

    header     = strc_write.make_pdb_headers()

    ### With the flavor options, one can control what is written
    ### to the pdb block.
    ###
    ### flavor: (optional)
    ### flavor & 1 : Write MODEL/ENDMDL lines around each record
    ### 
    ### flavor & 2 : Don’t write any CONECT records
    ### 
    ### flavor & 4 : Write CONECT records in both directions
    ### 
    ### flavor & 8 : Don’t use multiple CONECTs to encode bond order
    ### 
    ### flavor & 16 : Write MASTER record
    ### 
    ### flavor & 32 : Write TER record

    crds_block = Chem.MolToPDBBlock(mol_new, flavor=32)
    pdb_block  = header + crds_block

    return pdb_block


def get_pdb_str(
    replicated_mol_list,
    strc,
    a_min_max,
    b_min_max,
    c_min_max):

    """
    Get full pdb file as str.
    """

    import gemmi

    ### Write pdb file
    ### ==============
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

    pdb_block = get_pdb_block(replicated_mol_list, strc_write)

    return pdb_block
    

def parse_cif(cif_path):

    """
    Parse cif file as gemmis structure object.
    """

    import gemmi

    doc  = gemmi.cif.read(cif_path)[0]
    strc = gemmi.make_small_structure_from_block(doc)

    return strc


def main():

    """
    Run the workflow.
    """

    args = parse_arguments()
    strc = parse_cif(args.input)

    ### Build the supercell as a set of rdkit molecule objects
    ### ======================================================
    replicated_mol_list = generate_replicated_mol_list(
        strc,
        args.a_min_max,
        args.b_min_max,
        args.c_min_max
        )

    ### Write pdb file
    ### ==============
    pdb_str = get_pdb_str(
        replicated_mol_list, 
        strc,
        args.a_min_max,
        args.b_min_max,
        args.c_min_max
        )
    with open(args.output, "w") as fopen:
        fopen.write(pdb_str)

    ### Generate list of unique smiles for unique
    ### molecules in UC
    ### =========================================
    atom_crds_ortho, atom_num = make_P1(strc)
    unitcell_mol_list = make_supercell(
        strc,
        atom_crds_ortho, 
        atom_num,
        0,0,
        0,0,
        0,0,
        )

    smiles_list = [Chem.MolToSmiles(mol, isomericSmiles=True) for mol in unitcell_mol_list]
    smiles_list = set(smiles_list)
    smiles_list = list(smiles_list)

    ### Output final summary
    ### ====================
    a_len = np.max(args.a_min_max) - np.min(args.a_min_max) + 1.
    b_len = np.max(args.b_min_max) - np.min(args.b_min_max) + 1.
    c_len = np.max(args.c_min_max) - np.min(args.c_min_max) + 1.

    import gemmi
    doc           = gemmi.cif.read(args.input)[0]
    cif_info_dict = {
        "temperature"  : "Not found",
        "cell_setting" : "Not found",
        "space_group"  : "Not found",
        "density"      : "Not found"
    }
    for item in doc:
        if item.pair == None:
            continue
        key, value = item.pair
        if "_diffrn_ambient_temperature".lower() == key.lower():
            cif_info_dict["temperature"] = value
        elif "_symmetry_cell_setting".lower() == key.lower():
            cif_info_dict["cell_setting"] = value
        elif "_symmetry_space_group_name_H-M".lower() == key.lower():
            cif_info_dict["space_group"] = value
        elif "_exptl_crystal_density_diffrn".lower() == key.lower():
            cif_info_dict["density"] = value

    print(f"""
Summary:
========
Temperature [K]               : {cif_info_dict['temperature']},
Cell Setting                  : {cif_info_dict['cell_setting']},
Space Group H-M               : {cif_info_dict['space_group']},
Density [g/cm3]               : {cif_info_dict['density']},

Total number of molecules     : {len(replicated_mol_list)},
Total Length edge a [Ang]     : {strc.cell.a * a_len:4.2f},
Total Length edge b [Ang]     : {strc.cell.b * b_len:4.2f},
Total Length edge c [Ang]     : {strc.cell.c * c_len:4.2f},
Cell angle alpha [deg]        : {strc.cell.alpha:4.2f},
Cell angle beta  [deg]        : {strc.cell.beta:4.2f},
Cell angle gamma [deg]        : {strc.cell.gamma:4.2f},
Total Volume supercell [Ang3] : {strc.cell.volume * a_len * b_len * c_len:4.2f}
SMILES for molecules in UC    : {" ".join(smiles_list)}
"""
)

def entry_point():

    main()

if __name__ == "__main__":

    entry_point()