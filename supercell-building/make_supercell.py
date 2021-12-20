#!/usr/bin/env python3

import gemmi
import copy
import numpy as np
from scipy.spatial import distance
import string
from rdkit import Chem
from rdkit.Geometry import Point3D

from xyz2mol import xyz2mol

import argparse

def parse_arguments():

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

    strc.setup_cell_images()
    atom_crds_ortho = list()
    atom_num  = list()

    for site in strc.get_all_unit_cell_sites():
        pos = strc.cell.orthogonalize(site.fract)
        atom_crds_ortho.append([pos.x, pos.y, pos.z])
        atom_num.append(site.element.atomic_number)

    mol_list = Chem.GetMolFrags(
        xyz2mol(
            atom_num,
            atom_crds_ortho)[0], 
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

    a_replicate = np.arange(a_min,a_max+1, dtype=int)
    b_replicate = np.arange(b_min,b_max+1, dtype=int)
    c_replicate = np.arange(c_min,c_max+1, dtype=int)

    mol_list = Chem.GetMolFrags(
        xyz2mol(
            atom_num.tolist(), 
            atom_crds_ortho)[0], 
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
        mi.SetResidueName('MOL')
        mi.SetResidueNumber(mol_identifies[mol_idx])
        mi.SetOccupancy(1.0)
        mi.SetTempFactor(0.0)
        mi.SetChainId(f"{string.ascii_uppercase[mol_identifies[mol_idx]]}")
        
        atom_counts_dict = dict()
        for atom in mol.GetAtoms():
            atomic_num = atom.GetAtomicNum()
            atomic_ele = atom.GetSymbol()
            if not atomic_num in atom_counts_dict:
                atom_counts_dict[atomic_num] = 1
            else:
                atom_counts_dict[atomic_num] += 1
            mi.SetName(f"{atomic_ele[0]}{atom_counts_dict[atomic_num]}".ljust(4))
            atom.SetMonomerInfo(mi)

        replicated_mol_list[mol_idx] = mol

    return replicated_mol_list

def main():

    args = parse_arguments()

    a_len = np.max(args.a_min_max) - np.min(args.a_min_max) + 1.
    b_len = np.max(args.b_min_max) - np.min(args.b_min_max) + 1.
    c_len = np.max(args.c_min_max) - np.min(args.c_min_max) + 1.

    doc  = gemmi.cif.read(args.input)[0]
    strc = gemmi.make_small_structure_from_block(doc)

    atom_crds_ortho, atom_num = make_P1(strc)
    replicated_mol_list = make_supercell(
        strc,
        atom_crds_ortho, 
        atom_num,
        args.a_min_max[0], args.a_min_max[1],
        args.b_min_max[0], args.b_min_max[1],
        args.c_min_max[0], args.c_min_max[1],
        )

    N_mol = len(replicated_mol_list)
    for mol_idx in range(N_mol):
        mol = replicated_mol_list[mol_idx]
        if mol_idx == 0:
            mol_new = mol
        else:
            mol_new = Chem.CombineMols(mol_new, mol)

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

    header     = strc_write.make_pdb_headers()
    crds_block = Chem.MolToPDBBlock(mol_new)
    with open(args.output, "w") as fopen:
        fopen.write(header)
        fopen.write(crds_block)

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

    print(f"""
Summary:
========

Total number of molecules : {len(replicated_mol_list)},
Total Length edge a       : {strc.cell.a * a_len:4.2f},
Total Length edge b       : {strc.cell.b * b_len:4.2f},
Total Length edge c       : {strc.cell.c * c_len:4.2f},
Cell angle alpha          : {strc.cell.alpha:4.2f},
Cell angle beta           : {strc.cell.beta:4.2f},
Cell angle gamma          : {strc.cell.gamma:4.2f},
Total Volume supercell    : {strc.cell.volume * a_len * b_len * c_len:4.2f}
SMILES for molecules in UC: {" ".join(smiles_list)}
"""
)

if __name__ == "__main__":

    main()