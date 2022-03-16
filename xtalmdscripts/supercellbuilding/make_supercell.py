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
        '--prefix', 
        "-pre", 
        type=str, 
        help="Prefix for pdb and csv file.", 
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
        default=0,\
        )

    return parser.parse_args()


def get_nonoverlapping_atoms(atom_crds_ortho, filter_overlapping=False):

    """
    Retrieve indices of rows in `atom_crds_ortho` (N,3) list
    that correspond to non-overlapping atoms.
    If `filter_overlapping==True`, then one of the overlapping (not clear 
    which one though) atoms will be retained.

    """

    from scipy.spatial import distance

    atom_crds_ortho_cp = np.copy(atom_crds_ortho)
    dists  = distance.cdist(atom_crds_ortho_cp, atom_crds_ortho_cp)
    np.fill_diagonal(dists, np.inf)
    if filter_overlapping:
        tril_idxs = np.tril_indices(atom_crds_ortho.shape[0])
        dists[tril_idxs] = np.inf
    invalids = np.where(dists < 0.01)[0]
    invalids = np.unique(invalids)
    valids   = np.arange(atom_crds_ortho_cp.shape[0], dtype=int)
    valids   = np.delete(valids, invalids)

    return valids


def random_fill(
    strc, 
    mol_list, 
    N_per_unitcell, 
    radius=1.4, 
    smiles="[H]O[H]"):

    """
    Randomly fill unit cell with molecules.
    """

    from rdkit.Chem import AllChem as Chem
    from rdkit.Geometry import Point3D
    import numpy as np
    from scipy.spatial import distance
    import gemmi
    import copy

    vdw_dict = {
        1  : 1.09,  #H
        6  : 1.7,  #C
        7  : 1.55, #N
        8  : 1.52, #O
        9  : 1.47, #F
        15 : 1.8,  #P
        16 : 1.8,  #S
        17 : 1.75  #Cl
    }

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    Chem.EmbedMolecule(mol)
    Chem.UFFOptimizeMolecule(mol)

    conformer      = mol.GetConformer()
    atom_crds_mol  = conformer.GetPositions()
    atom_crds_mol  = np.array(atom_crds_mol)

    replicated_mol_list, _, _ = make_supercell(
        strc,
        mol_list,
        -1, 1,
        -1, 1,
        -1, 1)

    atom_crds_xtal = list()
    atom_radii_xtal = list()
    for mol_r in replicated_mol_list:
        conf     = mol_r.GetConformer(0)
        conf_pos = conf.GetPositions()
        for atom in mol_r.GetAtoms():
            if atom.GetAtomicNum() == 1:
                continue
            atom_radii_xtal.append(
                vdw_dict[atom.GetAtomicNum()]
            )
            atom_crds_xtal.append(
                conf_pos[atom.GetIdx()]
                )
    atom_crds_xtal = np.array(atom_crds_xtal)
    atom_radii_xtal = np.array(atom_radii_xtal)
    atom_radii_xtal += radius

    grid = list()
    for a in np.linspace(0., 1., 50, True):
        for b in np.linspace(0., 1., 50, True):
            for c in np.linspace(0., 1., 50, True):
                frac  = np.array([a,b,c], dtype=float)
                ortho = strc.cell.orthogonalize(
                    gemmi.Fractional(*frac)
                            ).tolist()
                dists      = distance.cdist(atom_crds_xtal, [ortho])
                is_outside = np.all(dists[:,0] > atom_radii_xtal)
                if is_outside:
                    grid.append(ortho)

    overlap = True
    grid    = np.array(grid)
    print(f"Found {grid.shape[0]} / {50**3} valid grid points.")
    print("Scanning overlap...")
    while overlap:
        mask           = np.arange(grid.shape[0], dtype=int)
        mask_selection = np.random.choice(
            mask, 
            size=N_per_unitcell, 
            replace=False
            )
        grid_selection = grid[mask_selection]
        grid_selection_query = np.copy(grid_selection).tolist()
        for crd in grid_selection:
            frac = strc.cell.fractionalize(
                gemmi.Position(
                    *crd
                    )
                ).tolist()
            frac  = np.array(frac)
            for a in [-1.,0.,1.]:
                for b in [-1.,0.,1.]:
                    for c in [-1.,0.,1.]:
                        if (a==0) & (b==0) & (c==0):
                            continue
                        frac += [a,b,c]
                        ortho = strc.cell.orthogonalize(
                            gemmi.Fractional(*frac)
                                    ).tolist()
                        grid_selection_query.append(ortho)
                        frac -= [a,b,c]

        grid_selection_query = np.array(grid_selection_query)
        dists = distance.cdist(grid_selection_query, grid_selection_query)
        np.fill_diagonal(dists, np.inf)
        min_dist = np.min(dists)
        if min_dist > 2.*radius:
            overlap = False

    for crds in grid_selection:
        mol_cp         = copy.deepcopy(mol)
        conformer      = mol_cp.GetConformer()
        atom_crds_mol  = conformer.GetPositions()
        trans          = crds - np.mean(atom_crds_mol, axis=0)
        atom_crds_mol += trans
        for atm_idx in range(mol_cp.GetNumAtoms()):
            conformer.SetAtomPosition(
                atm_idx,
                Point3D(*atom_crds_mol[atm_idx])
            )
        mol_list.append(mol_cp)

    return 1


def make_P1(strc, addhs=False):

    """
    Generate the P1 cell. Return tuple with atomic coordinates (in Ang) and
    atomic numbers of all atoms in P1 cell.
    """

    import networkx as nx
    from rdkit.Chem import AllChem as Chem
    from rdkit.Geometry import Point3D

    #strc.setup_cell_images()
    atom_crds_ortho = list()
    atom_num = list()

    for site in strc.get_all_unit_cell_sites():
        pos = strc.cell.orthogonalize(site.fract)
        if addhs:
            if int(site.element.atomic_number) != 1:
                atom_crds_ortho.append([pos.x, pos.y, pos.z])
                atom_num.append(site.element.atomic_number)
        else:
            atom_crds_ortho.append([pos.x, pos.y, pos.z])
            atom_num.append(site.element.atomic_number)
    atom_crds_ortho = np.array(atom_crds_ortho, dtype=float)
    atom_num = np.array(atom_num, dtype=int)
    nonoverlapping_idxs = get_nonoverlapping_atoms(atom_crds_ortho, filter_overlapping=True)
    atom_crds_ortho = atom_crds_ortho[nonoverlapping_idxs].tolist()
    atom_num = atom_num[nonoverlapping_idxs].tolist()
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
        acmatrix_new, _ = xyz2mol.xyz2AC(
            atom_num,
            atom_crds_ortho,
            0,
            )
        acmatrix_best = acmatrix_new

        ### Find the disconnected graphs from the adjacency matrix
        G           = nx.convert_matrix.from_numpy_matrix(acmatrix_new)
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
                            if not ortho in atom_crds_ortho_cp:
                                atom_crds_ortho_cp[atm_idx] = ortho
                        acmatrix_new, _ = xyz2mol.xyz2AC(
                            atom_num,
                            atom_crds_ortho,
                            0,
                        )
                        if not np.all(acmatrix_new == acmatrix_best) and np.sum(acmatrix_new) >= np.sum(acmatrix_best):
                            nonoverlapping_idxs = get_nonoverlapping_atoms(atom_crds_ortho_cp)
                            if nonoverlapping_idxs.size == N_atoms:
                                atom_crds_ortho = copy.deepcopy(atom_crds_ortho_cp)
                                acmatrix_best = acmatrix_new
                                found_bond = True

    acmatrix, _ = xyz2mol.xyz2AC(
        atom_num,
        atom_crds_ortho,
        0,
        )
    G           = nx.convert_matrix.from_numpy_matrix(acmatrix)
    G_node_list = list(nx.connected_components(G))
    atom_num    = np.array(atom_num, dtype=int)
    atom_crds_ortho = np.array(atom_crds_ortho)
    mol_list    = list()
    for g in G_node_list:
        g = list(g)
        _, mol = xyz2mol.xyz2AC(
                atom_num[g].tolist(), 
                atom_crds_ortho[g].tolist(),
                0)
        mol_list.append(mol)

    N_mol           = len(mol_list)
    atom_crds_ortho = list()
    atom_num        = list()
    mol_list_new    = list()
    for mol_idx in range(N_mol):
        mol  = mol_list[mol_idx]
        conf = mol.GetConformer()
        conf_pos = conf.GetPositions()
        frac_crds = list()
        for pos in conf_pos:
            frac_crds.append(
                strc.cell.fractionalize(
                    gemmi.Position(
                        *pos
                    )
                ).tolist()
            )
        frac_crds       = np.array(frac_crds)
        valid_atoms     = np.where(
            (frac_crds[:,0] > 0.) * (frac_crds[:,0] < 1.) *\
            (frac_crds[:,1] > 0.) * (frac_crds[:,1] < 1.) *\
            (frac_crds[:,2] > 0.) * (frac_crds[:,2] < 1.)
        )[0]
        ### If no atom is in uc, bring the molecule into unit cell
        if valid_atoms.size == 0:
            is_inside = False
            for a in [0,-1,1]:
                for b in [0,-1,1]:
                    for c in [0,-1,1]:
                        frac_crds += [a,b,c]
                        valid_atoms = np.where(
                            (frac_crds[:,0] > 0.) * (frac_crds[:,0] < 1.) *\
                            (frac_crds[:,1] > 0.) * (frac_crds[:,1] < 1.) *\
                            (frac_crds[:,2] > 0.) * (frac_crds[:,2] < 1.)
                        )[0]
                        if valid_atoms.size != 0:
                            is_inside = True
                        else:
                            frac_crds -= [a,b,c]
                        if is_inside:
                            break
                    if is_inside:
                        break
                if is_inside:
                    break
            if not is_inside:
                continue

        frac_crds_query = np.copy(frac_crds[valid_atoms[0]])
        ### Check for overlap
        overlap = False
        for a in [0,-1,1]:
            for b in [0,-1,1]:
                for c in [0,-1,1]:
                    frac_crds_query += [a,b,c]
                    ortho = strc.cell.orthogonalize(
                        gemmi.Fractional(*frac_crds_query)
                    ).tolist()
                    if len(atom_crds_ortho) > 0:
                        dists = distance.cdist([ortho], atom_crds_ortho)
                        valids = np.where(dists < 0.01)[0]
                        if valids.size > 0:
                            overlap = True
                    frac_crds_query -= [a,b,c]
        if not overlap:
            for atm_idx, frac in enumerate(frac_crds):
                ortho = strc.cell.orthogonalize(
                    gemmi.Fractional(*frac)
                ).tolist()
                atom_crds_ortho.append(ortho)
                rd_atom = mol.GetAtomWithIdx(atm_idx)
                atom_num.append(rd_atom.GetAtomicNum())
                conf.SetAtomPosition(
                    atm_idx,
                    Point3D(*ortho)
                )
            mol_list_new.append(mol)

    mol_list = list()
    if addhs:
        from openeye import oechem
        from openeye import oequacpac
        from xtalmdscripts.supercellbuilding.utils import rdmol_from_oemol
        from xtalmdscripts.supercellbuilding.utils import oemol_from_rdmol

        count = 0
        for mol in mol_list_new:
            oemol = oechem.OEMol()
            oemol.SetDimension(3)
            conf_pos = mol.GetConformer(0).GetPositions()
            crds = list()
            for atm_idx in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(atm_idx)
                oemol.NewAtom(int(atom.GetAtomicNum()))
                crds.extend(conf_pos[atm_idx])
            oemol.SetCoords(crds)

            oechem.OEDetermineConnectivity(oemol)
            oechem.OEFindRingAtomsAndBonds(oemol)
            oechem.OEAssignAromaticFlags(oemol)
            oechem.OEPerceiveBondOrders(oemol)
            oechem.OE3DToInternalStereo(oemol)
            oechem.OEPerceiveChiral(oemol)
            oechem.OEAssignImplicitHydrogens(oemol)
            oechem.OEAssignFormalCharges(oemol)

            oequacpac.OEGetReasonableProtomer(oemol)
            oechem.OEAddExplicitHydrogens(oemol)

            mol = rdmol_from_oemol(oemol)
            Chem.AssignStereochemistryFrom3D(mol)
            mol_list.append(mol)

            #with open(f"./test_{count}.pdb", "w") as fopen:
            #    fopen.write(Chem.MolToPDBBlock(mol))
            count += 1

    else:
        acmatrix, _ = xyz2mol.xyz2AC(
            atom_num,
            atom_crds_ortho,
            0,
            )
        G           = nx.convert_matrix.from_numpy_matrix(acmatrix)
        G_node_list = list(nx.connected_components(G))
        atom_num    = np.array(atom_num, dtype=int)
        atom_crds_ortho = np.array(atom_crds_ortho)
        for g in G_node_list:
            g = list(g)
            mol = Chem.GetMolFrags(
                xyz2mol.xyz2mol(
                    atom_num[g].tolist(), 
                    atom_crds_ortho[g].tolist(),
                    charge=0)[0], 
                asMols=True
            )[0]
            mol_list.append(mol)

    strc_write               = gemmi.Structure()
    strc_write.spacegroup_hm = strc.spacegroup_hm
    strc_write.cell          = strc.cell

    #with open("./make_p1_test.pdb", "w") as fopen:
    #    fopen.write(get_pdb_block(mol_list, strc_write))

    return mol_list


def make_supercell(
    strc,
    mol_list, 
    a_min, a_max,
    b_min, b_max,
    c_min, c_max):

    """
    Generate supercell based specified parameters. Assuming that `atom_crds_ortho`
    and `atom_num` are for P1 cell. See method `make_P1`. Returns list of rdkit mol
    objects with all molecules in supercell, list containing an int that is unique 
    among the molecules in a unit cell, list containing the frac coordinates of the
    unit cell origins in the basis of the supercell.
    """

    a_replicate = np.arange(a_min,a_max+1, dtype=int)
    b_replicate = np.arange(b_min,b_max+1, dtype=int)
    c_replicate = np.arange(c_min,c_max+1, dtype=int)

    N_mol               = len(mol_list)
    replicated_mol_list = list()
    mol_identifies      = list()
    unitcell_in_supercell_fracs = list()
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
                    mol_identifies.append(mol_idx)
                    unitcell_in_supercell_fracs.append([a,b,c])

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

    return replicated_mol_list, mol_identifies, unitcell_in_supercell_fracs


def get_unique_mapping(mol_list, stereochemistry=True):

    """
    Get unique mapping dict and list of unique rdkit mol objects in list of rdkit mol objects.
    if `stereochemistry=True`, the mapping will honor stereochemistry.
    """

    N_mol = len(mol_list)

    smiles_list = [Chem.MolToSmiles(mol, isomericSmiles=stereochemistry) for mol in mol_list]
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
                    unique_mapping[mol_idx] = smiles_unique_idx
                else:
                    unique_mapping[mol_idx] = smiles_unique_idx

    return unique_mapping, rdmol_list_unique


def equalize_rdmols(mol_list, stereochemistry=True):

    """
    Get list of rdkit mol objects in which all chemically idential mol objects
    have identical topology and pdb monomer info. Only difference are coordinates.
    If `stereochemistry=True` it will honor stereochemistry.
    """

    import copy

    unique_mapping, rdmol_list_unique = get_unique_mapping(mol_list, stereochemistry)

    for mol_idx in unique_mapping:
        ### This is the molecule that holds the correct coordinates
        ### and pdb monomer info.
        mol_1 = mol_list[mol_idx]
        ### This is the molecule that holds the correct names, ordering, etc...
        mol_2 = copy.deepcopy(rdmol_list_unique[unique_mapping[mol_idx]])
        match = mol_1.GetSubstructMatch(mol_2, useChirality=stereochemistry)
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
    c_min_max,
    addhs=False,
    addwater=0):

    """
    Generate rdkit mol object list for molecules in supercell. supercell is generated
    according input parameters.
    """

    mol_list = make_P1(strc, addhs)
    if addwater > 0:
        random_fill(
            strc,
            mol_list,
            N_per_unitcell=addwater,
            radius=0.5,
            smiles="O"
            )
    replicated_mol_list, mol_identifies, unitcell_in_supercell_fracs = make_supercell(
        strc,
        mol_list, 
        a_min_max[0], a_min_max[1],
        b_min_max[0], b_min_max[1],
        c_min_max[0], c_min_max[1],
        )

    replicated_mol_list = equalize_rdmols(replicated_mol_list)

    return replicated_mol_list, mol_identifies, unitcell_in_supercell_fracs


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
    header = strc_write.make_pdb_headers()

    ### With the flavor options, one can control what is written
    ### to the pdb block.
    ###
    ### flavor: (optional)
    ### flavor & 1 : Write MODEL/ENDMDL lines around each record
    ### flavor & 2 : Don’t write any CONECT records
    ### flavor & 4 : Write CONECT records in both directions
    ### flavor & 8 : Don’t use multiple CONECTs to encode bond order
    ### flavor & 16 : Write MASTER record
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


def get_supercell_info_str(mol_identifies, unitcell_in_supercell_fracs):

    """
    Takes list of in unitcell molecule identifiers and unitcell in supercell
    frac coordinates. See output generated by method `generate_replicated_mol_list`.
    Returns csv formatted info string.
    """

    info_str  = "#Mol_idx,"
    info_str += "mol_in_unitcell,"
    info_str += "unitcell_in_supercell_a,"
    info_str += "unitcell_in_supercell_b,"
    info_str += "unitcell_in_supercell_c\n"

    N_mols = len(mol_identifies)
    for mol_idx in range(N_mols):
        info_str += f"{mol_idx:d},"
        info_str += f"{mol_identifies[mol_idx]:d},"
        info_str += f"{unitcell_in_supercell_fracs[mol_idx][0]:d},"
        info_str += f"{unitcell_in_supercell_fracs[mol_idx][1]:d},"
        info_str += f"{unitcell_in_supercell_fracs[mol_idx][2]:d}\n"

    return info_str


def get_replicated_mol_list_json(replicated_mol_list):

    """
    Returns rdkit json string of collapsed replicated mol_list.
    """

    from rdkit import Chem

    mol_combo = Chem.Mol()
    for mol in replicated_mol_list:
        mol_combo = Chem.CombineMols(mol_combo, mol)
    return Chem.MolToJSON(mol_combo)


def main():

    """
    Run the workflow.
    """

    import gemmi

    args = parse_arguments()
    strc = parse_cif(args.input)

    ### Build the supercell as a set of rdkit molecule objects
    ### ======================================================
    replicated_mol_list, mol_identifies, unitcell_in_supercell_fracs = generate_replicated_mol_list(
        strc,
        args.a_min_max,
        args.b_min_max,
        args.c_min_max,
        args.addhs,
        args.addwater
        )

    ### Write pdb file
    ### ==============
    strc_write               = gemmi.Structure()
    strc_write.spacegroup_hm = "P1"
    strc_write.cell          = strc.cell
    pdb_str = get_pdb_str(
        replicated_mol_list, 
        strc_write,
        args.a_min_max,
        args.b_min_max,
        args.c_min_max
        )
    with open(f"{args.prefix}.pdb", "w") as fopen:
        fopen.write(pdb_str)
    with open(f"{args.prefix}.csv", "w") as fopen:
        info_str = get_supercell_info_str(
            mol_identifies, 
            unitcell_in_supercell_fracs
            )
        fopen.write(info_str)

    with open(f"{args.prefix}.json", "w") as fopen:
        json_str = get_replicated_mol_list_json(replicated_mol_list)
        fopen.write(json_str)

    ### Generate list of unique smiles for unique
    ### molecules in UC
    ### =========================================
    mol_list = make_P1(strc, args.addhs)
    if args.addwater > 0:
        random_fill(
            strc,
            mol_list,
            N_per_unitcell=args.addwater,
            radius=0.5,
            smiles="O"
            )

    unitcell_mol_list, _, _ = make_supercell(
        strc,
        mol_list,
        0,0,
        0,0,
        0,0,
        )
    from rdkit.Chem import Descriptors
    unitcell_weight = 0.
    smiles_list     = list()
    for mol in unitcell_mol_list:
        unitcell_weight += Descriptors.MolWt(mol)
        smiles_list.append(Chem.MolToSmiles(mol, isomericSmiles=True))

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
Expt:
-----
Temperature [K]               : {cif_info_dict['temperature']},
Cell Setting                  : {cif_info_dict['cell_setting']},
Space Group H-M               : {cif_info_dict['space_group']},
Density [g/cm3]               : {cif_info_dict['density']},

Supercell:
----------
Total number of molecules     : {len(replicated_mol_list)},
Total Length edge a [Ang]     : {strc.cell.a * a_len:4.2f},
Total Length edge b [Ang]     : {strc.cell.b * b_len:4.2f},
Total Length edge c [Ang]     : {strc.cell.c * c_len:4.2f},
Cell angle alpha [deg]        : {strc.cell.alpha:4.2f},
Cell angle beta  [deg]        : {strc.cell.beta:4.2f},
Cell angle gamma [deg]        : {strc.cell.gamma:4.2f},
Total Volume supercell [Ang3] : {strc.cell.volume * a_len * b_len * c_len:4.2f}
Density [g/cm3]               : {unitcell_weight / strc.cell.volume * 1.6605:4.2f}
SMILES for molecules in UC    : {" ".join(smiles_list)}
"""

### 1.6605 conversion g/mol/Ang^3 to g/cm^3
)

def entry_point():

    main()

if __name__ == "__main__":

    entry_point()