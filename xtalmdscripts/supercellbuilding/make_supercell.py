#!/usr/bin/env python3

import gemmi
import copy
import numpy as np
from scipy.spatial import distance
import string
from rdkit import Chem
from rdkit.Geometry import Point3D
import warnings

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
        type=float, 
        help="If this is two arguments: Minimum and maximum unit cell replicates along direction `a`. \
If this is a single argument: Minimum length of axis `a` in the supercell. Units [nm]", 
        required=False,
        default=[-1,1],
        nargs='+',
        )

    parser.add_argument(
        '--b_min_max', 
        "-b", 
        type=float, 
        help="If this is two arguments: Minimum and maximum unit cell replicates along direction `b`. \
If this is a single argument: Minimum length of axis `b` in the supercell. Units [nm]", 
        required=False,
        default=[-1,1],
        nargs='+',
        )

    parser.add_argument(
        '--c_min_max', 
        "-c", 
        type=float, 
        help="If this is two arguments: Minimum and maximum unit cell replicates along direction `c`. \
If this is a single argument: Minimum length of axis `c` in the supercell. Units [nm]", 
        required=False,
        default=[-1,1],
        nargs='+',
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

    return parser.parse_args()


def combine_mols(mol_list):

    mol_list  = copy.deepcopy(mol_list)
    N_mol_per_unitcell = len(mol_list)
    mol_combo = Chem.Mol()
    for mol in mol_list:
        mol_combo = Chem.CombineMols(mol_combo, mol)
    return mol_combo


class FfWrapper(object):

    def __init__(self, mol, only_hydrogen=True):

        from rdkit.Chem import AllChem as Chem

        self._mol = copy.deepcopy(mol)

        mp  = Chem.MMFFGetMoleculeProperties(mol)
        self._ffm = Chem.MMFFGetMoleculeForceField(mol, mp)

        self._H_atom_list   = list()
        self._all_atom_list = list()
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                self._H_atom_list.append(atom.GetIdx())
            self._all_atom_list.append(atom.GetIdx())

        self._H_atom_list = np.array(self._H_atom_list)
        self._H_idxs = np.arange(mol.GetNumAtoms()*3, dtype=int)
        self._H_idxs = self._H_idxs.reshape((mol.GetNumAtoms(), 3))
        self._H_idxs = self._H_idxs[self._H_atom_list]
        self._H_idxs_flat = self._H_idxs.flatten()

        self._all_atom_list = np.array(self._all_atom_list)
        self._all_idxs = np.arange(mol.GetNumAtoms()*3, dtype=int)
        self._all_idxs = self._all_idxs.reshape((mol.GetNumAtoms(), 3))
        self._all_idxs_flat = self._H_idxs.flatten()

        self._num_all_atoms = len(self._all_atom_list)
        self._num_H_atoms   = len(self._H_atom_list)

        self.only_hydrogen = only_hydrogen

    @property
    def num_all_atoms(self):
        return self._num_all_atoms

    @property
    def num_H_atoms(self):
        return self._num_H_atoms

    @property
    def H_atom_list(self):
        return self._H_atom_list

    @property
    def H_idxs(self):
        return self._H_idxs

    @property
    def H_idxs_flat(self):
        return self._H_idxs_flat

    @property
    def all_atom_list(self):
        return self._all_atom_list

    @property
    def all_idxs(self):
        return self._all_idxs

    @property
    def all_idxs_flat(self):
        return self._all_idxs_flat

    @property
    def ffm(self):
        return self._ffm

    @property
    def mol(self):
        return self._mol

    @property
    def pos(self):
        pos = np.array(self.ffm.Positions())
        if self.only_hydrogen:
            pos = pos[self.H_idxs_flat].tolist()
        return pos

    @pos.setter
    def pos(self, pos):

        from rdkit.Chem import AllChem as Chem

        conformer      = self._mol.GetConformer()

        if self.only_hydrogen:
            pos_stack = np.array(pos).reshape((self.num_H_atoms, 3))
            for H_atm_idx, atm_idx in enumerate(self.H_atom_list):
                conformer.SetAtomPosition(
                    int(atm_idx),
                    Point3D(
                        *pos_stack[H_atm_idx]
                        )
                    )
        else:
            pos_stack = np.array(pos).reshape((self.num_all_atoms, 3))
            for atm_idx in range(self._mol.GetNumAtoms()):
                conformer.SetAtomPosition(
                    int(atm_idx),
                    Point3D(
                        *pos_stack[atm_idx]
                        )
                    )

        mp  = Chem.MMFFGetMoleculeProperties(self._mol)
        self._ffm = Chem.MMFFGetMoleculeForceField(self._mol, mp)

    def ene(self,x):
        pos = self.ffm.Positions()
        pos = np.array(pos)
        if self.only_hydrogen:
            pos[self.H_idxs_flat] = x
        else:
            pos[:] = x
        pos = pos.tolist()
        return self.ffm.CalcEnergy(pos)

    def grad(self,x):
        pos = self.ffm.Positions()
        pos = np.array(pos)
        if self.only_hydrogen:
            pos[self.H_idxs_flat] = x
        else:
            pos[:] = x
        pos = pos.tolist()
        grad_vals = self.ffm.CalcGrad(pos)
        grad_vals = np.array(grad_vals)[self.H_idxs_flat]
        grad_vals = grad_vals.tolist()
        return grad_vals


def apply_delta_hydrogen(
    mol_list,
    delta_hydrogen_dict
    ):

    import numpy as np
    from rdkit.Chem import AllChem as Chem
    import copy

    mol_combo   = combine_mols(mol_list)
    emol        = Chem.EditableMol(mol_combo)
    offset_list = np.zeros(mol_combo.GetNumAtoms(), dtype=int)
    for atomidx, switch in delta_hydrogen_dict.items():
        atom    = mol_combo.GetAtomWithIdx(atomidx)
        charge  = atom.GetFormalCharge()
        if switch == 1:
            Hatom = Chem.Atom(1)
            idx1  = emol.AddAtom(Hatom)
            Patom = Chem.AtomFromSmarts(atom.GetSmarts())
            Patom.SetFormalCharge(charge+1)
            emol.ReplaceAtom(
                atomidx+int(offset_list[atomidx]), 
                Patom
                )
            emol.AddBond(
                idx1,
                atomidx+int(offset_list[atomidx]),
                Chem.BondType.SINGLE
            )
            
        elif switch == -1:
            for n_atom in atom.GetNeighbors():
                if n_atom.GetAtomicNum() == 1:
                    idx1  = n_atom.GetIdx()
                    Patom = Chem.AtomFromSmarts(atom.GetSmarts())
                    Patom.SetFormalCharge(charge-1)
                    emol.ReplaceAtom(
                        atomidx+int(offset_list[atomidx]),
                        Patom
                        )
                    emol.RemoveAtom(
                        idx1+int(offset_list[idx1])
                        )
                    offset_list[idx1:] -= 1
                    break
        else:
            ### Nix.
            pass
    mol_combo_prot = emol.GetMol()
    Chem.SanitizeMol(mol_combo_prot)
    return Chem.GetMolFrags(mol_combo_prot, asMols=True)
    

def minimize_H(
    mol_list,
    sanitize=True,
    ):

    import numpy as np
    from rdkit.Chem import AllChem as Chem
    from scipy import optimize
    import copy

    mol_combo = combine_mols(mol_list)
    if sanitize:
        Chem.SanitizeMol(mol_combo)
    ff        = FfWrapper(mol_combo, only_hydrogen=True)
    x0        = ff.pos
    
    result = optimize.minimize(
        fun=ff.ene,
        x0=x0,
        jac=ff.grad,
        method="L-BFGS-B"
        )
    ff.pos = result.x
    energy = ff.ene(result.x)

    return energy, Chem.GetMolFrags(ff.mol, asMols=True)


def assign_protonation_states(
    cell,
    mol_list,
    N_iterations=10):

    """
    Find the energetically best protonation state of the unitcell.
    """

    import numpy as np
    from rdkit.Chem import AllChem as Chem
    import copy

    mol_combo = combine_mols(mol_list)
    polar_heavy_atom = Chem.MolFromSmarts("[#7,#8,#16]")
    matches = mol_combo.GetSubstructMatches(
        polar_heavy_atom, 
        uniquify=False
        )
    N_matches = len(matches)
    delta_hydrogen_dict_best = dict()
    for atm_idx in matches:
        atm_idx = int(atm_idx[0])
        delta_hydrogen_dict_best[atm_idx] = 0
    energy_best = 999999999999999999999.
    mol_list_best = copy.deepcopy(mol_list)
    N_mol_per_unitcell = len(mol_list)

    charge0 = Chem.GetFormalCharge(mol_combo)

    #with open(f"./best_init.pdb", "w") as fopen:
    #    fopen.write(
    #        Chem.MolToPDBBlock(
    #            combine_mols(
    #                mol_list_best
    #                )
    #            )
    #        )
    #count = 0
    for _ in range(N_iterations):
        delta_hydrogen_dict = copy.deepcopy(
            delta_hydrogen_dict_best
            )
        for i in range(N_matches):
            atm_idx_i = int(matches[i][0])
            delta0_i  = delta_hydrogen_dict[atm_idx_i]
            for j in range(i+1, N_matches):
                atm_idx_j = int(matches[j][0])
                delta0_j  = delta_hydrogen_dict[atm_idx_j]
                for dij in [[-1,1],[1,-1]]:
                    ###  0: Do nothing
                    ### +1: Add hydrogen
                    ### -1: Remove hydrogen
                    delta_hydrogen_dict[atm_idx_i] = dij[0]
                    delta_hydrogen_dict[atm_idx_j] = dij[1]
                    mol_list_prot = apply_delta_hydrogen(
                        copy.deepcopy(mol_list),
                        delta_hydrogen_dict
                        )
                    charge = Chem.GetFormalCharge(
                        combine_mols(mol_list_prot)
                        )
                    if charge == charge0:
                        mol_list_prot, _, _ = make_supercell(
                            cell,
                            mol_list_prot,
                            0, 2,
                            0, 2,
                            0, 2)
                        energy, mol_list_prot = minimize_H(mol_list_prot)
                        if energy < energy_best:
                            energy_best = energy
                            delta_hydrogen_dict_best = copy.deepcopy(
                                delta_hydrogen_dict
                                )
                            mol_list_best = copy.deepcopy(
                                mol_list_prot[:N_mol_per_unitcell]
                                )
                delta_hydrogen_dict[atm_idx_j] = delta0_j
            delta_hydrogen_dict[atm_idx_i] = delta0_i
        #print(f"Best energy {energy_best:4.2f}")
        #with open(f"./best_{count}.pdb", "w") as fopen:
        #    fopen.write(
        #        Chem.MolToPDBBlock(
        #            combine_mols(
        #                mol_list_best
        #                )
        #            )
        #        )
        #count += 1

    return mol_list_best


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
    cell, 
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
        1  : 1.09, #H
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
        cell,
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
                ortho = cell.orthogonalize(
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
            frac = cell.fractionalize(
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
                        ortho = cell.orthogonalize(
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

    import copy
    mol_list_new = copy.deepcopy(mol_list)
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
        mol_list_new.append(mol_cp)

    return mol_list_new


def make_P1(
    cell, 
    atom_crds_ortho, 
    atom_num, 
    addhs=False, 
    use_openeye=False):

    """
    Generate the P1 cell. Return tuple with atomic coordinates (in Ang) and
    atomic numbers of all atoms in P1 cell.
    """

    import networkx as nx
    from rdkit.Chem import AllChem as Chem
    from rdkit.Geometry import Point3D
    import copy

    _atom_crds_ortho = list()
    _atom_num = list()
    N_atoms = len(atom_num)
    for i in range(N_atoms):
        if addhs:
            if atom_num[i] != 1:
                _atom_crds_ortho.append(copy.copy(atom_crds_ortho[i]))
                _atom_num.append(copy.copy(atom_num[i]))
        else:
            _atom_crds_ortho.append(copy.copy(atom_crds_ortho[i]))
            _atom_num.append(copy.copy(atom_num[i]))

    atom_crds_ortho = _atom_crds_ortho
    atom_num = _atom_num

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
            frac = cell.fractionalize(
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
                            ortho = cell.orthogonalize(
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
                cell.fractionalize(
                    gemmi.Position(
                        *pos
                    )
                ).tolist()
            )
        frac_crds   = np.array(frac_crds)
        valid_atoms = np.where(
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
                    ortho = cell.orthogonalize(
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
                ortho = cell.orthogonalize(
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
        if use_openeye:
            import warnings
            warnings.warn("With addhs=True, we automatically set use_openeye=True.")
        from openeye import oechem
        from openeye import oequacpac
        from xtalmdscripts.supercellbuilding.oe_utils import rdmol_from_oemol
        from xtalmdscripts.supercellbuilding.oe_utils import oemol_from_rdmol

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

            oechem.OEAssignAromaticFlags(oemol)
            oechem.OEDetermineConnectivity(oemol)
            oechem.OEFindRingAtomsAndBonds(oemol)
            oechem.OEPerceiveBondOrders(oemol)
            oechem.OE3DToInternalStereo(oemol)
            oechem.OEPerceiveChiral(oemol)
            oechem.OEAssignImplicitHydrogens(oemol)
            oechem.OEAssignFormalCharges(oemol)

            oechem.OEAddExplicitHydrogens(oemol)

            oechem.OEAssignAromaticFlags(oemol)

            mol = rdmol_from_oemol(oemol)
            Chem.AssignStereochemistryFrom3D(mol)
            mol_list.append(mol)

            #with open(f"./test_{count}.pdb", "w") as fopen:
            #    fopen.write(Chem.MolToPDBBlock(mol))
            count += 1

    else:
        if use_openeye:
            from openeye import oechem
            from openeye import oequacpac
            from xtalmdscripts.supercellbuilding.oe_utils import rdmol_from_oemol
            from xtalmdscripts.supercellbuilding.oe_utils import oemol_from_rdmol

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
                oechem.OEPerceiveBondOrders(oemol)
                oechem.OE3DToInternalStereo(oemol)
                oechem.OEPerceiveChiral(oemol)
                oechem.OEAssignImplicitHydrogens(oemol)
                oechem.OEAssignFormalCharges(oemol)

                oechem.OEAssignAromaticFlags(oemol)
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

    #strc_write               = gemmi.Structure()
    #strc_write.spacegroup_hm = "P1"
    #strc_write.cell          = cell
    #with open("./make_p1_test.pdb", "w") as fopen:
    #    fopen.write(get_pdb_block(mol_list, strc_write))

    return mol_list


def make_supercell(
    cell,
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
                        frac_pos = cell.fractionalize(
                            gemmi.Position(*pos)
                        )
                        frac_pos.x += a
                        frac_pos.y += b
                        frac_pos.z += c
                        
                        conf_pos_abc = cell.orthogonalize(frac_pos).tolist()
                            
                        conf.SetAtomPosition(
                            atom_idx,
                            Point3D(*conf_pos_abc)
                        )

                    replicated_mol_list.append(mol)
                    conf = mol.GetConformer(0)
                    mol_identifies.append(mol_idx)
                    unitcell_in_supercell_fracs.append([a,b,c])

    return replicated_mol_list, mol_identifies, unitcell_in_supercell_fracs


def clean_names(mol_list):

    """
    Uniquify atom and residue names. Equalize residue names for
    chemically equal residues.
    """

    import copy

    mol_list_new = list()
    N_mol = len(mol_list)
    for mol_idx in range(N_mol):
        mol = copy.deepcopy(mol_list[mol_idx])
        
        atom_counts_dict = dict()
        for atom in mol.GetAtoms():
            mi  = Chem.AtomPDBResidueInfo()
            mi.SetIsHeteroAtom(True)
            mi.SetResidueName(f'M{mol_idx}'.ljust(3))
            mi.SetResidueNumber(mol_idx + 1)
            mi.SetOccupancy(1.0)
            mi.SetTempFactor(0.0)
            atomic_num = atom.GetAtomicNum()
            atomic_ele = atom.GetSymbol()
            if not atomic_num in atom_counts_dict:
                atom_counts_dict[atomic_num] = 1
            else:
                atom_counts_dict[atomic_num] += 1
            mi.SetName(
                f"{atomic_ele}{atom_counts_dict[atomic_num]}".ljust(4)
                )
            atom.SetMonomerInfo(mi)

        mol_list_new.append(mol)

    return mol_list_new


def get_unique_mapping(
    mol_list, 
    stereochemistry=True
    ):

    import copy

    """
    Get unique mapping dict and list of unique rdkit mol objects in list of rdkit mol objects.
    if `stereochemistry=True`, the mapping will honor stereochemistry.
    """

    N_mol = len(mol_list)

    smiles_list = [Chem.MolToSmiles(mol, isomericSmiles=stereochemistry, allHsExplicit=True) for mol in mol_list]
    smiles_list_unique = list(set(smiles_list))
    ### Sort smiles according to built-in
    ### python `sorted` method.
    smiles_list_unique_argsort = sorted(
        range(
            len(smiles_list_unique)
            ), 
        key=smiles_list_unique.__getitem__
        )
    _smiles_list_unique = [smiles_list_unique[i] for i in smiles_list_unique_argsort]
    smiles_list_unique = _smiles_list_unique

    rdmol_list_unique  = list()
    unique_mapping     = dict()
    for smiles_unique_idx, smiles_unique in enumerate(smiles_list_unique):
        found_unique = False
        for mol_idx in range(N_mol):
            mol    = mol_list[mol_idx]
            smiles = smiles_list[mol_idx]
            if smiles == smiles_unique:
                if not found_unique:
                    rdmol_list_unique.append(
                        copy.deepcopy(mol)
                        )
                    found_unique = True
                    unique_mapping[mol_idx] = smiles_unique_idx
                else:
                    unique_mapping[mol_idx] = smiles_unique_idx

    assert len(unique_mapping) == N_mol

    return unique_mapping, rdmol_list_unique


def equalize_rdmols(
    mol_list, 
    stereochemistry=True
    ):

    """
    Get list of rdkit mol objects in which all chemically idential mol objects
    have identical topology and pdb monomer info. Only difference are coordinates.
    If `stereochemistry=True` it will honor stereochemistry.
    """

    import copy

    unique_mapping, rdmol_list_unique = get_unique_mapping(mol_list, stereochemistry)
    N_unique_mols = len(rdmol_list_unique)

    mol_list_new = copy.deepcopy(mol_list)
    for mol_idx in unique_mapping:
        ### This is the molecule that holds the correct coordinates
        ### and pdb monomer info.
        mol_crds = copy.deepcopy(mol_list[mol_idx])
        ### This is the molecule that holds the correct names, ordering, etc...
        mol_info = copy.deepcopy(rdmol_list_unique[unique_mapping[mol_idx]])
        match    = mol_crds.GetSubstructMatch(mol_info, useChirality=stereochemistry)

        conf_pos_crds = mol_crds.GetConformer(0).GetPositions()
        conf_info     = mol_info.GetConformer(0)
        for mol_info_atm_idx, mol_crds_atm_idx in enumerate(match):
            pos = conf_pos_crds[mol_crds_atm_idx]
            conf_info.SetAtomPosition(
                mol_info_atm_idx,
                Point3D(*pos)
            )
            ### Note, we cannot `copy.copy(mi_original)`
            ### or `copy.copy(mol_target.GetAtomWithIdx(atm_idx))`
            mi = mol_info.GetAtomWithIdx(mol_info_atm_idx).GetMonomerInfo()
            mi.SetResidueName(f'M{unique_mapping[mol_idx]}'.ljust(3))
            mi.SetResidueNumber(mol_idx + 1)
            mol_info.GetAtomWithIdx(mol_info_atm_idx).SetMonomerInfo(mi)

        mol_list_new[mol_idx] = mol_info

    return mol_list_new


def generate_replicated_mol_list(
    cell,
    atom_crds_ortho,
    atom_num, 
    a_min_max,
    b_min_max,
    c_min_max,
    addhs=False,
    protonate_unitcell=True,
    addwater=0,
    N_iterations_protonation=0,
    use_openeye=False,
    ):

    """
    Generate rdkit mol object list for molecules in supercell. supercell is generated
    according input parameters.
    """

    mol_list = make_P1(cell, atom_crds_ortho, atom_num, addhs, use_openeye)
    if N_iterations_protonation > 0: 
        mol_list = assign_protonation_states(
            cell=cell, 
            mol_list=mol_list, 
            N_iterations=N_iterations_protonation
            )
    if addwater > 0:
        mol_list = random_fill(
            cell,
            mol_list,
            N_per_unitcell=addwater,
            radius=0.5,
            smiles="O"
            )
    mol_list = clean_names(mol_list)

    replicated_mol_list, mol_identifies, unitcell_in_supercell_fracs = make_supercell(
        cell,
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

    crds_block = Chem.MolToPDBBlock(mol_new, flavor=8|32)
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
    

def parse_cif(
    cif_path, 
    use_symmetry_operations=False
    ):

    """
    Parse cif file as gemmis structure object.
    """

    import gemmi

    doc  = gemmi.cif.read(cif_path)[0]
    strc = gemmi.make_small_structure_from_block(doc)
    
    ### finding it by number is much better
    ### Sometimes the HM name cannot be found by gemmi.
    table_number = -1
    for item in doc:
        if item.pair == None:
            continue
        key, value = item.pair
        if "_symmetry_Int_Tables_number".lower() == key.lower():
            table_number=int(value)
            break
    if table_number > -1:
        strc.spacegroup_hm = gemmi.find_spacegroup_by_number(table_number).hm

    atom_crds_ortho = list()
    atom_num = list()
    if use_symmetry_operations:
        op_list  = doc.find_values('_symmetry_equiv_pos_as_xyz')
        gops     = gemmi.GroupOps([gemmi.Op(o) for o in op_list])
        for site in strc.sites:
            for op in gops:
                pos_frac = op.apply_to_xyz(site.fract.tolist())
                pos = strc.cell.orthogonalize(
                    gemmi.Fractional(
                        *pos_frac
                        )
                    ).tolist()
                atom_crds_ortho.append(pos)
                atom_num.append(site.element.atomic_number)
    else:
        for site in strc.get_all_unit_cell_sites():
            pos = strc.cell.orthogonalize(site.fract)
            atom_crds_ortho.append([pos.x, pos.y, pos.z])
            atom_num.append(site.element.atomic_number)

    return strc, atom_crds_ortho, atom_num


def get_supercell_info_str(
    mol_identifies, 
    unitcell_in_supercell_fracs
    ):

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
    strc, atom_crds_ortho, atom_num = parse_cif(args.input, args.use_symmetry_operations)

    if len(args.a_min_max) == 2:
        a_min_max = [int(args.a_min_max[0]), int(args.a_min_max[1])]
    elif len(args.a_min_max) == 1:
        uc_length_a = strc.cell.a * 0.1
        a_min_max = [0, np.ceil(args.a_min_max[0] / uc_length_a)]
    else:
        raise ValueError(
            "Argument a_min_max must be either pair (x,y) or single value (z)"
            )

    if len(args.b_min_max) == 2:
        b_min_max = [int(args.b_min_max[0]), int(args.b_min_max[1])]
    elif len(args.b_min_max) == 1:
        uc_length_b = strc.cell.b * 0.1
        b_min_max = [0, np.ceil(args.b_min_max[0] / uc_length_b)]
    else:
        raise ValueError(
            "Argument b_min_max must be either pair (x,y) or single value (z)"
            )

    if len(args.c_min_max) == 2:
        c_min_max = [int(args.c_min_max[0]), int(args.c_min_max[1])]
    elif len(args.c_min_max) == 1:
        uc_length_c = strc.cell.c * 0.1
        c_min_max = [0, np.ceil(args.c_min_max[0] / uc_length_c)]
    else:
        raise ValueError(
            "Argument c_min_max must be either pair (x,y) or single value (z)"
            )

    ### Build the supercell as a set of rdkit molecule objects
    ### ======================================================
    replicated_mol_list, mol_identifies, unitcell_in_supercell_fracs = generate_replicated_mol_list(
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
    pdb_str = get_pdb_str(
        replicated_mol_list, 
        strc_write,
        a_min_max,
        b_min_max,
        c_min_max
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
    mol_list = make_P1(strc.cell, atom_crds_ortho, atom_num, args.addhs, args.use_openeye)

    unitcell_mol_list, _, _ = make_supercell(
        strc.cell,
        mol_list,
        0,0,
        0,0,
        0,0,
        )
    from rdkit.Chem import Descriptors
    unitcell_weight = 0.
    if args.addwater > 0:
        unitcell_weight += (args.addwater * 18.01528)
    smiles_list     = list()
    for mol in unitcell_mol_list:
        unitcell_weight += Descriptors.MolWt(mol)
        smiles_list.append(Chem.MolToSmiles(mol, isomericSmiles=True))

    smiles_list = set(smiles_list)
    smiles_list = list(smiles_list)

    ### Output final summary
    ### ====================
    a_len = np.max(a_min_max) - np.min(a_min_max) + 1.
    b_len = np.max(b_min_max) - np.min(b_min_max) + 1.
    c_len = np.max(c_min_max) - np.min(c_min_max) + 1.

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
Density [g/cm3]               : {unitcell_weight / strc.cell.volume * 1.6605387823355087:4.2f}
SMILES for molecules in UC    : {" ".join(smiles_list)}
"""
### 1.6605 conversion g/mol/Ang^3 to g/cm^3
)

    try:
        diff = (unitcell_weight / strc.cell.volume * 1.6605) - float(cif_info_dict["density"])
        if abs(diff) > 0.1:
            warnings.warn(f"Density difference {diff}. Check structure.")
    except:
        pass

def entry_point():

    main()

if __name__ == "__main__":

    entry_point()