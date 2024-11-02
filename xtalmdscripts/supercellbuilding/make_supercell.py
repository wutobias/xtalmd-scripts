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
        '--label_residues', 
        "-lr", 
        action='store_true',
        help="This labels common amino acid residues to their commonly used name in PDB files.", 
        required=False,
        default=False,
        )

    return parser.parse_args()


def label_rdmol(rdmol):

    pass


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
    try:
        mol_combo_prot = emol.GetMol()
        Chem.SanitizeMol(mol_combo_prot)
        return True, Chem.GetMolFrags(mol_combo_prot, asMols=True)
    except:
        return False, copy.deepcopy(mol_list)

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
    ff = FfWrapper(mol_combo, only_hydrogen=True)
    x0 = ff.pos
    
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
    matches = mol_combo.GetSubstructMatches(
        Chem.MolFromSmarts("[#7,#8,#16]"),
        uniquify=True
        )
    N_matches = len(matches)
    delta_hydrogen_dict_best = dict()
    for atm_idx in matches:
        atm_idx = int(atm_idx[0])
        delta_hydrogen_dict_best[atm_idx] = 0
    mol_list_best = copy.deepcopy(mol_list)
    N_mol_per_unitcell = len(mol_list)
    while True:
        charge = Chem.GetFormalCharge(
            combine_mols(mol_list_best)
            )
        if charge == 0:
            break
        delta = -1 * charge/abs(charge)
        i = np.random.randint(0, N_matches)
        atm_idx_i = int(matches[i][0])
        delta_hydrogen_dict_best[atm_idx_i] = delta
        success, mol_list_prot = apply_delta_hydrogen(
            copy.deepcopy(mol_list_best),
            delta_hydrogen_dict_best
            )

        try:
            mol = combine_mols(mol_list_prot)
            mp  = Chem.MMFFGetMoleculeProperties(mol)
            r   = Chem.MMFFGetMoleculeForceField(mol, mp)
            if isinstance(r, type(None)):
                success = False
        except:
            success = False
        delta_hydrogen_dict_best[atm_idx_i] = 0
        if success:
            energy, mol_list_best = minimize_H(mol_list_prot)
            print(charge, energy)
        
    #mol_list_prot, _, _ = make_supercell(
    #    cell,
    #    mol_list_best,
    #    0, 2,
    #    0, 2,
    #    0, 2)
    mol_list_prot = mol_list_best
    energy_best, mol_list_prot = minimize_H(mol_list_prot)
    mol_list_best = copy.deepcopy(
        mol_list_prot[:N_mol_per_unitcell]
        )
    #with open(f"./best_init.pdb", "w") as fopen:
    #    fopen.write(
    #        Chem.MolToPDBBlock(
    #            combine_mols(
    #                mol_list_best
    #                )
    #            )
    #        )
    delta_hydrogen_dict = copy.deepcopy(
            delta_hydrogen_dict_best
            )
    for _ in range(N_iterations):
        matches_plus = mol_combo.GetSubstructMatches(
            Chem.MolFromSmarts("[#7,#8,#16&+0,+1,+2,+3]"),
            uniquify=True
            )
        matches_minus = mol_combo.GetSubstructMatches(
            Chem.MolFromSmarts("[#7,#8,#16&-0,-1,-2,-3]"),
            uniquify=True
            )
        matches_neutral = mol_combo.GetSubstructMatches(
            Chem.MolFromSmarts("[#7,#8,#16&-0]"),
            uniquify=True
            )
            
        matches_plus = [m[0] for m in matches_plus]
        matches_minus = [m[0] for m in matches_minus]
        matches_neutral = [m[0] for m in matches_neutral]
        
        N_matches_plus  = len(matches_plus)
        N_matches_minus = len(matches_minus)
        i = np.random.randint(0, N_matches_plus)
        j = np.random.randint(0, N_matches_minus)
        atm_idx_i = int(matches_plus[i])
        atm_idx_j = int(matches_minus[j])
        if atm_idx_i == atm_idx_j:
            continue
        if atm_idx_i in matches_neutral and atm_idx_j not in matches_neutral:
            continue
        if atm_idx_j in matches_neutral and atm_idx_i not in matches_neutral:
            continue
            
        ###  0: Do nothing
        ### +1: Add hydrogen
        ### -1: Remove hydrogen
        delta_hydrogen_dict[atm_idx_i] = -1
        delta_hydrogen_dict[atm_idx_j] = +1
        success, mol_list_prot = apply_delta_hydrogen(
            copy.deepcopy(mol_list_best),
            delta_hydrogen_dict
            )
        
        try:
            mol = combine_mols(mol_list_prot)
            mp  = Chem.MMFFGetMoleculeProperties(mol)
            r   = Chem.MMFFGetMoleculeForceField(mol, mp)
            if isinstance(r, type(None)):
                success = False
        except:
            success = False
        
        print(success,energy_best, Chem.GetFormalCharge(combine_mols(mol_list_best)))

        if success:
            #mol_list_prot, _, _ = make_supercell(
            #    cell,
            #    mol_list_prot,
            #    0, 2,
            #    0, 2,
            #    0, 2
            #    )
            try:
                energy, mol_list_prot = minimize_H(mol_list_prot)
            except:
                print("min failed")
                energy = energy_best+9999.
            print("new energy", energy)
            if energy < energy_best:
                energy_best = energy
                delta_hydrogen_dict[atm_idx_i] = 0
                delta_hydrogen_dict[atm_idx_j] = 0
                delta_hydrogen_dict_best = copy.deepcopy(
                    delta_hydrogen_dict
                    )
                mol_list_best = copy.deepcopy(
                    mol_list_prot[:N_mol_per_unitcell]
                    )
            
        delta_hydrogen_dict[atm_idx_i] = 0
        delta_hydrogen_dict[atm_idx_j] = 0
                    
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
    radius=1.7, 
    smiles="[H]O[H]",
    bin_width_real=0.5):

    """
    Randomly fill unit cell with molecules.
    """

    from rdkit.Chem import AllChem as Chem
    from rdkit.Geometry import Point3D
    import numpy as np
    from scipy.spatial import distance
    import gemmi
    import copy

    vdw_scaling    = 1.2

    solute_solvent_epsilon  = 1.
    solvent_solvent_epsilon = 1.

    vdw_dict = {
        1  : vdw_scaling*1.09, #H
        6  : vdw_scaling*1.70, #C
        7  : vdw_scaling*1.55, #N
        8  : vdw_scaling*1.52, #O
        9  : vdw_scaling*1.47, #F
        15 : vdw_scaling*1.80, #P
        16 : vdw_scaling*1.80, #S
        17 : vdw_scaling*1.75, #Cl
    }
    radius *= vdw_scaling

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    Chem.EmbedMolecule(mol)
    Chem.UFFOptimizeMolecule(mol)

    conformer      = mol.GetConformer()
    atom_crds_mol  = conformer.GetPositions()
    atom_crds_mol  = np.array(atom_crds_mol)

    atom_crds_xtal = list()
    atom_radii_xtal = list()
    frac = np.zeros(3, dtype=float)

    real2frac = np.array(
        cell.fractionalization_matrix.tolist()
        )
    frac2real = np.array(
        cell.orthogonalization_matrix.tolist()
        )
    
    all_rad_list = list()
    for mol_r in mol_list:
        all_rad_list.append(list())
        for atom in mol_r.GetAtoms():
            n = atom.GetAtomicNum()
            if n in vdw_dict:
                all_rad_list[-1].append(
                    vdw_dict[n]
                )
            else:
                raise ValueError(
                    f"Atomic number {n} not found."
                    )

    for mol_idx, mol_r in enumerate(mol_list):
        conf     = mol_r.GetConformer(0)
        sele     = [atom.GetIdx() for atom in mol_r.GetAtoms() if atom.GetAtomicNum() != 1]
        conf_pos = conf.GetPositions()[sele]
        rad_list = np.array(all_rad_list[mol_idx])[sele]
        frac_pos = np.matmul(real2frac, conf_pos.T).T
        for a in [-1.,0.,1.]:
            for b in [-1.,0.,1.]:
                for c in [-1.,0.,1.]:
                    frac_pos += [a,b,c]
                    conf_pos[:] = np.matmul(frac2real, frac_pos.T).T
                    atom_crds_xtal.extend(
                        conf_pos.tolist()
                        )
                    atom_radii_xtal.extend(
                        rad_list.tolist()
                        )
                    frac_pos -= [a,b,c]
    atom_crds_xtal   = np.array(atom_crds_xtal)
    atom_radii_xtal  = np.array(atom_radii_xtal)

    _x_extend = np.zeros((27,N_per_unitcell,3), dtype=float)
    _x = np.zeros((N_per_unitcell,3), dtype=float)
    radius_comb = np.sqrt(atom_radii_xtal*radius)[:,np.newaxis]
    r_min = 1.1225*radius_comb
    def func(x):

        _x[:] = 0.
        _x_extend[:] = 0.

        _x[:] = x.reshape((-1,3))
        frac  = np.matmul(real2frac, _x.T).T
        frac -= np.floor(frac)
        counts = 0
        for a in [-1.,0.,1.]:
            for b in [-1.,0.,1.]:
                for c in [-1.,0.,1.]:
                    frac[:] += [a,b,c]
                    _x_extend[counts] = np.matmul(frac2real, frac.T).T
                    frac[:] -= [a,b,c]
                    counts += 1

        x_extend = _x_extend.reshape((-1,3))
        dists_self = distance.pdist(
            x_extend,
            )
        valids_self = np.less(dists_self, 1.1225*radius)
        r_sigma_6 = (radius/dists_self[valids_self])**6
        score_selfclash = np.sum(
            4.*solvent_solvent_epsilon*(r_sigma_6*r_sigma_6-r_sigma_6)+solvent_solvent_epsilon
            )

        dists_atoms = distance.cdist(
            atom_crds_xtal,
            x_extend,
            )
        valids_atoms = np.less(dists_atoms, r_min)
        r_over_sigma = radius_comb/dists_atoms
        r_sigma_6 = (r_over_sigma[valids_atoms])**6
        score_atomclash = np.sum(
            4.*solute_solvent_epsilon*(r_sigma_6*r_sigma_6-r_sigma_6)+solute_solvent_epsilon
            )

        score_total = score_selfclash + score_atomclash

        return score_total

    o = cell.orthogonalize(
        gemmi.Fractional(0.,0.,0.)
        ).tolist()
    a1 = cell.orthogonalize(
        gemmi.Fractional(1.,0.,0.)
        ).tolist()
    b1 = cell.orthogonalize(
        gemmi.Fractional(0.,1.,0.)
        ).tolist()
    c1 = cell.orthogonalize(
        gemmi.Fractional(0.,0.,1.)
        ).tolist()

    Na = int(abs(a1[0]-o[0])/bin_width_real) + 1
    Nb = int(abs(b1[1]-o[1])/bin_width_real) + 1
    Nc = int(abs(c1[2]-o[2])/bin_width_real) + 1

    print(
        f"Generating {Na}x{Nb}x{Nc} vdw grid..."
        )
    grid = list()
    frac = np.zeros(3, dtype=float)
    for a in np.linspace(0., 1., Na, True):
        for b in np.linspace(0., 1., Nb, True):
            for c in np.linspace(0., 1., Nc, True):
                frac[:] = [a,b,c]
                ortho = cell.orthogonalize(
                    gemmi.Fractional(*frac)
                    ).tolist()
                dists = distance.cdist(
                    atom_crds_xtal, 
                    [ortho],
                    )
                is_outside = np.all(dists[:,0] > atom_radii_xtal)
                if is_outside:
                    grid.append(ortho)
    grid = np.array(grid)
    mask = np.arange(grid.shape[0], dtype=int)
    best_f = 999999999999999999999.
    best_x = np.zeros(3*N_per_unitcell, dtype=float)
    print(
        "Placing and optimizing positions in unit cell."
        )
    from scipy import optimize
    basinhopping = False
    iter_max = 5
    for i in range(iter_max):
        mask_selection = np.random.choice(
            mask, 
            size=N_per_unitcell, 
            replace=False
            )
        x0 = grid[mask_selection].flatten()
        print(
            f"Iteration {i+1}/{iter_max}",
            f"Initial func {func(x0)}"
            )
        if basinhopping:
            result = optimize.basinhopping(
                func, 
                x0, 
                minimizer_kwargs={
                    "method" : "BFGS"
                    },
                )
        else:
            result = optimize.minimize(
                    func, 
                    x0, 
                    method="BFGS",
                    )
        print(
            f"Final func {func(result.x)}"
            )
        if result.fun < best_f:
            best_f = result.fun
            best_x = result.x
        if best_f < 1.e-8:
            break

    print(
        f"Inserted {N_per_unitcell} molecules {Chem.MolToSmiles(mol)}."
        )

    import copy
    mol_list_new = copy.deepcopy(mol_list)
    for crds in best_x.reshape((-1,3)):
        frac  = np.matmul(real2frac, crds.T).T
        frac -= np.floor(frac)
        crds  = np.matmul(frac2real, frac.T).T
        mol_cp         = copy.deepcopy(mol)
        atom_crds_mol  = conformer.GetPositions()
        trans          = crds - np.mean(atom_crds_mol, axis=0)
        atom_crds_mol += trans
        conformer      = mol_cp.GetConformer()
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
    use_openeye=False,
    removewater=False,
    acmatrix=None,
    ):

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
    ucmatrix = np.array(cell.orth.mat.tolist())

    acmatrix_new, mol_new = xyz2mol.xyz2AC(
        atom_num,
        atom_crds_ortho,
        0,
        use_huckel=False,
        acmatrix=acmatrix,
        ucmatrix=ucmatrix,
        )
    mol_new = Chem.RWMol(mol_new)
    for i in range(N_atoms):
        for j in range(i+1,N_atoms):
            if acmatrix_new[i,j]:
                mol_new.AddBond(
                    i,j,Chem.BondType.UNSPECIFIED)
    mol_new = mol_new.GetMol()

    mol_list_new = Chem.GetMolFrags(
        mol_new, asMols=True)

    if removewater:
        mol_list_new = remove_water(mol_list_new)
    if addhs or use_openeye:
        if addhs and not use_openeye:
            import warnings
            warnings.warn("With addhs=True, we automatically set use_openeye=True.")
        from openeye import oechem
        from xtalmdscripts.supercellbuilding.oe_utils import rdmol_from_oemol

        _mol  = combine_mols(mol_list_new)
        oemol = oechem.OEMol()
        oemol.SetDimension(3)
        conf_pos = _mol.GetConformer(0).GetPositions()
        crds = list()
        map_atoms = dict()
        for atm_idx in range(N_atoms):
            atom = _mol.GetAtomWithIdx(atm_idx)
            a = oemol.NewAtom(
                atom.GetAtomicNum())
            map_atoms[atm_idx] = a
            crds.extend(conf_pos[atm_idx])
        oemol.SetCoords(crds)

        oechem.OEDetermineConnectivity(oemol)
        oechem.OEAssignAromaticFlags(oemol)
        oechem.OEPerceiveBondOrders(oemol)
        oechem.OEAssignAromaticFlags(oemol)
        oechem.OEFindRingAtomsAndBonds(oemol)
        oechem.OE3DToInternalStereo(oemol)
        oechem.OEPerceiveChiral(oemol)
        if addhs:
            oechem.OEAssignImplicitHydrogens(oemol)
        oechem.OEAssignFormalCharges(oemol)

        if addhs:
            oechem.OEAddExplicitHydrogens(oemol)

        oechem.OEAssignAromaticFlags(oemol)
        _mol = rdmol_from_oemol(oemol)

        #with open(f"./test_oechem.pdb", "w") as fopen:
        #    fopen.write(Chem.MolToPDBBlock(_mol))

        mol_list = list(Chem.GetMolFrags(_mol, asMols=True))
    else:
        from rdkit.Chem import rdDetermineBonds
        import itertools
        mol_list = list()
        N_mol = len(mol_list_new)
        Q = [-1,0,1]
        for q_list in itertools.product(Q, repeat=N_mol):
            if sum(q_list) == 0:
                mol_list = list()
                for i in range(N_mol):
                    q = q_list[i]
                    _mol = mol_list_new[i]
                    try:
                        rdDetermineBonds.DetermineBondOrders(
                            _mol, charge=q, embedChiral=False, allowChargedFragments=False)
                    except:
                        break
                    mol_list.append(Chem.Mol(_mol))
                if len(mol_list) == N_mol:
                    break
        
        #_mol = combine_mols(mol_list_new)
        #mol_list = Chem.GetMolFrags(_mol, asMols=True)
        #[_mol] = xyz2mol.xyz2mol(
        #    atom_num,
        #    atom_crds_ortho,
        #    charge=0,
        #    ucmatrix=ucmatrix)
        #mol_list = Chem.GetMolFrags(_mol, asMols=True)

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
            mi = atom.GetMonomerInfo()
            if mi == None:
                mi  = Chem.AtomPDBResidueInfo()
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
    stereochemistry=True,    
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
    stereochemistry=True,
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


def remove_water(mol_list):

    from rdkit import Chem
    import copy

    rdsma_list = [
        Chem.MolFromSmarts("[#1]-[#8]-[#1]"),
        Chem.MolFromSmarts("[#8]"),
        Chem.MolFromSmarts("[#8+1]"),
        Chem.MolFromSmarts("[#1]"),
        Chem.MolFromSmarts("[#1+1]"),
        Chem.MolFromSmarts("[#1]-[#8+1](-[#1])-[#1]"),
    ]
    mol_list_new = list()
    for mol in mol_list:
        _mol = copy.deepcopy(mol)
        nowat = 0
        for rdsma in rdsma_list:
            nowat += int(rdsma.HasSubstructMatch(_mol))
            if _mol.HasSubstructMatch(rdsma):
                matches = _mol.GetSubstructMatches(
                    rdsma, 
                    uniquify=False
                    )
                match_set = set()
                for match in matches:
                    match_set.update(set(match))
                if len(match_set) == _mol.GetNumAtoms():
                    nowat += 1
        if nowat == 0:
            mol_list_new.append(_mol)
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
    removewater=False,
    N_iterations_protonation=0,
    use_openeye=False,
    label_residues=False,
    acmatrix=None,
    ):

    """
    Generate rdkit mol object list for molecules in supercell. supercell is generated
    according input parameters.
    """

    mol_list = make_P1(cell, atom_crds_ortho, atom_num, addhs, use_openeye, removewater, acmatrix)

    if removewater:
        mol_list = remove_water(mol_list)

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
            radius=1.7,
            smiles="[H]O[H]",
            )
    if label_residues:
        N_mol = len(mol_list)
        for mol_idx in range(N_mol):
            rdmol = mol_list[mol_idx]
            success, _rdmol, atomorder = label_amino_acid_residues(rdmol, reorder=True)
            if success:
                mol_list[mol_idx] = _rdmol
    else:
        mol_list = clean_names(mol_list)

    replicated_mol_list, mol_identifies, unitcell_in_supercell_fracs = make_supercell(
        cell,
        mol_list, 
        a_min_max[0], a_min_max[1],
        b_min_max[0], b_min_max[1],
        c_min_max[0], c_min_max[1],
        )

    if not label_residues:
        replicated_mol_list = equalize_rdmols(replicated_mol_list)

    return replicated_mol_list, mol_identifies, unitcell_in_supercell_fracs


def get_pdb_block(
    replicated_mol_list, 
    strc_write=None,
    header=None):

    """
    Get pdb block as str. strc_write is gemmi structure object and must reflect
    the dimensions of the supercell.
    """

    import gemmi

    if isinstance(strc_write, (gemmi.Structure, gemmi.SmallStructure)):
        header = strc_write.make_pdb_headers()
    elif not isinstance(header, str):
        raise ValueError(
            "Must either provide header or Structure to provide header."
            )

    ### Combine all rdmols in a single big rdmol
    N_mol   = len(replicated_mol_list)
    mol_new = Chem.Mol()
    for mol_idx in range(N_mol):
        mol = copy.deepcopy(replicated_mol_list[mol_idx])
        mol_new = Chem.CombineMols(mol_new, mol)

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
    crds_block_new = []
    atom_ref = 0
    atom_count = 0
    for line in crds_block.split("\n"):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if atom_ref != atom_count:
                path = Chem.GetShortestPath(
                    mol_new, 
                    atom_ref, 
                    atom_count
                    )
                if len(path) == 0:
                    atom_ref = atom_count
                    crds_block_new.append("TER")
            atom_count += 1
        crds_block_new.append(line)

    pdb_block  = header + "\n" + "\n".join(crds_block_new)

    return pdb_block


def label_amino_acid_residues(rdmol, reorder=True):

    """
    Label amino acid residues according to standard pdb
    names. Does not rename atom names. This will also reorder,
    if `reorder=True`, the atoms in the molecule so that
    sequentiel residues are obtained.
    """
    
    import rdkit
    from rdkit import Chem
    
    generic = lambda x: f"[$([NX3H2,NX4H3+]),$([NX3H](C)(C)):1][CX4H:2]({x})[CX3:3](=[OX1:5])[OX2H,OX1-,N:4]"
    
    aa_dict = {
        "NME" : "[NX3H:1][CX4H3:4]",
        "ACE" : "[CX4H3:1][CX3:3](=[OX1:5])[N:4]",
        "ALA" : "[CH3X4:6]",
        "ARG" : "[CH2X4:6][CH2X4:7][CH2X4:8][NHX3:9][CH0X3:10](=[NH2X3+,NHX2+0:11])[NH2X3:12]",
        "ASN" : "[CH2X4:6][CX3:7](=[OX1:8])[NX3H2:9]",
        "ASP" : "[CH2X4:6][CX3:7](=[OX1:8])[OH0-,OH:9]",
        "CYS" : "[CH2X4:6][SX2H,SX1H0-,SX2H0:7]",
        "GLU" : "[CH2X4:6][CH2X4:7][CX3:8](=[OX1:9])[OH0-,OH:10]",
        "GLN" : "[CH2X4:6][CH2X4:7][CX3:8](=[OX1:9])[NH2X3:10]",
        #"GLY" : "[$([$([NX3H2,NX4H3+]),$([NX3H](C)(C)):1][CX4H2:2][CX3:3](=[OX1:5])[OX2H,OX1-,N:4])]",
        "GLY" : "[$([NX3H2,NX4H3+]),$([NX3H](C)(C)):1][CX4H2:2][CX3:3](=[OX1:5])[OX2H,OX1-,N:4]",
        "HIS" : "[CH2X4:6][#6X3:7]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H]):8]:[#6X3H:10]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H]):11]:[#6X3H:9]1",
        "ILE" : "[CHX4:6]([CH3X4:8])[CH2X4:7][CH3X4:9]",
        "LEU" : "[CH2X4:6][CHX4:7]([CH3X4:8])[CH3X4:9]",
        "LYS" : "[CH2X4:6][CH2X4:7][CH2X4:8][CH2X4:9][NX4+,NX3+0:10]",
        "MET" : "[CH2X4:6][CH2X4:7][SX2:8][CH3X4:9]",
        "PHE" : "[CH2X4:6][cX3:7]1[cX3H:8][cX3H:10][cX3H:12][cX3H:11][cX3H:9]1",
        "PRO" : "[$([NX3H,NX4H2+]),$([NX3](C)(C)(C)):1]1[CX4H:2]([CH2:6][CH2:7][CH2:8]1)[CX3:3](=[OX1:5])[OX2H,OX1-,N:4]",
        "SER" : "[CH2X4:6][OX2H:7]",
        "THR" : "[CHX4:6]([CH3X4:8])[OX2H:7]",
        "TRP" : "[CH2X4:6][cX3:7]1[cX3H:8][nX3H:10][cX3:11]2[cX3H:13][cX3H:15][cX3H:14][cX3H:12][cX3:9]12",
        "TYR" : "[CH2X4:6][cX3:7]1[cX3H:8][cX3H:10][cX3:12]([OHX2,OH0X1-:13])[cX3H:11][cX3H:9]1",
        "VAL" : "[CHX4:6]([CH3X4:7])[CH3X4:8]",
        "HOH" : "[#1:2]-[#8:1]-[#1:3]"
        }

    naming_order_dict = {
        "NME" : {1:"N",   2:"CA", 3:"C", 4:"CH3", 5:"O", },
        "ACE" : {1:"CH3", 2:"CA", 3:"C", 4:"N",   5:"O", },
        "ALA" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG"},
        "ARG" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG", 8:"CD",  9:"NE",   10:"CZ",  11:"NH1", 12:"NH2",},
        "ASN" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG", 8:"OD1", 9:"ND2",},
        "ASP" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG", 8:"OD1", 9:"OD2",},
        "CYS" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"SG"},
        "GLU" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG", 8:"CD", 9:"OE1",   10:"OE2"},
        "GLN" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG", 8:"CD", 9:"OE1",   10:"NE2"},
        "GLY" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", },
        "HIS" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG",  8:"ND1", 9:"CD2", 10:"CE1", 11:"NE2"},
        "ILE" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG1", 8:"CG2", 9:"CD1"},
        "LEU" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG",  8:"CD1", 9:"CD2"},
        "LYS" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG",  8:"CD",  9:"CE",  10:"NZ"},
        "MET" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG",  8:"SD",  9:"CE"},
        "PHE" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG",  8:"CD1", 9:"CD2", 10:"CE1", 11:"CE2", 12:"CZ"},
        "PRO" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG",  8:"CD"},
        "SER" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"OG"},
        "THR" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"OG1", 8:"CG2"},
        "TRP" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG",  8:"CD1", 9:"CD2", 10:"NE1", 11:"CE2", 12:"CE3", 13:"CZ2", 14:"CZ3", 15:"CH2"},
        "TYR" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG",  8:"CD1", 9:"CD2", 10:"CE1", 11:"CE2", 12:"CZ",  13:"OH"},
        "VAL" : {1:"N",   2:"CA", 3:"C", 4:"OXT", 5:"O", 6:"CB", 7:"CG1", 8:"CG2"},
        "HOH" : {1:"O",   2:"H1", 3:"H2" }
        }
    
    N_atoms    = rdmol.GetNumAtoms()
    ###             n  c  other
    names_list = [[-1,-1,-1] for _ in range(N_atoms)]
    resid_list = [[-1,-1,-1] for _ in range(N_atoms)]
    ###             n  c
    tag_list   = [[-1,-1,-1] for _ in range(N_atoms)]
    hname_list = [None for _ in range(N_atoms)]
    match_idx  = 0
    for aa in aa_dict:
        if aa in ["PRO", "GLY", "NME", "ACE", "HOH"]:
            smarts = aa_dict[aa]
        else:
            smarts = generic(aa_dict[aa])
        smarts_mol = Chem.MolFromSmarts(smarts)
        matches = rdmol.GetSubstructMatches(
            smarts_mol
        )
        for match in matches:
            for aa_idx, idx in enumerate(match):
                matchatom = smarts_mol.GetAtomWithIdx(aa_idx)
                tag       = matchatom.GetAtomMapNum()
                rdmolatom = rdmol.GetAtomWithIdx(idx)
                symbol    = rdmolatom.GetSymbol()
                ### n
                if tag == 1:
                    tag_list[idx][0]   = match_idx
                    names_list[idx][0] = aa
                    resid_list[idx][0] = match_idx
                ### c
                elif tag == 4:
                    tag_list[idx][1]   = match_idx
                    names_list[idx][1] = aa
                    resid_list[idx][1] = match_idx
                else:
                    tag_list[idx][2]   = tag
                    names_list[idx][2] = aa
                    resid_list[idx][2] = match_idx
                n_neighbors = 0
                for n in rdmolatom.GetNeighbors():
                    if n.GetAtomicNum() == 1:
                        n_neighbors += 1
                if tag in [1,4]:
                    hname = f"H{symbol}{n_neighbors}"
                else:
                    hname = "H"
                for n in rdmolatom.GetNeighbors():
                    if n.GetAtomicNum() == 1:
                        idx = n.GetIdx()
                        hname_list[idx]    = hname
                        tag_list[idx][2]   = tag
                        names_list[idx][2] = aa
                        resid_list[idx][2] = match_idx
            match_idx += 1

    ### Check if we have any unlabeled atoms
    for name in names_list:
        if name == [-1,-1,-1]:
            return False, rdmol, []

    N_resids = match_idx
    reordered_resid_list = [None for _ in range(N_resids)]
    ### First figure out n terminus
    found_n = False
    for idx, resid_pair in enumerate(tag_list):
        if resid_pair[0] > -1 and resid_pair[1] == -1:
            residq = resid_pair[0]
            reordered_resid_list[residq] = 0
            found_n = True
            break
    ### If not found, we are cyclic
    if not found_n:
        residq = 0
        reordered_resid_list[0] = 0
    ### Next, figure out everything else
    ### up to c terminus
    resid_count = 1
    found_c = False
    iter_count = 0
    while resid_count < N_resids and iter_count < 2*N_resids:
        for idx, resid_pair in enumerate(tag_list):
            if resid_pair[0] > -1 and resid_pair[1] == residq:
                residq = resid_pair[0]
                reordered_resid_list[residq] = resid_count
                resid_count += 1
                break
        iter_count += 1

    new_atom_order_dict = {i:[] for i in range(N_resids)}
    for idx in range(N_atoms):
        atom  = rdmol.GetAtomWithIdx(idx)
        ### C term
        if resid_list[idx][0] == -1 and resid_list[idx][1] > -1:
            name  = names_list[idx][1]
            resid = resid_list[idx][1]
            atom_name = "OXT"
        ### N term
        elif resid_list[idx][0] > -1:
            name = names_list[idx][0]
            resid = resid_list[idx][0]
            if name == "HOH":
                atom_name = "O"
            elif name == "ACE":
                atom_name = "CH3"
            else:
                atom_name = "N"
        ### C term
        elif resid_list[idx][1] > -1:
            name  = names_list[idx][1]
            resid = resid_list[idx][1]
            atom_name = "O"
        else:
            name = names_list[idx][2]
            resid = resid_list[idx][2]
            tag = tag_list[idx][2]
            symbol = atom.GetSymbol()
            if tag == -1:
                atom_name = symbol
            elif symbol == "H":
                if (tag in [1,4]) and (name not in ["ACE","NME"]):
                    atom_name = hname_list[idx]
                else:
                    atom_name = symbol
            else:
                atom_name = naming_order_dict[name][tag]
        resid = reordered_resid_list[resid]
        mi = Chem.AtomPDBResidueInfo()
        mi.SetIsHeteroAtom(False)
        mi.SetResidueName(name.ljust(3))
        mi.SetResidueNumber(resid+1)
        mi.SetName(atom_name.ljust(4))
        atom.SetMonomerInfo(mi)
        new_atom_order_dict[resid].append(idx)

    new_atom_order_list = []
    if reorder:
        for resid in range(N_resids):
            new_atom_order_list.extend(
                new_atom_order_dict[resid]
            )
        _rdmol = Chem.RenumberAtoms(rdmol, new_atom_order_list)
        import copy
        for idx_old, idx_new in enumerate(new_atom_order_list):
            atomold = rdmol.GetAtomWithIdx(idx_new)
            atomnew = _rdmol.GetAtomWithIdx(idx_old)
            miold = atomold.GetMonomerInfo()
            atomnew.SetMonomerInfo(miold)
    else:
        _rdmol = rdmol

    return True, _rdmol, new_atom_order_list


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
        removewater=args.removewater,
        N_iterations_protonation=args.n_protonation_attempts,
        use_openeye=args.use_openeye,
        label_residues=args.label_residues
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
    mol_list = make_P1(
        strc.cell, 
        atom_crds_ortho, 
        atom_num, 
        args.addhs, 
        args.use_openeye
        )
    if args.removewater:
        mol_list = remove_water(mol_list)

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
