"""
Collection of methods useful for comparing computed xtal structures
with xtal structures from experiment.
"""

def atoms_that_are_n_bonds_apart_dict(
    traj, 
    n=4, 
    exclude_hydrogen=True):
    
    """
    Find all atom pairs that are less then `n` bonds apart.
    Returns result as a dict where each key is an atom index and the
    corresponding value is the list of atom indices that are less than
    `n` bonds apart.
    """

    import networkx as nx
    import numpy as np
    import mdtraj as md
    
    G = traj.topology.to_bondgraph()

    n_bonds_apart_dict = dict()
    atom_count = traj.topology.n_atoms
    for atm_idx in range(atom_count):
        n_bonds_apart_dict[atm_idx] = list()
        atom = traj.topology.atom(atm_idx)
        if exclude_hydrogen:
            if atom.element.number == 1:
                continue
        path_dict = nx.shortest_path_length(
            G, 
            source=traj.topology.atom(atm_idx)
        )
        for atom in path_dict:
            if exclude_hydrogen:
                if atom.element.number == 1:
                    continue
            if path_dict[atom] < n:
                n_bonds_apart_dict[atm_idx].append(atom.index)

    return n_bonds_apart_dict


def build_pair_list(
    traj,
    distance_cutoff = 0.4, 
    bond_cutoff = 4, 
    exclude_hydrogen=True):

    """
    Find all atom pairs that are within `distance_cutoff` and at least
    `bond_cutoff` bonds apart. Return result as a list of atom_pair lists.
    """

    import numpy as np
    import mdtraj as md

    atom_count     = traj.topology.n_atoms

    ### Build exclusion dict first
    exclusion_dict = atoms_that_are_n_bonds_apart_dict(
        traj,
        bond_cutoff,
        exclude_hydrogen)
    ### Loop over all atoms and find pairs
    pair_list = list()
    for atm_idx in range(atom_count):
        atom = traj.topology.atom(atm_idx)
        if exclude_hydrogen:
            if atom.element.number == 1:
                continue
        haystack = np.arange(atom_count, dtype=int)
        haystack = np.delete(haystack, exclusion_dict[atm_idx])
        neighbor_list = md.compute_neighbors(
            traj,
            cutoff = 0.4, # cutoff in nm
            query_indices = [atm_idx],
            haystack_indices = haystack,
            periodic=True
        )[0]
        pair_list_for_atom = np.zeros((neighbor_list.size,2), dtype=int)
        pair_list_for_atom[:,0] = neighbor_list
        pair_list_for_atom[:,1] = atm_idx
        if atm_idx == 0:
            pair_list = pair_list_for_atom
        else:
            _pair_list = np.vstack((pair_list, pair_list_for_atom))
            pair_list  = _pair_list

    return pair_list


def compute_pairwise_distances(
    traj,
    pair_list
    ):

    """
    Compute pair-wise distances for atoms in `pair_list` in traj.
    """

    import mdtraj as md

    distances = md.compute_distances(
        traj, 
        atom_pairs=pair_list,
        periodic=True
    )

    return distances


def get_hbond_indices(
    traj,
    rdmol):

    """
    Get pair-wise atom indices of (h-atom-index, acceptor-atom-index) pairs.
    Returns three lists: 1.) 
    """

    import mdtraj as md
    from rdkit import Chem
    import warnings
    import numpy as np

    acc_O_single = rdmol.GetSubstructMatches(
        Chem.MolFromSmarts(
            "[*:1]-[#8:2]"
            )
        )
    acc_O_double = rdmol.GetSubstructMatches(
        Chem.MolFromSmarts(
            "[*:1]=[#8:2]"
            )
        )
    acc_N_single = rdmol.GetSubstructMatches(
        Chem.MolFromSmarts(
            "[*:1]-[#7:2]"
            )
        )
    acc_N_double = rdmol.GetSubstructMatches(
        Chem.MolFromSmarts(
            "[*:1]=[#7:2]"
            )
        )

    acc_O_single_atom_indices_set = set()
    for match in acc_O_single:
        for m in match:
            acc_O_single_atom_indices_set.add(m)

    acc_O_double_atom_indices_set = set()
    for match in acc_O_double:
        for m in match:
            acc_O_double_atom_indices_set.add(m)

    acc_N_single_atom_indices_set = set()
    for match in acc_N_single:
        for m in match:
            acc_N_single_atom_indices_set.add(m)

    acc_N_double_atom_indices_set = set()
    for match in acc_N_double:
        for m in match:
            acc_N_double_atom_indices_set.add(m)

    ### hbond_idxs: d_i, h_i, a_i
    hbond_idxs = md.baker_hubbard(
        traj,
        freq=0.0
    )

    acc_O_single_pair_list = list()
    acc_O_double_pair_list = list()
    acc_N_single_pair_list = list()
    acc_N_double_pair_list = list()

    for hbond_idx in hbond_idxs:
        d_i, h_i, a_i = hbond_idx
        if a_i in acc_O_single_atom_indices_set:
            acc_O_single_pair_list.append([h_i, a_i])
        elif a_i in acc_O_double_atom_indices_set:
            acc_O_double_pair_list.append([h_i, a_i])
        elif a_i in acc_N_single_atom_indices_set:
            acc_N_single_pair_list.append([h_i, a_i])
        elif a_i in acc_N_double_atom_indices_set:
            acc_N_double_pair_list.append([h_i, a_i])
        else:
            warnings.warn(
                f"Could not put hbond with d_i:{d_i}, h_i:{h_i}, a_i:{a_i} in any category."
                )

    acc_O_single_pair_list = np.array(acc_O_single_pair_list)
    acc_O_double_pair_list = np.array(acc_O_double_pair_list)
    acc_N_single_pair_list = np.array(acc_N_single_pair_list)
    acc_N_double_pair_list = np.array(acc_N_double_pair_list)

    return acc_O_single_pair_list, acc_O_double_pair_list, acc_N_single_pair_list, acc_N_double_pair_list


def compute_com_diff_per_residue(
    query_traj, 
    ref_strc):

    """
    Compute com difference for each residue in each frame of `query_traj`
    against the first frame in `ref_strc`.
    Returns list with com diffs (in nm) for each residue in each frame of
    `query_traj`.
    """

    import mdtraj as md
    import numpy as np

    assert ref_strc.n_residues == query_traj.n_residues
    assert ref_strc.n_frames == 1

    diff_distances = np.zeros((query_traj.n_residues, query_traj.n_frames), dtype=float)
    for res in ref_strc.topology.residues:
        atm_idxs_per_residue = list()
        for atom in res.atoms:
            atm_idxs_per_residue.append(atom.index)
       
        com_ref = md.compute_center_of_mass(
            traj=ref_strc.atom_slice(atm_idxs_per_residue)
        )
        com_query = md.compute_center_of_mass(
            traj=query_traj.atom_slice(atm_idxs_per_residue)
        )

        diff = com_ref-com_query
        dists = np.linalg.norm(diff, axis=1)
        diff_distances[res.index, :] = dists

    return diff_distances


def compute_pc_diff_per_residue(
    query_traj, 
    ref_strc):

    """
    Compute difference in principal axes for each residue in each frame of `query_traj`
    against the first frame in `ref_strc`. The difference is reported as angle (in degree) between
    corresponding principle axes.
    Returns list with principal axes angle diffs (in degree) for each residue in each frame of
    `query_traj`.
    """

    import mdtraj as md
    import numpy as np

    assert ref_strc.n_residues == query_traj.n_residues
    assert ref_strc.n_frames == 1

    diff_pcs = np.zeros((query_traj.n_residues, query_traj.n_frames, 3), dtype=float)
    for res in ref_strc.topology.residues:
        atm_idxs_per_residue = list()
        for atom in res.atoms:
            atm_idxs_per_residue.append(atom.index)

        pa_ref = md.compute_inertia_tensor(
            traj=ref_strc.atom_slice(atm_idxs_per_residue),
        )
        pa_query = md.compute_inertia_tensor(
            traj=query_traj.atom_slice(atm_idxs_per_residue),
        )
        ### Note: Eigenvectors are in columns
        eigvals_ref, eigvecs_ref = np.linalg.eig(pa_ref)
        eigvals_ref_sorted_idxs  = np.argsort(
            np.expand_dims(
                eigvals_ref,
                axis=1),
            axis=2
        )
        eigvecs_ref_sorted = np.take_along_axis(
            eigvecs_ref, 
            eigvals_ref_sorted_idxs,
            axis=2
        )

        eigvals_query, eigvecs_query = np.linalg.eig(pa_query)
        eigvals_query_sorted_idxs    = np.argsort(
            np.expand_dims(
                eigvals_query,
                axis=1),
            axis=2
        )
        eigvecs_query_sorted = np.take_along_axis(
            eigvecs_query, 
            eigvals_query_sorted_idxs,
            axis=2
        )

        ### Compare with directors. These must match. Directors just
        ### outputs the largest principal component vector
        #print("Eigen:", eigvals_query_sorted[0])
        #print("Directors:",
        #    md.compute_directors(
        #        traj=query_traj,
        #        indices=[atm_idxs_per_residue]
        #    )[0]
        #)
        ### Normalization.
        eigvecs_query_sorted = np.einsum(
            'ljk,lk->ljk',
            eigvecs_query_sorted,
            1./np.linalg.norm(eigvecs_query_sorted, axis=2)
        )
        eigvecs_ref_sorted = np.einsum(
            'ljk,lk->ljk',
            eigvecs_ref_sorted,
            1./np.linalg.norm(eigvecs_ref_sorted, axis=2)
        )
        dot_products  = np.einsum(
            'ljk,ajk->lk', 
            eigvecs_query_sorted, 
            eigvecs_ref_sorted
        )
        angular_diffs = np.arccos(
            dot_products
        ) * 180./np.pi
        angular_diffs_2 = np.zeros(
            angular_diffs.shape + (2,)
        )
        angular_diffs_2[...,-2] = angular_diffs
        angular_diffs_2[...,-1] = 180. - angular_diffs
        angular_diffs_min = np.min(angular_diffs_2, axis=-1)

        diff_pcs[res.index, :] = angular_diffs_min

        ### Serial way to compute dot products
        ### Should give same result as np.einsum
        #print(
        #   np.arccos(
        #       np.dot(
        #            eigvecs_query_sorted[0][:,1],
        #            eigvecs_ref_sorted[0][:,1]
        #        )
        #   )  * 180./np.pi,
        #   np.arccos(
        #       dot_products[0][1]
        #   ) * 180./np.pi
        #)

    return diff_pcs