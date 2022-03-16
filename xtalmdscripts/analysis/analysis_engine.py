"""
Collection of methods useful for comparing computed xtal structures
with xtal structures from experiment.
"""

def atoms_that_are_n_bonds_apart_dict(
    topology, 
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
    
    G = topology.to_bondgraph()

    n_bonds_apart_dict = dict()
    atom_count = topology.n_atoms
    for atm_idx in range(atom_count):
        n_bonds_apart_dict[atm_idx] = list()
        atom = topology.atom(atm_idx)
        if exclude_hydrogen:
            if atom.element.number == 1:
                continue
        path_dict = nx.shortest_path_length(
            G, 
            source=topology.atom(atm_idx)
        )
        for atom in path_dict:
            if exclude_hydrogen:
                if atom.element.number == 1:
                    continue
            if path_dict[atom] < n:
                n_bonds_apart_dict[atm_idx].append(atom.index)

    return n_bonds_apart_dict


def build_pair_ranks(
    topology,
    pair_list,
    residue_classes_list=None):
    
    """
    For a given list of atom pairs in a topology object, find the
    rank of each pair, i.e. classify each pair according to its 
    uniqueness in the crystal.    
    """

    import numpy as np
    
    if isinstance(residue_classes_list, type(None)):
        residue_classes_list = np.arange(
            topology.n_residues, 
            dtype=int
        )

    ### Build a signature for each atom
    atom_signature_list = np.zeros((topology.n_atoms,2), dtype=int)
    for residue in topology.residues:
        residue_class = residue_classes_list[residue.index]
        atom_index_in_residue = 0
        for atom in residue.atoms:
            atom_signature_list[atom.index] = [residue_class, atom_index_in_residue]
            atom_index_in_residue += 1
    unique_signatures_list = np.unique(atom_signature_list, axis=0)
    ### Convert the signatures into simple
    ### ranks for each atom.
    atom_rank_list = np.zeros(topology.n_atoms, dtype=int)
    for atom_idx in range(topology.n_atoms):
        atom_signatur = atom_signature_list[atom_idx]
        valids = np.where(
            (unique_signatures_list[:,0] == atom_signatur[0]) *\
            (unique_signatures_list[:,1] == atom_signatur[1])
            )
        atom_rank_list[atom_idx] = valids[0]
        
    pair_list_sorted = np.sort(pair_list, axis=1)
    
    ### From the atom ranks, get a signature for each pair
    pair_signature_list = np.zeros_like(pair_list_sorted, dtype=int)
    pair_signature_list[:,0] = atom_rank_list[pair_list_sorted[:,0]]
    pair_signature_list[:,1] = atom_rank_list[pair_list_sorted[:,1]]
    
    N_pairs = pair_list_sorted.shape[0]
    unique_signatures_list = np.unique(pair_signature_list, axis=0)
    pair_rank_list = np.zeros(N_pairs, dtype=int)
    for pair_idx in range(N_pairs):
        pair_signatur = pair_signature_list[pair_idx]
        valids        = np.where(
            (unique_signatures_list[:,0] == pair_signatur[0]) *\
            (unique_signatures_list[:,1] == pair_signatur[1])
        )
        pair_rank_list[pair_idx] = valids[0]
        
    return pair_rank_list


def build_pair_list(
    traj,
    rdmol,
    distance_cutoff = 0.4, # cutoff in nm
    bond_cutoff = 4, 
    exclude_hydrogen=True,
    exclude_water=True):

    """
    Find all atom pairs that are within `distance_cutoff` and at least
    `bond_cutoff` bonds apart. Return result as a list of atom_pair lists.
    """

    import numpy as np
    import mdtraj as md
    from rdkit import Chem

    water_list = np.array([], dtype=int)
    if exclude_water:
        params = Chem.SmilesParserParams()
        params.removeHs = False
        wat_rdmol = Chem.MolFromSmiles("[H]O[H]", params)
        water_list = get_all_matches(rdmol, wat_rdmol)
    
    ### Build exclusion dict first
    exclusion_dict = atoms_that_are_n_bonds_apart_dict(
        traj.topology,
        bond_cutoff,
        exclude_hydrogen
        )
    hydrogen_idxs = list()
    if exclude_hydrogen:
        for atom in traj.topology.atoms:
            if atom.element.number == 1:
                hydrogen_idxs.append(atom.index)

    ### This list contains all atom pairs
    pair_list = list()
    ### Loop over all atoms and find pairs
    for atom in traj.topology.atoms:
        if exclude_hydrogen:
            if atom.element.number == 1:
                continue
        haystack = np.arange(traj.topology.n_atoms, dtype=int)
        delete_list = exclusion_dict[atom.index]
        if exclude_hydrogen:
            delete_list.extend(hydrogen_idxs)
        if exclude_water:
            delete_list.extend(water_list)
        haystack = np.delete(haystack, delete_list)
        neighbor_list = md.compute_neighbors(
            traj,
            cutoff = 0.4, # cutoff in nm
            query_indices = [atom.index],
            haystack_indices = haystack,
            periodic=True
        )[0]
        pair_list_for_atom = np.zeros((neighbor_list.size,2), dtype=int)
        pair_list_for_atom[:,0] = neighbor_list
        pair_list_for_atom[:,1] = atom.index

        if atom.index == 0:
            pair_list           = pair_list_for_atom
        else:
            _pair_list = np.vstack(
                (
                    pair_list, 
                    pair_list_for_atom
                    )
                )
            pair_list  = _pair_list

    ### Sort them first, so that any identical i/j pair
    ### will appear identical.
    pair_list = np.sort(pair_list, axis=1)
    ### Then remove redundancies
    pair_list = np.unique(pair_list, axis=0)

    return pair_list


def compute_pairwise_distances(
    traj,
    pair_list,
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


def compute_tuplewise_angles(
    traj,
    tuple_list,
    ):

    """
    Compute pair-wise angles for atoms in `tuple_list` in traj.
    """

    import mdtraj as md

    angles = md.compute_angles(
        traj, 
        angle_indices=tuple_list,
        periodic=True
    )

    return angles


def get_all_matches(rdmol, query_mol):

    """
    Get all matching atom idxs in `rdmol` with `query_mol`
    """

    from rdkit import Chem
    import numpy as np

    match_list = list()
    frag_mols  = Chem.GetMolFrags(rdmol, asMols=True)
    frag_idxs  = Chem.GetMolFrags(rdmol, asMols=False)
    for frag, idxs in zip(frag_mols, frag_idxs):
        idxs = np.array(idxs)
        for match in frag.GetSubstructMatches(query_mol):
            match = np.array(match)
            match_list.extend(idxs[match].tolist())
    match_list = np.array(match_list, dtype=int)
    return match_list


def get_hbond_indices(
    traj,
    rdmol,
    exclude_water=True):

    """
    Get pair-wise atom indices of (h-atom-index, acceptor-atom-index) pairs.
    Returns three lists: 1.) 
    """

    import mdtraj as md
    from rdkit import Chem
    import warnings
    import numpy as np

    water_list = np.array([], dtype=int)
    if exclude_water:
        params = Chem.SmilesParserParams()
        params.removeHs = False
        wat_rdmol = Chem.MolFromSmiles("[H]O[H]", params)
        water_list = get_all_matches(rdmol, wat_rdmol)

    acc_O_single = get_all_matches(
        rdmol, 
        Chem.MolFromSmarts(
            "[*:1]-[#8:2]"
            )
        )
    acc_O_double = get_all_matches(
        rdmol, 
        Chem.MolFromSmarts(
            "[*:1]=[#8:2]"
            )
        )
    acc_N_single = get_all_matches(
        rdmol, 
        Chem.MolFromSmarts(
            "[*:1]-[#7:2]"
            )
        )
    acc_N_double = get_all_matches(
        rdmol, 
        Chem.MolFromSmarts(
            "[*:1]=[#7:2]"
            )
        )

    acc_O_single_atom_indices_set = set(acc_O_single)
    acc_O_double_atom_indices_set = set(acc_O_double)
    acc_N_single_atom_indices_set = set(acc_N_single)
    acc_N_double_atom_indices_set = set(acc_N_double)

    ### hbond_idxs: d_i, h_i, a_i
    hbond_idxs = md.baker_hubbard(
        traj,
        freq=0.0,
        exclude_water=False, # We exclude water using method below
        periodic=True,
    )

    acc_O_single_list = list()
    acc_O_double_list = list()
    acc_N_single_list = list()
    acc_N_double_list = list()

    for hbond_idx in hbond_idxs:
        if exclude_water and np.any(np.isin(water_list, hbond_idx)):
            continue
        d_i, h_i, a_i = hbond_idx
        if a_i in acc_O_single_atom_indices_set:
            acc_O_single_list.append(hbond_idx)
        elif a_i in acc_O_double_atom_indices_set:
            acc_O_double_list.append(hbond_idx)
        elif a_i in acc_N_single_atom_indices_set:
            acc_N_single_list.append(hbond_idx)
        elif a_i in acc_N_double_atom_indices_set:
            acc_N_double_list.append(hbond_idx)
        else:
            smiles_set = set([Chem.MolToSmiles(mol) for mol in Chem.GetMolFrags(rdmol, asMols=True)])
            warnings.warn(
                f"For {smiles_set}:Could not put hbond with d_i:{d_i}, h_i:{h_i}, a_i:{a_i} in any category."
                )

    acc_O_single_list = np.array(acc_O_single_list)
    acc_O_double_list = np.array(acc_O_double_list)
    acc_N_single_list = np.array(acc_N_single_list)
    acc_N_double_list = np.array(acc_N_double_list)

    ### Remove redundancies
    acc_O_single_list = np.unique(acc_O_single_list, axis=0)
    acc_O_double_list = np.unique(acc_O_double_list, axis=0)
    acc_N_single_list = np.unique(acc_N_single_list, axis=0)
    acc_N_double_list = np.unique(acc_N_double_list, axis=0)

    return acc_O_single_list, acc_O_double_list, acc_N_single_list, acc_N_double_list


def compute_com_diff_per_residue(
    query_traj, 
    ref_strc,
    rdmol,
    residue_classes_list,
    exclude_water=True):

    """
    Compute com difference for each residue in each frame of `query_traj`
    against the first frame in `ref_strc`.
    Returns list with com diffs (in nm) for each residue in each frame of
    `query_traj`.
    """

    import mdtraj as md
    import numpy as np
    from rdkit import Chem

    assert ref_strc.n_residues == query_traj.n_residues
    assert ref_strc.n_frames == 1

    water_list = np.array([], dtype=int)
    if exclude_water:
        params = Chem.SmilesParserParams()
        params.removeHs = False
        wat_rdmol = Chem.MolFromSmiles("[H]O[H]", params)
        water_list = get_all_matches(rdmol, wat_rdmol)
    
    a_scan = np.arange(-1,2, dtype=float)
    b_scan = np.arange(-1,2, dtype=float)
    c_scan = np.arange(-1,2, dtype=float)

    grid = np.array(
        np.meshgrid(
            *[a_scan, b_scan, c_scan]
        )
    )
    N_scans = 27
    grid    = grid.T.reshape((N_scans,-1))
    com_ref_big = np.zeros((ref_strc.n_residues,27, 3))
    diff_distances = np.zeros((query_traj.n_residues, query_traj.n_frames), dtype=float)
    count_start = 0
    count_stop  = 0
    for res in ref_strc.topology.residues:
        atm_idxs_per_residue = list()
        for atom in res.atoms:
            atm_idxs_per_residue.append(atom.index)
        if exclude_water and np.any(np.isin(water_list, atm_idxs_per_residue)):
            continue
        count_start = count_stop
        count_stop  = count_start + 27
        ### Compute com for given residue for each frame
        com_ref = md.compute_center_of_mass(
            traj=ref_strc.atom_slice(atm_idxs_per_residue)
        )
        box       = ref_strc.unitcell_vectors
        box_inv   = np.linalg.inv(box)
        com_frac  = np.einsum('ljk,lk->lj', box_inv, com_ref)

        com_frac  = np.expand_dims(com_frac, axis=0) + grid
        com_ref   = np.einsum('ljk,lhk->lhj', box, com_frac)
        com_ref_big[res.index] = com_ref
        
    for res in ref_strc.topology.residues:
        atm_idxs_per_residue = list()
        for atom in res.atoms:
            atm_idxs_per_residue.append(atom.index)

        ### Compute com for given residue for each frame
        com_traj = md.compute_center_of_mass(
            traj=query_traj.atom_slice(atm_idxs_per_residue)
        )
        
        valids = np.where(
            residue_classes_list == residue_classes_list[res.index]
        )[0]

        com_ref  = com_ref_big[valids].reshape(valids.size*27,3)
        diff     = com_traj - np.expand_dims(com_ref, axis=1)
        dists    = np.linalg.norm(diff, axis=2)
        diff_distances[res.index,:] = np.min(dists, axis=0)

    return diff_distances


def compute_com_diff_per_residue_DEPRECATED(
    query_traj, 
    ref_strc,
    N_unitcells_a,
    N_unitcells_b,
    N_unitcells_c):

    """
    Compute com difference for each residue in each frame of `query_traj`
    against the first frame in `ref_strc`.
    Returns list with com diffs (in nm) for each residue in each frame of
    `query_traj`.
    This works by translating each residue by the fraction of the number of
    unit cell repitions in each dimension until a minimal distance to the 
    same residue in the reference structure is found. This does not work 
    optimial and output trajectory do not look right...

    Use `compute_com_diff_per_residue` instead.
    """

    import mdtraj as md
    import numpy as np

    assert ref_strc.n_residues == query_traj.n_residues
    assert ref_strc.n_frames == 1
    
    a_scan = np.linspace(-1, 1, 2 * int(N_unitcells_a) + 1, True)
    b_scan = np.linspace(-1, 1, 2 * int(N_unitcells_b) + 1, True)
    c_scan = np.linspace(-1, 1, 2 * int(N_unitcells_c) + 1, True)

    grid = np.array(
        np.meshgrid(
            *[a_scan, b_scan, c_scan]
        )
    )
    N_scans  = 1
    N_scans *= 2 * int(N_unitcells_a) + 1
    N_scans *= 2 * int(N_unitcells_b) + 1
    N_scans *= 2 * int(N_unitcells_c) + 1
    grid = grid.T.reshape((N_scans,-1))

    ### Here we store the final min dist values
    diff_distances = np.zeros((query_traj.n_residues, query_traj.n_frames), dtype=float)
    for res in ref_strc.topology.residues:
        atm_idxs_per_residue = list()
        for atom in res.atoms:
            atm_idxs_per_residue.append(atom.index)
            
        ### Compute com for given residue for each frame
        com_ref = md.compute_center_of_mass(
            traj=ref_strc.atom_slice(atm_idxs_per_residue)
        )
        com_query = md.compute_center_of_mass(
            traj=query_traj.atom_slice(atm_idxs_per_residue)
        )

        ### Translate residue around box and find com minimum
        ### distance between ref and query
        box       = query_traj.unitcell_vectors
        box_inv   = np.linalg.inv(box)
        com_frac  = np.einsum('ljk,lk->lj', box_inv, com_query)

        com_frac_test = np.expand_dims(com_frac, axis=1) + grid
        com_test      = np.einsum('ljk,lhk->lhj', box, com_frac_test)
        diff          = com_ref - com_test
        dists         = np.linalg.norm(diff, axis=2)
        diff_distances[res.index, :] = np.min(dists, axis=1)
        
        ### The code below explictely loops over all axis and
        ### is here for testing purpose only. It is about 2 times 
        ### slower than the grid approach above.
        #
        #min_dists = np.zeros(query_traj.n_frames, dtype=float)
        #min_dists[:] = np.inf
        #for a in a_scan:
        #    for b in b_scan:
        #        for c in c_scan:
        #            com_frac_test = com_frac + [a,b,c]
        #            com_test      = np.einsum('ljk,lk->lj', box, com_frac_test)
        #            ### Sanity check
        #            #print(
        #            #    com_test[0],
        #            #    np.matmul(box[0], com_frac_test[0].T).T
        #            #)
        #            
        #            diff = com_ref-com_test
        #            dists = np.linalg.norm(diff, axis=1)
        #            valids = np.where(dists < min_dists)
        #            min_dists[valids] = np.copy(dists[valids])
        #diff_distances[res.index, :] = np.copy(min_dists)

    return diff_distances


def compute_pc_diff_per_residue(
    query_traj, 
    ref_strc,
    rdmol,
    exclude_water=True):

    """
    Compute difference in principal axes for each residue in each frame of `query_traj`
    against the first frame in `ref_strc`. The difference is reported as angle (in degree) between
    corresponding principle axes.
    Returns list with principal axes angle diffs (in degree) for each residue in each frame of
    `query_traj`.
    """

    import mdtraj as md
    import numpy as np
    from rdkit import Chem

    assert ref_strc.n_residues == query_traj.n_residues
    assert ref_strc.n_frames == 1

    water_list = np.array([], dtype=int)
    if exclude_water:
        params = Chem.SmilesParserParams()
        params.removeHs = False
        wat_rdmol = Chem.MolFromSmiles("[H]O[H]", params)
        water_list = get_all_matches(rdmol, wat_rdmol)

    diff_pcs = np.zeros((query_traj.n_residues, query_traj.n_frames, 3), dtype=float)
    for res in ref_strc.topology.residues:
        atm_idxs_per_residue = list()
        for atom in res.atoms:
            atm_idxs_per_residue.append(atom.index)
        if exclude_water and np.any(np.isin(water_list, atm_idxs_per_residue)):
            continue

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