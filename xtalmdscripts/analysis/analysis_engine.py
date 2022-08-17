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
    Distances computed in nanometer.
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
    Angles computed in radians.
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
    ON_bonds=False,
    exclude_water=True):

    """
    Get pair-wise atom indices of (h-atom-index, acceptor-atom-index) pairs.
    Returns list of list:
        [d_i, h_i, a_i]

        d_i: list of donor atom idxs
        h_i: list of hydrogen atom idxs
        a_i: list of acceptor atom idxs

    If `ON_bonds=True`, will return four lists with h-bond indices for O single bonds,
    O double bonds, N single bonds, N double bonds.

        acc_O_single_list, acc_O_double_list, acc_N_single_list, acc_N_double_list

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

    if ON_bonds:
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

        acc_O_single_list = list()
        acc_O_double_list = list()
        acc_N_single_list = list()
        acc_N_double_list = list()

    else:
        acc_list = list()

    ### hbond_idxs: d_i, h_i, a_i
    hbond_idxs = md.baker_hubbard(
        traj,
        freq=0.0,
        exclude_water=False, # We exclude water using method below
        periodic=True,
    )

    for hbond_idx in hbond_idxs:
        if exclude_water and np.any(np.isin(water_list, hbond_idx)):
            continue
        if ON_bonds:
            if a_i in acc_O_single_atom_indices_set:
                acc_O_single_list.append(hbond_idx)
            elif a_i in acc_O_double_atom_indices_set:
                acc_O_double_list.append(hbond_idx)
            elif a_i in acc_N_single_atom_indices_set:
                acc_N_single_list.append(hbond_idx)
            elif a_i in acc_N_double_atom_indices_set:
                acc_N_double_list.append(hbond_idx)
            else:
                d_i, h_i, a_i = hbond_idx
                smiles_set = set([Chem.MolToSmiles(mol) for mol in Chem.GetMolFrags(rdmol, asMols=True)])
                warnings.warn(
                    f"For {smiles_set}:Could not put hbond with d_i:{d_i}, h_i:{h_i}, a_i:{a_i} in any category."
                    )
        else:
            acc_list.append(hbond_idx)

    if ON_bonds:
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

    else:
        acc_list = np.array(acc_list)

        ### Remove redundancies
        acc_list = np.unique(acc_list, axis=0)

        return acc_list    


def compute_com_diff_per_residue(
    query_traj, 
    ref_strc,
    rdmol,
    N_closest_molecules,
    exclude_water=True):

    __doc__ = """
    Find com distances between each molecule and its `N_closest_molecules` 
    closest neighboring molecules in trajectory `query_traj`. The closest 
    molecules are determined and indexed as found in the reference structure 
    `ref_strc`.
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

    atm_idxs_per_residue_dict = dict()
    for res in ref_strc.topology.residues:
        atm_idxs_per_residue_dict[res.index] = list()
        for atom in res.atoms:
            atm_idxs_per_residue_dict[res.index].append(atom.index)
    
    a_scan = np.array([0,-1,1], dtype=float)
    b_scan = np.array([0,-1,1], dtype=float)
    c_scan = np.array([0,-1,1], dtype=float)

    grid = np.array(
        np.meshgrid(
            *[a_scan, b_scan, c_scan]
        )
    )
    N_scans = 27
    grid    = grid.T.reshape((N_scans,-1))
    com_ref_big = np.zeros((ref_strc.n_residues, 27, 3))
    com_query_big = np.zeros((query_traj.n_frames, ref_strc.n_residues, 27, 3))
    diff_distances = np.zeros(
        (query_traj.n_residues, query_traj.n_frames), 
        dtype=float
    )
    
    #top = md.Topology()
    #chain = top.add_chain()
    #for _ in range(N_closest_molecules):
    #    _res = top.add_residue("BLA", chain)
    #    top.add_atom("C", md.element.carbon, _res)
    #query_traj.save_pdb("traj.pdb")
    #ref_strc.save_pdb("ref_strc.pdb")
    
    for res in ref_strc.topology.residues:
        if exclude_water and np.any(np.isin(water_list, atm_idxs_per_residue_dict[res.index])):
            continue
        ### Compute com for given residue in reference structure
        com_ref = md.compute_center_of_mass(
            traj = ref_strc.atom_slice(atm_idxs_per_residue_dict[res.index])
        )
        box      = ref_strc.unitcell_vectors
        box_inv  = np.linalg.inv(box)
        com_frac = np.einsum('lkj,lk->lj', box_inv, com_ref)

        com_frac = np.expand_dims(com_frac, axis=0) + grid
        com_ref  = np.einsum('lkj,lhk->lhj', box, com_frac)
        com_ref_big[res.index] = com_ref
        
        ### Compute com for given residue in each frame
        ### in query structure.
        com_query = md.compute_center_of_mass(
            traj = query_traj.atom_slice(
                atm_idxs_per_residue_dict[res.index]
                )
            )
        box      = query_traj.unitcell_vectors
        box_inv  = np.linalg.inv(box)
        com_frac = np.einsum('lkj,lk->lj', box_inv, com_query)

        com_frac  = np.expand_dims(com_frac, axis=1) + grid
        com_query = np.einsum('lkj,lhk->lhj', box, com_frac)
        com_query_big[:,res.index] = com_query

    com_diff_result = np.zeros(
        (
            query_traj.n_frames, 
            query_traj.n_residues, 
            N_closest_molecules
        ),
        dtype=float    
    )
    for res in ref_strc.topology.residues:
        if exclude_water and np.any(np.isin(water_list, atm_idxs_per_residue_dict[res.index])):
            continue
        res_com = com_ref_big[res.index,0]
        diff    = com_ref_big - np.expand_dims(res_com, axis=0)
        dists   = np.linalg.norm(diff, axis=2)
        ### Get the index of the closest image for each residue
        ### to the reference residue.
        min_res_idxs  = np.argmin(dists, axis=1)
        ### Get the distance of the closest image for each residue
        ### to the reference residue.
        min_res_dists = np.choose(min_res_idxs, dists.T)
        ### Order them  by distance and only keep the `N_closest_molecules` ones.
        ### The closest one should always be residue itself at the image [0,0,0]
        ### We can skip that.
        min_res_idxs_sorted = np.argsort(min_res_dists)[1:N_closest_molecules+1]
        ### These are the images for each of the closest residues found.
        min_image_idxs_sorted = min_res_idxs[min_res_idxs_sorted]
        
        min_dists_ref = min_res_dists[min_res_idxs_sorted]
        
        #traj = md.Trajectory(
        #    com_ref_big[min_res_idxs_sorted, min_image_idxs_sorted],
        #    topology=top
        #)
        #traj.save_pdb(
        #    f"closest_neighbor_res{res.index}.pdb"
        #)
        
        res_com = com_query_big[:,res.index,0]
        res_com = np.expand_dims(res_com, axis=(1,2))
        diff    = res_com - com_query_big[:,min_res_idxs_sorted]
        dists   = np.linalg.norm(diff, axis=3)
        min_dists = np.min(dists, axis=2)

        com_diff_result[:,res.index] = min_dists

    return com_diff_result


def compute_pc_diff_per_residue(
    query_traj, 
    ref_strc,
    rdmol,
    N_closest_molecules,
    exclude_water=True):

    """
    Compute difference in principal axes for each residue in each frame of `query_traj`
    against the first frame in `ref_strc`. The difference is reported as angle (in radians) between
    corresponding principle axes.
    Returns list with principal axes angle diffs (in radians) for each residue in each frame of
    `query_traj`. Also return difference in principal axis angle (in radians) for each
    residue with respect to each of its `N_closest_molecules` neighboring molecules.
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
        
    atm_idxs_per_residue_dict = dict()
    for res in ref_strc.topology.residues:
        atm_idxs_per_residue_dict[res.index] = list()
        for atom in res.atoms:
            atm_idxs_per_residue_dict[res.index].append(atom.index)
            
    a_scan = np.array([0,-1,1], dtype=float)
    b_scan = np.array([0,-1,1], dtype=float)
    c_scan = np.array([0,-1,1], dtype=float)

    grid = np.array(
        np.meshgrid(
            *[a_scan, b_scan, c_scan]
        )
    )
    N_scans = 27
    grid    = grid.T.reshape((N_scans,-1))
    com_ref_big = np.zeros((ref_strc.n_residues, 27, 3))
    com_query_big = np.zeros((query_traj.n_frames, ref_strc.n_residues, 27, 3))
    pc_ref_big = np.zeros((ref_strc.n_residues, 3, 3))
    pc_query_big = np.zeros((query_traj.n_frames, ref_strc.n_residues, 3, 3))
    diff_distances = np.zeros(
        (query_traj.n_residues, query_traj.n_frames), 
        dtype=float
    )
    
    for res in ref_strc.topology.residues:
        if exclude_water and np.any(np.isin(water_list, atm_idxs_per_residue_dict[res.index])):
            continue
        ### Compute com for given residue in reference structure
        com_ref = md.compute_center_of_mass(
            traj = ref_strc.atom_slice(atm_idxs_per_residue_dict[res.index])
        )
        box      = ref_strc.unitcell_vectors
        box_inv  = np.linalg.inv(box)
        com_frac = np.einsum('lkj,lk->lj', box_inv, com_ref)

        com_frac = np.expand_dims(com_frac, axis=0) + grid
        com_ref  = np.einsum('lkj,lhk->lhj', box, com_frac)
        com_ref_big[res.index] = com_ref
        
        ### Compute com for given residue in each frame
        ### in query structure.
        com_query = md.compute_center_of_mass(
            traj = query_traj.atom_slice(atm_idxs_per_residue_dict[res.index])
        )
        box      = query_traj.unitcell_vectors
        box_inv  = np.linalg.inv(box)
        com_frac = np.einsum('lkj,lk->lj', box_inv, com_query)

        com_frac  = np.expand_dims(com_frac, axis=1) + grid
        com_query = np.einsum('lkj,lhk->lhj', box, com_frac)
        com_query_big[:,res.index] = com_query

        ### Compute principal axis.
        pa_ref = md.compute_inertia_tensor(
            traj=ref_strc.atom_slice(
                atm_idxs_per_residue_dict[res.index]
            ),
        )
        pa_query = md.compute_inertia_tensor(
            traj=query_traj.atom_slice(
                atm_idxs_per_residue_dict[res.index]
            ),
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
        
        pc_ref_big[res.index] = eigvecs_ref_sorted
        pc_query_big[:,res.index] = eigvecs_query_sorted
            
    diff_pcs_neighbors = np.zeros(
        (
            query_traj.n_frames,
            query_traj.n_residues, 
            N_closest_molecules,            
            3
        ), 
        dtype=float
    )
    diff_pcs_self = np.zeros(
        (
            query_traj.n_frames,
            query_traj.n_residues, 
            3
        ), 
        dtype=float
    )
    for res in ref_strc.topology.residues:
        if exclude_water and np.any(np.isin(water_list, atm_idxs_per_residue_dict[res.index])):
            continue

        ### =================================================== ###
        ### CALCULATE PC ANGLES BETWEEN QUERY AND REF STRUCTURE ###
        ### =================================================== ###
        ### Dot products query structure and reference structure
        dot_products_self  = np.einsum(
            'njk,jk->nj', #n: number of frames, #j,k: matrix indices
            pc_query_big[:,res.index],
            pc_ref_big[res.index],
        )

        ### Some dot products will be slightly larger than 1.
        ### Make them exactly one.
        valids = np.where(dot_products_self > 1.)
        dot_products_self[valids] = 1.
        valids = np.where(dot_products_self < -1.)
        dot_products_self[valids] = -1.

        angular_diffs_self = np.arccos(
            dot_products_self
        )
        
        angular_diffs_self_2 = np.zeros(
            (1,) + angular_diffs_self.shape + (2,)
        )

        angular_diffs_self_2[...,-2] = angular_diffs_self
        angular_diffs_self_2[...,-1] = np.pi - angular_diffs_self
        angular_diffs_self_min = np.min(angular_diffs_self_2, axis=-1)

        diff_pcs_self[:,res.index] = angular_diffs_self_min

        ### ============================================================================ ###
        ### CALCULATE PC ANGLES BETWEEN EACH MOLECULE AND ITS CLOSEST NEIGHBORS IN QUERY ###
        ### ============================================================================ ###
        res_com = com_ref_big[res.index,0]
        diff    = com_ref_big - np.expand_dims(res_com, axis=0)
        dists   = np.linalg.norm(diff, axis=2)
        ### Get the index of the closest image for each residue
        ### to the reference residue.
        min_res_idxs  = np.argmin(dists, axis=1)
        ### Get the distance of the closest image for each residue
        ### to the reference residue.
        min_res_dists = np.choose(min_res_idxs, dists.T)
        ### Order them  by distance and only keep the `N_closest_molecules` ones.
        ### The closest one should always be residue itself at the image [0,0,0]
        ### We can skip that.
        min_res_idxs_sorted = np.argsort(min_res_dists)[1:N_closest_molecules+1]

        ### Dot products with nearest neighbors in ref structure
        dot_products_ref  = np.einsum(
            'ljk,jk->lj', #l: N_closest_molecules #j,k: matrix indices
            pc_ref_big[min_res_idxs_sorted],
            pc_ref_big[res.index],
        )
        
        ### Dot products with nearest neighbors in query structure
        dot_products_query  = np.einsum(
            'nljk,njk->nlj', #n: number of frames #l: N_closest_molecules #j,k: matrix indices
            pc_query_big[:,min_res_idxs_sorted],
            pc_query_big[:,res.index],
        )
        
        ### Some dot products will be slightly larger than 1.
        ### Make them exactly one.
        valids = np.where(dot_products_ref > 1.)
        dot_products_ref[valids] = 1.
        valids = np.where(dot_products_ref < -1.)
        dot_products_ref[valids] = -1.
        
        valids = np.where(dot_products_query > 1.)
        dot_products_query[valids] = 1.
        valids = np.where(dot_products_query < -1.)
        dot_products_query[valids] = -1.

        angular_diffs_ref = np.arccos(
            dot_products_ref
        )
        angular_diffs_query = np.arccos(
            dot_products_query
        )
        
        angular_diffs_ref_2 = np.zeros(
            (1,) + angular_diffs_ref.shape + (2,)
        )
        angular_diffs_query_2 = np.zeros(
            angular_diffs_query.shape + (2,)
        )
        
        angular_diffs_ref_2[...,-2] = angular_diffs_ref
        angular_diffs_ref_2[...,-1] = np.pi - angular_diffs_ref
        angular_diffs_ref_min = np.min(angular_diffs_ref_2, axis=-1)
        
        angular_diffs_query_2[...,-2] = angular_diffs_query
        angular_diffs_query_2[...,-1] = np.pi - angular_diffs_query
        angular_diffs_query_min = np.min(angular_diffs_query_2, axis=-1)

        diff_pcs_neighbors[:,res.index] = angular_diffs_query_min

#        ### Serial way to compute differences
#        ### Should give same result as np.einsum
#        for frame_idx in range(query_traj.n_frames):
#            for res_j in range(N_closest_molecules):
#                for pc_idx in [0,1,2]:
#                    dot_ref = np.dot(
#                        pc_ref_big[min_res_idxs_sorted[res_j],pc_idx], 
#                        pc_ref_big[res.index,pc_idx]
#                    )
#                    dot_query = np.dot(
#                        pc_query_big[frame_idx,min_res_idxs_sorted[res_j],pc_idx], 
#                        pc_query_big[frame_idx,res.index,pc_idx]
#                    )
#
#                    if dot_ref > 1.:
#                        dot_ref =1.
#                    if dot_query > 1.:
#                        dot_query =1.
#
#                    angular_ref = np.arccos(
#                        dot_ref
#                    )
#                    angular_query = np.arccos(
#                        dot_query
#                    )
#
#                    angular_ref = np.min(
#                        [angular_ref,
#                         np.pi - angular_ref
#                        ]
#                    )
#                    angular_query = np.min(
#                        [angular_query,
#                         np.pi - angular_query
#                        ]
#                    )
#
#                    diff_pcs_neighbors[frame_idx, res.index, res_j, pc_idx] = angular_query

    return diff_pcs_neighbors, diff_pcs_self