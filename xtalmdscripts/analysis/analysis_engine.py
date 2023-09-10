"""
Collection of methods useful for comparing computed xtal structures
with xtal structures from experiment.
"""

def center_on_ref(ref_strc, query_traj, frame_idx=None):

    """
    Center query_traj on ref_strc in-place.
    """

    import numpy as np

    if frame_idx == None:
        frame_list = np.arange(query_traj.n_frames, dtype=int)
    else:
        frame_list = np.array([frame_idx], dtype=int)

    ref_cog  = np.mean(ref_strc.xyz[0], axis=0)
    for i in frame_list:
        
        x = query_traj.xyz[i]
        traj_cog = np.mean(x, axis=0)
        
        displ = ref_cog - traj_cog
        query_traj.xyz[i] += displ

    return query_traj


def center_frac(traj, frame_idx=None, mol_idx=None):

    """
    Shifts the traj (with frame `frame_idx` and molecule `mol_idx`)
    *in-place* so that the supercell origin is in the origin of the 
    labframe.
    """

    import numpy as np

    molecule_list = list(traj.topology.find_molecules())
    N_molecules   = len(molecule_list)
    
    if frame_idx == None:
        frame_list = np.arange(traj.n_frames, dtype=int)
    else:
        frame_list = np.array([frame_idx], dtype=int)
    if mol_idx == None:
        mol_list = np.arange(N_molecules, dtype=int)
    else:
        mol_list = np.array([mol_idx], dtype=int)
        
    _N_frames    = frame_list.size
    _N_molecules = mol_list.size
    
    r_fract_com_traj = np.zeros((_N_frames, _N_molecules, 3), dtype=float)
    r_com_traj = np.zeros((_N_frames, _N_molecules, 3), dtype=float)

    r_com       = np.zeros((_N_molecules, 3), dtype=float)
    r_fract_com = np.zeros((_N_molecules, 3), dtype=float)
    for i, frame_i in enumerate(frame_list):
        
        r_com[:] = 0.
        r_fract_com[:] = 0.

        uc = traj.unitcell_vectors[frame_i].T
        r  = traj.xyz[frame_i]

        ucinv   = np.linalg.inv(uc)
        r_fract = np.matmul(ucinv, r.T).T

        for j, mol_j in enumerate(mol_list):
            molecule = molecule_list[mol_j]
            total_mass = 0.
            for atom in molecule:
                total_mass += atom.element.mass
                r_com[j] += r[atom.index] * atom.element.mass
            r_com[j] /= total_mass
            r_fract_com[j]  = np.matmul(ucinv, r_com[j].T).T
            r_fract_com[j] -= np.floor(r_fract_com[j])
            _r_com = np.matmul(uc, r_fract_com[j].T).T
            shift = _r_com - r_com[j]
            r_com[j] = _r_com
            for atom in molecule:
                r[atom.index] += shift
        traj.xyz[frame_i] = r[:]
        r_fract_com_traj[i] = r_fract_com[:]
        r_com_traj[i]       = r_com[:]

    return r_fract_com_traj, r_com_traj
    

def unwrap_trajectory(query_traj, ref_strc, min_real=True):

    """
    Unwraps the molecules in a xtal trajectory `query_traj` based on the COM positions
    in the reference structure `ref_strc`. Unwrapping is based on the difference
    in fractional or real (if `min_real=True`) coordinate w.r.t. the reference structure.
    """

    import numpy as np
    import copy

    translation_length = 1

    query_traj_cp = copy.deepcopy(query_traj)
    ref_strc_cp   = copy.deepcopy(ref_strc)

    uc_ref = ref_strc_cp.unitcell_vectors[0].T

    molecule_list = list(ref_strc_cp.topology.find_molecules())
    N_molecules   = len(molecule_list)
    
    r_fract_com_traj, r_com_traj = center_frac(ref_strc_cp, frame_idx=0)
    r_fract_com_ref = r_fract_com_traj[0]
    r_com_ref = r_com_traj[0]

    r_com = np.zeros(3, dtype=float)
    r_fract = np.zeros(3, dtype=float)
    r  = np.zeros((query_traj_cp.n_atoms, 3), dtype=float)
    uc = np.zeros((3,3), dtype=float)
    ucinv = np.zeros((3,3), dtype=float)
    for i in range(query_traj_cp.n_frames):

        uc[:]  = np.abs(query_traj_cp.unitcell_vectors[i].T)
        uc[:] *= np.sign(uc_ref)
        ucinv[:] = np.linalg.inv(uc)

        for mol_idx in range(N_molecules):
            molecule = molecule_list[mol_idx]

            center_on_ref(ref_strc_cp, query_traj_cp, i)
            r[:] = np.copy(query_traj_cp.xyz[i])

            r_com[:]   = 0.
            total_mass = 0.
            for atom in molecule:
                total_mass += atom.element.mass
                r_com += r[atom.index] * atom.element.mass
            r_com  /= total_mass
            r_fract[:] = np.matmul(ucinv, r_com.T).T

            min_diff = np.inf
            min_r_fract = r_fract
            for a in np.arange(-translation_length,translation_length+1,1):
                for b in np.arange(-translation_length,translation_length+1,1):
                    for c in np.arange(-translation_length,translation_length+1,1):
                        _r_fract  = r_fract + np.array([a,b,c])
                        if min_real:
                            diff = np.linalg.norm(
                                np.matmul(uc, _r_fract.T).T-
                                r_com_ref[mol_idx]
                            )
                        else:
                            diff = np.linalg.norm(
                                _r_fract-
                                r_fract_com_ref[mol_idx]
                            )
                        if diff < min_diff:
                            min_diff  = diff
                            min_r_fract = _r_fract

            _r_com = np.matmul(uc, min_r_fract.T).T
            r_com_shift = _r_com - r_com
            for atom in molecule:
                query_traj_cp.xyz[i,atom.index] += r_com_shift

    return query_traj_cp



def get_dihedral_indices(
    rdmol,
    exclude_hydrogen=True
    ):

    __doc__ = """
    Retrieve list of dihedrals and list of dihedral ranks.
    """

    from rdkit import Chem
    import numpy as np

    if exclude_hydrogen:
        smarts = Chem.MolFromSmarts("[!#1:1]~[!#1:2]~[!#1:3]~[!#1:4]")
    else:
        smarts = Chem.MolFromSmarts("[*:1]~[*:2]~[*:3]~[*:4]")
    matches = list(
        rdmol.GetSubstructMatches(smarts)
        )
    ranks = list(
        Chem.CanonicalRankAtoms(rdmol, breakTies=False)
        )
    ### First, assign each dihedral a signature using
    ### the rank of the atoms
    signature_list = list()
    index_list = list()
    rank_list = list()
    for match in matches:
        s1 = ranks[match[0]]
        s2 = ranks[match[1]]
        s3 = ranks[match[2]]
        s4 = ranks[match[3]]
        sig1 = (s1,s2,s3,s4)
        sig2 = (s4,s3,s2,s1)
        if sig1 in signature_list:
            idx = signature_list.index(sig1)
            rank_list.append(idx)
        elif sig2 in signature_list:
            idx = signature_list.index(sig2)
            rank_list.append(idx)
        else:
            signature_list.append(sig1)
            idx = signature_list.index(sig1)
            rank_list.append(idx)

        index_list.append(list(match))

    index_list = np.array(index_list, dtype=int)
    rank_list  = np.array(rank_list, dtype=int)

    return index_list, rank_list


def get_angle_indices(
    rdmol,
    exclude_hydrogen=True
    ):

    __doc__ = """
    Retrieve list of angles and list of angle ranks.
    """

    from rdkit import Chem
    import numpy as np

    if exclude_hydrogen:
        smarts = Chem.MolFromSmarts("[!#1:1]~[!#1:2]~[!#1:3]")
    else:
        smarts = Chem.MolFromSmarts("[*:1]~[*:2]~[*:3]")
    matches = list(
        rdmol.GetSubstructMatches(smarts)
        )
    ranks = list(
        Chem.CanonicalRankAtoms(rdmol, breakTies=False)
        )
    ### First, assign each angle a signature using
    ### the rank of the atoms
    signature_list = list()
    index_list = list()
    rank_list = list()
    for match in matches:
        s1 = ranks[match[0]]
        s2 = ranks[match[1]]
        s3 = ranks[match[2]]
        sig1 = (s1,s2,s3)
        sig2 = (s3,s2,s1)
        if sig1 in signature_list:
            idx = signature_list.index(sig1)
            rank_list.append(idx)
        elif sig2 in signature_list:
            idx = signature_list.index(sig2)
            rank_list.append(idx)
        else:
            signature_list.append(sig1)
            idx = signature_list.index(sig1)
            rank_list.append(idx)

        index_list.append(list(match))

    index_list = np.array(index_list, dtype=int)
    rank_list  = np.array(rank_list, dtype=int)

    return index_list, rank_list


def get_bond_indices(
    rdmol,
    exclude_hydrogen=True
    ):

    __doc__ = """
    Retrieve list of bonds and list of bond ranks.
    """

    from rdkit import Chem
    import numpy as np

    if exclude_hydrogen:
        smarts = Chem.MolFromSmarts("[!#1:1]~[!#1:2]")
    else:
        smarts = Chem.MolFromSmarts("[*:1]~[*:2]")
    matches = list(
        rdmol.GetSubstructMatches(smarts)
        )
    ranks = list(
        Chem.CanonicalRankAtoms(rdmol, breakTies=False)
        )
    ### First, assign each bond a signature using
    ### the rank of the atoms
    signature_list = list()
    index_list = list()
    rank_list = list()
    for match in matches:
        s1 = ranks[match[0]]
        s2 = ranks[match[1]]
        sig1 = (s1,s2)
        sig2 = (s2,s1)
        if sig1 in signature_list:
            idx = signature_list.index(sig1)
            rank_list.append(idx)
        elif sig2 in signature_list:
            idx = signature_list.index(sig2)
            rank_list.append(idx)
        else:
            signature_list.append(sig1)
            idx = signature_list.index(sig1)
            rank_list.append(idx)

        index_list.append(list(match))

    index_list = np.array(index_list, dtype=int)
    rank_list  = np.array(rank_list, dtype=int)

    return index_list, rank_list


def compute_dihedrals(
    traj,
    dihedral_indices
    ):

    __doc__ = """
    Compute dihedrals for all atoms in dihedral_indices.
    """

    import mdtraj as md
    import numpy as np

    dihedrals  = md.compute_dihedrals(traj, dihedral_indices)
    dihedrals *= 180./np.pi
    ### Shift dihedrals to range [0,360]
    dihedrals += 180.

    return dihedrals


def compute_angles(
    traj,
    angle_indices
    ):

    __doc__ = """
    Compute angles for all atoms in angle_indices.
    """

    import mdtraj as md
    import numpy as np

    angles  = md.compute_angles(traj, angle_indices)
    angles *= 180./np.pi

    return angles


def compute_bonds(
    traj,
    bond_indices
    ):

    __doc__ = """
    Compute bond lengths for all atoms in bond_indices.
    """

    import mdtraj as md
    import numpy as np

    bond_lengths  = md.compute_distances(traj, bond_indices)

    return bond_lengths


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
    exclude_water=True):

    __doc__ = """
    Compute displacement of com positions for each molecule in query_traj.
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

    com_diff_result = np.zeros(
        (
            query_traj.n_frames, 
            query_traj.n_residues,
        ),
        dtype=float    
    )

    for res in ref_strc.topology.residues:
        if exclude_water and np.any(np.isin(water_list, atm_idxs_per_residue_dict[res.index])):
            continue
        com_query = md.compute_center_of_mass(
            traj = query_traj.atom_slice(atm_idxs_per_residue_dict[res.index])
        )
        com_ref   = md.compute_center_of_mass(
            traj = ref_strc.atom_slice(atm_idxs_per_residue_dict[res.index])
        )
        com_diff_result[:,res.index] = np.linalg.norm(
            com_query - com_ref, 
            axis=1
            )
    return com_diff_result


def compute_rmsd(
    query_traj,
    ref_strc,
    rdmol,
    exclude_hydrogen=True,
    exclude_water=True,
    ):

    __doc__="""
    Compute rmsd for every frame in query_traj w.r.t. ref_strc (frame 0).
    Returns overall rmsd and per residue rmsd.
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
            if exclude_hydrogen and atom.element.atomic_number > 1:
                atm_idxs_per_residue_dict[res.index].append(atom.index)

    N_residues = len(atm_idxs_per_residue_dict)

    rmsd_result = np.zeros(
        (
            query_traj.n_frames, 
        ),
        dtype=float    
    )

    rmsd_per_residue_result = np.zeros(
        (
            query_traj.n_frames, 
            N_residues
        ),
        dtype=float    
    )
    res_count = 0
    atom_count = 0
    for res in ref_strc.topology.residues:
        if exclude_water and np.any(np.isin(water_list, atm_idxs_per_residue_dict[res.index])):
            continue
        xyz_query = query_traj.atom_slice(atm_idxs_per_residue_dict[res.index]).xyz
        xyz_ref   = ref_strc.atom_slice(atm_idxs_per_residue_dict[res.index]).xyz
        diff      = xyz_query - xyz_ref
        diff2_sum = np.sum(diff**2, axis=(1,2))
        rmsd_result[:] += diff2_sum
        rmsd_per_residue_result[:,res_count] = np.sqrt(diff2_sum/res.n_atoms)
        res_count += 1
        atom_count += res.n_atoms

    rmsd_result = rmsd_result/float(atom_count)
    rmsd_result = np.sqrt(rmsd_result)

    return rmsd_result, rmsd_per_residue_result


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
            
    N_residues = len(atm_idxs_per_residue_dict)

    diff_pcs_neighbors = np.zeros(
        (
            query_traj.n_frames,
            N_residues, 
            N_closest_molecules,            
            3
        ), 
        dtype=float
    )
    diff_pcs_self = np.zeros(
        (
            query_traj.n_frames,
            N_residues, 
            3
        ), 
        dtype=float
    )
    res_count = 0
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

        diff_pcs_self[:,res_count] = angular_diffs_self_min

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

        diff_pcs_neighbors[:,res_count] = angular_diffs_query_min

#        ### Serial way to compute differences
#        ### Should give same result as np.einsum
#        for frame_idx in range(query_traj.n_frames):
#            for res_j in range(N_closest_molecules):
#                for pc_idx in [0,1,2]:
#                    dot_ref = np.dot(
#                        pc_ref_big[min_res_idxs_sorted[res_j],pc_idx], 
#                        pc_ref_big[res_count,pc_idx]
#                    )
#                    dot_query = np.dot(
#                        pc_query_big[frame_idx,min_res_idxs_sorted[res_j],pc_idx], 
#                        pc_query_big[frame_idx,res_count,pc_idx]
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
#                    diff_pcs_neighbors[frame_idx, res_count, res_j, pc_idx] = angular_query

        res_count += 1

    return diff_pcs_neighbors, diff_pcs_self