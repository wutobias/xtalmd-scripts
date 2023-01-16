import numpy as np
from simtk.openmm import app
import os
from pyxtal.operations import angle, create_matrix
from pyxtal.constants import deg, rad, ltype_keywords


# This script will calculate RMSD with different options
# Calculate the number of atoms in the unitcell for RMSD calculation
def cal_number_of_atom(pdb_path_1):
    monomer_path = ".." + pdb_path_1.strip(".pdb") + "_monomer0.pdb"
    cif_file = open(monomer_path, 'r')
    count = 0
    for line in cif_file.readlines():
        if line.startswith("HETATM"):
            count += 1
    return count


# Calculate the number of atoms in the supercell for RMSD calculation
def cal_number_of_molecule_in_supercell(pdb_path_1):
    cif_file = open(pdb_path_1, 'r')
    count = 0
    for line in cif_file.readlines():
        if line.startswith("HETATM"):
            count += 1
    number_of_atom = cal_number_of_atom(pdb_path_1)
    number_of_molecule_in_supercell = count/number_of_atom

    return number_of_molecule_in_supercell


# Remove periodicity to avoid periodic boundary conditions misleading RMSD calculation
def rmv_periodicity(pos1Array, pos2Array, box_vector_new):
    """Calculate Remove Periodicity RMSD

    keyword arguments:
    pos1Array(Array)   -- the array of atom position before minimization
    pos2Array(Array)   -- the array of atom position after minimization

    return
    pos2Array_r(Array)   -- the array that had been removed Periodicity
    """

    box_parameter_new = matrix2para(box_vector_new, radians=False)
    pos2Array_r = pos2Array.copy()
    for i in range(len(pos1Array)):
        for j in range(3):
            diff = abs((pos1Array[i][j] - pos2Array_r[i][j]))
            if diff > 0.5 * box_parameter_new[j]:
                if pos2Array_r[i][j] > pos1Array[i][j]:
                    tmp = pos2Array_r[i][j] - box_parameter_new[j]
                elif pos2Array_r[i][j] < pos1Array[i][j]:
                    tmp = pos2Array_r[i][j] + box_parameter_new[j]
                else:
                    tmp = pos2Array_r[i][j]
                pos2Array_r[i][j] = tmp

    return pos2Array_r


# Transform box vectors(matrix) to box parameter(matrix)
def matrix2para(matrix, radians=False):
    """
    Given a 3x3 matrix representing a unit cell, outputs a list of lattice
    parameters.
    Args:
        matrix: a 3x3 array or list, where the first, second, and third rows
            represent the a, b, and c vectors respectively
        radians: if True, outputs angles in radians. If False, outputs in
            degrees
    Returns:
        a 1x6 list of lattice parameters [a, b, c, alpha, beta, gamma]. a, b,
        and c are the length of the lattice vectos, and alpha, beta, and gamma
        are the angles between these vectors (in radians by default)
    """
    cell_para = np.zeros(6)
    # a
    cell_para[0] = np.linalg.norm(matrix[0])
    # b
    cell_para[1] = np.linalg.norm(matrix[1])
    # c
    cell_para[2] = np.linalg.norm(matrix[2])
    # alpha
    cell_para[3] = angle(matrix[1], matrix[2])
    # beta
    cell_para[4] = angle(matrix[0], matrix[2])
    # gamma
    cell_para[5] = angle(matrix[0], matrix[1])

    if not radians:
        # convert radians to degrees
        deg = 180.0 / np.pi
        cell_para[3] *= deg
        cell_para[4] *= deg
        cell_para[5] *= deg
    return cell_para


# RMSD 20 calculation options
def get_rmsd_option(pdb_path_1, pdb_path_2, box_vector_old, box_vector_new, option):
    """Calculate Remove Periodicity RMSD

    keyword arguments:
    pdb_path_1(list)   -- the list of atom position before minimization
    pdb_path_2(list)   -- the list of atom position after minimization

    return:
    option = 0
    rmsd(float)        -- RMSD
    option = r
    rmsd(float)        -- RMSD that remove periodicity
    option = int
    rmsd(list)         -- RMSD for option molecule

    """

    if not os.path.exists(pdb_path_1):
        return None
    if not os.path.exists(pdb_path_2):
        return None

    box_parameter = matrix2para(box_vector_old, radians=False)
    box_parameter_new = matrix2para(box_vector_new, radians=False)
    number_of_atom = cal_number_of_atom(pdb_path_1)
    number_of_molecule_in_supercell = cal_number_of_molecule_in_supercell(pdb_path_1)

    # pos1Array (exp) 2D Array - 1. atom 2. x y z
    # pos2Array (simulation) 2D Array - 1. atom 2. x y z
    pos1Array = np.array(app.PDBFile(pdb_path_1).getPositions(asNumpy=True))
    pos2Array = np.array(app.PDBFile(pdb_path_2).getPositions(asNumpy=True))

    # Record non hydrogen and remove hydrogen for RMSD calculation
    topology = app.PDBFile(pdb_path_1).getTopology()
    non_H_idxs = list()
    for atom_idx, atom in enumerate(topology.atoms()):
        if atom.element.atomic_number != 1:
            non_H_idxs.append(atom_idx)

    # Tobias RMSD (option = o)
    if option == "o":
        diff = np.linalg.norm(
            pos1Array[non_H_idxs] - pos2Array[non_H_idxs],
            axis=1
        )
        rmsd = np.sqrt(
            np.mean(
                diff ** 2
            )
        )


    # Remove Periodicity RMSD (option = r)
    elif option == "r":

        # remove periodicity
        pos2Array_r = rmv_periodicity(pos1Array, pos2Array, box_vector_new)

        diff = np.linalg.norm(
            pos1Array[non_H_idxs] - pos2Array_r[non_H_idxs],
            axis=1
        )

        rmsd = np.sqrt(
            np.mean(
                diff ** 2
            )
        )


    # RMSD_20 Calculation (option = 20)
    elif isinstance(option, int):

        n = option
        rmsd = np.zeros(int(number_of_molecule_in_supercell))

        # remove periodicity
        pos2Array_r = rmv_periodicity(pos1Array, pos2Array, box_vector_new)

        # Remove periodicity for calculation of distance. From first molecule, second, ...
        for target in range(int(number_of_molecule_in_supercell)):
            pos1Array_mol = pos1Array.copy()
            pos2Array_mol = pos2Array_r.copy()
            pos1Array_mol = pos1Array_mol.reshape(int(number_of_molecule_in_supercell), number_of_atom, 3)
            pos2Array_mol = pos2Array_mol.reshape(int(number_of_molecule_in_supercell), number_of_atom, 3)

            for i in range(int(number_of_molecule_in_supercell)):
                for j in range(number_of_atom):
                    for k in range(3):
                        diff = abs((pos2Array_mol[target][j][k] - pos2Array_mol[i][j][k]))
                        if diff > 0.5 * box_parameter_new[k]:
                            if pos2Array_mol[target][j][k] < pos2Array_mol[i][j][k]:
                                tmp1 = pos1Array_mol[i][j][k] - box_parameter[k]
                                tmp2 = pos2Array_mol[i][j][k] - box_parameter_new[k]
                            elif pos2Array_mol[target][j][k] > pos2Array_mol[i][j][k]:
                                tmp1 = pos1Array_mol[i][j][k] + box_parameter[k]
                                tmp2 = pos2Array_mol[i][j][k] + box_parameter_new[k]
                            else:
                                tmp1 = pos1Array_mol[i][j][k]
                                tmp2 = pos2Array_mol[i][j][k]
                            pos1Array_mol[i][j][k] = tmp1
                            pos2Array_mol[i][j][k] = tmp2

            # Calculate average molecule position
            mol_pos = np.mean(pos2Array_mol, axis=1)

            # Distance calculation
            mol_dis = np.zeros(int(number_of_molecule_in_supercell))
            for i in range(int(number_of_molecule_in_supercell)):
                mol_dis_diff = np.array(mol_pos[i]) - np.array(mol_pos[target])
                mol_dis[i] = np.dot(mol_dis_diff, mol_dis_diff)

            # sort and find the index of n the closest molecule
            close_idx = np.argsort(mol_dis)[:n].tolist()

            # Calculate the diff for n closest molecule
            n_pos_diff = np.zeros((int(number_of_molecule_in_supercell), int(number_of_atom), 3))
            for i in range(len(pos1Array_mol)):
                if i in close_idx:
                    n_pos_diff[i] = pos1Array_mol[i] - pos2Array_mol[i]

            # reshape format
            n_pos_diff = n_pos_diff.flatten().reshape(int(number_of_molecule_in_supercell) * number_of_atom, 3)

            diff = np.linalg.norm(
                n_pos_diff[non_H_idxs],
                axis=1
            )

            rmsd_target = np.sqrt(
                np.mean(
                    diff ** 2
                )
            )
            rmsd[target] = rmsd_target
    return rmsd


