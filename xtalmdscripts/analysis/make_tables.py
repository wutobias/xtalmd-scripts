#!/usr/bin/env python3

import xlsxwriter
import numpy as np
import argparse
import gemmi
import yaml
import glob
import os
import copy
import mdtraj as md
from collections import OrderedDict
from rdkit import Chem
from openmm.app.internal.unitcell import computeLengthsAndAngles, reducePeriodicBoxVectors

from . import analysis_engine

_KJ_2_KCAL = 1./4.184
_NM_2_ANG  = 10.
_RAD_2_DEG = 180./np.pi
_DEG_2_RAD = np.pi/180.
_GASCONST_KCAL = 8.31446261815324 * _KJ_2_KCAL / 1000.
_ANG_FACTOR = 0.5

def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Python script for merging simulation data for xtal MD project."
        )

    parser.add_argument('--input',  "-i", type=str, help="Input yaml file.", required=True)
    parser.add_argument('--output', "-o", type=str, help="Output xlsx MS Excel file", required=True)

    return parser.parse_args()


def read_csv(csv_path):

    """
    Read csv file. Save each column as key:value pair in dict.
    """

    data_dict   = OrderedDict()
    header_dict = OrderedDict()
    with open(csv_path, "r") as fopen:
        read_first_line = True
        for line in fopen:
            if read_first_line:
                if line[0] == "#":
                    line = line[1:]
                line = line.lstrip().rstrip().split(",")
                for column, key in enumerate(line):
                    data_dict[key] = list()
                    header_dict[column] = key
                read_first_line = False
                continue
            else:
                line = line.lstrip().rstrip().split(",")
                for column in header_dict:
                    key   = header_dict[column]
                    value = line[column]
                    data_dict[key].append(float(value))
    if len(data_dict) > 0:
        for key in data_dict:
            data_dict[key] = np.array(data_dict[key])
        data_dict["N_rows"] = data_dict[key].size
        data_dict["N_columns"] = len(data_dict[key])
    else:
        data_dict["N_rows"]    = 0
        data_dict["N_columns"] = 0

    return data_dict


class WorkbookWrapper(object):

    def __init__(self, output):

        """
        Constructor for the whole thing. `output` is the xlsx output path on disk.
        """

        all_labels = [
            "Thermo data",
            # ========== #
            "Sublimation Energy",
            "Energy Lattice",
            "Energy Gas",
            "Density",

            "Cell parameters",
            # ============== #
            "a",
            "b",
            "c",
            "alpha",
            "beta",
            "gamma",
            "cell error",

            "RMSD noH",
            # ==================== #
            "<RMSD>",
            "Max <RMSD>",

            "Interatomic distances",
            # ==================== #
            "<(d < 4Å)>",
            "<RMSD [(d < 4Å)]>",
            "Max <RMSD [(d < 4Å)]>",
            "<[Δ(d < 4Å)]>",
            "Max <[Δ(d < 4Å)]>",

            "Bonded Geometry",
            # ==================== #
            "<RMSD τ>",
            "Max <RMSD τ>",
            "<RMSD θ>",
            "Max <RMSD θ>",
            "<RMSD r>",
            "Max <RMSD r>",

            "H-bond geometry abs",
            # ==================== #
            "<[d(D-H•••A)]>",
            "<[∠(D-H•••A)]>",

            "H-bond geometry delta",
            # ==================== #
            "<Δ[d(D-H•••A)]>",
            "Max <Δ[d(D-H•••A)]>",
            "<Δ[∠(D-H•••A)]>",
            "Max <Δ[∠(D-H•••A)]>",
            
            "Translation/Rotation",
            # =================== #
            "<Δ[d(COM)]>",
            "Max <Δ[d(COM)]>",
            "<{∠PA1-PA1(expt)}>",
            "<{∠PA2-PA2(expt)}>",
            "<{∠PA3-PA3(expt)}>",
            "Max <{∠PA1-PA1(expt)}>",
            "Max <{∠PA2-PA2(expt)}>",
            "Max <{∠PA3-PA3(expt)}>",
            "<{∠N_PA1-∠N_PA1(expt)}>",
            "<{∠N_PA2-∠N_PA2(expt)}>",
            "<{∠N_PA3-∠N_PA3(expt)}>",
            "Max <{∠N_PA1-∠N_PA1(expt)}>",
            "Max <{∠N_PA2-∠N_PA2(expt)}>",
            "Max <{∠N_PA3-∠N_PA3(expt)}>",
        ]

        self.workbook  = xlsxwriter.Workbook(
            output,
            {'strings_to_numbers': True}
            )

        ### Define some formats for later use
        ### =================================
        self.header_format_1 = self.workbook.add_format(
            {
                'bold'   : 1,
                'border' : 2,
                'align'  : 'center',
                'valign' : 'vcenter',
                'bottom' : True,
                'top'    : True,
                'left'   : False,
                'right'  : False,
            }
        )

        self.header_format_2 = self.workbook.add_format(
            {
                'bold'   : 1,
                'border' : 1,
                'align'  : 'center',
                'valign' : 'vcenter',
                'bottom' : True,
                'top'    : True,
                'left'   : False,
                'right'  : False,
            }
        )

        self.data_format_1 = self.workbook.add_format(
            {
                'align'  : 'center',
                'valign' : 'vcenter',
                'num_format': '#,##0.00', 
            }
        )

        self.data_format_2 = self.workbook.add_format(
            {
                'align'  : 'center',
                'valign' : 'vcenter',
                'italic' : 1,
            }
        )

        self.labels_dict_row = OrderedDict()
        for label_row, label_name in enumerate(all_labels):
            self.labels_dict_row[label_name] = label_row + 3

        self.sub_titles_row = [
            "Thermo data",
            "Cell parameters",
            "Interatomic distances",
            "RMSD noH",
            "Bonded Geometry",
            "H-bond geometry abs",
            "H-bond geometry delta",
            "Translation/Rotation",
        ]

        self.worksheet_dict   = OrderedDict()
        self.force_field_dict = OrderedDict()

    def add_xtal(self, crystal_name):

        """
        Add a new crystal with name `crystal_name` to the workbook. This will come as a new
        worksheet in the xlsx file.
        """

        if len(crystal_name) > 31:
            crystal_name = crystal_name[:31]

        worksheet = self.workbook.add_worksheet(crystal_name)

        worksheet.set_column(
            0,
            0,
            20.0
        )

        worksheet.set_row(
            0,
            30.0
        )

        ### Make header
        ### ===========
        worksheet.write(
            0,
            0,
            "Crystal/Property",
            self.header_format_1
        )

        ### Write row labels
        ### ================
        worksheet.write(
            1,
            0,
            crystal_name, 
            self.header_format_1
        )

        worksheet.write(
            2,
            0,
            "",
            self.header_format_2
        )

        for label, row_idx in self.labels_dict_row.items():
            if label in self.sub_titles_row:
                worksheet.write(
                    row_idx,
                    0,
                    label,
                    self.data_format_2
                )
            else:
                worksheet.write(
                    row_idx,
                    0,
                    label,
                    self.data_format_1
                )

        self.worksheet_dict[crystal_name] = worksheet
        self.force_field_dict[crystal_name] = list()


    def add_forcefield(self, forcefield_name, crystal_name):

        """
        Add a force field with name `forcefield_name` to worksheet
        of xtal with name `crystal_name`.
        """

        if len(crystal_name) > 31:
            crystal_name = crystal_name[:31]

        self.force_field_dict[crystal_name].append(forcefield_name)
        N_force_fields = len(self.force_field_dict[crystal_name])
        force_field_idx = N_force_fields - 1
        expt_col = 1 + N_force_fields * 3
        worksheet = self.worksheet_dict[crystal_name]

        worksheet.merge_range(
            first_row=1,
            first_col=1 + force_field_idx * 3, 
            last_row=1, 
            last_col=3 + force_field_idx * 3,
            data=forcefield_name,
            cell_format=self.header_format_1
        )

        worksheet.write(
            2,
            1 + force_field_idx * 3,
            "< >",
            self.header_format_2
        )

        worksheet.write(
            2,
            2 + force_field_idx * 3,
            "Std.",
            self.header_format_2
        )

        worksheet.write(
            2,
            3 + force_field_idx * 3,
            "% Dev",
            self.header_format_2
        )


        worksheet.write(
            2,
            expt_col,
            "< >",
            self.header_format_2
        )

        worksheet.write(
            2,
            1 + expt_col,
            "Std.",
            self.header_format_2
        )

        worksheet.write(
            2,
            2 + expt_col,
            "% Dev",
            self.header_format_2
        )

        worksheet.merge_range(
            first_row=0,
            first_col=1, 
            last_row=0, 
            last_col=3 * N_force_fields,
            data="Force Field",
            cell_format=self.header_format_1
        )

        worksheet.merge_range(
            first_row=1,
            first_col=expt_col, 
            last_row=1, 
            last_col=2 + expt_col,
            data="Expt.",
            cell_format=self.header_format_1
        )


    def add_data(
        self,
        data_value, 
        data_std,
        data_dev,
        data_name, 
        forcefield_name, 
        crystal_name):

        """
        Add data of value `data_value` and name `data_name` to worksheet
        of crystal `crystal_name` in column of force field `forcefield_name`.
        """

        if len(crystal_name) > 31:
            crystal_name = crystal_name[:31]

        if forcefield_name == "experiment":
            N_force_fields  = len(self.force_field_dict[crystal_name])
            force_field_idx = N_force_fields
        else:
            force_field_idx = self.force_field_dict[crystal_name].index(forcefield_name)

        worksheet = self.worksheet_dict[crystal_name]
        ### Add in data
        ### ===========
        if isinstance(data_value, str):
            worksheet.write(
                self.labels_dict_row[data_name],
                1 + force_field_idx * 3,
                data_value,
                self.data_format_1
            )
        elif isinstance(data_value, type(None)):
            worksheet.write(
                self.labels_dict_row[data_name],
                1 + force_field_idx * 3,
                "--",
                self.data_format_1
            )
        elif np.isnan(data_value) or np.isinf(data_value):
            worksheet.write(
                self.labels_dict_row[data_name],
                1 + force_field_idx * 3,
                "--",
                self.data_format_1
            )
        else:
            worksheet.write(
                self.labels_dict_row[data_name],
                1 + force_field_idx * 3,
                data_value,
                self.data_format_1
            )

        ### Add in std
        ### ==========
        if isinstance(data_std, str):
            worksheet.write(
                self.labels_dict_row[data_name],
                2 + force_field_idx * 3,
                data_std,
                self.data_format_1
            )
        elif isinstance(data_std, type(None)):
            worksheet.write(
                self.labels_dict_row[data_name],
                2 + force_field_idx * 3,
                "--",
                self.data_format_1
            )
        elif np.isnan(data_std) or np.isinf(data_std):
            worksheet.write(
                self.labels_dict_row[data_name],
                2 + force_field_idx * 3,
                "--",
                self.data_format_1
            )
        else:
            worksheet.write(
                self.labels_dict_row[data_name],
                2 + force_field_idx * 3,
                data_std,
                self.data_format_1
            )

        ### Add in dev
        ### ==========
        if isinstance(data_dev, str):
            worksheet.write(
                self.labels_dict_row[data_name],
                3 + force_field_idx * 3,
                data_dev,
                self.data_format_1
            )
        elif isinstance(data_dev, type(None)):
            worksheet.write(
                self.labels_dict_row[data_name],
                3 + force_field_idx * 3,
                "--",
                self.data_format_1
            )
        elif np.isnan(data_dev) or np.isinf(data_dev):
            worksheet.write(
                self.labels_dict_row[data_name],
                3 + force_field_idx * 3,
                "--",
                self.data_format_1
            )
        else:
            worksheet.write(
                self.labels_dict_row[data_name],
                3 + force_field_idx * 3,
                data_dev,
                self.data_format_1
            )


    def close(self):

        """
        Add deviations from expt. data for each force field
        and then close the workbook.
        """

        data_format_1_left_border = self.workbook.add_format(
            {
                'align'  : 'center',
                'valign' : 'vcenter',
                'num_format': '#,##0.00', 
                'border' : 1,
                'bottom' : False,
                'top'    : False,
                'left'   : True,
                'right'  : False,
            }
        )

        for crystal_name in self.worksheet_dict:
            if len(crystal_name) > 31:
                crystal_name = crystal_name[:31]
            worksheet = self.worksheet_dict[crystal_name]
            force_field_idx  = len(self.force_field_dict[crystal_name])
            expt_val_dict = dict()
            for data_name in self.labels_dict_row:
                row = self.labels_dict_row[data_name]
                row_dict = worksheet.table.get(row, None)
                data_val_col = 1 + force_field_idx * 3
                if data_val_col in row_dict:
                    val = row_dict.get(data_val_col, None)
                    if isinstance(val, xlsxwriter.worksheet.cell_number_tuple):
                        worksheet.write(
                            row,
                            data_val_col,
                            val.number,
                            data_format_1_left_border
                        )
                        expt_val_dict[data_name] = val.number
                    elif data_name in self.sub_titles_row:
                        worksheet.write(
                            row,
                            data_val_col,
                            "",
                            data_format_1_left_border
                        )
                    else:
                        worksheet.write(
                            row,
                            data_val_col,
                            "--",
                            data_format_1_left_border
                        )                        
                else:
                    worksheet.write(
                        row,
                        data_val_col,
                        "",
                        data_format_1_left_border
                    )

            for forcefield_name in self.force_field_dict[crystal_name]:
                force_field_idx = self.force_field_dict[crystal_name].index(forcefield_name)

                for data_name in self.labels_dict_row:
                    row = self.labels_dict_row[data_name]
                    row_dict = worksheet.table.get(row, None)
                    data_val_col = 1 + force_field_idx * 3
                    data_dev_col = 3 + force_field_idx * 3
                    val = row_dict.get(data_val_col, None)
                    if isinstance(val, xlsxwriter.worksheet.cell_number_tuple):
                        ### If the value is present,
                        ### just write it down.
                        worksheet.write(
                            row,
                            data_val_col,
                            val.number,
                            data_format_1_left_border
                        )
                    ### If we are in a title row,
                    ### Don't write anything
                    elif data_name in self.sub_titles_row:
                        worksheet.write(
                            row,
                            data_val_col,
                            "",
                            data_format_1_left_border
                        )
                    ### Else, Just write "--"
                    else:
                        worksheet.write(
                            row,
                            data_val_col,
                            "--",
                            data_format_1_left_border
                        )

                    if data_dev_col in row_dict and isinstance(row_dict.get(data_dev_col, None), xlsxwriter.worksheet.cell_number_tuple):
                        continue
                    elif not data_name in expt_val_dict:
                        worksheet.write(
                            row,
                            data_dev_col,
                            "--",
                            self.data_format_1,
                            )
                    elif data_val_col in row_dict:
                        if isinstance(row_dict.get(data_val_col, None), xlsxwriter.worksheet.cell_number_tuple):
                            ### Compute perc. deviation from experiment
                            dev = row_dict.get(data_val_col, None).number - expt_val_dict[data_name]
                            dev /= expt_val_dict[data_name]
                            dev *= 100.
                            dev  = abs(dev)
                            worksheet.write(
                                row,
                                data_dev_col,
                                dev,
                                self.data_format_1
                            )
                        else:
                            worksheet.write(
                                row,
                                data_dev_col,
                                "--",
                                self.data_format_1,
                                )
                    elif data_name in self.sub_titles_row:
                        worksheet.write(
                            row,
                            data_dev_col,
                            "",
                            self.data_format_1,
                        )                    
                    else:
                        worksheet.write(
                            row,
                            data_dev_col,
                            "--",
                            self.data_format_1
                        )

        self.workbook.close()


def main():

    """
    Main routine to run the analysis workflow.
    """

    args = parse_arguments()
    with open(args.input, "r") as fopen:
        input_dict = yaml.safe_load(fopen)

    workbook_wrap = WorkbookWrapper(args.output)

    ### Loop over each xtal and put in new workbook
    ### ===========================================
    for crystal_name in input_dict:
        print(f"Pre-Processing {crystal_name}")

        workbook_wrap.add_xtal(crystal_name)

        ref_distances = list()
        ref_dihedrals = list()
        ref_pc1_neighbors = list()
        ref_pc2_neighbors = list()
        ref_pc3_neighbors = list()
        ref_pc1_self = list()
        ref_pc2_self = list()
        ref_pc3_self = list()
        ref_hbond_distances = list()
        ref_hbond_angles = list()

        ref_strc   = md.load(
            input_dict[crystal_name]["experiment"]["supercell-pdb"],
            top=input_dict[crystal_name]["experiment"]["supercell-pdb"]
            )
        with open(input_dict[crystal_name]["experiment"]["supercell-rdkit"], "r") as fopen:
            if input_dict[crystal_name]["experiment"]["supercell-rdkit"].lower().endswith(".pdb"):
                rdmol = Chem.MolFromPDBBlock(fopen.read())
            else:
                rdmol = Chem.JSONToMols(fopen.read())[0]

        ucinfo = read_csv(input_dict[crystal_name]["experiment"]["supercell-ucinfo"])

        a_idxs = [ucinfo["unitcell_in_supercell_a"][mol_idx] for mol_idx in range(ucinfo["N_rows"])]
        b_idxs = [ucinfo["unitcell_in_supercell_b"][mol_idx] for mol_idx in range(ucinfo["N_rows"])]
        c_idxs = [ucinfo["unitcell_in_supercell_c"][mol_idx] for mol_idx in range(ucinfo["N_rows"])]

        N_unitcells_a = np.max(a_idxs) - np.min(a_idxs) + 1
        N_unitcells_b = np.max(b_idxs) - np.min(b_idxs) + 1
        N_unitcells_c = np.max(c_idxs) - np.min(c_idxs) + 1

        residue_classes_list = np.array(ucinfo["mol_in_unitcell"], dtype=int)

        ### N_rows is number of molecules in supercell info csv
        N_molecules = ucinfo["N_rows"]

        dist_pair_list = analysis_engine.build_pair_list(
            traj = ref_strc,
            rdmol = rdmol,
            distance_cutoff = 0.4, 
            bond_cutoff = 4, 
            exclude_hydrogen=True,
            )
        dist_pair_rank_list = analysis_engine.build_pair_ranks(
            topology = ref_strc.topology,
            pair_list = dist_pair_list,
            residue_classes_list = residue_classes_list
            )
        ref_distances = analysis_engine.compute_pairwise_distances(
            ref_strc, 
            dist_pair_list
            )
        dihedral_indices, dihedral_ranks = analysis_engine.get_dihedral_indices(
            rdmol
            )
        angle_indices, angle_ranks = analysis_engine.get_angle_indices(
            rdmol
            )
        bond_indices, bond_ranks = analysis_engine.get_bond_indices(
            rdmol
            )
        
        if dihedral_indices.size > 0:
            ref_dihedrals = analysis_engine.compute_dihedrals(
                ref_strc,
                dihedral_indices
                )
        if angle_indices.size > 0:
            ref_angles = analysis_engine.compute_angles(
                ref_strc,
                angle_indices
                )
        if bond_indices.size > 0:
            ref_bonds = analysis_engine.compute_bonds(
                ref_strc,
                bond_indices
                )
        ref_pc_neighbors, ref_pc_self = analysis_engine.compute_pc_diff_per_residue(
            ref_strc, 
            ref_strc,
            rdmol,
            N_closest_molecules=6,
            exclude_water=True
            )

        ref_pc1_neighbors = ref_pc_neighbors[...,0]
        ref_pc2_neighbors = ref_pc_neighbors[...,1]
        ref_pc3_neighbors = ref_pc_neighbors[...,2]

        ref_pc1_self = ref_pc_self[...,0]
        ref_pc2_self = ref_pc_self[...,1]
        ref_pc3_self = ref_pc_self[...,2]

        a_h_d_list = analysis_engine.get_hbond_indices(
            ref_strc,
            rdmol
            )

        if a_h_d_list.size > 0:
            ref_hbond_distances = analysis_engine.compute_pairwise_distances(
                ref_strc, 
                a_h_d_list[:,[1,2]]
                )

            ref_hbond_angles = analysis_engine.compute_tuplewise_angles(
                ref_strc, 
                a_h_d_list
                )

            hbond_pair_rank_list = analysis_engine.build_pair_ranks(
                topology = ref_strc.topology,
                pair_list = a_h_d_list[:,[1,2]],
                residue_classes_list = residue_classes_list
            )

        doc  = gemmi.cif.read(input_dict[crystal_name]["experiment"]["experiment-cif"])[0]
        strc = gemmi.make_small_structure_from_block(doc)
        ref_density = doc.find_value("_exptl_crystal_density_diffrn")
        if ref_density != None:
            try:
                ref_density = float(ref_density)
            except:
                ref_density = "--"
        elif "density" in input_dict[crystal_name]["experiment"]:
            ref_density = input_dict[crystal_name]["experiment"]["density"]
        else:
            ref_density = md.density(ref_strc)[0] / 1000.

        ### Make sure everything is np.ndarray
        ### Also convert to final units.
        #ref_a_len = strc.cell.a
        #ref_b_len = strc.cell.b
        #ref_c_len = strc.cell.c
        #ref_alpha = strc.cell.alpha
        #ref_beta  = strc.cell.beta
        #ref_gamma = strc.cell.gamma

        reduced_box   = reducePeriodicBoxVectors(ref_strc.unitcell_vectors[0])
        length_angles = computeLengthsAndAngles(reduced_box)
        ref_a_len = length_angles[0] / N_unitcells_a * _NM_2_ANG
        ref_b_len = length_angles[1] / N_unitcells_b * _NM_2_ANG
        ref_c_len = length_angles[2] / N_unitcells_c * _NM_2_ANG
        ref_alpha = length_angles[3] * _RAD_2_DEG
        ref_beta  = length_angles[4] * _RAD_2_DEG
        ref_gamma = length_angles[5] * _RAD_2_DEG
        ref_distances = np.array(ref_distances) * _NM_2_ANG
        ref_bonds = np.array(ref_bonds) * _NM_2_ANG
        ref_pc1_neighbors = np.array(ref_pc1_neighbors) * _RAD_2_DEG
        ref_pc2_neighbors = np.array(ref_pc2_neighbors) * _RAD_2_DEG
        ref_pc3_neighbors = np.array(ref_pc3_neighbors) * _RAD_2_DEG
        ref_pc1_self = np.array(ref_pc1_self) * _RAD_2_DEG
        ref_pc2_self = np.array(ref_pc2_self) * _RAD_2_DEG
        ref_pc3_self = np.array(ref_pc3_self) * _RAD_2_DEG
        ref_hbond_distances = np.array(ref_hbond_distances) * _NM_2_ANG
        ref_hbond_angles = np.array(ref_hbond_angles) * _RAD_2_DEG

        ### Carry out analysis for each forcefield
        ### ======================================
        for forcefield_name in input_dict[crystal_name]:
            if forcefield_name.lower() == "experiment":
                continue
            print(f"Processing {crystal_name} / {forcefield_name}")

            ### First check if we have any data
            xtal_trajectory = input_dict[crystal_name][forcefield_name]["xtal-trajectory"]
            xtal_output     = input_dict[crystal_name][forcefield_name]["xtal-output"]

            if len(list(glob.glob(xtal_output))) == 0:
                import warnings
                warnings.warn(
                    f"No xtal csv output found for {crystal_name} / {forcefield_name}. Skipping."
                    )
                continue
            if len(list(glob.glob(xtal_trajectory))) == 0:
                import warnings
                warnings.warn(
                    f"No xtal trajectory output found for {crystal_name} / {forcefield_name}. Skipping."
                    )
                continue

            workbook_wrap.add_forcefield(forcefield_name, crystal_name)

            gas_output_list = list()
            stoichiometry_list = list()
            for key in input_dict[crystal_name][forcefield_name].keys():
                if key.startswith("gas-output"):
                    gasidx = key.split("-")[-1]
                    gas_output = input_dict[crystal_name][forcefield_name][key]
                    stoi       = input_dict[crystal_name][forcefield_name][f"stoichiometry-{gasidx}"]
                    stoi       = float(stoi)
                    stoichiometry_list.append(stoi)
                    gas_output_list.append(gas_output)
            N_components = len(gas_output_list)
            stoichiometry_list = np.array(stoichiometry_list)
            stoichiometry_list /= np.sum(stoichiometry_list)

            a_len = list()
            b_len = list()
            c_len = list()

            alpha = list()
            beta  = list()
            gamma = list()

            ene_xtal = list()
            ene_gas  = list()

            density_xtal = list()

            distances = list()
            dihedrals = list()
            angles    = list()
            bonds     = list()
            com_diff  = list()
            rmsd      = list()
            rmsd_per_residue = list()
            pc1_neighbors = list()
            pc2_neighbors = list()
            pc3_neighbors = list()

            pc1_self = list()
            pc2_self = list()
            pc3_self = list()

            hbond_distances = list()
            hbond_angles = list()

            bad_xtal_list = list()
            bad_gas_list  = list()
            for output_csv in glob.glob(xtal_output):

                data = read_csv(output_csv)
                if data["N_rows"] == 0 or data["N_columns"] == 0:
                    basename, _ = os.path.splitext(output_csv)
                    bad_xtal_list.append(basename)
                    continue                    

                _potential = data["Potential"] # Potential energy in kJ/mol
                if np.isnan(_potential).any():
                    basename, _ = os.path.splitext(output_csv)
                    bad_xtal_list.append(basename)
                    continue
                ene_xtal.extend(_potential.tolist())

            for comp_idx in range(N_components):
                gas_output = gas_output_list[comp_idx]
                stoi       = stoichiometry_list[comp_idx]
                for output_csv in glob.glob(gas_output):

                    data = read_csv(output_csv)
                    if data["N_rows"] == 0 or data["N_columns"] == 0:
                        basename, _ = os.path.splitext(output_csv)
                        bad_gas_list.append(basename)
                        continue

                    _potential = data["Potential"] # Potential energy in kJ/mol
                    if np.isnan(_potential).any():
                        basename, _ = os.path.splitext(output_csv)
                        bad_gas_list.append(basename)
                        continue
                    _potential  = np.array(_potential) * stoi
                    ene_gas.extend(_potential)

            ### For each frame (idx) save the path
            frame_and_path_list = list()
            ### For each frame (idx) save the frame idx in the original traj
            frame_and_sub_frame_list = list()
            for output_traj in glob.glob(xtal_trajectory):

                basename, _ = os.path.splitext(output_traj)
                if basename in bad_xtal_list:
                    continue

                query_traj = md.load(
                    output_traj,
                    top=ref_strc.topology
                    )
                    
                query_traj_wrapped = analysis_engine.unwrap_trajectory(
                    query_traj, 
                    ref_strc, 
                    )

                frame_and_path_list.extend(
                        [output_traj for _ in range(query_traj.n_frames)]
                    )

                frame_and_sub_frame_list.extend(
                        [i for i in range(query_traj.n_frames)]
                    )

                _unitcell_angles = query_traj.unitcell_angles
                alpha.extend(_unitcell_angles[:,0].tolist())
                beta.extend(_unitcell_angles[:,1].tolist())
                gamma.extend(_unitcell_angles[:,2].tolist())

                ### Very important do make deepcopy here.
                ### Otherwise many other calcs on `query_traj` will 
                ### get messed up due the way mdtraj processes things
                ### internally.
                _unitcell_lengths = copy.deepcopy(query_traj.unitcell_lengths)
                ### Correct for number of unitcells along each direction
                _unitcell_lengths[:,0] /= N_unitcells_a
                _unitcell_lengths[:,1] /= N_unitcells_b
                _unitcell_lengths[:,2] /= N_unitcells_c
                a_len.extend(_unitcell_lengths[:,0].tolist())
                b_len.extend(_unitcell_lengths[:,1].tolist())
                c_len.extend(_unitcell_lengths[:,2].tolist())

                _density = md.density(query_traj)
                density_xtal.extend(_density)

                _com_diff = analysis_engine.compute_com_diff_per_residue(
                    query_traj_wrapped, 
                    ref_strc,
                    rdmol,
                    exclude_water=True
                    )
                _pc_neighbors, _pc_self = analysis_engine.compute_pc_diff_per_residue(
                    query_traj_wrapped, 
                    ref_strc,
                    rdmol,
                    N_closest_molecules=6,
                    exclude_water=True
                    )
                _distances = analysis_engine.compute_pairwise_distances(
                    query_traj, 
                    dist_pair_list
                    )
                _rmsd, _rmsd_per_residue = analysis_engine.compute_rmsd(
                        query_traj_wrapped,
                        ref_strc,
                        rdmol,
                        exclude_hydrogen=True,
                        exclude_water=True,
                    )
                if dihedral_indices.size > 0:
                    _dihedrals = analysis_engine.compute_dihedrals(
                        query_traj,
                        dihedral_indices
                        )
                if angle_indices.size > 0:
                    _angles = analysis_engine.compute_angles(
                        query_traj,
                        angle_indices
                        )
                if bond_indices.size > 0:
                    _bonds = analysis_engine.compute_bonds(
                        query_traj,
                        bond_indices
                        )

                if a_h_d_list.size > 0:
                    _hbond_distances = analysis_engine.compute_pairwise_distances(
                        query_traj, 
                        a_h_d_list[:,[1,2]]
                        )

                    _hbond_angles = analysis_engine.compute_tuplewise_angles(
                        query_traj, 
                        a_h_d_list,
                        )

                ### This means, we only do this the first iteration
                if len(com_diff) == 0:
                    com_diff = _com_diff
                    distances = _distances
                    rmsd = _rmsd
                    rmsd_per_residue = _rmsd_per_residue
                    if dihedral_indices.size > 0:
                        dihedrals = _dihedrals
                    if angle_indices.size > 0:
                        angles = _angles
                    if bond_indices.size > 0:
                        bonds = _bonds
                    pc1_neighbors = _pc_neighbors[...,0]
                    pc2_neighbors = _pc_neighbors[...,1]
                    pc3_neighbors = _pc_neighbors[...,2]
                    pc1_self = _pc_self[...,0]
                    pc2_self = _pc_self[...,1]
                    pc3_self = _pc_self[...,2]
                    if a_h_d_list.size > 0:
                        hbond_distances  = _hbond_distances
                        hbond_angles = _hbond_angles

                else:
                    com_diff = np.vstack((com_diff, _com_diff))
                    distances = np.vstack((distances, _distances))
                    rmsd = np.hstack((rmsd, _rmsd))
                    rmsd_per_residue = np.vstack((rmsd_per_residue, _rmsd_per_residue))
                    if dihedral_indices.size > 0:
                        dihedrals = np.vstack((dihedrals, _dihedrals))
                    if angle_indices.size > 0:
                        angles = np.vstack((angles, _angles))
                    if bond_indices.size > 0:
                        bonds = np.vstack((bonds, _bonds))
                    pc1_neighbors = np.vstack((pc1_neighbors, _pc_neighbors[...,0]))
                    pc2_neighbors = np.vstack((pc2_neighbors, _pc_neighbors[...,1]))
                    pc3_neighbors = np.vstack((pc3_neighbors, _pc_neighbors[...,2]))
                    pc1_self = np.vstack((pc1_self, _pc_self[...,0]))
                    pc2_self = np.vstack((pc2_self, _pc_self[...,1]))
                    pc3_self = np.vstack((pc3_self, _pc_self[...,2]))
                    if a_h_d_list.size > 0:
                        hbond_distances  = np.vstack((hbond_distances,  _hbond_distances))
                        hbond_angles = np.vstack((hbond_angles, _hbond_angles))

            ### Make sure everything is np.ndarray
            ### Also convert to final units.
            a_len = np.array(a_len) * _NM_2_ANG
            b_len = np.array(b_len) * _NM_2_ANG
            c_len = np.array(c_len) * _NM_2_ANG
            alpha = np.array(alpha)
            beta  = np.array(beta)
            gamma = np.array(gamma)
            ene_xtal = np.array(ene_xtal) * _KJ_2_KCAL
            ene_gas  = np.array(ene_gas) * _KJ_2_KCAL
            ### Devide by 1000 to get g/cm^3
            density_xtal = np.array(density_xtal) / 1000.
            distances = np.array(distances) * _NM_2_ANG
            bonds     = np.array(bonds) * _NM_2_ANG
            dihedrals = np.array(dihedrals)
            angles    = np.array(angles)
            rmsd = np.array(rmsd) * _NM_2_ANG
            rmsd_per_residue = np.array(rmsd_per_residue) * _NM_2_ANG
            com_diff  = np.array(com_diff) * _NM_2_ANG
            pc1_neighbors = np.array(pc1_neighbors) * _RAD_2_DEG
            pc2_neighbors = np.array(pc2_neighbors) * _RAD_2_DEG
            pc3_neighbors = np.array(pc3_neighbors) * _RAD_2_DEG
            pc1_self = np.array(pc1_self) * _RAD_2_DEG
            pc2_self = np.array(pc2_self) * _RAD_2_DEG
            pc3_self = np.array(pc3_self) * _RAD_2_DEG
            hbond_distances = np.array(hbond_distances) * _NM_2_ANG
            hbond_angles = np.array(hbond_angles) * _RAD_2_DEG

            ### If we don't have any data:
            if distances.size < 2:
                continue

            ### Write distance diff data ###
            ### ======================== ###
            avg   = np.mean(distances)
            std   = np.std(distances)
            dev   = "--"
            workbook_wrap.add_data(
                avg,
                std,
                dev,
                "<(d < 4Å)>",
                forcefield_name, 
                crystal_name
                )

            diffs2   = (distances - ref_distances)**2
            rms_dist = np.sqrt(
                np.mean(
                    diffs2, axis=1
                    )
                )
            avg   = np.mean(rms_dist)
            std   = np.std(rms_dist)
            dev   = "--"
            workbook_wrap.add_data(
                avg,
                std,
                dev,
                "<RMSD [(d < 4Å)]>",
                forcefield_name, 
                crystal_name
                )

            max_avg = 0.
            max_std = 0.
            max_dev = "--"
            for unique_rank in np.unique(dist_pair_rank_list):
                valids   = np.where(unique_rank == dist_pair_rank_list)[0]
                _max_avg = np.sqrt(
                    np.mean(diffs2[:,valids], axis=1)
                    )
                if np.mean(_max_avg) > max_avg:
                    max_avg = np.mean(_max_avg)
                    max_std = np.std(
                        np.sqrt(
                            np.mean(diffs2[:,valids], axis=1)
                            )
                        )

            workbook_wrap.add_data(
                max_avg,
                max_std,
                max_dev,
                "Max <RMSD [(d < 4Å)]>",
                forcefield_name, 
                crystal_name
                )

            diffs = np.abs(distances - ref_distances)
            avg   = np.mean(diffs)
            std   = np.std(diffs)
            dev   = np.mean(diffs/ref_distances) * 100.
            workbook_wrap.add_data(
                avg,
                std,
                dev,
                "<[Δ(d < 4Å)]>",
                forcefield_name, 
                crystal_name
                )

            max_avg = 0.
            max_std = 0.
            max_dev = 0.
            N_max_pairs = 50
            N_max_instances = 10
            max_tolerance = 0.1
            max_avg_list = list()
            max_valids_list = list()
            max_frame_list = list()
            max_pair_list = list()
            max_dist_list = list()
            for unique_rank in np.unique(dist_pair_rank_list):
                valids   = np.where(unique_rank == dist_pair_rank_list)[0]
                _max_avg = np.mean(diffs[:,valids])
                if _max_avg > max_tolerance:
                    max_avg_list.append(_max_avg)
                    max_valids_list.append(valids)
                if _max_avg > max_avg:
                    max_avg = _max_avg
                    max_std = np.std(diffs[:,valids])
                    max_dev = np.mean(diffs[:,valids]/ref_distances[:,valids]) * 100.
            max_idx_list = np.argsort(max_avg_list)[::-1]
            if max_idx_list.size > N_max_pairs:
                max_idx_list = max_idx_list[:N_max_pairs]
            else:
                N_max_pairs = max_idx_list.size
            for max_idx in max_idx_list:
                valids = max_valids_list[max_idx]
                mean_per_frame = np.mean(diffs[:,valids], axis=1)
                max_frame = np.argmax(mean_per_frame)
                ### Sort pairs from highest to lowest deviation
                max_pairs_idx = np.argsort(diffs[max_frame,valids])[::-1]
                if max_pairs_idx.size > N_max_instances:
                    max_pairs_idx = max_pairs_idx[:N_max_instances]
                max_pairs = dist_pair_list[valids[max_pairs_idx]]
                max_dists = diffs[max_frame,valids[max_pairs_idx]]

                max_frame_list.append(max_frame)
                max_pair_list.append(max_pairs)
                max_dist_list.append(max_dists)

            workbook_wrap.add_data(
                max_avg,
                max_std,
                max_dev,
                "Max <[Δ(d < 4Å)]>",
                forcefield_name, 
                crystal_name
                )

            tmp_pdbfile = f"/tmp/{crystal_name}-expt.pdb"
            ref_strc.save(
                tmp_pdbfile
                )
            with open(tmp_pdbfile, "r") as fopen:
                pdbstring = fopen.read()

            uc_ref = ref_strc.unitcell_vectors[0].T
            frac_pos = np.eye(3, dtype=float)
            _abc_ref = np.matmul(uc_ref, frac_pos.T).T * 10.

            pymolstr = f"""#!/usr/bin/env python3

from pymol.cgo import cyl_text, CYLINDER
from pymol import cmd
from pymol.vfont import plain

obj_expt = [\\
   CYLINDER, 0., 0., 0., {_abc_ref[0,0]},{_abc_ref[0,1]},{_abc_ref[0,2]}, 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_ref[1,0]},{_abc_ref[1,1]},{_abc_ref[1,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_ref[2,0]},{_abc_ref[2,1]},{_abc_ref[2,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,\\
   ]

cyl_text(obj_expt,plain,[{_abc_ref[0,0]},{_abc_ref[0,1]},{_abc_ref[0,2]}],'a',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_expt,plain,[{_abc_ref[1,0]},{_abc_ref[1,1]},{_abc_ref[1,2]}],'b',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_expt,plain,[{_abc_ref[2,0]},{_abc_ref[2,1]},{_abc_ref[2,2]}],'c',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])

cmd.load_cgo(obj_expt,'abc_expt')
cmd.cmd.read_pdbstr(
    \"\"\"{pdbstring}\"\"\",
    \"strc_expt\"
)
cmd.group(\"expt\", \"abc_expt\")
cmd.group(\"expt\", \"strc_expt\")
"""
            for max_pair_idx in range(N_max_pairs): 
                
                max_frame = max_frame_list[max_pair_idx]
                max_pairs = max_pair_list[max_pair_idx]
                max_dists = max_dist_list[max_pair_idx]

                query_traj = md.load(
                    frame_and_path_list[max_frame],
                    top=ref_strc.topology
                    )
                query_traj = analysis_engine.unwrap_trajectory(
                    query_traj, 
                    ref_strc, 
                    )

                tmp_pdbfile = f"/tmp/{crystal_name}-{forcefield_name}.pdb"
                query_traj[frame_and_sub_frame_list[max_frame]].save(
                    tmp_pdbfile
                    )
                with open(tmp_pdbfile, "r") as fopen:
                    pdbstring = fopen.read()

                uc_ref = ref_strc.unitcell_vectors[0].T
                uc_md  = query_traj.unitcell_vectors[frame_and_sub_frame_list[max_frame]].T
                frac_pos = np.eye(3, dtype=float)
                _abc_ref = np.matmul(uc_ref, frac_pos.T).T * 10.
                _abc_md  = np.matmul(uc_md, frac_pos.T).T  * 10.

                pymolstr += f"""
obj_md_{max_pair_idx} = [\\
   CYLINDER, 0., 0., 0., {_abc_md[0,0]},{_abc_md[0,1]},{_abc_md[0,2]}, 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_md[1,0]},{_abc_md[1,1]},{_abc_md[1,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_md[2,0]},{_abc_md[2,1]},{_abc_md[2,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,\\
   ]

cyl_text(obj_md_{max_pair_idx},plain,[{_abc_md[0,0]},{_abc_md[0,1]},{_abc_md[0,2]}],'a',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_md_{max_pair_idx},plain,[{_abc_md[1,0]},{_abc_md[1,1]},{_abc_md[1,2]}],'b',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_md_{max_pair_idx},plain,[{_abc_md[2,0]},{_abc_md[2,1]},{_abc_md[2,2]}],'c',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])

cmd.load_cgo(obj_md_{max_pair_idx},  'abc_md_{max_pair_idx}')

cmd.cmd.read_pdbstr(
    \"\"\"{pdbstring}\"\"\",
    \"strc_{max_pair_idx}\"
)
cmd.group(\"pair-{max_pair_idx}\", \"abc_md_{max_pair_idx}\")
cmd.group(\"pair-{max_pair_idx}\", \"strc_{max_pair_idx}\")
"""
                for pair_idx, pair in enumerate(max_pairs):
                    name = f"d_{max_pair_idx}_{pair_idx}-{max_dists[pair_idx]:4.2f}"
                    pymolstr += f"""
cmd.distance(
    name = \"{name}\",
    selection1 = \"rank {pair[0]} and strc_{max_pair_idx}\",
    selection2 = \"rank {pair[1]} and strc_{max_pair_idx}\",
    cutoff = 999.,
    mode = 0
)
cmd.group(\"pair-{max_pair_idx}\", \"{name}\")
"""
            with open(f"{crystal_name}-{forcefield_name}-d-lt4.py", "w") as fopen:
                fopen.write(pymolstr)

            ### Write dihedral data ###
            ### =================== ###

            if dihedral_indices.size > 0:
                diffs = np.abs(ref_dihedrals-dihedrals)
                shifted = np.where(
                    (360.-diffs) < diffs
                    )
                diffs[shifted] = 360.-diffs[shifted]
                diffs2 = diffs**2
                rms_dist = np.sqrt(
                    np.mean(
                        diffs2, axis=1
                        )
                    )
                avg = np.mean(rms_dist)
                std = np.std(rms_dist)
                dev = "--"
                workbook_wrap.add_data(
                    avg,
                    std,
                    dev,
                    "<RMSD τ>",
                    forcefield_name, 
                    crystal_name
                    )

                max_avg = 0.
                max_std = 0.
                max_dev = "--"
                N_max_pairs = 50
                N_max_instances = 10
                max_tolerance = 15.
                max_avg_list = list()
                max_valids_list = list()
                max_frame_list = list()
                max_pair_list = list()
                max_dist_list = list()
                for unique_rank in np.unique(dihedral_ranks):
                    valids   = np.where(unique_rank == dihedral_ranks)[0]
                    _max_avg = np.sqrt(
                        np.mean(diffs2[:,valids], axis=1)
                        )
                    if np.mean(_max_avg) > max_tolerance:
                        max_avg_list.append(np.mean(_max_avg))
                        max_valids_list.append(valids)
                    if np.mean(_max_avg) > max_avg:
                        max_avg = np.mean(_max_avg)
                        max_std = np.std(
                            np.sqrt(
                                np.mean(diffs2[:,valids], axis=1)
                                )
                            )

                max_idx_list = np.argsort(max_avg_list)[::-1]
                if max_idx_list.size > N_max_pairs:
                    max_idx_list = max_idx_list[:N_max_pairs]
                else:
                    N_max_pairs = max_idx_list.size
                for max_idx in max_idx_list:
                    valids = max_valids_list[max_idx]
                    mean_per_frame = np.mean(diffs2[:,valids], axis=1)
                    max_frame = np.argmax(mean_per_frame)
                    ### Sort pairs from highest to lowest deviation
                    max_pairs_idx = np.argsort(diffs2[max_frame,valids])[::-1]
                    if max_pairs_idx.size > N_max_instances:
                        max_pairs_idx = max_pairs_idx[:N_max_instances]
                    max_pairs = dihedral_indices[valids[max_pairs_idx]]
                    max_dists = np.sqrt(diffs2[max_frame,valids[max_pairs_idx]])

                    max_frame_list.append(max_frame)
                    max_pair_list.append(max_pairs)
                    max_dist_list.append(max_dists)

                workbook_wrap.add_data(
                    max_avg,
                    max_std,
                    max_dev,
                    "Max <RMSD τ>",
                    forcefield_name, 
                    crystal_name
                    )

                tmp_pdbfile = f"/tmp/{crystal_name}-expt.pdb"
                ref_strc.save(
                    tmp_pdbfile
                    )
                with open(tmp_pdbfile, "r") as fopen:
                    pdbstring = fopen.read()

                uc_ref = ref_strc.unitcell_vectors[0].T
                frac_pos = np.eye(3, dtype=float)
                _abc_ref = np.matmul(uc_ref, frac_pos.T).T * 10.

                pymolstr = f"""#!/usr/bin/env python3

from pymol.cgo import cyl_text, CYLINDER
from pymol import cmd
from pymol.vfont import plain

obj_expt = [\\
   CYLINDER, 0., 0., 0., {_abc_ref[0,0]},{_abc_ref[0,1]},{_abc_ref[0,2]}, 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_ref[1,0]},{_abc_ref[1,1]},{_abc_ref[1,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_ref[2,0]},{_abc_ref[2,1]},{_abc_ref[2,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,\\
   ]

cyl_text(obj_expt,plain,[{_abc_ref[0,0]},{_abc_ref[0,1]},{_abc_ref[0,2]}],'a',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_expt,plain,[{_abc_ref[1,0]},{_abc_ref[1,1]},{_abc_ref[1,2]}],'b',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_expt,plain,[{_abc_ref[2,0]},{_abc_ref[2,1]},{_abc_ref[2,2]}],'c',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])

cmd.load_cgo(obj_expt,'abc_expt')
cmd.cmd.read_pdbstr(
    \"\"\"{pdbstring}\"\"\",
    \"strc_expt\"
)
cmd.group(\"expt\", \"abc_expt\")
cmd.group(\"expt\", \"strc_expt\")
"""
                for max_pair_idx in range(N_max_pairs): 
                    
                    max_frame = max_frame_list[max_pair_idx]
                    max_pairs = max_pair_list[max_pair_idx]
                    max_dists = max_dist_list[max_pair_idx]

                    query_traj = md.load(
                        frame_and_path_list[max_frame],
                        top=ref_strc.topology
                        )
                    query_traj = analysis_engine.unwrap_trajectory(
                        query_traj, 
                        ref_strc, 
                        )

                    tmp_pdbfile = f"/tmp/{crystal_name}-{forcefield_name}.pdb"
                    query_traj[frame_and_sub_frame_list[max_frame]].save(
                        tmp_pdbfile
                        )
                    with open(tmp_pdbfile, "r") as fopen:
                        pdbstring = fopen.read()

                    uc_ref = ref_strc.unitcell_vectors[0].T
                    uc_md  = query_traj.unitcell_vectors[frame_and_sub_frame_list[max_frame]].T
                    frac_pos = np.eye(3, dtype=float)
                    _abc_ref = np.matmul(uc_ref, frac_pos.T).T * 10.
                    _abc_md  = np.matmul(uc_md, frac_pos.T).T  * 10.

                    pymolstr += f"""
obj_md_{max_pair_idx} = [\\
   CYLINDER, 0., 0., 0., {_abc_md[0,0]},{_abc_md[0,1]},{_abc_md[0,2]}, 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_md[1,0]},{_abc_md[1,1]},{_abc_md[1,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_md[2,0]},{_abc_md[2,1]},{_abc_md[2,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,\\
   ]

cyl_text(obj_md_{max_pair_idx},plain,[{_abc_md[0,0]},{_abc_md[0,1]},{_abc_md[0,2]}],'a',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_md_{max_pair_idx},plain,[{_abc_md[1,0]},{_abc_md[1,1]},{_abc_md[1,2]}],'b',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_md_{max_pair_idx},plain,[{_abc_md[2,0]},{_abc_md[2,1]},{_abc_md[2,2]}],'c',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])

cmd.load_cgo(obj_md_{max_pair_idx},  'abc_md_{max_pair_idx}')

cmd.cmd.read_pdbstr(
    \"\"\"{pdbstring}\"\"\",
    \"strc_{max_pair_idx}\"
)
cmd.group(\"dih-{max_pair_idx}\", \"abc_md_{max_pair_idx}\")
cmd.group(\"dih-{max_pair_idx}\", \"strc_{max_pair_idx}\")
"""
                    for pair_idx, pair in enumerate(max_pairs):
                        name = f"dih_{max_pair_idx}_{pair_idx}-{max_dists[pair_idx]:4.2f}"
                        pymolstr += f"""
cmd.dihedral(
    name = \"{name}\",
    selection1 = \"rank {pair[0]} and strc_{max_pair_idx}\",
    selection2 = \"rank {pair[1]} and strc_{max_pair_idx}\",
    selection3 = \"rank {pair[2]} and strc_{max_pair_idx}\",
    selection4 = \"rank {pair[3]} and strc_{max_pair_idx}\",
)
cmd.group(\"dih-{max_pair_idx}\", \"{name}\")
"""
                with open(f"{crystal_name}-{forcefield_name}-dih.py", "w") as fopen:
                    fopen.write(pymolstr)

            ### Write angle data ###
            ### ================ ###

            if angle_indices.size > 0:
                diffs = np.abs(ref_angles-angles)
                shifted = np.where(
                    (360.-diffs) < diffs
                    )
                diffs[shifted] = 360.-diffs[shifted]
                diffs2 = diffs**2
                rms_dist = np.sqrt(
                    np.mean(
                        diffs2, axis=1
                        )
                    )
                avg = np.mean(rms_dist)
                std = np.std(rms_dist)
                dev = "--"
                workbook_wrap.add_data(
                    avg,
                    std,
                    dev,
                    "<RMSD θ>",
                    forcefield_name, 
                    crystal_name
                    )

                max_avg = 0.
                max_std = 0.
                max_dev = "--"
                N_max_pairs = 50
                N_max_instances = 10
                max_tolerance = 10.
                max_avg_list = list()
                max_valids_list = list()
                max_frame_list = list()
                max_pair_list = list()
                max_dist_list = list()
                for unique_rank in np.unique(angle_ranks):
                    valids   = np.where(unique_rank == angle_ranks)[0]
                    _max_avg = np.sqrt(
                        np.mean(diffs2[:,valids], axis=1)
                        )
                    if np.mean(_max_avg) > max_tolerance:
                        max_avg_list.append(np.mean(_max_avg))
                        max_valids_list.append(valids)
                    if np.mean(_max_avg) > max_avg:
                        max_avg = np.mean(_max_avg)
                        max_std = np.std(
                            np.sqrt(
                                np.mean(diffs2[:,valids], axis=1)
                                )
                            )

                max_idx_list = np.argsort(max_avg_list)[::-1]
                if max_idx_list.size > N_max_pairs:
                    max_idx_list = max_idx_list[:N_max_pairs]
                else:
                    N_max_pairs = max_idx_list.size
                for max_idx in max_idx_list:
                    valids = max_valids_list[max_idx]
                    mean_per_frame = np.mean(diffs2[:,valids], axis=1)
                    max_frame = np.argmax(mean_per_frame)
                    ### Sort pairs from highest to lowest deviation
                    max_pairs_idx = np.argsort(diffs2[max_frame,valids])[::-1]
                    if max_pairs_idx.size > N_max_instances:
                        max_pairs_idx = max_pairs_idx[:N_max_instances]
                    max_pairs = angle_indices[valids[max_pairs_idx]]
                    max_dists = np.sqrt(diffs2[max_frame,valids[max_pairs_idx]])

                    max_frame_list.append(max_frame)
                    max_pair_list.append(max_pairs)
                    max_dist_list.append(max_dists)

                workbook_wrap.add_data(
                    max_avg,
                    max_std,
                    max_dev,
                    "Max <RMSD θ>",
                    forcefield_name, 
                    crystal_name
                    )

                tmp_pdbfile = f"/tmp/{crystal_name}-expt.pdb"
                ref_strc.save(
                    tmp_pdbfile
                    )
                with open(tmp_pdbfile, "r") as fopen:
                    pdbstring = fopen.read()

                uc_ref = ref_strc.unitcell_vectors[0].T
                frac_pos = np.eye(3, dtype=float)
                _abc_ref = np.matmul(uc_ref, frac_pos.T).T * 10.

                pymolstr = f"""#!/usr/bin/env python3

from pymol.cgo import cyl_text, CYLINDER
from pymol import cmd
from pymol.vfont import plain

obj_expt = [\\
   CYLINDER, 0., 0., 0., {_abc_ref[0,0]},{_abc_ref[0,1]},{_abc_ref[0,2]}, 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_ref[1,0]},{_abc_ref[1,1]},{_abc_ref[1,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_ref[2,0]},{_abc_ref[2,1]},{_abc_ref[2,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,\\
   ]

cyl_text(obj_expt,plain,[{_abc_ref[0,0]},{_abc_ref[0,1]},{_abc_ref[0,2]}],'a',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_expt,plain,[{_abc_ref[1,0]},{_abc_ref[1,1]},{_abc_ref[1,2]}],'b',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_expt,plain,[{_abc_ref[2,0]},{_abc_ref[2,1]},{_abc_ref[2,2]}],'c',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])

cmd.load_cgo(obj_expt,'abc_expt')
cmd.cmd.read_pdbstr(
    \"\"\"{pdbstring}\"\"\",
    \"strc_expt\"
)
cmd.group(\"expt\", \"abc_expt\")
cmd.group(\"expt\", \"strc_expt\")
"""
                for max_pair_idx in range(N_max_pairs): 
                    
                    max_frame = max_frame_list[max_pair_idx]
                    max_pairs = max_pair_list[max_pair_idx]
                    max_dists = max_dist_list[max_pair_idx]

                    query_traj = md.load(
                        frame_and_path_list[max_frame],
                        top=ref_strc.topology
                        )
                    query_traj = analysis_engine.unwrap_trajectory(
                        query_traj, 
                        ref_strc, 
                        )

                    tmp_pdbfile = f"/tmp/{crystal_name}-{forcefield_name}.pdb"
                    query_traj[frame_and_sub_frame_list[max_frame]].save(
                        tmp_pdbfile
                        )
                    with open(tmp_pdbfile, "r") as fopen:
                        pdbstring = fopen.read()

                    uc_ref = ref_strc.unitcell_vectors[0].T
                    uc_md  = query_traj.unitcell_vectors[frame_and_sub_frame_list[max_frame]].T
                    frac_pos = np.eye(3, dtype=float)
                    _abc_ref = np.matmul(uc_ref, frac_pos.T).T * 10.
                    _abc_md  = np.matmul(uc_md, frac_pos.T).T  * 10.

                    pymolstr += f"""
obj_md_{max_pair_idx} = [\\
   CYLINDER, 0., 0., 0., {_abc_md[0,0]},{_abc_md[0,1]},{_abc_md[0,2]}, 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_md[1,0]},{_abc_md[1,1]},{_abc_md[1,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_md[2,0]},{_abc_md[2,1]},{_abc_md[2,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,\\
   ]

cyl_text(obj_md_{max_pair_idx},plain,[{_abc_md[0,0]},{_abc_md[0,1]},{_abc_md[0,2]}],'a',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_md_{max_pair_idx},plain,[{_abc_md[1,0]},{_abc_md[1,1]},{_abc_md[1,2]}],'b',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_md_{max_pair_idx},plain,[{_abc_md[2,0]},{_abc_md[2,1]},{_abc_md[2,2]}],'c',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])

cmd.load_cgo(obj_md_{max_pair_idx},  'abc_md_{max_pair_idx}')

cmd.cmd.read_pdbstr(
    \"\"\"{pdbstring}\"\"\",
    \"strc_{max_pair_idx}\"
)
cmd.group(\"ang-{max_pair_idx}\", \"abc_md_{max_pair_idx}\")
cmd.group(\"ang-{max_pair_idx}\", \"strc_{max_pair_idx}\")
"""
                    for pair_idx, pair in enumerate(max_pairs):
                        name = f"ang_{max_pair_idx}_{pair_idx}-{max_dists[pair_idx]:4.2f}"
                        pymolstr += f"""
cmd.angle(
    name = \"{name}\",
    selection1 = \"rank {pair[0]} and strc_{max_pair_idx}\",
    selection2 = \"rank {pair[1]} and strc_{max_pair_idx}\",
    selection3 = \"rank {pair[2]} and strc_{max_pair_idx}\",
)
cmd.group(\"ang-{max_pair_idx}\", \"{name}\")
"""
                with open(f"{crystal_name}-{forcefield_name}-ang.py", "w") as fopen:
                    fopen.write(pymolstr)

            
            ### Write bond length data ###
            ### ====================== ###

            if bond_indices.size > 0:
                diffs2   = (bonds - ref_bonds)**2
                rms_dist = np.sqrt(
                    np.mean(
                        diffs2, axis=1
                        )
                    )
                avg   = np.mean(rms_dist)
                std   = np.std(rms_dist)
                dev   = "--"

                workbook_wrap.add_data(
                    avg,
                    std,
                    dev,
                    "<RMSD r>",
                    forcefield_name, 
                    crystal_name
                    )

                max_avg = 0.
                max_std = 0.
                max_dev = "--"
                N_max_pairs = 50
                N_max_instances = 10
                max_tolerance = 1.e-2
                max_avg_list = list()
                max_valids_list = list()
                max_frame_list = list()
                max_pair_list = list()
                max_dist_list = list()
                for unique_rank in np.unique(bond_ranks):
                    valids   = np.where(unique_rank == bond_ranks)[0]
                    _max_avg = np.sqrt(
                        np.mean(diffs2[:,valids], axis=1)
                        )
                    if np.mean(_max_avg) > max_tolerance:
                        max_avg_list.append(np.mean(_max_avg))
                        max_valids_list.append(valids)
                    if np.mean(_max_avg) > max_avg:
                        max_avg = np.mean(_max_avg)
                        max_std = np.std(
                            np.sqrt(
                                np.mean(diffs2[:,valids], axis=1)
                                )
                            )
                max_idx_list = np.argsort(max_avg_list)[::-1]
                if max_idx_list.size > N_max_pairs:
                    max_idx_list = max_idx_list[:N_max_pairs]
                else:
                    N_max_pairs = max_idx_list.size
                for max_idx in max_idx_list:
                    valids = max_valids_list[max_idx]
                    mean_per_frame = np.mean(diffs2[:,valids], axis=1)
                    max_frame = np.argmax(mean_per_frame)
                    ### Sort pairs from highest to lowest deviation
                    max_pairs_idx = np.argsort(diffs2[max_frame,valids])[::-1]
                    if max_pairs_idx.size > N_max_instances:
                        max_pairs_idx = max_pairs_idx[:N_max_instances]
                    max_pairs = bond_indices[valids[max_pairs_idx]]
                    max_dists = np.sqrt(diffs2[max_frame,valids[max_pairs_idx]])

                    max_frame_list.append(max_frame)
                    max_pair_list.append(max_pairs)
                    max_dist_list.append(max_dists)

                workbook_wrap.add_data(
                    max_avg,
                    max_std,
                    max_dev,
                    "Max <RMSD r>",
                    forcefield_name, 
                    crystal_name
                    )

                tmp_pdbfile = f"/tmp/{crystal_name}-expt.pdb"
                ref_strc.save(
                    tmp_pdbfile
                    )
                with open(tmp_pdbfile, "r") as fopen:
                    pdbstring = fopen.read()

                uc_ref = ref_strc.unitcell_vectors[0].T
                frac_pos = np.eye(3, dtype=float)
                _abc_ref = np.matmul(uc_ref, frac_pos.T).T * 10.

                pymolstr = f"""#!/usr/bin/env python3

from pymol.cgo import cyl_text, CYLINDER
from pymol import cmd
from pymol.vfont import plain

obj_expt = [\\
   CYLINDER, 0., 0., 0., {_abc_ref[0,0]},{_abc_ref[0,1]},{_abc_ref[0,2]}, 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_ref[1,0]},{_abc_ref[1,1]},{_abc_ref[1,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_ref[2,0]},{_abc_ref[2,1]},{_abc_ref[2,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,\\
   ]

cyl_text(obj_expt,plain,[{_abc_ref[0,0]},{_abc_ref[0,1]},{_abc_ref[0,2]}],'a',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_expt,plain,[{_abc_ref[1,0]},{_abc_ref[1,1]},{_abc_ref[1,2]}],'b',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_expt,plain,[{_abc_ref[2,0]},{_abc_ref[2,1]},{_abc_ref[2,2]}],'c',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])

cmd.load_cgo(obj_expt,'abc_expt')
cmd.cmd.read_pdbstr(
    \"\"\"{pdbstring}\"\"\",
    \"strc_expt\"
)
cmd.group(\"expt\", \"abc_expt\")
cmd.group(\"expt\", \"strc_expt\")
"""
                for max_pair_idx in range(N_max_pairs): 
                    
                    max_frame = max_frame_list[max_pair_idx]
                    max_pairs = max_pair_list[max_pair_idx]
                    max_dists = max_dist_list[max_pair_idx]

                    query_traj = md.load(
                        frame_and_path_list[max_frame],
                        top=ref_strc.topology
                        )
                    query_traj = analysis_engine.unwrap_trajectory(
                        query_traj, 
                        ref_strc, 
                        )

                    tmp_pdbfile = f"/tmp/{crystal_name}-{forcefield_name}.pdb"
                    query_traj[frame_and_sub_frame_list[max_frame]].save(
                        tmp_pdbfile
                        )
                    with open(tmp_pdbfile, "r") as fopen:
                        pdbstring = fopen.read()

                    uc_ref = ref_strc.unitcell_vectors[0].T
                    uc_md  = query_traj.unitcell_vectors[frame_and_sub_frame_list[max_frame]].T
                    frac_pos = np.eye(3, dtype=float)
                    _abc_ref = np.matmul(uc_ref, frac_pos.T).T * 10.
                    _abc_md  = np.matmul(uc_md, frac_pos.T).T  * 10.

                    pymolstr += f"""
obj_md_{max_pair_idx} = [\\
   CYLINDER, 0., 0., 0., {_abc_md[0,0]},{_abc_md[0,1]},{_abc_md[0,2]}, 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_md[1,0]},{_abc_md[1,1]},{_abc_md[1,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_md[2,0]},{_abc_md[2,1]},{_abc_md[2,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,\\
   ]

cyl_text(obj_md_{max_pair_idx},plain,[{_abc_md[0,0]},{_abc_md[0,1]},{_abc_md[0,2]}],'a',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_md_{max_pair_idx},plain,[{_abc_md[1,0]},{_abc_md[1,1]},{_abc_md[1,2]}],'b',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_md_{max_pair_idx},plain,[{_abc_md[2,0]},{_abc_md[2,1]},{_abc_md[2,2]}],'c',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])

cmd.load_cgo(obj_md_{max_pair_idx},  'abc_md_{max_pair_idx}')

cmd.cmd.read_pdbstr(
    \"\"\"{pdbstring}\"\"\",
    \"strc_{max_pair_idx}\"
)
cmd.group(\"bon-{max_pair_idx}\", \"abc_md_{max_pair_idx}\")
cmd.group(\"bon-{max_pair_idx}\", \"strc_{max_pair_idx}\")
"""
                    for pair_idx, pair in enumerate(max_pairs):
                        name = f"bon_{max_pair_idx}_{pair_idx}-{max_dists[pair_idx]:4.2f}"
                        pymolstr += f"""
cmd.distance(
    name = \"{name}\",
    selection1 = \"rank {pair[0]} and strc_{max_pair_idx}\",
    selection2 = \"rank {pair[1]} and strc_{max_pair_idx}\",
    cutoff = 999.,
    mode = 0
)
cmd.group(\"bon-{max_pair_idx}\", \"{name}\")
"""
                with open(f"{crystal_name}-{forcefield_name}-bon.py", "w") as fopen:
                    fopen.write(pymolstr)


            ### Write com diff data ###
            ### =================== ###
            diffs = com_diff
            avg   = np.mean(diffs)
            std   = np.std(diffs)
            dev   = "--"
            workbook_wrap.add_data(
                avg ,
                std,
                dev,
                "<Δ[d(COM)]>",
                forcefield_name, 
                crystal_name
                )

            res_avg = np.mean(diffs, axis=0)
            max_idx = np.unravel_index(
                np.argmax(res_avg),
                res_avg.shape
                )
            max_avg = res_avg[max_idx]
            max_std = np.std(diffs[:,max_idx])
            max_dev = "--"
            workbook_wrap.add_data(
                max_avg,
                max_std,
                max_dev,
                "Max <Δ[d(COM)]>",
                forcefield_name, 
                crystal_name
                )

            ### Write rmsd data ###
            ### =============== ###
            avg   = np.mean(rmsd)
            std   = np.std(rmsd)
            dev   = None
            workbook_wrap.add_data(
                avg ,
                std,
                dev,
                "<RMSD>",
                forcefield_name, 
                crystal_name
                )

            avg     = np.mean(rmsd_per_residue, axis=0)
            max_idx = np.argmax(avg)
            max_avg = avg[max_idx]
            max_std = np.std(rmsd_per_residue[:,max_idx])
            max_dev = None
            workbook_wrap.add_data(
                max_avg,
                max_std,
                max_dev,
                "Max <RMSD>",
                forcefield_name, 
                crystal_name
                )

            ### Write pc diff data ###
            ### ================== ###

            ### PC1 data
            ### --------
            diffs = np.abs(pc1_self)
            avg   = np.mean(diffs)
            std   = np.std(diffs)
            dev   = "--"
            workbook_wrap.add_data(
                avg,
                std,
                dev,
                "<{∠PA1-PA1(expt)}>", 
                forcefield_name, 
                crystal_name
                )

            res_avg = np.mean(diffs, axis=0)
            max_idx = np.unravel_index(
                np.argmax(res_avg),
                res_avg.shape
                )
            max_avg = res_avg[max_idx]
            max_std = np.std(diffs[:,max_idx])
            max_dev = "--"
            workbook_wrap.add_data(
                max_avg,
                max_std,
                max_dev,
                "Max <{∠PA1-PA1(expt)}>", 
                forcefield_name, 
                crystal_name
                )

            diffs = np.abs(ref_pc1_neighbors-pc1_neighbors)
            shifted = np.where(
                (360.-diffs) < diffs
                )
            diffs[shifted] = 360.-diffs[shifted]
            avg   = np.mean(diffs)
            std   = np.std(diffs)
            dev   = "--"
            workbook_wrap.add_data(
                avg,
                std,
                dev,
                "<{∠N_PA1-∠N_PA1(expt)}>",
                forcefield_name, 
                crystal_name
                )

            res_avg = np.mean(diffs, axis=0)
            max_idx = np.unravel_index(
                np.argmax(res_avg),
                res_avg.shape
                )
            max_avg = res_avg[max_idx]
            max_std = np.std(diffs[:,max_idx])
            max_dev = "--"
            workbook_wrap.add_data(
                max_avg,
                max_std,
                max_dev,
                "Max <{∠N_PA1-∠N_PA1(expt)}>",
                forcefield_name, 
                crystal_name
                )

            ### PC2 data
            ### --------
            diffs = np.abs(pc2_self)
            avg   = np.mean(diffs)
            std   = np.std(diffs)
            dev   = "--"
            workbook_wrap.add_data(
                avg,
                std,
                dev,
                "<{∠PA2-PA2(expt)}>", 
                forcefield_name, 
                crystal_name
                )

            res_avg = np.mean(diffs, axis=0)
            max_idx = np.unravel_index(
                np.argmax(res_avg),
                res_avg.shape
                )
            max_avg = res_avg[max_idx]
            max_std = np.std(diffs[:,max_idx])
            max_dev = "--"
            workbook_wrap.add_data(
                max_avg,
                max_std,
                max_dev,
                "Max <{∠PA2-PA2(expt)}>", 
                forcefield_name, 
                crystal_name
                )

            diffs = np.abs(ref_pc2_neighbors-pc2_neighbors)
            shifted = np.where(
                (360.-diffs) < diffs
                )
            diffs[shifted] = 360.-diffs[shifted]
            avg   = np.mean(diffs)
            std   = np.std(diffs)
            dev   = "--"
            workbook_wrap.add_data(
                avg,
                std,
                dev,
                "<{∠N_PA2-∠N_PA2(expt)}>",
                forcefield_name, 
                crystal_name
                )

            res_avg = np.mean(diffs, axis=0)
            max_idx = np.unravel_index(
                np.argmax(res_avg),
                res_avg.shape
                )
            max_avg = res_avg[max_idx]
            max_std = np.std(diffs[:,max_idx])
            max_dev = "--"
            workbook_wrap.add_data(
                max_avg,
                max_std,
                max_dev,
                "Max <{∠N_PA2-∠N_PA2(expt)}>",
                forcefield_name, 
                crystal_name
                )

            ### PC3 data
            ### --------
            diffs = np.abs(pc3_self)
            avg   = np.mean(diffs)
            std   = np.std(diffs)
            dev   = "--"
            workbook_wrap.add_data(
                avg,
                std,
                dev,
                "<{∠PA3-PA3(expt)}>", 
                forcefield_name, 
                crystal_name
                )

            res_avg = np.mean(diffs, axis=0)
            max_idx = np.unravel_index(
                np.argmax(res_avg),
                res_avg.shape
                )
            max_avg = res_avg[max_idx]
            max_std = np.std(diffs[:,max_idx])
            #max_dev = np.mean(diffs[:,max_idx]/ref_pc3[:,max_idx]) * 100.
            max_dev = "--"
            workbook_wrap.add_data(
                max_avg,
                max_std,
                max_dev,
                "Max <{∠PA3-PA3(expt)}>", 
                forcefield_name, 
                crystal_name
                )

            diffs = np.abs(ref_pc3_neighbors-pc3_neighbors)
            shifted = np.where(
                (360.-diffs) < diffs
                )
            diffs[shifted] = 360.-diffs[shifted]
            avg   = np.mean(diffs)
            std   = np.std(diffs)
            dev   = "--"
            workbook_wrap.add_data(
                avg,
                std,
                dev,
                "<{∠N_PA3-∠N_PA3(expt)}>",
                forcefield_name, 
                crystal_name
                )

            res_avg = np.mean(diffs, axis=0)
            max_idx = np.unravel_index(
                np.argmax(res_avg),
                res_avg.shape
                )
            max_avg = res_avg[max_idx]
            max_std = np.std(diffs[:,max_idx])
            max_dev = "--"
            workbook_wrap.add_data(
                max_avg,
                max_std,
                max_dev,
                "Max <{∠N_PA3-∠N_PA3(expt)}>",
                forcefield_name, 
                crystal_name
                )

            ### ================ ###
            ### Write hbond data ###
            ### ================ ###

            ### ... Absolute Distances ...

            if len(hbond_distances) > 0:
                avg   = np.mean(hbond_distances)
                std   = np.std(hbond_distances)
                dev   = "--"
            else:
                avg = "--"
                std = "--"
                dev = "--"
            workbook_wrap.add_data(
                avg,
                std,
                dev,
                "<[d(D-H•••A)]>",
                forcefield_name, 
                crystal_name
                )

            ### ... Absolute angles ...

            if len(hbond_angles) > 0:
                avg   = np.mean(np.abs(hbond_angles))
                std   = np.std(np.abs(hbond_angles))
                dev   = "--"
            else:
                avg = "--"
                std = "--"
                dev = "--"
            workbook_wrap.add_data(
                avg,
                std,
                dev,
                "<[∠(D-H•••A)]>",
                forcefield_name, 
                crystal_name
                )

            ### =================== ###
            ### Write distance data ###
            ### =================== ###

            ### ... Less then 4 Ang Distances ...

            avg   = np.mean(distances)
            std   = np.std(distances)
            dev   = "--"
            workbook_wrap.add_data(
                avg,
                std,
                dev,
                "<(d < 4Å)>",
                forcefield_name, 
                crystal_name
                )

            diffs = np.abs(distances - ref_distances)
            avg   = np.mean(diffs)
            std   = np.std(diffs)
            dev   = np.mean(diffs/ref_distances) * 100.
            workbook_wrap.add_data(
                avg,
                std,
                dev,
                "<[Δ(d < 4Å)]>",
                forcefield_name, 
                crystal_name
                )

            if len(hbond_distances) > 0:
                diffs = np.abs(
                    hbond_distances - ref_hbond_distances
                    )
                avg   = np.mean(diffs)
                std   = np.std(diffs)
                dev   = np.mean(diffs/ref_hbond_distances) * 100.

                max_avg = 0.
                max_std = 0.
                max_dev = 0.
                N_max_pairs = 50
                N_max_instances = 10
                max_tolerance = 0.1
                max_avg_list = list()
                max_valids_list = list()
                max_frame_list = list()
                max_pair_list = list()
                max_dist_list = list()
                for unique_rank in np.unique(hbond_pair_rank_list):
                    valids   = np.where(unique_rank == hbond_pair_rank_list)[0]
                    _max_avg = np.mean(diffs[:,valids])
                    if _max_avg > max_tolerance:
                        max_avg_list.append(_max_avg)
                        max_valids_list.append(valids)
                    if _max_avg > max_avg:
                        max_avg = _max_avg
                        max_std = np.std(diffs[:,valids])
                        max_dev = np.mean(diffs[:,valids]/ref_hbond_distances[:,valids]) * 100.
                max_idx_list = np.argsort(max_avg_list)[::-1]
                if max_idx_list.size > N_max_pairs:
                    max_idx_list = max_idx_list[:N_max_pairs]
                else:
                    N_max_pairs = max_idx_list.size
                for max_idx in max_idx_list:
                    valids = max_valids_list[max_idx]
                    mean_per_frame = np.mean(diffs[:,valids], axis=1)
                    max_frame = np.argmax(mean_per_frame)
                    ### Sort pairs from highest to lowest deviation
                    max_pairs_idx = np.argsort(diffs[max_frame,valids])[::-1]
                    if max_pairs_idx.size > N_max_instances:
                        max_pairs_idx = max_pairs_idx[:N_max_instances]
                    max_pairs = a_h_d_list[valids[max_pairs_idx]][:,[1,2]]
                    max_dists = diffs[max_frame,valids[max_pairs_idx]]

                    max_frame_list.append(max_frame)
                    max_pair_list.append(max_pairs)
                    max_dist_list.append(max_dists)

                workbook_wrap.add_data(
                    max_avg,
                    max_std,
                    max_dev,
                    "Max <Δ[d(D-H•••A)]>",
                    forcefield_name, 
                    crystal_name
                    )

                tmp_pdbfile = f"/tmp/{crystal_name}-expt.pdb"
                ref_strc.save(
                    tmp_pdbfile
                    )
                with open(tmp_pdbfile, "r") as fopen:
                    pdbstring = fopen.read()

                uc_ref = ref_strc.unitcell_vectors[0].T
                frac_pos = np.eye(3, dtype=float)
                _abc_ref = np.matmul(uc_ref, frac_pos.T).T * 10.

                pymolstr = f"""#!/usr/bin/env python3

from pymol.cgo import cyl_text, CYLINDER
from pymol import cmd
from pymol.vfont import plain

obj_expt = [\\
   CYLINDER, 0., 0., 0., {_abc_ref[0,0]},{_abc_ref[0,1]},{_abc_ref[0,2]}, 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_ref[1,0]},{_abc_ref[1,1]},{_abc_ref[1,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_ref[2,0]},{_abc_ref[2,1]},{_abc_ref[2,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,\\
   ]

cyl_text(obj_expt,plain,[{_abc_ref[0,0]},{_abc_ref[0,1]},{_abc_ref[0,2]}],'a',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_expt,plain,[{_abc_ref[1,0]},{_abc_ref[1,1]},{_abc_ref[1,2]}],'b',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_expt,plain,[{_abc_ref[2,0]},{_abc_ref[2,1]},{_abc_ref[2,2]}],'c',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])

cmd.load_cgo(obj_expt,'abc_expt')
cmd.cmd.read_pdbstr(
    \"\"\"{pdbstring}\"\"\",
    \"strc_expt\"
)
cmd.group(\"expt\", \"abc_expt\")
cmd.group(\"expt\", \"strc_expt\")
"""
                for max_pair_idx in range(N_max_pairs): 
                    
                    max_frame = max_frame_list[max_pair_idx]
                    max_pairs = max_pair_list[max_pair_idx]
                    max_dists = max_dist_list[max_pair_idx]

                    query_traj = md.load(
                        frame_and_path_list[max_frame],
                        top=ref_strc.topology
                        )
                    query_traj = analysis_engine.unwrap_trajectory(
                        query_traj, 
                        ref_strc,
                        )

                    tmp_pdbfile = f"/tmp/{crystal_name}-{forcefield_name}.pdb"
                    query_traj[frame_and_sub_frame_list[max_frame]].save(
                        tmp_pdbfile
                        )
                    with open(tmp_pdbfile, "r") as fopen:
                        pdbstring = fopen.read()

                    uc_ref = ref_strc.unitcell_vectors[0].T
                    uc_md  = query_traj.unitcell_vectors[frame_and_sub_frame_list[max_frame]].T
                    frac_pos = np.eye(3, dtype=float)
                    _abc_ref = np.matmul(uc_ref, frac_pos.T).T * 10.
                    _abc_md  = np.matmul(uc_md, frac_pos.T).T  * 10.

                    pymolstr += f"""
obj_md_{max_pair_idx} = [\\
   CYLINDER, 0., 0., 0., {_abc_md[0,0]},{_abc_md[0,1]},{_abc_md[0,2]}, 0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_md[1,0]},{_abc_md[1,1]},{_abc_md[1,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0,\\
   CYLINDER, 0., 0., 0., {_abc_md[2,0]},{_abc_md[2,1]},{_abc_md[2,2]}, 0.2, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,\\
   ]

cyl_text(obj_md_{max_pair_idx},plain,[{_abc_md[0,0]},{_abc_md[0,1]},{_abc_md[0,2]}],'a',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_md_{max_pair_idx},plain,[{_abc_md[1,0]},{_abc_md[1,1]},{_abc_md[1,2]}],'b',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])
cyl_text(obj_md_{max_pair_idx},plain,[{_abc_md[2,0]},{_abc_md[2,1]},{_abc_md[2,2]}],'c',0.20,axes=[[1,0,0],[0,1,0],[0,0,1]])

cmd.load_cgo(obj_md_{max_pair_idx},  'abc_md_{max_pair_idx}')

cmd.cmd.read_pdbstr(
    \"\"\"{pdbstring}\"\"\",
    \"strc_{max_pair_idx}\"
)
cmd.group(\"pair-{max_pair_idx}\", \"abc_md_{max_pair_idx}\")
cmd.group(\"pair-{max_pair_idx}\", \"strc_{max_pair_idx}\")
"""
                    for pair_idx, pair in enumerate(max_pairs):
                        name = f"d_{max_pair_idx}_{pair_idx}-{max_dists[pair_idx]:4.2f}"
                        pymolstr += f"""
cmd.distance(
    name = \"{name}\",
    selection1 = \"rank {pair[0]} and strc_{max_pair_idx}\",
    selection2 = \"rank {pair[1]} and strc_{max_pair_idx}\",
    cutoff = 999.,
    mode = 0
)
cmd.group(\"pair-{max_pair_idx}\", \"{name}\")
"""
                with open(f"{crystal_name}-{forcefield_name}-hbond.py", "w") as fopen:
                    fopen.write(pymolstr)

            else:
                avg = "--"
                std = "--"
                dev = "--"
                max_avg = "--"
                max_std = "--"
                max_dev = "--"

            workbook_wrap.add_data(
                avg,
                std,
                dev,
                "<Δ[d(D-H•••A)]>",
                forcefield_name, 
                crystal_name
                )
            workbook_wrap.add_data(
                max_avg,
                max_std,
                max_dev,
                "Max <Δ[d(D-H•••A)]>",
                forcefield_name, 
                crystal_name
                )


            ### ... Deviations Angles ...

            if len(hbond_angles) > 0:
                diffs = np.abs(ref_hbond_angles-hbond_angles)
                shifted = np.where(
                    (360.-diffs) < diffs
                    )
                diffs[shifted] = 360.-diffs[shifted]
                avg   = np.mean(diffs)
                std   = np.std(diffs)
                dev   = "--"
                max_avg = 0.
                max_std = 0.
                max_dev = "--"
                for unique_rank in np.unique(hbond_pair_rank_list):
                    valids   = np.where(unique_rank == hbond_pair_rank_list)[0]
                    _max_avg = np.mean(diffs[:,valids])
                    if _max_avg > max_avg:
                        max_avg = _max_avg
                        max_std = np.std(diffs[:,valids])
            else:
                avg = "--"
                std = "--"
                dev = "--"
                max_avg = "--"
                max_std = "--"
                max_dev = "--"
            workbook_wrap.add_data(
                avg,
                std,
                dev,
                "<Δ[∠(D-H•••A)]>",
                forcefield_name, 
                crystal_name
                )
            workbook_wrap.add_data(
                max_avg,
                max_std,
                max_dev,
                "Max <Δ[∠(D-H•••A)]>",
                forcefield_name, 
                crystal_name
                )


            ### ======================= ###
            ### Write box vector length ###
            ### ======================= ###

            avg   = np.mean(a_len)
            std   = np.std(a_len)
            workbook_wrap.add_data(
                avg,
                std,
                "--",
                "a", 
                forcefield_name, 
                crystal_name
                )

            avg  = np.mean(b_len)
            std  = np.std(b_len)
            workbook_wrap.add_data(
                avg,
                std,
                "--",
                "b", 
                forcefield_name, 
                crystal_name
                )

            avg  = np.mean(c_len)
            std  = np.std(c_len)
            workbook_wrap.add_data(
                avg,
                std,
                "--",
                "c", 
                forcefield_name, 
                crystal_name
                )
            
            ### Write box vector angles ###
            ### ======================= ###
            avg  = np.mean(alpha)
            std  = np.std(alpha)
            workbook_wrap.add_data(
                avg,
                std,
                "--",
                "alpha", 
                forcefield_name, 
                crystal_name
                )

            avg  = np.mean(beta)
            std  = np.std(beta)
            workbook_wrap.add_data(
                avg,
                std,
                "--",
                "beta", 
                forcefield_name, 
                crystal_name
                )

            avg  = np.mean(gamma)
            std  = np.std(gamma)
            workbook_wrap.add_data(
                avg,
                std,
                "--",
                "gamma", 
                forcefield_name, 
                crystal_name
                )
                
            ### Write combined cell error
            ### =========================
            delta_a2 = (a_len - ref_a_len)**2
            delta_b2 = (b_len - ref_b_len)**2
            delta_c2 = (c_len - ref_c_len)**2

            ### Alpha
            delta_alpha = np.abs(ref_alpha-alpha)
            shifted = np.where(
                (360.-delta_alpha) < delta_alpha
                )
            delta_alpha[shifted] = 360.-delta_alpha[shifted]

            ### Beta
            delta_beta = np.abs(ref_beta-beta)
            shifted = np.where(
                (360.-delta_beta) < delta_beta
                )
            delta_beta[shifted] = 360.-delta_beta[shifted]

            ### Gamma
            delta_gamma = np.abs(ref_gamma-gamma)
            shifted = np.where(
                (360.-delta_gamma) < delta_gamma
                )
            delta_gamma[shifted] = 360.-delta_gamma[shifted]

            delta_alpha2 = delta_alpha**2
            delta_beta2  = delta_beta**2
            delta_gamma2 = delta_gamma**2

            cell_error  = delta_a2
            cell_error += delta_b2
            cell_error += delta_c2
            cell_error += _ANG_FACTOR * delta_alpha2
            cell_error += _ANG_FACTOR * delta_beta2
            cell_error += _ANG_FACTOR * delta_gamma2
            cell_error  = np.sqrt(cell_error)

            avg = np.mean(cell_error)
            std = np.std(cell_error)
            workbook_wrap.add_data(
                avg,
                std,
                "--",
                "cell error", 
                forcefield_name, 
                crystal_name
                )

            ### Write sublimation enthalpy ###
            ### ========================== ###
            ene_xtal = np.array(ene_xtal)
            ene_gas  = np.array(ene_gas)
            if ene_xtal.size > 1 and ene_gas.size > 1:
                ene_xtal_avg     = ene_xtal / float(N_molecules)
                ### Lattice energy. This is by definition a negative quantity.
                ene_lattice      = np.mean(ene_xtal_avg) - np.mean(ene_gas)
                RT               = _GASCONST_KCAL * input_dict[crystal_name]["experiment"]["temperature"]
                sublimation_avg  = - ene_lattice - 2. * RT
                sublimation_std  = np.var(ene_xtal_avg) + np.var(ene_gas)
                sublimation_std  = np.sqrt(sublimation_std)
                workbook_wrap.add_data(
                    sublimation_avg,
                    sublimation_std,
                    "--",
                    "Sublimation Energy", 
                    forcefield_name, 
                    crystal_name
                    )
            else:
                workbook_wrap.add_data(
                    "--",
                    "--",
                    "--",
                    "Sublimation Energy", 
                    forcefield_name, 
                    crystal_name
                    )
            if ene_xtal.size:
                ene_xtal_avg = ene_xtal / float(N_molecules)
                ene_avg = np.mean(ene_xtal_avg)
                ene_std = np.std(ene_xtal_avg)
                workbook_wrap.add_data(
                    ene_avg,
                    ene_std,
                    "--",
                    "Energy Lattice", 
                    forcefield_name, 
                    crystal_name
                    )
            else:
                workbook_wrap.add_data(
                    "--",
                    "--",
                    "--",
                    "Energy Lattice", 
                    forcefield_name, 
                    crystal_name
                    )
            if ene_gas.size:
                ene_avg = np.mean(ene_gas)
                ene_std = np.std(ene_gas)
                workbook_wrap.add_data(
                    ene_avg,
                    ene_std,
                    "--",
                    "Energy Gas", 
                    forcefield_name, 
                    crystal_name
                    )
            else:
                workbook_wrap.add_data(
                    "--",
                    "--",
                    "--",
                    "Energy Gas", 
                    forcefield_name, 
                    crystal_name
                    )

            ### Write density ###
            ### ============= ###
            avg = np.mean(density_xtal)
            std = np.std(density_xtal)
            workbook_wrap.add_data(
                avg,
                std,
                "--",
                "Density", 
                forcefield_name, 
                crystal_name
                )

        ### Parse in experimental data ###
        ### ========================== ###
        workbook_wrap.add_data(
            ref_a_len,
            "--",
            "--",
            "a", 
            "experiment", 
            crystal_name
            )

        workbook_wrap.add_data(
            ref_b_len,
            "--",
            "--",
            "b", 
            "experiment", 
            crystal_name
            )

        workbook_wrap.add_data(
            ref_c_len,
            "--",
            "--",
            "c", 
            "experiment", 
            crystal_name
            )

        workbook_wrap.add_data(
            ref_alpha,
            "--",
            "--",
            "alpha", 
            "experiment", 
            crystal_name
            )

        workbook_wrap.add_data(
            ref_beta,
            "--",
            "--",
            "beta", 
            "experiment", 
            crystal_name
            )

        workbook_wrap.add_data(
            ref_gamma,
            "--",
            "--",
            "gamma", 
            "experiment", 
            crystal_name
            )

        if "sublimation-enthalpy" in input_dict[crystal_name]["experiment"]:
            workbook_wrap.add_data(
                input_dict[crystal_name]["experiment"]["sublimation-enthalpy"],
                input_dict[crystal_name]["experiment"]["sublimation-enthalpy-std"],
                "--",
                "Sublimation Energy",
                "experiment", 
                crystal_name
                )
        else:
            workbook_wrap.add_data(
                "--",
                "--",
                "--",
                "Sublimation Energy",
                "experiment", 
                crystal_name
                )

        workbook_wrap.add_data(
            ref_density,
            "--",
            "--",
            "Density", 
            "experiment", 
            crystal_name
            )

        ### Write distance diff data ###
        ### ======================== ###
        avg   = np.mean(ref_distances)
        std   = "--"
        dev   = "--"
        workbook_wrap.add_data(
            avg,
            std,
            dev,
            "<(d < 4Å)>",
            "experiment", 
            crystal_name
            )

        ### ================ ###
        ### Write hbond data ###
        ### ================ ###

        if len(ref_hbond_distances) > 0:
            avg   = np.mean(ref_hbond_distances)
            std   = "--"
            dev   = "--"
            max_avg = 0.
            max_std = "--"
            max_dev = "--"
            for unique_rank in np.unique(hbond_pair_rank_list):
                valids   = np.where(unique_rank == hbond_pair_rank_list)[0]
                _max_avg = np.mean(ref_hbond_distances[:,valids])
                if _max_avg > max_avg:
                    max_avg = _max_avg
        else:
            avg = "--"
            std = "--"
            dev = "--"
            max_avg = "--"
            max_std = "--"
            max_dev = "--"
        workbook_wrap.add_data(
            avg,
            std,
            dev,
            "<[d(D-H•••A)]>",
            "experiment",
            crystal_name
            )

        if len(ref_hbond_angles) > 0:
            avg   = np.mean(ref_hbond_angles)
            std   = "--"
            dev   = "--"
            max_avg = 0.
            max_std = "--"
            max_dev = "--"
            for unique_rank in np.unique(hbond_pair_rank_list):
                valids   = np.where(unique_rank == hbond_pair_rank_list)[0]
                _max_avg = np.mean(ref_hbond_angles[:,valids])
                if _max_avg > max_avg:
                    max_avg = _max_avg
        else:
            avg = "--"
            std = "--"
            dev = "--"
            max_avg = "--"
            max_std = "--"
            max_dev = "--"
        workbook_wrap.add_data(
            avg,
            std,
            dev,
            "<[∠(D-H•••A)]>", 
            "experiment",
            crystal_name
            )

    workbook_wrap.close()


def entry_point():

    main()

if __name__ == "__main__":

    entry_point()