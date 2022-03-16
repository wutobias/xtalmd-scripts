#!/usr/bin/env python3

import xlsxwriter
import numpy as np
import argparse
import gemmi
import yaml
import glob
import os
import mdtraj as md
from collections import OrderedDict
from rdkit import Chem

from . import analysis_engine

_KJ_2_KCAL = 1./4.184
_NM_2_ANG  = 10.
_RAD_2_DEG = 180./np.pi
_GASCONST_KCAL = 8.31446261815324 * _KJ_2_KCAL / 1000.


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

        self.workbook  = xlsxwriter.Workbook(output)

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

        self.worksheet_dict   = OrderedDict()
        self.force_field_dict = OrderedDict()

        self.labels_dict_row = OrderedDict()
        self.labels_dict_row["Sublimation Energy"]     = 3
        self.labels_dict_row["a"]                      = 4
        self.labels_dict_row["b"]                      = 5
        self.labels_dict_row["c"]                      = 6
        self.labels_dict_row["α"]                      = 7
        self.labels_dict_row["β"]                      = 8
        self.labels_dict_row["γ"]                      = 9
        self.labels_dict_row["alpha"]                  = 7
        self.labels_dict_row["beta"]                   = 8
        self.labels_dict_row["gamma"]                  = 9
        self.labels_dict_row["<[Δ(d < 4Å)]>"]          = 10
        self.labels_dict_row["Max <[Δ(d < 4Å)]>"]      = 11
        self.labels_dict_row["H-bond geometry"]        = 12
        self.labels_dict_row["<[d(X-H•••O)]>"]         = 13
        self.labels_dict_row["<[d(X-H•••O=C)]>"]       = 14
        self.labels_dict_row["<[d(X-H•••N)]>"]         = 15
        self.labels_dict_row["<[d(X-H•••N=C)]>"]       = 16
        self.labels_dict_row["Max <[d(X-H•••O)]>"]     = 17
        self.labels_dict_row["Max <[d(X-H•••O=C)]>"]   = 18
        self.labels_dict_row["Max <[d(X-H•••N)]>"]     = 19
        self.labels_dict_row["Max <[d(X-H•••N=C)]>"]   = 20
        self.labels_dict_row["<[∠(X-H•••O)]>"]         = 21
        self.labels_dict_row["<[∠(X-H•••O=C)]>"]       = 22
        self.labels_dict_row["<[∠(X-H•••N)]>"]         = 23
        self.labels_dict_row["<[∠(X-H•••N=C)]>"]       = 24
        self.labels_dict_row["Max <[∠(X-H•••O)]>"]     = 25
        self.labels_dict_row["Max <[∠(X-H•••O=C)]>"]   = 26
        self.labels_dict_row["Max <[∠(X-H•••N)]>"]     = 27
        self.labels_dict_row["Max <[∠(X-H•••N=C)]>"]   = 28
        self.labels_dict_row["Translation/Rotation"]   = 29
        self.labels_dict_row["<[Δ(dCOM)]>"]            = 30
        self.labels_dict_row["Max <[Δ(dCOM)]>"]        = 31
        self.labels_dict_row["<{∠PA1-PA1(expt)}>"]     = 32
        self.labels_dict_row["<{∠PA2-PA2(expt)}>"]     = 33
        self.labels_dict_row["<{∠PA3-PA3(expt)}>"]     = 34
        self.labels_dict_row["Max <{∠PA1-PA1(expt)}>"] = 35
        self.labels_dict_row["Max <{∠PA2-PA2(expt)}>"] = 36
        self.labels_dict_row["Max <{∠PA3-PA3(expt)}>"] = 37
        self.labels_dict_row["Density"]                = 38

        self.sub_titles_row = [
            "H-bond geometry",
            "Translation/Rotation",
        ]

    
    def add_xtal(self, crystal_name):

        """
        Add a new crystal with name `crystal_name` to the workbook. This will come as a new
        worksheet in the xlsx file.
        """

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

        self.force_field_dict[crystal_name].append(forcefield_name)
        N_force_fields = len(self.force_field_dict[crystal_name])
        force_field_idx = N_force_fields - 1
        expt_col = 1 + N_force_fields * 2
        worksheet = self.worksheet_dict[crystal_name]

        worksheet.merge_range(
            first_row=1,
            first_col=1 + force_field_idx * 2, 
            last_row=1, 
            last_col=2 + force_field_idx * 2,
            data=forcefield_name,
            cell_format=self.header_format_1
        )

        worksheet.write(
            2,
            1 + force_field_idx * 2,
            "< >",
            self.header_format_2
        )

        worksheet.write(
            2,
            2 + force_field_idx * 2,
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
            "% Dev",
            self.header_format_2
        )

        worksheet.merge_range(
            first_row=0,
            first_col=1, 
            last_row=0, 
            last_col=2 * N_force_fields,
            data="Force Field",
            cell_format=self.header_format_1
        )

        worksheet.merge_range(
            first_row=1,
            first_col=expt_col, 
            last_row=1, 
            last_col=1 + expt_col,
            data="Expt.",
            cell_format=self.header_format_1
        )


    def add_data(
        self,
        data_value, 
        data_std, 
        data_name, 
        forcefield_name, 
        crystal_name):

        """
        Add data of value `data_value` and name `data_name` to worksheet
        of crystal `crystal_name` in column of force field `forcefield_name`.
        """

        if forcefield_name == "experiment":
            N_force_fields  = len(self.force_field_dict[crystal_name])
            force_field_idx = N_force_fields
        else:
            force_field_idx = self.force_field_dict[crystal_name].index(forcefield_name)
        
        worksheet = self.worksheet_dict[crystal_name]
        if isinstance(data_value, str):
            worksheet.write(
                self.labels_dict_row[data_name],
                1 + force_field_idx * 2,
                data_value,
                self.data_format_1
            )
        elif isinstance(data_value, type(None)):
            worksheet.write(
                self.labels_dict_row[data_name],
                1 + force_field_idx * 2,
                "--",
                self.data_format_1
            )
        elif np.isnan(data_value) or np.isinf(data_value):
            worksheet.write(
                self.labels_dict_row[data_name],
                1 + force_field_idx * 2,
                "--",
                self.data_format_1
            )
        else:
            worksheet.write(
                self.labels_dict_row[data_name],
                1 + force_field_idx * 2,
                data_value,
                self.data_format_1
            )

        if isinstance(data_std, str):
            worksheet.write(
                self.labels_dict_row[data_name],
                2 + force_field_idx * 2,
                data_std,
                self.data_format_1
            )
        elif isinstance(data_std, type(None)):
            worksheet.write(
                self.labels_dict_row[data_name],
                2 + force_field_idx * 2,
                "--",
                self.data_format_1
            )
        elif np.isnan(data_std) or np.isinf(data_std):
            worksheet.write(
                self.labels_dict_row[data_name],
                2 + force_field_idx * 2,
                "--",
                self.data_format_1
            )
        else:
            worksheet.write(
                self.labels_dict_row[data_name],
                2 + force_field_idx * 2,
                data_std,
                self.data_format_1
            )


    def close(self):

        """
        Close the workbook.
        """

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
        workbook_wrap.add_xtal(crystal_name)

        ref_strc   = md.load(input_dict[crystal_name]["experiment"]["supercell-pdb"])
        with open(input_dict[crystal_name]["experiment"]["supercell-rdkit"], "r") as fopen:
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
        ref_distances = analysis_engine.compute_pairwise_distances(ref_strc, dist_pair_list)

        acc_O_single_list,\
        acc_O_double_list,\
        acc_N_single_list,\
        acc_N_double_list = analysis_engine.get_hbond_indices(
            ref_strc,
            rdmol
            )

        if acc_O_single_list.size > 0:
            ref_hbond_O_single_diffs = analysis_engine.compute_pairwise_distances(
                ref_strc, 
                acc_O_single_list[:,[0,2]]
                )
            ref_hbond_O_single_diffs *= _NM_2_ANG

            ref_hbond_O_single_angles = analysis_engine.compute_tuplewise_angles(
                ref_strc, 
                acc_O_single_list
                )
            ref_hbond_O_single_angles *= _RAD_2_DEG

            acc_O_single_pair_rank_list = analysis_engine.build_pair_ranks(
                topology = ref_strc.topology,
                pair_list = acc_O_single_list[:,[0,2]],
                residue_classes_list = residue_classes_list
            )

        if acc_O_double_list.size > 0:
            ref_hbond_O_double_diffs = analysis_engine.compute_pairwise_distances(
                ref_strc, 
                acc_O_double_list[:,[0,2]]
                )
            ref_hbond_O_double_diffs *= _NM_2_ANG

            ref_hbond_O_double_angles = analysis_engine.compute_tuplewise_angles(
                ref_strc, 
                acc_O_double_list
                )
            ref_hbond_O_double_angles *= _RAD_2_DEG

            acc_O_double_pair_rank_list = analysis_engine.build_pair_ranks(
                topology = ref_strc.topology,
                pair_list = acc_O_double_list[:,[0,2]],
                residue_classes_list = residue_classes_list
            )

        if acc_N_single_list.size > 0:
            ref_hbond_N_single_diffs = analysis_engine.compute_pairwise_distances(
                ref_strc, 
                acc_N_single_list[:,[0,2]]
                )
            ref_hbond_N_single_diffs *= _NM_2_ANG

            ref_hbond_N_single_angles = analysis_engine.compute_tuplewise_angles(
                ref_strc, 
                acc_N_single_list
                )
            ref_hbond_N_single_angles *= _RAD_2_DEG

            acc_N_single_pair_rank_list = analysis_engine.build_pair_ranks(
                topology = ref_strc.topology,
                pair_list = acc_N_single_list[:,[0,2]],
                residue_classes_list = residue_classes_list
            )

        if acc_N_double_list.size > 0:
            ref_hbond_N_double_diffs = analysis_engine.compute_pairwise_distances(
                ref_strc, 
                acc_N_double_list[:,[0,2]]
                )
            ref_hbond_N_double_diffs *= _NM_2_ANG

            ref_hbond_N_double_angles = analysis_engine.compute_tuplewise_angles(
                ref_strc, 
                acc_N_double_list
                )
            ref_hbond_N_double_angles *= _RAD_2_DEG

            acc_N_double_pair_rank_list = analysis_engine.build_pair_ranks(
                topology = ref_strc.topology,
                pair_list = acc_N_double_list[:,[0,2]],
                residue_classes_list = residue_classes_list
            )

        ### Carry out analysis for each forcefield
        ### ======================================
        for forcefield_name in input_dict[crystal_name]:
            if forcefield_name.lower() == "experiment":
                continue
            print(f"Processing {crystal_name} / {forcefield_name}")

            workbook_wrap.add_forcefield(forcefield_name, crystal_name)

            xtal_topology   = input_dict[crystal_name][forcefield_name]["xtal-topology"]
            xtal_trajectory = input_dict[crystal_name][forcefield_name]["xtal-trajectory"]
            xtal_output     = input_dict[crystal_name][forcefield_name]["xtal-output"]

            gas_output = input_dict[crystal_name][forcefield_name]["gas-output"]

            a_len = list()
            b_len = list()
            c_len = list()

            alpha = list()
            beta  = list()
            gamma = list()

            ene_xtal = list()
            ene_gas  = list()

            density_xtal = list()

            distance_diffs = list()
            com_diffs      = list()
            pc1_diffs      = list()
            pc2_diffs      = list()
            pc3_diffs      = list()

            hbond_O_single_diffs = list()
            hbond_O_double_diffs = list()
            hbond_N_single_diffs = list()
            hbond_N_double_diffs = list()

            hbond_O_single_angles = list()
            hbond_O_double_angles = list()
            hbond_N_single_angles = list()
            hbond_N_double_angles = list()

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
                ene_gas.extend(_potential)

            for output_traj in glob.glob(xtal_trajectory):

                basename, _ = os.path.splitext(output_traj)
                if basename in bad_xtal_list:
                    continue

                query_traj = md.load(
                    output_traj,
                    top=ref_strc.topology
                    )

                _unitcell_angles = query_traj.unitcell_angles
                alpha.extend(_unitcell_angles[:,0].tolist())
                beta.extend(_unitcell_angles[:,1].tolist())
                gamma.extend(_unitcell_angles[:,2].tolist())

                ### Multiply by 10 to get Ang
                _unitcell_lengths = query_traj.unitcell_lengths * 10.
                ### Correct for number of unitcells along each direction
                _unitcell_lengths[:,0] /= N_unitcells_a
                _unitcell_lengths[:,1] /= N_unitcells_b
                _unitcell_lengths[:,2] /= N_unitcells_c
                a_len.extend(_unitcell_lengths[:,0].tolist())
                b_len.extend(_unitcell_lengths[:,1].tolist())
                c_len.extend(_unitcell_lengths[:,2].tolist())

                ### Devide by 1000 to get g/cm^3
                _density = md.density(query_traj) / 1000.
                density_xtal.extend(_density)

                _com_diffs      = analysis_engine.compute_com_diff_per_residue(
                    query_traj, 
                    ref_strc,
                    rdmol,
                    residue_classes_list,
                    )
                _pc_diffs       = analysis_engine.compute_pc_diff_per_residue(
                    query_traj, 
                    ref_strc,
                    rdmol
                    )
                _distance_diffs = ref_distances - analysis_engine.compute_pairwise_distances(
                    query_traj, 
                    dist_pair_list
                    )
                _distance_diffs = np.abs(_distance_diffs) * _NM_2_ANG
                _com_diffs      = np.abs(_com_diffs) * _NM_2_ANG

                if acc_O_single_list.size > 0:
                    _hbond_O_single_diffs = analysis_engine.compute_pairwise_distances(
                        query_traj, 
                        acc_O_single_list[:,[0,2]]
                        )
                    _hbond_O_single_diffs *= _NM_2_ANG
                    _hbond_O_single_diffs  = ref_hbond_O_single_diffs - _hbond_O_single_diffs

                    _hbond_O_single_angles = analysis_engine.compute_tuplewise_angles(
                        query_traj, 
                        acc_O_single_list,
                        )
                    _hbond_O_single_angles *= _RAD_2_DEG
                    _hbond_O_single_angles  = ref_hbond_O_single_angles - _hbond_O_single_angles

                if acc_O_double_list.size > 0:
                    _hbond_O_double_diffs = analysis_engine.compute_pairwise_distances(
                        query_traj, 
                        acc_O_double_list[:,[0,2]]
                        )
                    _hbond_O_double_diffs *= _NM_2_ANG
                    _hbond_O_double_diffs  = ref_hbond_O_double_diffs - _hbond_O_double_diffs
                    
                    _hbond_O_double_angles = analysis_engine.compute_tuplewise_angles(
                        query_traj, 
                        acc_O_double_list,
                        )
                    _hbond_O_double_angles *= _RAD_2_DEG
                    _hbond_O_double_angles  = ref_hbond_O_double_angles - _hbond_O_double_angles

                if acc_N_single_list.size > 0:
                    _hbond_N_single_diffs = analysis_engine.compute_pairwise_distances(
                        query_traj, 
                        acc_N_single_list[:,[0,2]]
                        )
                    _hbond_N_single_diffs *= _NM_2_ANG
                    _hbond_N_single_diffs  = ref_hbond_N_single_diffs - _hbond_N_single_diffs
                    
                    _hbond_N_single_angles = analysis_engine.compute_tuplewise_angles(
                        query_traj, 
                        acc_N_single_list,
                        )
                    _hbond_N_single_angles *= _RAD_2_DEG
                    _hbond_N_single_angles  = ref_hbond_N_single_angles - _hbond_N_single_angles

                if acc_N_double_list.size > 0:
                    _hbond_N_double_diffs = analysis_engine.compute_pairwise_distances(
                        query_traj, 
                        acc_N_double_list[:,[0,2]]
                        )
                    _hbond_N_double_diffs *= _NM_2_ANG
                    _hbond_N_double_diffs  = ref_hbond_N_double_diffs - _hbond_N_double_diffs
                    
                    _hbond_N_double_angles = analysis_engine.compute_tuplewise_angles(
                        query_traj, 
                        acc_N_double_list,
                        )
                    _hbond_N_double_angles *= _RAD_2_DEG
                    _hbond_N_double_angles  = ref_hbond_N_double_angles - _hbond_N_double_angles

                ### This means, we only do this the first iteration
                if len(com_diffs) == 0:
                    com_diffs      = _com_diffs
                    distance_diffs = _distance_diffs
                    pc1_diffs = _pc_diffs[:,:,0]
                    pc2_diffs = _pc_diffs[:,:,1]
                    pc3_diffs = _pc_diffs[:,:,2]
                    if acc_O_single_list.size > 0:
                        hbond_O_single_diffs  = _hbond_O_single_diffs
                        hbond_O_single_angles = _hbond_O_single_angles
                    if acc_O_double_list.size > 0:
                        hbond_O_double_diffs  = _hbond_O_double_diffs
                        hbond_O_double_angles = _hbond_O_double_angles
                    if acc_N_single_list.size > 0:
                        hbond_N_single_diffs  = _hbond_N_single_diffs
                        hbond_N_single_angles = _hbond_N_single_angles
                    if acc_N_double_list.size > 0:
                        hbond_N_double_diffs  = _hbond_N_double_diffs
                        hbond_N_double_angles = _hbond_N_double_angles

                else:
                    com_diffs      = np.vstack((com_diffs, _com_diffs))
                    distance_diffs = np.vstack((distance_diffs, _distance_diffs))
                    pc1_diffs = np.vstack((pc1_diffs, _pc_diffs[:,:,0]))
                    pc2_diffs = np.vstack((pc2_diffs, _pc_diffs[:,:,1]))
                    pc3_diffs = np.vstack((pc3_diffs, _pc_diffs[:,:,2]))
                    if acc_O_single_list.size > 0:
                        hbond_O_single_diffs  = np.vstack((hbond_O_single_diffs,  _hbond_O_single_diffs))
                        hbond_O_single_angles = np.vstack((hbond_O_single_angles, _hbond_O_single_angles))
                    if acc_O_double_list.size > 0:
                        hbond_O_double_diffs  = np.vstack((hbond_O_double_diffs, _hbond_O_double_diffs))
                        hbond_O_double_angles = np.vstack((hbond_O_double_angles, _hbond_O_double_angles))
                    if acc_N_single_list.size > 0:
                        hbond_N_single_diffs  = np.vstack((hbond_N_single_diffs, _hbond_N_single_diffs))
                        hbond_N_single_angles = np.vstack((hbond_N_single_angles, _hbond_N_single_angles))
                    if acc_N_double_list.size > 0:
                        hbond_N_double_diffs  = np.vstack((hbond_N_double_diffs, _hbond_N_double_diffs))
                        hbond_N_double_angles = np.vstack((hbond_N_double_angles, _hbond_N_double_angles))

            ### Make sure everything is np.ndarray
            a_len = np.array(a_len)
            b_len = np.array(b_len)
            c_len = np.array(c_len)
            alpha = np.array(alpha)
            beta  = np.array(beta)
            gamma = np.array(gamma)
            ene_xtal = np.array(ene_xtal)
            ene_gas  = np.array(ene_gas)
            density_xtal = np.array(density_xtal)
            distance_diffs = np.array(distance_diffs)
            com_diffs      = np.array(com_diffs)
            pc1_diffs      = np.array(pc1_diffs)
            pc2_diffs      = np.array(pc2_diffs)
            pc3_diffs      = np.array(pc3_diffs)
            hbond_O_single_diffs = np.array(hbond_O_single_diffs)
            hbond_O_double_diffs = np.array(hbond_O_double_diffs)
            hbond_N_single_diffs = np.array(hbond_N_single_diffs)
            hbond_N_double_diffs = np.array(hbond_N_double_diffs)
            hbond_O_single_angles = np.array(hbond_O_single_angles)
            hbond_O_double_angles = np.array(hbond_O_double_angles)
            hbond_N_single_angles = np.array(hbond_N_single_angles)
            hbond_N_double_angles = np.array(hbond_N_double_angles)

            ### If we don't have any data:
            if distance_diffs.size < 2:
                continue

            ### Write distance diff data ###
            ### ======================== ###
            avg  = np.mean(distance_diffs)
            std  = np.std(distance_diffs)
            workbook_wrap.add_data(
                avg,
                std/avg*100.,
                "<[Δ(d < 4Å)]>", 
                forcefield_name, 
                crystal_name
                )

            max_avg = 0.
            max_std = 0.
            for unique_rank in np.unique(dist_pair_rank_list):
                valids   = np.where(unique_rank == dist_pair_rank_list)[0]
                _max_avg = np.mean(distance_diffs[:,valids])
                if _max_avg > max_avg:
                    max_avg = _max_avg
                    max_std = np.std(distance_diffs[:,valids])
            workbook_wrap.add_data(
                max_avg,
                max_std/max_avg*100.,
                "Max <[Δ(d < 4Å)]>", 
                forcefield_name, 
                crystal_name
                )

            ### Write com diff data ###
            ### =================== ###
            avg  = np.mean(com_diffs)
            std  = np.std(com_diffs)
            workbook_wrap.add_data(
                avg,
                std/avg*100.,
                "<[Δ(dCOM)]>",
                forcefield_name, 
                crystal_name
                )

            res_avg = np.mean(com_diffs, axis=0)
            max_idx = np.argmax(res_avg)
            max_avg = res_avg[max_idx]
            max_std = np.std(com_diffs[:,max_idx])
            workbook_wrap.add_data(
                max_avg,
                max_std/max_avg*100.,
                "Max <[Δ(dCOM)]>",
                forcefield_name, 
                crystal_name
                )

            ### Write pc diff data ###
            ### ================== ###
            avg  = np.mean(pc1_diffs)
            std  = np.std(pc1_diffs)
            workbook_wrap.add_data(
                avg,
                std/avg*100.,
                "<{∠PA1-PA1(expt)}>", 
                forcefield_name, 
                crystal_name
                )

            res_avg = np.mean(pc1_diffs, axis=0)
            max_idx = np.argmax(res_avg)
            max_avg = res_avg[max_idx]
            max_std = np.std(pc1_diffs[:,max_idx])
            workbook_wrap.add_data(
                max_avg,
                max_std/max_avg*100.,
                "Max <{∠PA1-PA1(expt)}>", 
                forcefield_name, 
                crystal_name
                )

            avg  = np.mean(pc2_diffs)
            std  = np.std(pc2_diffs)
            workbook_wrap.add_data(
                avg,
                std/avg*100.,
                "<{∠PA2-PA2(expt)}>", 
                forcefield_name, 
                crystal_name
                )

            res_avg = np.mean(pc2_diffs, axis=0)
            max_idx = np.argmax(res_avg)
            max_avg = res_avg[max_idx]
            max_std = np.std(pc2_diffs[:,max_idx])
            workbook_wrap.add_data(
                max_avg,
                max_std/max_avg*100.,
                "Max <{∠PA2-PA2(expt)}>", 
                forcefield_name, 
                crystal_name
                )

            avg  = np.mean(pc3_diffs)
            std  = np.std(pc3_diffs)
            workbook_wrap.add_data(
                avg,
                std/avg*100.,
                "<{∠PA3-PA3(expt)}>", 
                forcefield_name, 
                crystal_name
                )

            res_avg = np.mean(pc3_diffs, axis=0)
            max_idx = np.argmax(res_avg)
            max_avg = res_avg[max_idx]
            max_std = np.std(pc3_diffs[:,max_idx])
            workbook_wrap.add_data(
                max_avg,
                max_std/max_avg*100.,
                "Max <{∠PA3-PA3(expt)}>", 
                forcefield_name, 
                crystal_name
                )

            ### Write hbond data ###
            ### ================ ###

            ### X-H•••O
            ### -------
            if len(hbond_O_single_diffs) > 0:
                avg  = np.mean(hbond_O_single_diffs)
                std  = np.std(hbond_O_single_diffs)
                std  = std/avg*100.
                max_avg = 0.
                max_std = 0.
                for unique_rank in np.unique(acc_O_single_pair_rank_list):
                    valids   = np.where(unique_rank == acc_O_single_pair_rank_list)[0]
                    _max_avg = np.mean(hbond_O_single_diffs[:,valids])
                    if _max_avg > max_avg:
                        max_avg = _max_avg
                        max_std = np.std(hbond_O_single_diffs[:,valids])
            else:
                avg = "--"
                std = "--"
                max_avg = "--"
                max_std = "--"
            workbook_wrap.add_data(
                avg,
                std,
                "<[d(X-H•••O)]>", 
                forcefield_name, 
                crystal_name
                )
            workbook_wrap.add_data(
                max_avg,
                max_std,
                "Max <[d(X-H•••O)]>", 
                forcefield_name, 
                crystal_name
                )

            if len(hbond_O_single_angles) > 0:
                avg  = np.mean(hbond_O_single_angles)
                std  = np.std(hbond_O_single_angles)
                std  = std/avg*100.
                max_avg = -np.inf
                max_std = 0.
                min_avg = np.inf
                min_std = 0.
                for unique_rank in np.unique(acc_O_single_pair_rank_list):
                    valids   = np.where(unique_rank == acc_O_single_pair_rank_list)[0]
                    rank_avg = np.mean(hbond_O_single_angles[:,valids])
                    if rank_avg > max_avg:
                        max_avg = rank_avg
                        max_std = np.std(hbond_O_single_angles[:,valids])
                    if rank_avg < min_avg:
                        min_avg = rank_avg
                        min_std = np.std(hbond_O_single_angles[:,valids])
            else:
                avg = "--"
                std = "--"
                max_avg = "--"
                max_std = "--"
            workbook_wrap.add_data(
                avg,
                std,
                "<[∠(X-H•••O)]>", 
                forcefield_name, 
                crystal_name
                )
            workbook_wrap.add_data(
                max_avg,
                max_std,
                "Max <[∠(X-H•••O)]>", 
                forcefield_name, 
                crystal_name
                )

            ### X-H•••O=C
            ### ---------
            if len(hbond_O_double_diffs) > 0:
                avg  = np.mean(hbond_O_double_diffs)
                std  = np.std(hbond_O_double_diffs)
                std  = std/avg*100.
                max_avg = 0.
                max_std = 0.
                for unique_rank in np.unique(acc_O_double_pair_rank_list):
                    valids   = np.where(unique_rank == acc_O_double_pair_rank_list)[0]
                    _max_avg = np.mean(hbond_O_double_diffs[:,valids])
                    if _max_avg > max_avg:
                        max_avg = _max_avg
                        max_std = np.std(hbond_O_double_diffs[:,valids])
            else:
                avg = "--"
                std = "--"
                max_avg = "--"
                max_std = "--"
            workbook_wrap.add_data(
                avg,
                std,
                "<[d(X-H•••O=C)]>", 
                forcefield_name, 
                crystal_name
                )
            workbook_wrap.add_data(
                max_avg,
                max_std,
                "Max <[d(X-H•••O=C)]>", 
                forcefield_name, 
                crystal_name
                )

            if len(hbond_O_double_angles) > 0:
                avg  = np.mean(hbond_O_double_angles)
                std  = np.std(hbond_O_double_angles)
                std  = std/avg*100.
                max_avg = -np.inf
                max_std = 0.
                min_avg = np.inf
                min_std = 0.
                for unique_rank in np.unique(acc_O_double_pair_rank_list):
                    valids   = np.where(unique_rank == acc_O_double_pair_rank_list)[0]
                    rank_avg = np.mean(hbond_O_double_angles[:,valids])
                    if rank_avg > max_avg:
                        max_avg = rank_avg
                        max_std = np.std(hbond_O_double_angles[:,valids])
                    if rank_avg < min_avg:
                        min_avg = rank_avg
                        min_std = np.std(hbond_O_double_angles[:,valids])
            else:
                avg = "--"
                std = "--"
                max_avg = "--"
                max_std = "--"
            workbook_wrap.add_data(
                avg,
                std,
                "<[∠(X-H•••O=C)]>", 
                forcefield_name, 
                crystal_name
                )
            workbook_wrap.add_data(
                max_avg,
                max_std,
                "Max <[∠(X-H•••O=C)]>", 
                forcefield_name, 
                crystal_name
                )

            ### X-H•••N
            ### -------
            if len(hbond_N_single_diffs) > 0:
                avg  = np.mean(hbond_N_single_diffs)
                std  = np.std(hbond_N_single_diffs)
                std  = std/avg*100.
                max_avg = 0.
                max_std = 0.
                for unique_rank in np.unique(acc_N_single_pair_rank_list):
                    valids   = np.where(unique_rank == acc_N_single_pair_rank_list)[0]
                    _max_avg = np.mean(hbond_N_single_diffs[:,valids])
                    if _max_avg > max_avg:
                        max_avg = _max_avg
                        max_std = np.std(hbond_N_single_diffs[:,valids])
            else:
                avg = "--"
                std = "--"
                max_avg = "--"
                max_std = "--"
            workbook_wrap.add_data(
                avg,
                std,
                "<[d(X-H•••N)]>", 
                forcefield_name, 
                crystal_name
                )
            workbook_wrap.add_data(
                max_avg,
                max_std,
                "Max <[d(X-H•••N)]>", 
                forcefield_name, 
                crystal_name
                )

            if len(hbond_N_single_angles) > 0:
                avg  = np.mean(hbond_N_single_angles)
                std  = np.std(hbond_N_single_angles)
                std  = std/avg*100.
                max_avg = -np.inf
                max_std = 0.
                min_avg = np.inf
                min_std = 0.
                for unique_rank in np.unique(acc_N_single_pair_rank_list):
                    valids   = np.where(unique_rank == acc_N_single_pair_rank_list)[0]
                    rank_avg = np.mean(hbond_N_single_angles[:,valids])
                    if rank_avg > max_avg:
                        max_avg = rank_avg
                        max_std = np.std(hbond_N_single_angles[:,valids])
                    if rank_avg < min_avg:
                        min_avg = rank_avg
                        min_std = np.std(hbond_N_single_angles[:,valids])
            else:
                avg = "--"
                std = "--"
                max_avg = "--"
                max_std = "--"
            workbook_wrap.add_data(
                avg,
                std,
                "<[∠(X-H•••N)]>", 
                forcefield_name, 
                crystal_name
                )
            workbook_wrap.add_data(
                max_avg,
                max_std,
                "Max <[∠(X-H•••N)]>", 
                forcefield_name, 
                crystal_name
                )

            ### X-H•••N=C
            ### ---------
            if len(hbond_N_double_diffs) > 0:
                avg  = np.mean(hbond_N_double_diffs)
                std  = np.std(hbond_N_double_diffs)
                std  = std/avg*100.
                max_avg = 0.
                max_std = 0.
                for unique_rank in np.unique(acc_N_double_pair_rank_list):
                    valids   = np.where(unique_rank == acc_N_double_pair_rank_list)[0]
                    _max_avg = np.mean(hbond_N_double_diffs[:,valids])
                    if _max_avg > max_avg:
                        max_avg = _max_avg
                        max_std = np.std(hbond_N_double_diffs[:,valids])
            else:
                avg = "--"
                std = "--"
                max_avg = "--"
                max_std = "--"
            workbook_wrap.add_data(
                avg,
                std,
                "<[d(X-H•••N=C)]>", 
                forcefield_name, 
                crystal_name
                )
            workbook_wrap.add_data(
                max_avg,
                max_std,
                "Max <[d(X-H•••N=C)]>", 
                forcefield_name, 
                crystal_name
                )

            if len(hbond_N_double_angles) > 0:
                avg  = np.mean(hbond_N_double_angles)
                std  = np.std(hbond_N_double_angles)
                std  = std/avg*100.
                max_avg = -np.inf
                max_std = 0.
                min_avg = np.inf
                min_std = 0.
                for unique_rank in np.unique(acc_N_double_pair_rank_list):
                    valids   = np.where(unique_rank == acc_N_double_pair_rank_list)[0]
                    rank_avg = np.mean(hbond_N_double_angles[:,valids])
                    if rank_avg > max_avg:
                        max_avg = rank_avg
                        max_std = np.std(hbond_N_double_angles[:,valids])
                    if rank_avg < min_avg:
                        min_avg = rank_avg
                        min_std = np.std(hbond_N_double_angles[:,valids])
            else:
                avg = "--"
                std = "--"
                max_avg = "--"
                max_std = "--"
            workbook_wrap.add_data(
                avg,
                std,
                "<[∠(X-H•••N=C)]>", 
                forcefield_name, 
                crystal_name
                )
            workbook_wrap.add_data(
                max_avg,
                max_std,
                "Max <[∠(X-H•••N=C)]>", 
                forcefield_name, 
                crystal_name
                )

            ### Write box vector length ###
            ### ======================= ###
            avg   = np.mean(a_len)
            std   = np.std(a_len)
            workbook_wrap.add_data(
                avg,
                std/avg*100.,
                "a", 
                forcefield_name, 
                crystal_name
                )

            avg  = np.mean(b_len)
            std  = np.std(b_len)
            workbook_wrap.add_data(
                avg,
                std/avg*100.,
                "b", 
                forcefield_name, 
                crystal_name
                )

            avg  = np.mean(c_len)
            std  = np.std(c_len)
            workbook_wrap.add_data(
                avg,
                std/avg*100.,
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
                std/avg*100.,
                "alpha", 
                forcefield_name, 
                crystal_name
                )

            avg  = np.mean(beta)
            std  = np.std(beta)
            workbook_wrap.add_data(
                avg,
                std/avg*100.,
                "beta", 
                forcefield_name, 
                crystal_name
                )

            avg  = np.mean(gamma)
            std  = np.std(gamma)
            workbook_wrap.add_data(
                avg,
                std/avg*100.,
                "gamma", 
                forcefield_name, 
                crystal_name
                )

            ### Write sublimation enthalpy ###
            ### ========================== ###
            ene_xtal = np.array(ene_xtal) * _KJ_2_KCAL
            ene_gas  = np.array(ene_gas) * _KJ_2_KCAL
            if ene_xtal.size > 1 and ene_gas.size > 1:
                ene_xtal        /= float(N_molecules)
                sublimation_avg  = np.mean(ene_gas) - np.mean(ene_xtal)/float(N_molecules)
                sublimation_avg += (_GASCONST_KCAL * input_dict[crystal_name]["experiment"]["temperature"])
                sublimation_std  = np.var(ene_xtal) + np.var(ene_gas)
                sublimation_std  = np.sqrt(sublimation_std)
                workbook_wrap.add_data(
                    sublimation_avg,
                    sublimation_std/sublimation_avg * 100.,
                    "Sublimation Energy", 
                    forcefield_name, 
                    crystal_name
                    )
            else:
                workbook_wrap.add_data(
                    "--",
                    "--",
                    "Sublimation Energy", 
                    forcefield_name, 
                    crystal_name
                    )

            ### Write density ###
            ### ============= ###
            avg = np.mean(density_xtal)
            std = np.std(density_xtal)
            workbook_wrap.add_data(
                avg,
                std/avg*100.,
                "Density", 
                forcefield_name, 
                crystal_name
                )

        ### Parse in experimental data ###
        ### ========================== ###
        doc     = gemmi.cif.read(input_dict[crystal_name]["experiment"]["experiment-cif"])[0]
        strc    = gemmi.make_small_structure_from_block(doc)
        density = doc.find_value("_exptl_crystal_density_diffrn")
        if density != None:
            try:
                density = float(density)
            except:
                density = "--"
        elif "density" in input_dict[crystal_name]["experiment"]:
            density = input_dict[crystal_name]["experiment"]["density"]
        else:
            ### Then we don't have density
            density = "--"

        workbook_wrap.add_data(
            strc.cell.a,
            "--",
            "a", 
            "experiment", 
            crystal_name
            )

        workbook_wrap.add_data(
            strc.cell.b,
            "--",
            "b", 
            "experiment", 
            crystal_name
            )

        workbook_wrap.add_data(
            strc.cell.c,
            "--",
            "c", 
            "experiment", 
            crystal_name
            )

        workbook_wrap.add_data(
            strc.cell.alpha,
            "--",
            "alpha", 
            "experiment", 
            crystal_name
            )

        workbook_wrap.add_data(
            strc.cell.beta,
            "--",
            "beta", 
            "experiment", 
            crystal_name
            )

        workbook_wrap.add_data(
            strc.cell.gamma,
            "--",
            "gamma", 
            "experiment", 
            crystal_name
            )

        ### Write hbond data ###
        ### ================ ###

#        ### X-H•••O
#        ### -------
#        if acc_O_single_list.size > 0:
#            avg  = np.mean(ref_hbond_O_single_diffs)
#            std  = np.std(ref_hbond_O_single_diffs)
#            std  = std/avg*100.
#            max_avg = 0.
#            max_std = "--"
#            for unique_rank in np.unique(acc_O_single_pair_rank_list):
#                valids   = np.where(unique_rank == acc_O_single_pair_rank_list)[0]
#                _max_avg = np.mean(ref_hbond_O_single_diffs[:,valids])
#                if _max_avg > max_avg:
#                    max_avg = _max_avg
#        else:
#            avg = "--"
#            std = "--"
#            max_avg = "--"
#            max_std = "--"
#        workbook_wrap.add_data(
#            avg,
#            std,
#            "<[d(X-H•••O)]>", 
#            "experiment", 
#            crystal_name
#            )
#        workbook_wrap.add_data(
#            max_avg,
#            max_std,
#            "Max <[d(X-H•••O)]>", 
#            "experiment", 
#            crystal_name
#            )
#
#        if acc_O_single_list.size > 0:
#            avg  = np.mean(ref_hbond_O_single_angles)
#            std  = np.std(ref_hbond_O_single_angles)
#            std  = std/avg*100.
#            max_avg = -np.inf
#            max_std = 0.
#            min_avg = np.inf
#            min_std = 0.
#            for unique_rank in np.unique(acc_O_single_pair_rank_list):
#                valids   = np.where(unique_rank == acc_O_single_pair_rank_list)[0]
#                rank_avg = np.mean(ref_hbond_O_single_angles[:,valids])
#                if rank_avg > max_avg:
#                    max_avg = rank_avg
#                    max_std = np.std(ref_hbond_O_single_angles[:,valids])
#                if rank_avg < min_avg:
#                    min_avg = rank_avg
#                    min_std = np.std(ref_hbond_O_single_angles[:,valids])
#        else:
#            avg = "--"
#            std = "--"
#            max_avg = "--"
#            max_std = "--"
#        workbook_wrap.add_data(
#            avg,
#            std,
#            "<[∠(X-H•••O)]>", 
#            "experiment",
#            crystal_name
#            )
#        workbook_wrap.add_data(
#            max_avg,
#            max_std,
#            "Max <[∠(X-H•••O)]>", 
#            "experiment",
#            crystal_name
#            )
#
#        ### X-H•••O=C
#        ### ---------
#        if acc_O_double_list.size > 0:
#            avg  = np.mean(ref_hbond_O_double_diffs)
#            std  = np.std(ref_hbond_O_double_diffs)
#            std  = std/avg*100.
#            max_avg = 0.
#            max_std = "--"
#            for unique_rank in np.unique(acc_O_double_pair_rank_list):
#                valids   = np.where(unique_rank == acc_O_double_pair_rank_list)[0]
#                _max_avg = np.mean(ref_hbond_O_double_diffs[:,valids])
#                if _max_avg > max_avg:
#                    max_avg = _max_avg
#        else:
#            avg = "--"
#            std = "--"
#            max_avg = "--"
#            max_std = "--"
#        workbook_wrap.add_data(
#            avg,
#            std,
#            "<[d(X-H•••O=C)]>", 
#            "experiment", 
#            crystal_name
#            )
#        workbook_wrap.add_data(
#            max_avg,
#            max_std,
#            "Max <[d(X-H•••O=C)]>", 
#            "experiment", 
#            crystal_name
#            )
#
#        if acc_O_double_list.size > 0:
#            avg  = np.mean(ref_hbond_O_double_angles)
#            std  = np.std(ref_hbond_O_double_angles)
#            std  = std/avg*100.
#            max_avg = -np.inf
#            max_std = 0.
#            min_avg = np.inf
#            min_std = 0.
#            for unique_rank in np.unique(acc_O_double_pair_rank_list):
#                valids   = np.where(unique_rank == acc_O_double_pair_rank_list)[0]
#                rank_avg = np.mean(ref_hbond_O_double_angles[:,valids])
#                if rank_avg > max_avg:
#                    max_avg = rank_avg
#                    max_std = np.std(ref_hbond_O_double_angles[:,valids])
#                if rank_avg < min_avg:
#                    min_avg = rank_avg
#                    min_std = np.std(ref_hbond_O_double_angles[:,valids])
#        else:
#            avg = "--"
#            std = "--"
#            max_avg = "--"
#            max_std = "--"
#        workbook_wrap.add_data(
#            avg,
#            std,
#            "<[∠(X-H•••O=C)]>", 
#            "experiment",
#            crystal_name
#            )
#        workbook_wrap.add_data(
#            max_avg,
#            max_std,
#            "Max <[∠(X-H•••O=C)]>", 
#            "experiment",
#            crystal_name
#            )
#
#        ### X-H•••N
#        ### -------
#        if acc_N_single_list.size > 0:
#            avg  = np.mean(ref_hbond_N_single_diffs)
#            std  = np.std(ref_hbond_N_single_diffs)
#            std  = std/avg*100.
#            max_avg = 0.
#            max_std = "--"
#            for unique_rank in np.unique(acc_N_single_pair_rank_list):
#                valids   = np.where(unique_rank == acc_N_single_pair_rank_list)[0]
#                _max_avg = np.mean(ref_hbond_N_single_diffs[:,valids])
#                if _max_avg > max_avg:
#                    max_avg = _max_avg
#        else:
#            avg = "--"
#            std = "--"
#            max_avg = "--"
#            max_std = "--"
#        workbook_wrap.add_data(
#            avg,
#            std,
#            "<[d(X-H•••N)]>", 
#            "experiment", 
#            crystal_name
#            )
#        workbook_wrap.add_data(
#            max_avg,
#            max_std,
#            "Max <[d(X-H•••N)]>", 
#            "experiment", 
#            crystal_name
#            )
#
#        if acc_N_single_list.size > 0:
#            avg  = np.mean(ref_hbond_N_single_angles)
#            std  = np.std(ref_hbond_N_single_angles)
#            std  = std/avg*100.
#            max_avg = -np.inf
#            max_std = 0.
#            min_avg = np.inf
#            min_std = 0.
#            for unique_rank in np.unique(acc_N_single_pair_rank_list):
#                valids   = np.where(unique_rank == acc_N_single_pair_rank_list)[0]
#                rank_avg = np.mean(ref_hbond_N_single_angles[:,valids])
#                if rank_avg > max_avg:
#                    max_avg = rank_avg
#                    max_std = np.std(ref_hbond_N_single_angles[:,valids])
#                if rank_avg < min_avg:
#                    min_avg = rank_avg
#                    min_std = np.std(ref_hbond_N_single_angles[:,valids])
#        else:
#            avg = "--"
#            std = "--"
#            max_avg = "--"
#            max_std = "--"
#        workbook_wrap.add_data(
#            avg,
#            std,
#            "<[∠(X-H•••N)]>", 
#            "experiment",
#            crystal_name
#            )
#        workbook_wrap.add_data(
#            max_avg,
#            max_std,
#            "Max <[∠(X-H•••N)]>", 
#            "experiment",
#            crystal_name
#            )
#
#        ### X-H•••N=C
#        ### ---------
#        if acc_N_double_list.size > 0:
#            avg  = np.mean(ref_hbond_N_double_diffs)
#            std  = np.std(ref_hbond_N_double_diffs)
#            std  = std/avg*100.
#            max_avg = 0.
#            max_std = "--"
#            for unique_rank in np.unique(acc_N_double_pair_rank_list):
#                valids   = np.where(unique_rank == acc_N_double_pair_rank_list)[0]
#                _max_avg = np.mean(ref_hbond_N_double_diffs[:,valids])
#                if _max_avg > max_avg:
#                    max_avg = _max_avg
#        else:
#            avg = "--"
#            std = "--"
#            max_avg = "--"
#            max_std = "--"
#        workbook_wrap.add_data(
#            avg,
#            std,
#            "<[d(X-H•••N=C)]>", 
#            "experiment", 
#            crystal_name
#            )
#        workbook_wrap.add_data(
#            max_avg,
#            max_std,
#            "Max <[d(X-H•••N=C)]>", 
#            "experiment", 
#            crystal_name
#            )
#
#        if acc_N_double_list.size > 0:
#            avg  = np.mean(ref_hbond_N_double_angles)
#            std  = np.std(ref_hbond_N_double_angles)
#            std  = std/avg*100.
#            max_avg = -np.inf
#            max_std = 0.
#            min_avg = np.inf
#            min_std = 0.
#            for unique_rank in np.unique(acc_N_double_pair_rank_list):
#                valids   = np.where(unique_rank == acc_N_double_pair_rank_list)[0]
#                rank_avg = np.mean(ref_hbond_N_double_angles[:,valids])
#                if rank_avg > max_avg:
#                    max_avg = rank_avg
#                    max_std = np.std(ref_hbond_N_double_angles[:,valids])
#                if rank_avg < min_avg:
#                    min_avg = rank_avg
#                    min_std = np.std(ref_hbond_N_double_angles[:,valids])
#        else:
#            avg = "--"
#            std = "--"
#            max_avg = "--"
#            max_std = "--"
#        workbook_wrap.add_data(
#            avg,
#            std,
#            "<[∠(X-H•••N=C)]>", 
#            "experiment",
#            crystal_name
#            )
#        workbook_wrap.add_data(
#            max_avg,
#            max_std,
#            "Max <[∠(X-H•••N=C)]>", 
#            "experiment",
#            crystal_name
#            )

        ### Sublimation Enthalpy ###
        ### -------------------- ###

        if "sublimation-enthalpy" in input_dict[crystal_name]["experiment"]:
            workbook_wrap.add_data(
                input_dict[crystal_name]["experiment"]["sublimation-enthalpy"],
                input_dict[crystal_name]["experiment"]["sublimation-enthalpy-std"],
                "Sublimation Energy",
                "experiment", 
                crystal_name
                )

        workbook_wrap.add_data(
            density,
            "--",
            "Density", 
            "experiment", 
            crystal_name
            )

    workbook_wrap.close()


def entry_point():

    main()

if __name__ == "__main__":

    entry_point()