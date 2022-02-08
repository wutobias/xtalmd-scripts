#!/usr/bin/env python3

import xlsxwriter
import numpy as np
import argparse
import gemmi
import yaml
import glob
import mdtraj as md
from collections import OrderedDict
from rdkit import Chem

from . import analysis_engine

_KJ_2_KCAL = 1./4.184
_NM_2_ANG  = 10.
_GASCONST_KCAL = 8.31446261815324 * _KJ_2_KCAL

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
    for key in data_dict:
        data_dict[key] = np.array(data_dict[key])
    data_dict["N_rows"] = data_dict[key].size
    data_dict["N_columns"] = len(data_dict[key])

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
            }
        )

        self.worksheet_dict   = OrderedDict()
        self.force_field_dict = OrderedDict()

        self.labels_dict_row = OrderedDict()
        self.labels_dict_row["Sublimation Energy"]   = 3
        self.labels_dict_row["a"]                    = 4
        self.labels_dict_row["b"]                    = 5
        self.labels_dict_row["c"]                    = 6
        self.labels_dict_row["α"]                    = 7
        self.labels_dict_row["β"]                    = 8
        self.labels_dict_row["γ"]                    = 9
        self.labels_dict_row["alpha"]                = 7
        self.labels_dict_row["beta"]                 = 8
        self.labels_dict_row["gamma"]                = 9
        self.labels_dict_row["<[Δ(d < 4Å)]>"]        = 10
        self.labels_dict_row["Max Δ(d < 4Å)"]        = 11
        self.labels_dict_row["H-bond geometry"]      = 12
        self.labels_dict_row["X-H•••O"]              = 13
        self.labels_dict_row["X-H•••O=C"]            = 14
        self.labels_dict_row["X-H•••N"]              = 15
        self.labels_dict_row["X-H•••N=C"]            = 16
        self.labels_dict_row["Translation/Rotation"] = 17
        self.labels_dict_row["<[D(dCOM)]>"]          = 18
        self.labels_dict_row["Max D(dCOM)"]          = 19
        self.labels_dict_row["<{∠PA1-PA1(expt)}>"]   = 20
        self.labels_dict_row["<{∠PA2-PA2(expt)}>"]   = 21
        self.labels_dict_row["<{∠PA3-PA3(expt)}>"]   = 22
        self.labels_dict_row["Max ∠PA1-PA1(expt)"]   = 23
        self.labels_dict_row["Max ∠PA2-PA2(expt)"]   = 24
        self.labels_dict_row["Max ∠PA3-PA3(expt)"]   = 25
        self.labels_dict_row["Density"]              = 26

    
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
        worksheet.write(
            self.labels_dict_row[data_name],
            1 + force_field_idx * 2,
            data_value,
            self.data_format_1
        )
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

    args = parse_arguments()
    with open(args.input, "r") as fopen:
        input_dict = yaml.safe_load(fopen)

    workbook_wrap  = WorkbookWrapper(args.output)

    ### Loop over each xtal and put in new workbook
    ### ===========================================
    for crystal_name in input_dict:
        workbook_wrap.add_xtal(crystal_name)

        ref_strc = md.load(input_dict[crystal_name]["experiment"]["supercell-pdb"])
        with open(input_dict[crystal_name]["experiment"]["supercell-rdkit"], "r") as fopen:
            rdmol = Chem.JSONToMols(fopen.read())[0]
        dist_pair_list = analysis_engine.build_pair_list(
            traj = ref_strc,
            distance_cutoff = 0.4, 
            bond_cutoff = 4, 
            exclude_hydrogen=False
            )
        acc_O_single_pair_list,\
        acc_O_double_pair_list,\
        acc_N_single_pair_list,\
        acc_N_double_pair_list = analysis_engine.get_hbond_indices(
            ref_strc,
            rdmol
            )

        ref_distances = analysis_engine.compute_pairwise_distances(ref_strc, dist_pair_list)

        if acc_O_single_pair_list.size > 0:
            ref_hbond_O_single_diffs = analysis_engine.compute_pairwise_distances(
                ref_strc, 
                acc_O_single_pair_list
                )

        if acc_O_double_pair_list.size > 0:
            ref_hbond_O_double_diffs = analysis_engine.compute_pairwise_distances(
                ref_strc, 
                acc_O_double_pair_list
                )

        if acc_N_single_pair_list.size > 0:
            ref_hbond_N_single_diffs = analysis_engine.compute_pairwise_distances(
                ref_strc, 
                acc_N_single_pair_list
                )

        if acc_N_double_pair_list.size > 0:
            ref_hbond_N_double_diffs = analysis_engine.compute_pairwise_distances(
                ref_strc, 
                acc_N_double_pair_list
                )

        ### Carry out analysis for each forcefield
        ### ======================================
        for forcefield_name in input_dict[crystal_name]:
            if forcefield_name.lower() == "experiment":
                continue
            workbook_wrap.add_forcefield(forcefield_name, crystal_name)

            xtal_ucinfo     = input_dict[crystal_name][forcefield_name]["xtal-ucinfo"]
            xtal_topology   = input_dict[crystal_name][forcefield_name]["xtal-topology"]
            xtal_trajectory = input_dict[crystal_name][forcefield_name]["xtal-trajectory"]
            xtal_output     = input_dict[crystal_name][forcefield_name]["xtal-output"]

            gas_output = input_dict[crystal_name][forcefield_name]["gas-output"]

            ucinfo = read_csv(xtal_ucinfo)

            a_idxs = [ucinfo["unitcell_in_supercell_a"][mol_idx] for mol_idx in range(ucinfo["N_rows"])]
            b_idxs = [ucinfo["unitcell_in_supercell_b"][mol_idx] for mol_idx in range(ucinfo["N_rows"])]
            c_idxs = [ucinfo["unitcell_in_supercell_c"][mol_idx] for mol_idx in range(ucinfo["N_rows"])]

            N_unitcells_a = np.max(a_idxs) - np.min(a_idxs)
            N_unitcells_b = np.max(b_idxs) - np.min(b_idxs)
            N_unitcells_c = np.max(c_idxs) - np.min(c_idxs)

            ### N_rows is number of molecules in supercell info csv
            N_molecules = ucinfo["N_rows"]

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
            pc_diffs       = [list(), list(), list()]

            hbond_O_single_diffs = list()
            hbond_O_double_diffs = list()
            hbond_N_single_diffs = list()
            hbond_N_double_diffs = list()

            for output_csv in glob.glob(xtal_output):
                data = read_csv(output_csv)

                ### N_rows is number of frames in mdtraj output csv
                box_a  = np.zeros((data["N_rows"], 3), dtype=float)
                box_b  = np.zeros((data["N_rows"], 3), dtype=float)
                box_c  = np.zeros((data["N_rows"], 3), dtype=float)

                box_a[:,0] = data["Box-XX"]
                box_b[:,0] = data["Box-YX"]
                box_b[:,1] = data["Box-YY"]
                box_c[:,0] = data["Box-ZX"]
                box_c[:,1] = data["Box-ZY"]
                box_c[:,2] = data["Box-ZZ"]

                ### Calc vector norm
                norm_a = np.linalg.norm(box_a, axis=1)
                norm_b = np.linalg.norm(box_b, axis=1)
                norm_c = np.linalg.norm(box_c, axis=1)

                ### Calc vector norm and convert to Ang
                _a_len = norm_a / float(N_unitcells_a) * _NM_2_ANG
                _b_len = norm_b / float(N_unitcells_b) * _NM_2_ANG
                _c_len = norm_c / float(N_unitcells_c) * _NM_2_ANG

                a_len.extend(_a_len.tolist())
                b_len.extend(_b_len.tolist())
                c_len.extend(_c_len.tolist())

                ### Normalize boxvectors and calc box angles.
                _a_len = np.linalg.norm(box_a, axis=1)
                _b_len = np.linalg.norm(box_b, axis=1)
                _c_len = np.linalg.norm(box_c, axis=1)

                box_a = (box_a.T / norm_a).T
                box_b = (box_b.T / norm_b).T
                box_c = (box_c.T / norm_c).T

                _alpha = np.arccos(np.einsum('ij,ij->i', box_b, box_c)) * 180. / np.pi
                _beta  = np.arccos(np.einsum('ij,ij->i', box_a, box_c)) * 180. / np.pi
                _gamma = np.arccos(np.einsum('ij,ij->i', box_a, box_b)) * 180. / np.pi

                alpha.extend(_alpha.tolist())
                beta.extend(_beta.tolist())
                gamma.extend(_gamma.tolist())

                _potential = data["Potential"] # Potential energy in kJ/mol
                ene_xtal.extend(_potential.tolist())

                _density = data["Density"]
                density_xtal.extend(_density.tolist())

            for output_csv in glob.glob(gas_output):

                data = read_csv(output_csv)
                _potential = data["Potential"] # Potential energy in kJ/mol
                ene_gas.extend(_potential)

            for output_dcd in glob.glob(xtal_trajectory):

                query_traj = md.load(
                    output_dcd,
                    top=ref_strc.topology
                    )

                _com_diffs      = analysis_engine.compute_com_diff_per_residue(query_traj, ref_strc)
                _pc_diffs       = analysis_engine.compute_pc_diff_per_residue(query_traj, ref_strc)
                _distance_diffs = ref_distances - analysis_engine.compute_pairwise_distances(query_traj, dist_pair_list)
                _distance_diffs = np.abs(_distance_diffs) * _NM_2_ANG
                _com_diffs      = np.abs(_com_diffs) * _NM_2_ANG

                com_diffs.extend(_com_diffs.tolist())
                pc_diffs[0].extend(_pc_diffs[:,:,0].tolist())
                pc_diffs[1].extend(_pc_diffs[:,:,1].tolist())
                pc_diffs[2].extend(_pc_diffs[:,:,2].tolist())
                distance_diffs.extend(_distance_diffs.tolist())

                if acc_O_single_pair_list.size > 0:
                    _hbond_O_single_diffs = analysis_engine.compute_pairwise_distances(
                        query_traj, 
                        acc_O_single_pair_list
                        )
                    diffs = ref_hbond_O_double_diffs - _hbond_O_single_diffs
                    hbond_O_single_diffs.extend(diffs.tolist())

                if acc_O_double_pair_list.size > 0:
                    _hbond_O_double_diffs = analysis_engine.compute_pairwise_distances(
                        query_traj, 
                        acc_O_double_pair_list
                        )
                    diffs = ref_hbond_O_double_diffs - _hbond_O_double_diffs
                    hbond_O_double_diffs.extend(diffs.tolist())

                if acc_N_single_pair_list.size > 0:
                    _hbond_N_single_diffs = analysis_engine.compute_pairwise_distances(
                        query_traj, 
                        acc_N_single_pair_list
                        )
                    diffs = ref_hbond_N_single_diffs - _hbond_N_single_diffs
                    hbond_N_single_diffs.extend(diffs.tolist())

                if acc_N_double_pair_list.size > 0:
                    _hbond_N_double_diffs = analysis_engine.compute_pairwise_distances(
                        query_traj, 
                        acc_N_double_pair_list
                        )
                    diffs = ref_hbond_N_double_diffs - _hbond_N_double_diffs
                    hbond_N_double_diffs.extend(diffs.tolist())

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

            ### Write com diff data ###
            ### =================== ###
            avg  = np.mean(com_diffs)
            std  = np.std(com_diffs)
            workbook_wrap.add_data(
                avg,
                std/avg*100.,
                "<[D(dCOM)]>",
                forcefield_name, 
                crystal_name
                )

            ### Write pc diff data ###
            ### ================== ###
            avg  = np.mean(pc_diffs[0])
            std  = np.std(pc_diffs[0])
            workbook_wrap.add_data(
                avg,
                std/avg*100.,
                "<{∠PA1-PA1(expt)}>", 
                forcefield_name, 
                crystal_name
                )

            avg  = np.mean(pc_diffs[1])
            std  = np.std(pc_diffs[1])
            workbook_wrap.add_data(
                avg,
                std/avg*100.,
                "<{∠PA2-PA2(expt)}>", 
                forcefield_name, 
                crystal_name
                )

            avg  = np.mean(pc_diffs[2])
            std  = np.std(pc_diffs[2])
            workbook_wrap.add_data(
                avg,
                std/avg*100.,
                "<{∠PA3-PA3(expt)}>", 
                forcefield_name, 
                crystal_name
                )

            ### Write hbond data ###
            ### ================ ###

            if len(hbond_O_single_diffs) > 0:
                avg  = np.mean(hbond_O_single_diffs)
                std  = np.std(hbond_O_single_diffs)
                std  = std/avg*100.
            else:
                avg = "--"
                std = "--"
            workbook_wrap.add_data(
                avg,
                std,
                "X-H•••O", 
                forcefield_name, 
                crystal_name
                )

            if len(hbond_O_double_diffs) > 0:
                avg  = np.mean(hbond_O_double_diffs)
                std  = np.std(hbond_O_double_diffs)
                std  = std/avg*100.
            else:
                avg = "--"
                std = "--"
            workbook_wrap.add_data(
                avg,
                std,
                "X-H•••O=C", 
                forcefield_name, 
                crystal_name
                )

            if len(hbond_N_single_diffs) > 0:
                avg  = np.mean(hbond_N_single_diffs)
                std  = np.std(hbond_N_single_diffs)
                std  = std/avg*100.
            else:
                avg = "--"
                std = "--"
            workbook_wrap.add_data(
                avg,
                std,
                "X-H•••N", 
                forcefield_name, 
                crystal_name
                )

            if len(hbond_N_double_diffs) > 0:
                avg  = np.mean(hbond_N_double_diffs)
                std  = np.std(hbond_N_double_diffs)
                std  = std/avg*100.
            else:
                avg = "--"
                std = "--"
            workbook_wrap.add_data(
                avg,
                std,
                "X-H•••N=C", 
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
            ene_xtal = np.array(ene_xtal)
            ene_gas = np.array(ene_gas)
            sublimation_avg  = np.mean(ene_xtal)/float(N_molecules) - np.mean(ene_gas)
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
            density = float(density)
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