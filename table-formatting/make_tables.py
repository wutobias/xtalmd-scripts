#!/usr/bin/env python3

import xlsxwriter
import numpy as np
import os
import argparse
import json
from collections import OrderedDict
import gemmi

_KJ_2_KCAL = 1./4.184
_NM_2_ANG  = 10.

def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Python script for merging simulation data for xtal MD project."
        )

    parser.add_argument('--input',  "-i", type=str, help="Input file", required=True)
    parser.add_argument('--output', "-o", type=str, help="Output xlsx MS Excel file", required=True)

    return parser.parse_args()


def parse_data(file_path, file_format=None):

    if file_format == None:
        filename, file_format = os.path.splitext(file_path)
        file_format = file_format.replace(".","")

    if file_format.lower() == "xvg":

        data_dict   = dict()
        legend_list = list()
        with open(file_path, "r") as fopen:
            for line in fopen:
                line = line.lstrip().rstrip().split()
                if line[0].startswith("#"):
                    continue
                elif line[0] == "gmx" and\
                   line[1] == "energy" and\
                   "-nmol" in " ".join(line):
                   nmol_index = line.index("-nmol")
                   data_dict["nmol"] = int(line[nmol_index + 1])
                   
                elif line[0].startswith("@"):
                    if len(line) < 4:
                        continue
                    if line[2] == "legend":
                        legend = line[3].replace('"',"")
                        legend_list.append(legend)
                        data_dict[legend] = list()
                    continue
                else:
                    for legend_idx, value in enumerate(line[1:]):
                        legend = legend_list[legend_idx]
                        data_dict[legend].append(float(value))
                        
        for legend in legend_list:
            data_dict[legend] = np.array(data_dict[legend])
    else:
        raise ValueError(
            f"File format {file_format} not known."
            )

    return data_dict


def main():

    args       = parse_arguments()
    with open(args.input, "r") as fopen:
        input_list = json.load(fopen)

    workbook  = xlsxwriter.Workbook(args.output)

    ### Define some formats for later use
    ### =================================
    header_format_1 = workbook.add_format(
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

    header_format_2 = workbook.add_format(
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

    data_format_1 = workbook.add_format(
        {
            'align'  : 'center',
            'valign' : 'vcenter',
        }
    )

    ### Loop over each xtal and put in new notebook
    ### ===========================================
    for crystal_name in input_list:
        worksheet      = workbook.add_worksheet(crystal_name)
        N_force_fields = len(input_list[crystal_name]["FORCEFIELD"])

        ### Set column width and row height for header
        ### ==========================================
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
            header_format_1
        )

        worksheet.merge_range(
            first_row=0,
            first_col=1, 
            last_row=0, 
            last_col=2 * N_force_fields,
            data="Force Field",
            cell_format=header_format_1
        )

        ### Write row labels
        ### ================
        worksheet.write(
            1,
            0,
            crystal_name, 
            header_format_1
        )

        worksheet.write(
            2,
            0,
            "",
            header_format_2
        )

        labels_dict_row = OrderedDict()
        labels_dict_row["Sublimation Energy"]   = 3
        labels_dict_row["a"]                    = 4
        labels_dict_row["b"]                    = 5
        labels_dict_row["c"]                    = 6
        labels_dict_row["α"]                    = 7
        labels_dict_row["β"]                    = 8
        labels_dict_row["γ"]                    = 9
        labels_dict_row["<[Δ(d < 4Å)]>"]        = 10
        labels_dict_row["Max Δ(d < 4Å)"]        = 11
        labels_dict_row["H-bond geometry"]      = 12
        labels_dict_row["X-H•••O1"]             = 13
        labels_dict_row["H•••O1=C"]             = 14
        labels_dict_row["X-H•••O2"]             = 15
        labels_dict_row["H•••O2=C"]             = 16
        labels_dict_row["Translation/Rotation"] = 17
        labels_dict_row["<[D(dCOM)]>"]          = 18
        labels_dict_row["Max D(dCOM)"]          = 19
        labels_dict_row["<{∠PA1-PA1(expt)}>"]   = 20
        labels_dict_row["<{∠PA2-PA2(expt)}>"]   = 21
        labels_dict_row["<{∠PA3-PA3(expt)}>"]   = 22
        labels_dict_row["Max ∠PA1-PA1(expt)"]   = 23
        labels_dict_row["Max ∠PA2-PA2(expt)"]   = 24
        labels_dict_row["Max ∠PA3-PA3(expt)"]   = 25

        for label, row_idx in labels_dict_row.items():
            worksheet.write(
                row_idx,
                0,
                label,
                data_format_1
            )

        ### Write data for each FF
        ### ======================
        force_field_dict_list = input_list[crystal_name]["FORCEFIELD"]
        for force_field_idx, force_field_dicts in enumerate(force_field_dict_list):
            xtal_energy_data = parse_data(force_field_dicts["XTAL_ENERGY"])
            xtal_cell_data   = parse_data(force_field_dicts["XTAL_CELL"])
            gas_energy_data  = parse_data(force_field_dicts["GAS_ENERGY"])
            ff_name          = force_field_dicts["NAME"]

            worksheet.merge_range(
                first_row=1,
                first_col=1 + force_field_idx * 2, 
                last_row=1, 
                last_col=2 + force_field_idx * 2,
                data=ff_name,
                cell_format=header_format_1
            )

            worksheet.write(
                2,
                1 + force_field_idx * 2,
                "< >",
                header_format_2
            )

            worksheet.write(
                2,
                2 + force_field_idx * 2,
                "% Dev",
                header_format_2
            )

            ### Compute Cell data
            ### =================
            a = np.stack((
                xtal_cell_data["Box-XX"],
                xtal_cell_data["Box-YX"],
                xtal_cell_data["Box-ZX"]
                )).T * _NM_2_ANG
            a_len = np.linalg.norm(a, axis=1)

            b = np.stack((
                xtal_cell_data["Box-YX"],
                xtal_cell_data["Box-YY"],
                xtal_cell_data["Box-ZY"]
                )).T * _NM_2_ANG
            b_len = np.linalg.norm(b, axis=1)

            c = np.stack((
                xtal_cell_data["Box-ZX"],
                xtal_cell_data["Box-ZY"],
                xtal_cell_data["Box-ZZ"]
                )).T * _NM_2_ANG
            c_len = np.linalg.norm(c, axis=1)
            
            a = (a.T / a_len).T
            b = (b.T / b_len).T
            c = (c.T / c_len).T

            ### Calculate cell angles
            ### =====================
            alpha = np.arccos(np.einsum('ij,ij->i', b, c)) * 180. / np.pi
            beta  = np.arccos(np.einsum('ij,ij->i', a, c)) * 180. / np.pi
            gamma = np.arccos(np.einsum('ij,ij->i', a, b)) * 180. / np.pi

            ### Correct unit cell lengths
            ### =========================
            a_len /= float(force_field_dicts["UNITCELLS_A"])
            b_len /= float(force_field_dicts["UNITCELLS_B"])
            c_len /= float(force_field_dicts["UNITCELLS_C"])

            worksheet.write(
                labels_dict_row["a"],
                1 + force_field_idx * 2,
                np.mean(a_len),
                data_format_1
            )
            worksheet.write(
                labels_dict_row["a"],
                2 + force_field_idx * 2,
                np.std(a_len)/np.mean(a_len)*100.,
                data_format_1
            )

            worksheet.write(
                labels_dict_row["b"],
                1 + force_field_idx * 2,
                np.mean(b_len),
                data_format_1
            )
            worksheet.write(
                labels_dict_row["b"],
                2 + force_field_idx * 2,
                np.std(b_len)/np.mean(b_len)*100.,
                data_format_1
            )

            worksheet.write(
                labels_dict_row["c"],
                1 + force_field_idx * 2,
                np.mean(c_len),
                data_format_1
            )
            worksheet.write(
                labels_dict_row["c"],
                2 + force_field_idx * 2,
                np.std(c_len)/np.mean(c_len)*100.,
                data_format_1
            )

            worksheet.write(
                labels_dict_row["α"],
                1 + force_field_idx * 2,
                np.mean(alpha),
                data_format_1
            )
            worksheet.write(
                labels_dict_row["α"],
                2 + force_field_idx * 2,
                np.std(alpha)/np.mean(alpha)*100.,
                data_format_1
            )

            worksheet.write(
                labels_dict_row["β"],
                1 + force_field_idx * 2,
                np.mean(beta),
                data_format_1
            )
            worksheet.write(
                labels_dict_row["β"],
                2 + force_field_idx * 2,
                np.std(beta)/np.mean(beta)*100.,
                data_format_1
            )

            worksheet.write(
                labels_dict_row["γ"],
                1 + force_field_idx * 2,
                np.mean(gamma),
                data_format_1
            )
            worksheet.write(
                labels_dict_row["γ"],
                2 + force_field_idx * 2,
                np.std(gamma)/np.mean(gamma)*100.,
                data_format_1
            )

            ### Write sublimation enthalpy
            ### ==========================
            
            ### If we have "nmol" in xtal_energy_data, it means we
            ### the field "Potential" has been corrected for the
            ### number of molecules.
            if "nmol" in xtal_energy_data:
                num_molecules_total = 1. 
            
            num_molecules_total      = force_field_dicts["UNITCELLS_A"]
            num_molecules_total     *= force_field_dicts["UNITCELLS_B"]
            num_molecules_total     *= force_field_dicts["UNITCELLS_C"]
            num_molecules_total     *= force_field_dicts["MOLS_PER_CELL"]
            sublimation_energy_mean  = np.mean(xtal_energy_data["Potential"]/num_molecules_total)
            sublimation_energy_mean -= np.mean(gas_energy_data["Potential"])
            sublimation_energy_mean *= _KJ_2_KCAL
            ### Error propgation
            sublimation_energy_std   = np.var(xtal_energy_data["Potential"]/num_molecules_total)
            sublimation_energy_std  += np.var(gas_energy_data["Potential"])
            sublimation_energy_std   = np.sqrt(sublimation_energy_std)
            sublimation_energy_std  *= _KJ_2_KCAL
            worksheet.write(
                labels_dict_row["Sublimation Energy"],
                1 + force_field_idx * 2,
                sublimation_energy_mean,
                data_format_1
            )
            if sublimation_energy_mean > 0.:
                worksheet.write(
                    labels_dict_row["Sublimation Energy"],
                    2 + force_field_idx * 2,
                    sublimation_energy_std/sublimation_energy_mean*100.,
                    data_format_1
                )

        ### Write data for expt
        ### ===================
        expt_col = 1 + N_force_fields * 2

        worksheet.merge_range(
            first_row=1,
            first_col=expt_col, 
            last_row=1, 
            last_col=1 + expt_col,
            data="Expt.",
            cell_format=header_format_1
        )

        worksheet.write(
            2,
            expt_col,
            "< >",
            header_format_2
        )

        worksheet.write(
            2,
            1 + expt_col,
            "% Dev",
            header_format_2
        )

        expt_dict = input_list[crystal_name]["EXPT"]
        doc       = gemmi.cif.read(expt_dict["CIF"])[0]
        strc      = gemmi.make_small_structure_from_block(doc)

        worksheet.write(
            labels_dict_row["a"],
            expt_col,
            ### Ang to nm
            strc.cell.a * 0.1,
            data_format_1
        )

        worksheet.write(
            labels_dict_row["b"],
            expt_col,
            ### Ang to nm
            strc.cell.b * 0.1,
            data_format_1
        )

        worksheet.write(
            labels_dict_row["c"],
            expt_col,
            ### Ang to nm
            strc.cell.c * 0.1,
            data_format_1
        )

        worksheet.write(
            labels_dict_row["α"],
            expt_col,
            strc.cell.alpha,
            data_format_1
        )

        worksheet.write(
            labels_dict_row["β"],
            expt_col,
            strc.cell.beta,
            data_format_1
        )

        worksheet.write(
            labels_dict_row["γ"],
            expt_col,
            strc.cell.gamma,
            data_format_1
        )

        worksheet.write(
            labels_dict_row["Sublimation Energy"],
            expt_col,
            expt_dict["SUBLIMATION_ENTHALPY"],
            data_format_1
        )

        worksheet.write(
            labels_dict_row["Sublimation Energy"],
            expt_col + 1,
            expt_dict["SUBLIMATION_ENTHALPY_ERROR"],
            data_format_1
        )

    workbook.close()

if __name__ == "__main__":

    main()