
def get_params_path_list_charmm36m(toppar_dir_path):

    params_path_list = [
        f"{toppar_dir_path}/top_all36_prot.rtf",
        f"{toppar_dir_path}/par_all36m_prot.prm",
        f"{toppar_dir_path}/top_all36_na.rtf",
        f"{toppar_dir_path}/par_all36_na.prm",
        f"{toppar_dir_path}/top_all36_carb.rtf",
        f"{toppar_dir_path}/par_all36_carb.prm",
        f"{toppar_dir_path}/top_all36_lipid.rtf",
        f"{toppar_dir_path}/par_all36_lipid.prm",
        f"{toppar_dir_path}/top_all36_cgenff.rtf",
        f"{toppar_dir_path}/par_all36_cgenff.prm",
        f"{toppar_dir_path}/toppar_water_ions.str"
        ]

    return params_path_list


def get_params_path_list_charmm36(toppar_dir_path):

    params_path_list = [
        f"{toppar_dir_path}/top_all36_prot.rtf",
        f"{toppar_dir_path}/par_all36_prot.prm",
        f"{toppar_dir_path}/top_all36_na.rtf",
        f"{toppar_dir_path}/par_all36_na.prm",
        f"{toppar_dir_path}/top_all36_carb.rtf",
        f"{toppar_dir_path}/par_all36_carb.prm",
        f"{toppar_dir_path}/top_all36_lipid.rtf",
        f"{toppar_dir_path}/par_all36_lipid.prm",
        f"{toppar_dir_path}/top_all36_cgenff.rtf",
        f"{toppar_dir_path}/par_all36_cgenff.prm",
        f"{toppar_dir_path}/toppar_water_ions.str",
        ]

    return params_path_list


def get_params_path_list_opls(toppar_dir_path):

    params_path_list = [
        f"{toppar_dir_path}/top_opls_aam.inp",
        f"{toppar_dir_path}/par_opls_aam.inp",
        ]

    return params_path_list


def get_charmm_inp_header_charmm36(toppar_dir_path):

    charmm_inp_header = f"""
! protein topology and parameter
open read card unit 10 name {toppar_dir_path}/top_all36_prot.rtf
read  rtf card unit 10

open read card unit 20 name {toppar_dir_path}/par_all36_prot.prm
read para card unit 20 flex

! nucleic acids
open read card unit 10 name {toppar_dir_path}/top_all36_na.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir_path}/par_all36_na.prm
read para card unit 20 append flex

! carbohydrates
open read card unit 10 name {toppar_dir_path}/top_all36_carb.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir_path}/par_all36_carb.prm
read para card unit 20 append flex

! lipids
open read card unit 10 name {toppar_dir_path}/top_all36_lipid.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir_path}/par_all36_lipid.prm
read para card unit 20 append flex

! CGENFF
open read card unit 10 name {toppar_dir_path}/top_all36_cgenff.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir_path}/par_all36_cgenff.prm
read para card unit 20 append flex

stream {toppar_dir_path}/toppar_water_ions.str
"""

    return charmm_inp_header


def get_charmm_inp_header_charmm36m(toppar_dir_path):

    charmm_inp_header = f"""
! protein topology and parameter
open read card unit 10 name {toppar_dir_path}/top_all36_prot.rtf
read  rtf card unit 10

open read card unit 20 name {toppar_dir_path}/par_all36m_prot.prm
read para card unit 20 flex

! nucleic acids
open read card unit 10 name {toppar_dir_path}/top_all36_na.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir_path}/par_all36_na.prm
read para card unit 20 append flex

! carbohydrates
open read card unit 10 name {toppar_dir_path}/top_all36_carb.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir_path}/par_all36_carb.prm
read para card unit 20 append flex

! lipids
open read card unit 10 name {toppar_dir_path}/top_all36_lipid.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir_path}/par_all36_lipid.prm
read para card unit 20 append flex

! CGENFF
open read card unit 10 name {toppar_dir_path}/top_all36_cgenff.rtf
read  rtf card unit 10 append

open read card unit 20 name {toppar_dir_path}/par_all36_cgenff.prm
read para card unit 20 append flex

stream {toppar_dir_path}/toppar_water_ions.str
"""

    return charmm_inp_header


def get_charmm_inp_header_opls(toppar_dir_path):

    charmm_inp_header = f"""
! protein topology and parameter for opls forcefield
open read card unit 10 name {toppar_dir_path}/top_opls_aam.inp
read  rtf card unit 10

open read card unit 20 name {toppar_dir_path}/par_opls_aam.inp
read para card unit 20
"""

    return charmm_inp_header