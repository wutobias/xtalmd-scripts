

def get_params_path_list(toppar_dir_path):

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
###    f"{toppar_dir_path}/toppar_all36_moreions.str",
###    f"{toppar_dir_path}/top_interface.rtf",
###    f"{toppar_dir_path}/par_interface.prm",
###    f"{toppar_dir_path}/toppar_all36_nano_lig.str",
###    f"{toppar_dir_path}/toppar_all36_nano_lig_patch.str",
###    f"{toppar_dir_path}/toppar_all36_synthetic_polymer.str",
###    f"{toppar_dir_path}/toppar_all36_synthetic_polymer_patch.str",
###    f"{toppar_dir_path}/toppar_all36_polymer_solvent.str",
###    f"{toppar_dir_path}/toppar_water_ions.str",
###    f"{toppar_dir_path}/toppar_dum_noble_gases.str",
###    f"{toppar_dir_path}/toppar_ions_won.str",
###    f"{toppar_dir_path}/toppar_all36_prot_arg0.str",
###    f"{toppar_dir_path}/toppar_all36_prot_c36m_d_aminoacids.str",
###    f"{toppar_dir_path}/toppar_all36_prot_fluoro_alkanes.str",
###    f"{toppar_dir_path}/toppar_all36_prot_heme.str",
###    f"{toppar_dir_path}/toppar_all36_prot_na_combined.str",
###    f"{toppar_dir_path}/toppar_all36_prot_retinol.str",
###    f"{toppar_dir_path}/toppar_all36_prot_model.str",
###    f"{toppar_dir_path}/toppar_all36_prot_modify_res.str",
###    f"{toppar_dir_path}/toppar_all36_na_nad_ppi.str",
###    f"{toppar_dir_path}/toppar_all36_na_rna_modified.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_sphingo.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_archaeal.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_bacterial.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_cardiolipin.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_cholesterol.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_dag.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_inositol.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_lnp.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_lps.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_mycobacterial.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_miscellaneous.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_model.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_prot.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_tag.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_yeast.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_hmmm.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_detergent.str",
###    f"{toppar_dir_path}/toppar_all36_lipid_ether.str",
###    f"{toppar_dir_path}/toppar_all36_carb_glycolipid.str",
###    f"{toppar_dir_path}/toppar_all36_carb_glycopeptide.str",
###    f"{toppar_dir_path}/toppar_all36_carb_imlab.str",
###    f"{toppar_dir_path}/toppar_all36_label_spin.str",
###    f"{toppar_dir_path}/toppar_all36_label_fluorophore.str",
    ]

    return params_path_list


def get_charmm_inp_header(toppar_dir_path):

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

!!! 
!!! 
!!! ! Phosphates and Sulfates
!!! !stream {toppar_dir_path}/toppar_all36_moreions.str
!!! 
!!! ! Interface FF
!!! open read card unit 10 name {toppar_dir_path}/top_interface.rtf
!!! read  rtf card unit 10 append
!!! 
!!! open read card unit 10 name {toppar_dir_path}/par_interface.prm
!!! read para card unit 10 append flex
!!! 
!!! stream {toppar_dir_path}/toppar_all36_nano_lig.str
!!! stream {toppar_dir_path}/toppar_all36_nano_lig_patch.str
!!! 
!!! ! Additional topologies and parameters for synthetic polymer
!!! stream {toppar_dir_path}/toppar_all36_synthetic_polymer.str
!!! stream {toppar_dir_path}/toppar_all36_synthetic_polymer_patch.str
!!! stream {toppar_dir_path}/toppar_all36_polymer_solvent.str
!!! 
!!! ! Additional topologies and parameters for water and ions
!!! stream {toppar_dir_path}/toppar_water_ions.str
!!! stream {toppar_dir_path}/toppar_dum_noble_gases.str
!!! stream {toppar_dir_path}/toppar_ions_won.str
!!! 
!!! ! Additional topologies and parameters for protein
!!! stream {toppar_dir_path}/toppar_all36_prot_arg0.str
!!! stream {toppar_dir_path}/toppar_all36_prot_c36m_d_aminoacids.str
!!! stream {toppar_dir_path}/toppar_all36_prot_fluoro_alkanes.str
!!! stream {toppar_dir_path}/toppar_all36_prot_heme.str
!!! stream {toppar_dir_path}/toppar_all36_prot_na_combined.str
!!! stream {toppar_dir_path}/toppar_all36_prot_retinol.str
!!! stream {toppar_dir_path}/toppar_all36_prot_model.str
!!! stream {toppar_dir_path}/toppar_all36_prot_modify_res.str
!!! 
!!! ! Additional topologies and parameters for nucleic acids
!!! stream {toppar_dir_path}/toppar_all36_na_nad_ppi.str
!!! stream {toppar_dir_path}/toppar_all36_na_rna_modified.str
!!! 
!!! ! Additional topologies and parameters for lipids
!!! stream {toppar_dir_path}/toppar_all36_lipid_sphingo.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_archaeal.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_bacterial.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_cardiolipin.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_cholesterol.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_dag.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_inositol.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_lnp.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_lps.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_mycobacterial.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_miscellaneous.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_model.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_prot.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_tag.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_yeast.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_hmmm.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_detergent.str
!!! stream {toppar_dir_path}/toppar_all36_lipid_ether.str
!!! 
!!! ! Additional topologies and parameters for carbohydrates
!!! stream {toppar_dir_path}/toppar_all36_carb_glycolipid.str
!!! stream {toppar_dir_path}/toppar_all36_carb_glycopeptide.str
!!! stream {toppar_dir_path}/toppar_all36_carb_imlab.str
!!! 
!!! ! Additional topologies and parameters for spin/fluorophore labels
!!! stream {toppar_dir_path}/toppar_all36_label_spin.str
!!! stream {toppar_dir_path}/toppar_all36_label_fluorophore.str

"""

    return charmm_inp_header