#!/bin/bash -f

mkdir -p sage_min
mkdir -p sage_min_alternating
mkdir -p sage_min_combo
mkdir -p parsley_min
mkdir -p parsley_min_alternating
mkdir -p parsley_min_combo
#for cif in `cat ./success_xml.txt`; do
for cif in "2100147"; do

    prefix=`basename -s .cif ${cif}`
    
    if [ ! "${prefix}" = "2100147" ]; then
        continue
    fi

    run_xtal_min xml -i sage_xml/${prefix}.xml \
                     -p sage_xml/${prefix}.pdb \
                     -pre sage_min/${prefix} \
                     --method L-BFGS-B \
                     --epsilon 1.e-6 \
                     --steps 100

    run_xtal_min xml -i parsley_xml/${prefix}.xml \
                     -p parsley_xml/${prefix}.pdb \
                     -pre parsley_min/${prefix} \
                     --method L-BFGS-B \
                     --epsilon 1.e-6 \
                     --steps 100

    run_xtal_min xml -i sage_xml/${prefix}.xml \
                     -p sage_xml/${prefix}.pdb \
                     -pre sage_min_alternating/${prefix} \
                     --alternating \
                     --method L-BFGS-B \
                     --epsilon 1.e-6 \
                     --steps 100

    run_xtal_min xml -i parsley_xml/${prefix}.xml \
                    -p parsley_xml/${prefix}.pdb \
                    -pre parsley_min_alternating/${prefix} \
                     --alternating \
                     --method L-BFGS-B \
                     --epsilon 1.e-6 \
                     --steps 100

    run_xtal_min xml -i sage_xml/${prefix}.xml \
                     -p sage_min_alternating/${prefix}.pdb \
                     -pre sage_min_combo/${prefix} \
                     --method L-BFGS-B \
                     --epsilon 1.e-6 \
                     --steps 100

    run_xtal_min xml -i parsley_xml/${prefix}.xml \
                    -p parsley_min_alternating/${prefix}.pdb \
                     -pre parsley_min_combo/${prefix} \
                     --method L-BFGS-B \
                     --epsilon 1.e-6 \
                     --steps 100

    exit

done

