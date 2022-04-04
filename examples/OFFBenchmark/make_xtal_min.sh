#!/bin/bash -f

mkdir -p sage_min
mkdir -p sage_min_alternating
mkdir -p sage_min_combo
mkdir -p parsley_min
mkdir -p parsley_min_alternating
mkdir -p parsley_min_combo
for cif in `cat ./success_xml.txt`; do

    prefix=`basename -s .cif ${cif}`

    run_xtal_min xml -i sage_xml/${prefix}.xml \
                     -p sage_xml/${prefix}.pdb \
                     -pre sage_min/${prefix} \
                     --method BFGS \
                     --epsilon 1.e-8 \
                     --steps 100

    run_xtal_min xml -i parsley_xml/${prefix}.xml \
                     -p parsley_xml/${prefix}.pdb \
                     -pre parsley_min/${prefix} \
                     --method BFGS \
                     --epsilon 1.e-8 \
                     --steps 100

    run_xtal_min xml -i sage_xml/${prefix}.xml \
                     -p sage_xml/${prefix}.pdb \
                     -pre sage_min_alternating/${prefix} \
                     --alternating \
                     --method BFGS \
                     --epsilon 1.e-8 \
                     --steps 100

    run_xtal_min xml -i parsley_xml/${prefix}.xml \
                     -p parsley_xml/${prefix}.pdb \
                     -pre parsley_min_alternating/${prefix} \
                     --alternating \
                     --method BFGS \
                     --epsilon 1.e-8 \
                     --steps 100

    run_xtal_min xml -i sage_xml/${prefix}.xml \
                     -p sage_min_alternating/${prefix}.pdb \
                     -pre sage_min_combo/${prefix} \
                     --method BFGS \
                     --epsilon 1.e-8 \
                     --steps 100

    run_xtal_min xml -i parsley_xml/${prefix}.xml \
                     -p parsley_min_alternating/${prefix}.pdb \
                     -pre parsley_min_combo/${prefix} \
                     --method BFGS \
                     --epsilon 1.e-8 \
                     --steps 100


done

