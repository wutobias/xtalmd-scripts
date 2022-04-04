#!/bin/bash -f

echo "" > success_xml.txt
echo "" > failed_xml.txt
echo "" > success_pdb.txt
echo "" > failed_pdb.txt

mkdir -p sage_xml
mkdir -p parsley_xml
for cif in CIF/*.cif; do

    prefix=`basename -s .cif ${cif}`
    timeout 10s build_system -i CIF/${prefix}.cif \
                 -pre sage_xml/${prefix} \
                 -ff sage \
                 -ax 3.5
    timeout 10s build_system -i CIF/${prefix}.cif \
                 -pre parsley_xml/${prefix} \
                 -ff parsley \
                 -ax 3.5
    ### PDB success/fail
    if [ -f sage_xml/${prefix}.pdb ]; then
        echo "CIF/${prefix}.cif" >> success_pdb.txt
    else
        echo "CIF/${prefix}.cif" >> failed_pdb.txt
    fi

    ### XML success/fail
    if [ -f sage_xml/${prefix}.xml ]; then
        echo "CIF/${prefix}.cif" >> success_xml.txt
    else
        echo "CIF/${prefix}.cif" >> failed_xml.txt
    fi
done