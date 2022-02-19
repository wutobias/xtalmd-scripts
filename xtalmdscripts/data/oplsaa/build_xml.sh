#!/bin/bash

if [ ! $# -eq 4 ]; then
    echo "Usage: build_xml.sh <molfile> <resname> <charge> <1/0 for LBCC>"
    exit
fi

infile=$1
resname=$2
charge=$3
lbcc=$4

echo "BOSSdir $BOSSdir ..."

TMPDIR="/tmp/ligpargen${RANDOM}"

if [ $lbcc_flag ]; then
    ligpargen --ifile ${infile} \
              --path ${TMPDIR} \
              --molname ${resname} \
              --resname ${resname} \
              --charge ${charge} \
              --cgen CM1A-LBCC
else
    ligpargen --ifile ${infile} \
              --path ${TMPDIR} \
              --molname ${resname} \
              --resname ${resname} \
              --charge ${charge} \
              --cgen CM1A
fi

molfilename=`basename -s .pdb ${infile}`
molfilename=`basename -s .mol ${infile}`


cp -f ${TMPDIR}/${resname}.openmm.xml .
cp -f ${TMPDIR}/${resname}.openmm.pdb .
rm -rf ${TMPDIR}

echo "Finished."
