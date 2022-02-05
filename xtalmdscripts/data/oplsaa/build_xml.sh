#!/bin/bash

if [ ! $# -eq 4 ]; then
	echo "Usage: build_xml.sh <pdbfile> <resname> <charge> <1/0 for LBCC>"
	exit
fi

molfile=$1
resname=$2
charge=$3
lbcc=$4

echo "Activate conda env ..."
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/tobias/progs/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/tobias/progs/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/tobias/progs/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/tobias/progs/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
conda activate ligpargen

echo "BOSSdir $BOSSdir ..."

if [ $lbcc_flag ]; then
	LigParGen --mol ${molfile} --resname ${resname} --charge ${charge} --lbcc
else
	LigParGen --mol ${molfile} --resname ${resname} --charge ${charge}
fi

molfilename=`basename -s .mol ${molfile}`

cp -f /tmp/${resname}.xml .
cp -f /tmp/${resname}.pdb .
rm -f /tmp/${resname}.*
rm -f /tmp/${molfilename}.*
rm -f /tmp/LL
rm -f /tmp/slvzmat
rm -f /tmp/clu.pdb
rm -f /tmp/optzmat

sleep 2s

echo "Finished."
