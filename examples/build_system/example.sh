#!/bin/bash

toppar="/nfs/data_onsager/projects/xtalmd/toppar_c36_jul21"

build_system -i ../data/oxamide/1473522.cif   -pre oxamide/1473522_gaff1   -ff gaff1
build_system -i ../data/oxamide/1473522.cif   -pre oxamide/1473522_gaff2   -ff gaff2
build_system -i ../data/oxamide/1473522.cif   -pre oxamide/1473522_parsley -ff parsley
build_system -i ../data/oxamide/1473522.cif   -pre oxamide/1473522_sage    -ff sage
build_system -i ../data/oxamide/1473522.cif   -pre oxamide/1473522_cgenff  -ff cgenff --toppar ${toppar}
build_system -i ../data/oxamide/1473522.cif   -pre oxamide/1473522_oplsaa  -ff oplsaa

build_system -i ../data/adipamide/1101307.cif   -pre adipamide/1101307_gaff1   -ff gaff1
build_system -i ../data/adipamide/1101307.cif   -pre adipamide/1101307_gaff2   -ff gaff2
build_system -i ../data/adipamide/1101307.cif   -pre adipamide/1101307_parsley -ff parsley
build_system -i ../data/adipamide/1101307.cif   -pre adipamide/1101307_sage    -ff sage
build_system -i ../data/adipamide/1101307.cif   -pre adipamide/1101307_cgenff  -ff cgenff --toppar ${toppar}
build_system -i ../data/adipamide/1101307.cif   -pre adipamide/1101307_oplsaa  -ff oplsaa
