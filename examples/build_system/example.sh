#!/bin/bash

build_system -i ../data/oxamide/1473522.cif -o oxamide/1473522_gaff1.xml -ff gaff1
build_system -i ../data/oxamide/1473522.cif -o oxamide/1473522_gaff2.xml -ff gaff2
build_system -i ../data/adipamide/1101307.cif -o adipamide/1101307.xml -ff gaff1
build_system -i ../data/adipamide/1101307.cif -o adipamide/1101307.xml -ff gaff2