#!/bin/bash

run_xtal_md xml -i ../build_system/oxamide/1473522_gaff1.xml \
                -p ../xtal_min/oxamide-gaff1/1473522_gaff1_min.pdb \
                --temperature=298.15 \
                --nanoseconds 2 \
                --prefix oxamide_xtal
run_xtal_md xml -i ../build_system/oxamide/1473522_gaff1_monomer0.xml \
                -p ../build_system/oxamide/1473522_gaff1_monomer0.pdb \
                --nvt \
                --temperature=298.15 \
                --nanoseconds 2 \
                --prefix oxamide_gas
