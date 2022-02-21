#!/bin/bash

run_xtal_min xml -i ../build_system/oxamide/1473522_oplsaa.xml \
                 -p ../build_system/oxamide/1473522_oplsaa.pdb \
                 -pre oxamide_oplsaa_min \
                 --method BFGS \
                 --steps 100
