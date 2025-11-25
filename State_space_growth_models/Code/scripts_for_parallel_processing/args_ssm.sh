#!/bin/bash
R --vanilla -f State_space_growth_models/Code/run_ssm.R --args ${1} ${2} ${3} "p1c1" 'c1' > ${1}_${2}_${3}_p1c1.Rout