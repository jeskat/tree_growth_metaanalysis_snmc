#!/bin/bash


#SBATCH --job-name=chain1
#SBATCH --account=fc_lmklab
#SBATCH --partition=savio3_htc
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#
## Command(s) to run
module load gcc r 

parallel -j 20 -a State_space_growth_models/Code/scripts_for_parallel_processing/pft.lst -a State_space_growth_models/Code/scripts_for_parallel_processing/site.lst -a State_space_growth_models/Code/scripts_for_parallel_processing/model.lst bash State_space_growth_models/Code/scripts_for_parallel_processing/args_ssm.sh {1} {2} {3} 