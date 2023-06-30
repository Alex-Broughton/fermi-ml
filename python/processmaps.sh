#!/bin/bash

#SBATCH -p standard                      ## run on the standard partition
#SBATCH -A SMURGIA_LAB		      ## Specify lab account
#SBATCH -N 5                          ## run on a single node
#SBATCH -n 5                          ## request 1 task (1 CPU)
#SBATCH --mem-per-cpu=64G             ## Mem per cpu
#SBATCH -o /pub/abrought/GALPROPModels/process.out ##STDOUT
#SBATCH -e /pub/abrought/GALPROPModels/process.err ##STDERR
##SBATCH --mail-type=end,time_limit,fail,invalid_depend
##SBATCH --mail-user=abrought@uci.edu  ## use this email address
#SBATCH -t 10:00:00

#python reprojectmaps.py
#python combinemaps.py
python maskmaps.py

