#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=50:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1  # number of nodes
#SBATCH --ntasks-per-node=36 # 36 processor core(s) per node
#SBATCH --mem=100G   # maximum memory per node
#SBATCH --job-name="GRN-code"
#SBATCH --mail-user=muskan@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="12Aug'24-genie3.out" # job standard output file (%j replaced by job id)
#SBATCH --error="Error_genie3-12Aug'24" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load miniconda3
module load python

#pip install numpy scipy

# activating virtual environemnt

source activate cis

# path to working directory
cd /work/ABG/MUSKAN/Terra-Test

# Run python script
python new-genie3.py
