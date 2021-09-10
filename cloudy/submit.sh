#!/bin/bash -login
# Propogate environment variables to compute node
#SBATCH --export=ALL

# Select partition
#SBATCH --partition=iist-all

# set the number of nodes and processes per node
#SBATCH --nodes=7

# set the number of tasks (processes) per node.
#SBATCH --ntasks-per-node=28

# set name of job
#SBATCH --job-name=abhi

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL

# send mail to this address
#SBATCH --mail-user=vikramkhaire@iist.ac.in


echo $SLURM_JOB_ID
echo $SLURM_NPROCS


export PATH="/home/vikram/anaconda3/bin:$PATH"

python grid_cloudy.py

