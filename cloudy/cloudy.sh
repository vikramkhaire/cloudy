#!/bin/bash
#PBS -l select=2:ncpus=23:mpiprocs=23
#PBS -l walltime=400:00:00
#PBS -N test
#PBS -q express
#PBS -o /home/abhisek/soft/garbage/$PBS_JOBNAME.$PBS_ARRAYID
#PBS -e /home/abhisek/soft/garbage/err_file.txt 
cd $PBS_O_WORKDIR
source /home/abhisek/.bashrc
python /home/abhisek/soft/cld/grid_cloudy.py
