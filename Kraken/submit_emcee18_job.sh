#!/bin/bash

#PBS -A TG-AST080002N
#PBS -l size=384,walltime=24:00:00
#PBS -N emcee_Mr18_Consuelo
#PBS -o /lustre/scratch/jenpi/Outputs/Consuelo/Main18/emcee_Consuelo_Mr18.o$PBS_JOBID
#PBS -e /lustre/scratch/jenpi/Outputs/Consuelo/Main18/emcee_Consuelo_Mr18.e$PBS_JOBID
#PBS -q small 
#PBS -m abe
#PBS -M j.piscionere@vanderbilt.edu

cd $PBS_O_WORKDIR
module swap PrgEnv-pgi PrgEnv-gnu
module swap python/2.7.1 python/2.7.1-cnl
aprun -n 252 -S 4 python emcee_consuelo18_run.py 
