#!/bin/bash

#SBATCH -J emcee_Mr18_Consuelo
#SBATCH -o /scratch/02497/jenpi/Outputs/Consuelo/Main18/emcee_Consuelo_Mr18.o%j 
#SBATCH -e /scratch/02497/jenpi/Outputs/Consuelo/Main18/emcee_Consuelo_Mr18.e%j
#SBATCH -n 256
#SBATCH -p normal 
#SBATCH -t 24:00:00                              
#SBATCH  --mail-type=ALL
#SBATCH --mail-user=j.piscionere@vanderbilt.edu
ibrun /opt/apps/intel13/mvapich2_1_9/mpi4py/1.3/lib/python/mpi4py/bin/python-mpi /home1/02497/jenpi/Wtheta/emcee_consuelo18_run.py 
