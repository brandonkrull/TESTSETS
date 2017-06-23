#!/bin/bash -l
#SBATCH -J napthalene-opt
#SBATCH -p regular
#SBATCH -N 1
#SBATCH -t 50:00
#SBATCH -o test_nwchem.o%j
#SBATCH -C haswell 

module load nwchem 
srun -n 4 nwchem test.nw >& test1.out
