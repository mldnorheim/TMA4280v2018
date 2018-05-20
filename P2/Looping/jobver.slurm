#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:15:00
#SBATCH --nodes=2
# Number of MPI processes per node
#SBATCH --ntasks-per-node=2
# Number of OpenMP threads per process
#SBATCH -c 8
#SBATCH --mem=12000
# Memory per node specification is in MB. It is optional. 
# The default limit is 3000MB per core.
#SBATCH --job-name="test"
#SBATCH --output=testPoissonverP4t8_4.out
#SBATCH --mail-user=mlnorhei@ntnu.no
#SBATCH --mail-type=ALL

module load GCC OpenMPI
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

mpirun ./poissonver
