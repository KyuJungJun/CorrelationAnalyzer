#!/bin/bash
# FILENAME:  myjobsubmissionfile

#SBATCH -A dmr970008  #yallocation  # Allocation name
#SBATCH --nodes=1        # Total # of nodes 
#SBATCH --ntasks=128     # Total # of MPI tasks
#SBATCH --time=04:00:00   # Total run time limit (hh:mm:ss)
#SBATCH -J vasptest    # Job name
#SBATCH -o myjob.o%j     # Name of stdout output file
#SBATCH -e myjob.e%j     # Name of stderr error file
#SBATCH -p wholenode     # Queue (partition) name
#SBATCH --mail-user=computation.management@gmail.com #useremailaddress
#SBATCH--mail-type=all   # Send email to above address at begin and end of job

# Manage processing environment, load compilers and applications.
module purge
#module load compilername
#module load mpilibrary
#module load applicationname
module load intel
module load intel-mkl
module list

export I_MPI_FABRICS=shm
export OMP_NUM_THREADS=1
conda activate diffusion

python myhops_2023.py 600;
