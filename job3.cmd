#!/bin/bash
#SBATCH -J 08_run #extraction          # Job name
#SBATCH -o out.%j       # Name of stdout output file
#SBATCH -e error.%j       # Name of stderr error file
#SBATCH -p skx-dev #icx-normal          # Queue name
#SBATCH -N 1               # Total # of nodes (now required)
#SBATCH -n 48              # Total # of mpi tasks, 320 = 5 node * 64 cores/node
#SBATCH -t 02:00:00        # Run time (hh:mm:ss)
#SBATCH -A TG-DMR970008S       # Allocation name (req'd if more than 1)
#SBATCH --mail-user=computation.management@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end    # email me when the job finishes

export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export LD_PRELOAD=/home1/apps/tacc-patches/getcwd-patch.so:$LD_PRELOAD


conda activate base
python count.py
