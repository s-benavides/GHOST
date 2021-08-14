#!/bin/sh
#SBATCH --nodes=1          # Number of nodes (1 node = 16 cores)
#SBATCH --tasks-per-node=1 # 128^3: 16 procs, 256^3: [32-64],  512^2:[128-512]
#SBATCH --time=0-00:05:00   # 1 day and 3 hours
#SBATCH --constraint=centos6
##SBATCH -p newnodes      # partition name
#SBATCH -p sched_mit_hill
##SBATCH -p sched_any_quicktest
#SBATCH -J Compile  # sensible name for the job

## load up the correct modules, if required
. /etc/profile.d/modules.sh
module load engaging/openmpi/2.0.3

echo $SLURM_JOB_NODELIST

## launch the code
make dist main

echo FINISHED
