#!/bin/bash

#SBATCH --qos=priority
#SBATCH --partition=priority
#SBATCH --job-name=Roessler
#SBATCH --account=synet
#SBATCH --output=name-%j.out
#SBATCH --error=name-%j.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8


echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "------------------------------------------------------------"

module load julia/1.7.0
julia comm_roessler_cl.jl $SLURM_NTASKS
