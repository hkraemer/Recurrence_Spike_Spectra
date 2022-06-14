#!/bin/bash

#SBATCH --qos=priority
#SBATCH --partition=priority
#SBATCH --job-name=kuramoto
#SBATCH --account=synet
#SBATCH --output=name-%j.out
#SBATCH --error=name-%j.err
#SBATCH --ntasks-per-node=32

module load julia/1.7.0
julia analyze_system_with_transition_cl.jl
