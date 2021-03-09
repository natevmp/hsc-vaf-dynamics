#!/bin/bash
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request core
#$ -l h_rt=240:0:0  # Request hour runtime
#$ -l h_vmem=32G   # Request 1GB RAM
module load julia
julia runBigvafSim.jl
