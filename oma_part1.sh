#!/bin/bash
#SBATCH --partition=<your_partition>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --job-name=oma1
#SBATCH --output=logs/oma1-%J.log
#SBATCH --export=None
#SBATCH --error=logs/oma1-%J.err

# This is the database conversion part
cd <full_path_to_OMA.latest> ./bin/oma -c
