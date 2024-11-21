#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=2GB
#SBATCH --job-name=oma1
#SBATCH --output=oma1-%J.log
#SBATCH --error=oma1-%J.err

## This script converts your OMA database ready for the AllAll allignment phase
./bin/oma -c
