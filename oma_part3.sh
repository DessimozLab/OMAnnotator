#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50GB
#SBATCH --time=72:00:00
#SBATCH --job-name=oma3
#SBATCH --output=logs/oma3-%J.log
#SBATCH --export=None
#SBATCH --error=logs/oma3-%J.err

## This step the orthology inference phase. HOGs are inferred from the AllAll alignments computed in step 2. 
./bin/oma
