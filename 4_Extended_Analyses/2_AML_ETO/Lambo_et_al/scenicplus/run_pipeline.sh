#!/bin/bash
#SBATCH --job-name=scplus
#SBATCH --partition=hmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=2000G
#SBATCH --output=scplus.%j.out
#SBATCH --error=scplus.%j.err
#SBATCH --time=2-00:00:00
#SBATCH --dependency=afterok:5173601

source /data/stemcell/jwhittle/mambaforge/etc/profile.d/conda.sh

conda activate scenic-plus

snakemake --cores 48

conda deactivate
