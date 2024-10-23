#!/bin/bash
#SBATCH --job-name=aucell
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=220G
#SBATCH --output=py_aucell.%j.out
#SBATCH --error=py_aucell.%j.err
#SBATCH --dependency=afterok:3968640

##1=loom file

ml apps/pyscenic

pyscenic aucell ${1} outs/reg.csv --output "outs/$(basename -s .loom "${1}")_AUCell.loom" --num_workers 48
