#!/bin/bash
#SBATCH --job-name=arboerto
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=244G
#SBATCH --time=7-00:00:00
#SBATCH --output=py_arboreto.%j.out
#SBATCH --error=py_arboreto.%j.err
#SBATCH --dependency=afterok:3873581

##1=loom

source /data/stemcell/jwhittle/mambaforge/etc/profile.d/conda.sh

conda activate scenic

mkdir -p outs

arboreto_with_multiprocessing.py ${1} /data/stemcell/jwhittle/ref/scenic/allTFs_hg38.txt \
--num_workers 48 --method grnboost2 --output outs/adj.tsv  --seed 777

conda deactivate
