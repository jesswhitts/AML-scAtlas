#!/bin/bash
#SBATCH --job-name=cis_topic_model
#SBATCH --partition=hmem
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=2000G
#SBATCH --output=cis_topic_model.%j.out
#SBATCH --error=cis_topic_model.%j.err
#SBATCH --time=4-00:00:00
#SBATCH --dependency=afterok:3951498

source /data/stemcell/jwhittle/mambaforge/etc/profile.d/conda.sh

conda activate scenic-plus

pycistopic topic_modeling mallet \
	--input scATAC/cistopic_obj.pkl \
	--output scATAC/models.pkl \
	--temp_dir /scratch/wsspaces/jwhittle-tmp2/ \
	--topics 5 10 15 20 25 30 35 40 45 50 55 60 70 80 90 100 150 200 250 300 \
	--iterations 300 \
	--alpha 50 \
	--parallel 48 \
	--keep True \
	--seed 555 \
	--mallet_path /data/stemcell/jwhittle/scenicplus/Mallet-202108/bin/mallet

conda deactivate
