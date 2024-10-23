#!/bin/bash
#SBATCH --job-name=ctx
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=220G
#SBATCH --output=py_ctx.%j.out
#SBATCH --error=py_ctx.%j.err
#SBATCH --dependency=afterok:3873582

##1=loom

ml apps/pyscenic

pyscenic ctx outs/adj.tsv \
/data/stemcell/jwhittle/ref/scenic/hg38__refseq-r80__10kb_up_and_down_tss.genes_vs_motifs.rankings.feather \
/data/stemcell/jwhittle/ref/scenic/encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.genes_vs_tracks.rankings.feather \
--annotations_fname /data/stemcell/jwhittle/ref/scenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname ${1} --output outs/reg.csv --mask_dropouts --num_workers 48
