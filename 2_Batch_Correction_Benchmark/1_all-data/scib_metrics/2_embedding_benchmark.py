#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
################################################## AML ATLAS BATCH CORRECTION ###################################################
######################################## Batch Correction Benchmarking with scIB Metrics ########################################
#################################################################################################################################


import os
import warnings
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
from scib_metrics.nearest_neighbors import NeighborsOutput
from scib_metrics.benchmark import Benchmarker, BioConservation
import time
from rich import print

warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=DeprecationWarning)
warnings.simplefilter(action="ignore", category=UserWarning)
sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(4, 4))

# Set Batch Key
benchmarker_batch = sys.argv[1]

# Set Label Key
benchmarker_label = sys.argv[2]

# Select Embedding
embeddings = ["X_uncorrected", "X_harmony", "X_scvi_study", "X_scvi_sample",
              "X_scanvi_study", "X_scanvi_sample"]
emb = int(sys.argv[3])
embedding_key = embeddings[emb]

benchmark_id = benchmarker_batch + "_" + benchmarker_label + "_" + embedding_key
print(benchmark_id)

# Load Anndata
adata = sc.read_h5ad("outs/combined_embeddings.h5ad")



#################################################################################################################################
##################################################### Run scIB Benchmarker ######################################################
#################################################################################################################################

biocons = BioConservation(isolated_labels=False, nmi_ari_cluster_labels_leiden=False, nmi_ari_cluster_labels_kmeans=False)

start = time.time()
bm = Benchmarker(
    adata,
    batch_key=benchmarker_batch,
    label_key=benchmarker_label,
    embedding_obsm_keys=[embedding_key],
    pre_integrated_embedding_obsm_key="X_pca",
    bio_conservation_metrics=biocons,
    n_jobs=48,
)

bm.benchmark()
end = time.time()
print(f"Time: {int((end - start) / 60)} min {int((end - start) % 60)} sec")

df = bm.get_results(min_max_scale=False)
df.transpose()
df.to_csv("outs/" + benchmark_id + ".csv")

df = df.rename(index={embedding_key:benchmark_id})
df.to_csv("outs/" + benchmark_id + ".csv")
