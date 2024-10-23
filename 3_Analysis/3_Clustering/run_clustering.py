#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
##################################################### AML ATLAS: CLUSTERING #####################################################
########################################### Run Some Selected Clustering Resolutions ############################################
#################################################################################################################################

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.metrics import davies_bouldin_score
import warnings

#configure scanpy settings
sc.settings.verbosity = 0   ### suppress output of scanpy
warnings.filterwarnings('ignore') # suppress Warnings
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

# create output directories
if not os.path.exists("outs"):
    os.makedirs("outs")

#set results output file
results_file = 'outs/scvi_clustered.h5ad'

#read in anndata object
adata = sc.read_h5ad('../dimensionality-reduction/outs/scvi_dimsred.h5ad')



#################################################################################################################################
###################################################### Run Leiden Clustering ####################################################
#################################################################################################################################

print("Running default clusterings")

# Compute neighbourhood graph
sc.pp.neighbors(adata, use_rep="X_scVI", key_added="scvi_neighbors")

#3 - Starndard clustering with default settings
sc.tl.leiden(adata, neighbors_key="scvi_neighbors", key_added = "leiden_nn15_1.0")
sc.pl.umap(adata, color=["leiden_nn15_1.0"], title="Leiden Resolution 1.0", save="_leiden_1.0.png")

#2 - Low resolution, with default nn (15)
sc.tl.leiden(adata, resolution=0.5, neighbors_key="scvi_neighbors", key_added = "leiden_nn15_0.5")
sc.pl.umap(adata, color=["leiden_nn15_0.5"], title="Leiden Resolution 0.5", save="_leiden_0.5.png")

#3 - High resolution, with default nn (15)
sc.tl.leiden(adata, resolution=2.0, neighbors_key="scvi_neighbors", key_added = "leiden_nn15_2.0")
sc.pl.umap(adata, color=["leiden_nn15_2.0"], title="Leiden Resolution 2.0", save="_leiden_2.0.png")


# # Run Optimal Clustering
sc.tl.leiden(adata, resolution=1.4, neighbors_key="scvi_neighbors", key_added = "leiden_nn15_1.4")
adata.obs["cluster"] = adata.obs["leiden_nn15_1.4"]
sc.pl.umap(adata, color=["cluster"], title="Leiden Resolution 1.4, KNN 15", save="_leiden_clusters.png")

sc.tl.dendrogram(adata, groupby = "cluster")
sc.pl.dendrogram(adata, groupby = "cluster", save="_leiden_dendogram.png")

print("Saving file")

adata.write_h5ad(results_file)

print(adata.obs.cluster.value_counts())
