#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
################################################ AML ATLAS: OPTIMISE CLUSTERING #################################################
###################### Iterate Through Various Leiden KNN Weights and Resolutions - Run As Job Array (n=18) #####################
#################################################################################################################################

from __future__ import print_function
from time import sleep
import warnings
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.metrics import davies_bouldin_score

print("Packages imported")

# Configure scanpy settings
sc.settings.verbosity = 0   ### suppress output of scanpy
warnings.filterwarnings('ignore') # suppress Warnings
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')



#################################################################################################################################
################################################## Set Environment Variables ####################################################
#################################################################################################################################

# Set KNN
knn = int(sys.argv[1])

# Set Resolution Parameter
res = list(range(20, 310, 10))
r = int(sys.argv[2])
resolution = res[r]

print("Reading anndata...")
adata = sc.read_h5ad('../dimensionality-reduction/outs/scvi_dimsred.h5ad')



#################################################################################################################################
####################################################### Run Clustering ##########################################################
#################################################################################################################################

print("Calculating neighbour graph...")
sc.pp.neighbors(adata, n_neighbors=knn, use_rep="X_scVI", key_added="scvi_neighbors")

print("Running UMAP...")
sc.tl.umap(adata, neighbors_key="scvi_neighbors")

j = resolution/100
sc.tl.leiden(adata, resolution=j, neighbors_key="scvi_neighbors")

silhouette_avg = silhouette_score(adata.obsm['X_scVI'], adata.obs['leiden'])

davies_bouldin_avg = davies_bouldin_score(adata.obsm['X_scVI'], adata.obs['leiden'])

df = pd.DataFrame(adata.obs['leiden'].value_counts())
print(df)
value = len(df)-1
min_cell = df["leiden"][value]

result = pd.DataFrame([['leiden', knn, j, max(adata.obs["leiden"].cat.codes), silhouette_avg, davies_bouldin_avg, 
                        min_cell]],columns=['method','KNN', 'resolution','number_of_clusters','sil',
                                            'davie_bould','min_cell_n'])

result.to_csv("outs/" + str(knn) + "_knn_clusters_" + str(resolution) + ".csv")
print(j)

result['number_of_clusters'] = result['number_of_clusters'].astype(float)+1
result.to_csv("outs/" + str(knn) + "_knn_clusters_" + str(resolution) + ".csv")

print("Process complete")



