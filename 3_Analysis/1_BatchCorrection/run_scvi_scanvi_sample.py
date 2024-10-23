#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
################################################## AML ATLAS BATCH CORRECTION ###################################################
################################################## Batch Correction with scVI ###################################################
#################################################################################################################################


import os
import sys
import warnings
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import anndata as ad
from rich import print
from scvi.model.utils import mde

#configure settings
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
warnings.simplefilter(action="ignore", category=FutureWarning)
sc.set_figure_params(figsize=(4, 4))
scvi.settings.seed = 666 

#create output directories
if not os.path.exists("outs"):
    os.makedirs("outs")

#set results output files
results_file = 'outs/scvi_sample_corrected.h5ad'
embedding_file = "outs/scvi_sample_embedding.h5ad"
batch = "scvi_sample"

#read in data
adata_ref = sc.read_h5ad(
    "../outs/filtered_atlasData_uncorrected.h5ad")


#################################################################################################################################
######################################################## Run scVI Model #########################################################
#################################################################################################################################

print(adata_ref.n_obs, adata_ref.n_vars)

scvi.model.SCVI.setup_anndata(adata_ref, layer="counts", batch_key="sample")

vae = scvi.model.SCVI(adata_ref, n_layers=2, n_latent=30, gene_likelihood="nb")

vae.train()

model_dir = os.path.join("./" + batch + "_model")
vae.save(model_dir, overwrite=True)

adata_ref.obsm["X_scVI"] = vae.get_latent_representation()


#################################################################################################################################
####################################################### Visualise Results #######################################################
#################################################################################################################################

# Run UMAP
sc.pp.neighbors(adata_ref, use_rep="X_scVI", key_added="scVI_neighbors")
sc.tl.umap(adata_ref, neighbors_key="scVI_neighbors")
adata_ref.obsm[f"X_scVI_umap"] = adata_ref.obsm['X_umap']

# Plot UMAP
sc.pl.embedding(adata_ref, basis="X_scVI_umap", color=["study"], save=batch + "_study.png",
               palette=list(matplotlib.colors.CSS4_COLORS.values()))
sc.pl.embedding(adata_ref, basis="X_scVI_umap", color=["main_original_celltype"], save=batch + "_main_celltype.png")
sc.pl.embedding(adata_ref, basis="X_scVI_umap", color=["author_original_celltype"], save=batch +"_author_celltype.png",
               palette=list(matplotlib.colors.CSS4_COLORS.values()))

# Plot UMAP Embedding Density
adata_ref.obsm[f"X_scvi_umap"] = adata_ref.obsm['X_scVI_umap']
sc.tl.embedding_density(adata_ref, basis='scvi_umap', groupby='study')
sc.pl.embedding_density(adata_ref, basis='scvi_umap', key='scvi_umap_density_study',
                        save=batch + "_density_map_study.png")
del(adata_ref.obsm[f"X_scvi_umap"])



#################################################################################################################################
####################################################### Cell Type Check #########################################################
#################################################################################################################################

subset = adata_ref[adata_ref.obs["main_original_celltype"]!="Unknown"]
subset = subset[subset.obs["main_original_celltype"]!="CD34+ Blasts"]
sc.pl.umap(subset, color="main_original_celltype", palette="tab10", save=batch + "_celltype1.png")

subset.obs["main_original_celltype"] = subset.obs["main_original_celltype"].replace(to_replace="CD4+ T Cells", value="T Cells")
subset.obs["main_original_celltype"] = subset.obs["main_original_celltype"].replace(to_replace="CD8+ T Cells", value="T Cells")

celltypes = ["Monocytes","T Cells","B Cells","HSC/MPPs","Granulocytes","NK Cells","Erythroid Cells"]
s = np.array([s in celltypes for s in subset.obs.main_original_celltype])
subset2 = subset[s].copy()
sc.pl.umap(subset2, color="main_original_celltype", palette='tab10', save=batch  + "_celltype2.png")

subset3 = subset2[subset2.obs["clinical"]=='normal']
sc.pl.umap(subset3, color="main_original_celltype", palette="tab10", save=batch + "_celltype3.png")

sc.pp.neighbors(subset, use_rep="X_scVI", key_added="scvi_neighbors")
sc.tl.umap(subset, neighbors_key="scvi_neighbors")
sc.pl.umap(subset, color="main_original_celltype", palette="tab10", save=batch + "_celltype1b.png")

sc.pp.neighbors(subset2, use_rep="X_scVI", key_added="scvi_neighbors")
sc.tl.umap(subset2, neighbors_key="scvi_neighbors")
sc.pl.umap(subset2, color="main_original_celltype", palette="tab10", save=batch + "_celltype2b.png")

sc.pp.neighbors(subset3, use_rep="X_scVI", key_added="scvi_neighbors")
sc.tl.umap(subset3, neighbors_key="scvi_neighbors")
sc.pl.umap(subset3, color="main_original_celltype", palette="tab10", save=batch + "_celltype3b.png")



#################################################################################################################################
######################################################## scVI Complete ##########################################################
#################################################################################################################################

#save anndata
adata_ref.write_h5ad(results_file)

#save embedding
emb = ad.AnnData(adata_ref.obsm["X_scVI"])
emb.write_h5ad(embedding_file)