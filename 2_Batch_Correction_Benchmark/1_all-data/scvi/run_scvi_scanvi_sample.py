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



#################################################################################################################################
####################################################### Prepare scANVI Model ####################################################
#################################################################################################################################

#set results output files
results_file = 'outs/scanvi_sample_corrected.h5ad'
embedding_file = 'outs/scanvi_sample_embedding.h5ad'
scvi_model = "scvi_sample_model/"
batch = "scanvi_sample"

# Harmonise input celltypes
adata_ref.obs["celltype"] = adata_ref.obs["celltype"].replace({"Perivascular cell":"Stromal Cells","Megakaryocyte":"Megakaryocytes",
                                                       "Plasma cell":"Plasma Cells","Plasmablast":"Plasma Cells",
                                                       "cDC1":"Conventional Dendritic Cells",
                                                       "cDC2":"Conventional Dendritic Cells",
                                                       "DC precursor":"Dendritic Progenitor Cell","Ery":"Erythroid Cells",
                                                       "Dendritic Cells":"Conventional Dendritic Cells","MAIT":"T Cells",
                                                       "CD8+ T":"CD8+ T Cells","CD4+ T":"CD4+ T Cells","gd T":"T Cells",
                                                       "Pre-B":"B Cell Precursors","Pro-B":"B Cell Precursors",
                                                       "Myelocytes":"Granulocytes","Granulocyte":"Granulocytes",
                                                       "Promyelocytes":"Granulocytes","HLA-II+ monocyte":"Monocytes",
                                                       "HSC":"HSC/MPPs","MPP":"HSC/MPPs","CD14+ monocyte":"Monocytes",
                                                       "CD11c+":"Unknown","LymP":"Unknown"
}) 

# Only input normal celltypes into model
subset = adata_ref.obs[adata_ref.obs["clinical"]!="normal"]
s = subset.index
adata_ref.obs.loc[s, 'celltype'] = "Unknown"

vae = scvi.model.SCVI.load(scvi_model, adata_ref)


#################################################################################################################################
######################################################## Run scANVI Model #######################################################
#################################################################################################################################

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata_ref,
    labels_key="celltype",
    unlabeled_category="Unknown",
)

lvae.train(max_epochs=20, n_samples_per_label=100)

adata_ref.obsm["X_scANVI"] = lvae.get_latent_representation(adata_ref)

# Infer celltypes from scANVI
adata_ref.obs["scANVI_annotations"] = lvae.predict(adata_ref)

# Save model
model_dir = os.path.join("./" + batch + "_model")
lvae.save(model_dir, overwrite=True)


#################################################################################################################################
###################################################### scANVI Visualisation #####################################################
#################################################################################################################################

# Run UMAP
sc.pp.neighbors(adata_ref, use_rep="X_scANVI", key_added="scANVI_neighbors")
sc.tl.umap(adata_ref, neighbors_key="scANVI_neighbors")
adata_ref.obsm[f"X_scANVI_umap"] = adata_ref.obsm['X_umap']

# Plot UMAP
sc.pl.embedding(adata_ref, basis="X_scANVI_umap", color=["study"], save=batch+"_study.png",
               palette=list(matplotlib.colors.CSS4_COLORS.values()))
sc.pl.embedding(adata_ref, basis="X_scANVI_umap", color=["main_original_celltype"], save=batch+"_main_celltype.png")
sc.pl.embedding(adata_ref, basis="X_scANVI_umap", color=["author_original_celltype"], save=batch+"_author_celltype.png",
               palette=list(matplotlib.colors.CSS4_COLORS.values()))

# Plot UMAP Embedding Density
adata_ref.obsm[f"X_scanvi_umap"] = adata_ref.obsm['X_scANVI_umap']
sc.tl.embedding_density(adata_ref, basis='scanvi_umap', groupby='study')
sc.pl.embedding_density(adata_ref, basis='scanvi_umap', key='scanvi_umap_density_study',
                        save=batch+"_density_map_study.png")
del(adata_ref.obsm[f"X_scanvi_umap"])



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

sc.pp.neighbors(subset, use_rep="X_scANVI", key_added="scanvi_neighbors")
sc.tl.umap(subset, neighbors_key="scanvi_neighbors")
sc.pl.umap(subset, color="main_original_celltype", palette="tab10", save=batch + "_celltype1b.png")

sc.pp.neighbors(subset2, use_rep="X_scANVI", key_added="scanvi_neighbors")
sc.tl.umap(subset2, neighbors_key="scanvi_neighbors")
sc.pl.umap(subset2, color="main_original_celltype", palette="tab10", save=batch + "_celltype2b.png")

sc.pp.neighbors(subset3, use_rep="X_scANVI", key_added="scanvi_neighbors")
sc.tl.umap(subset3, neighbors_key="scanvi_neighbors")
sc.pl.umap(subset3, color="main_original_celltype", palette="tab10", save=batch + "_celltype3b.png")



#################################################################################################################################
######################################################## scANVI Complete ########################################################
#################################################################################################################################

adata_ref.write_h5ad(results_file)


emb = ad.AnnData(adata_ref.obsm["X_scANVI"])
emb.write_h5ad(embedding_file)
