#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
################################################ AML ATLAS UNCORRECTED ANALYSIS #################################################
######## Perform Quality Check and Final Filtering, Run Dimensionality Reduction and Clustering, Check for Batch Effects ########
#################################################################################################################################


# Import packages
import os
import sys
import warnings
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import pymde as mde
import torch
import matplotlib
import seaborn as sns
#from MulticoreTSNE import MulticoreTSNE as TSNE
from rich import print


#configure settings
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')
warnings.simplefilter(action="ignore", category=FutureWarning)
sc.set_figure_params(figsize=(4, 4))

# create output directories
if not os.path.exists("outs"):
    os.makedirs("outs")


#set results output files
original_celltypes_file = "outs/filtered_atlasData_original_celltypes.h5ad.gz"
uncorrected_file = 'outs/filtered_atlasData_uncorrected.h5ad'
pca_embedding = 'outs/uncorrected_pca_embedding.h5ad'
umap_embedding = 'outs/uncorrected_umap_embedding.h5ad'
tsne_embedding = 'outs/uncorrected_tsne_embedding.h5ad'



#################################################################################################################################
##################################################### Read Combined Dataset #####################################################
#################################################################################################################################

# Read full dataset
adata = sc.read_h5ad("../../filtered_atlasData_original_celltypes.h5ad.gz")
print(adata.n_obs, adata.n_vars)

# Additional Cell Filtering
sc.pp.filter_cells(adata, min_counts=1000)
sc.pp.filter_cells(adata, min_genes=300)

# Remove Samples with Low Cell Counts
cell_counts = adata.obs['sample'].value_counts()
keep = cell_counts.index[cell_counts >= 50]
adata = adata[adata.obs['sample'].isin(keep)].copy()

# Remove cells with Hb contamination
#remove = pd.read_csv("../contaminated_cells.csv", index_col=0)
#remove = list(remove["0"])

#subset = np.array([s in remove for s in adata.obs.index])
#adata = adata[~subset]

# Remove Zhai 2022
#adata = adata[adata.obs["study"]!="zhai_2022"]

# Remove Wu 2020
#adata = adata[adata.obs["study"]!="wu_2020"]

# Filter genes on remaining data
sc.pp.filter_genes(adata, min_cells=30)
print(adata.n_obs, adata.n_vars)

print(adata.obs["study"].value_counts())
pd.DataFrame(adata.obs["study"].value_counts()).to_csv("outs/study_counts.csv")
print(adata.obs["sample"].value_counts())
pd.DataFrame(adata.obs["sample"].value_counts()).to_csv("outs/sample_counts.csv")

# Harmonise Celltypes
adata.obs["main_original_celltype"] = adata.obs[
    "main_original_celltype"].replace({"Perivascular cell":"Stromal Cells","Megakaryocyte":"Megakaryocytes",
                                       "Plasma cell":"Plasma Cells","Plasmablast":"Plasma Cells",
                                       "cDC1":"Conventional Dendritic Cells","cDC2":"Conventional Dendritic Cells",
                                       "DC precursor":"Dendritic Progenitor Cell","Ery":"Erythroid Cells",
                                       "Dendritic Cells":"Conventional Dendritic Cells","MAIT":"T Cells",
                                       "CD8+ T":"CD8+ T Cells","CD4+ T":"CD4+ T Cells","gd T":"T Cells",
                                       "Pre-B":"B Cell Precursors","Pro-B":"B Cell Precursors",
                                       "Myelocytes":"Granulocytes","Granulocyte":"Granulocytes",
                                       "Promyelocytes":"Granulocytes","HLA-II+ monocyte":"Monocytes",
                                       "HSC":"HSC/MPPs","MPP":"HSC/MPPs","CD14+ monocyte":"Monocytes",
                                       "CD11c+":"Unknown","LymP":"Unknown"
})


#################################################################################################################################
##################################################### Check Quality Control #####################################################
#################################################################################################################################

#annotate mitochondrial/ribosomal/hb genes
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL"))
adata.var['hb'] = adata.var_names.str.contains(("^HB[^(P)]"))

#run qc function
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], 
                           inplace=True, percent_top=[20], log1p=True)

#find proportion of mito genes
mito_genes = adata.var_names.str.startswith('MT-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mt2'] = np.sum(
    adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)

# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1)

#plot qc
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo', 'pct_counts_hb'],
             jitter=0.4, groupby = 'sample', rotation= 45, save="_full_qc.png")
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save="_mt_by_counts.png")
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color="pct_counts_mt", 
              save="_nfeatures_by_counts_mito.png")

# histograms 
fig = sns.displot(adata.obs["total_counts"], bins=100, kde=False).fig
fig.savefig("figures/hist_prefilt_total_counts.png") 

fig = sns.displot(adata.obs["n_genes_by_counts"], bins=100, kde=False).fig
fig.savefig("figures/hist_prefilt_ngenes_counts.png") 

#improve visibility for lower counts
sub = adata[adata.obs["total_counts"]<5000]
fig = sns.displot(sub.obs["total_counts"], bins=100, kde=False).fig
fig.savefig("figures/hist_prefilt_total_counts_lowcounts.png")

fig = sns.displot(sub.obs["n_genes_by_counts"], bins=100, kde=False).fig
fig.savefig("figures/hist_prefilt_ngenes_counts_lowcounts.png")



#################################################################################################################################
############################################## Normalisation and Feature Selection ##############################################
#################################################################################################################################

print("Log Normalisation")
#log normalise - depth 10 000
sc.pp.normalize_total(adata, target_sum=1e4)
# logaritmise
sc.pp.log1p(adata)


#save normalised counts in raw slot.
adata.raw = adata
adata.layers["normalised_counts"] = adata.X

# Save full anndata
adata.write_h5ad("core_atlasData.h5ad")

#get highly variable genes
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=4000,
    flavor="seurat_v3",
    layer="counts",
    batch_key="sample",
    subset=False,
    span=0.8
)



#################################################################################################################################
################################################### Dimensionality Reduction ####################################################
#################################################################################################################################

print("PCA")
#compute pca
sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=True)
sc.pl.pca_variance_ratio(adata, log=True, save="_principal_components.png")

#Save PCA
pca = adata.obsm[f"X_pca"]
pca = ad.AnnData(pca)
pca.write_h5ad(pca_embedding)


# Run UMAP
print("UMAP")
#calculate neighbours graph - save as unique key to retain uncorrected parameters
sc.pp.neighbors(adata, n_pcs = 30, key_added="uncorrected_neighbors")
sc.tl.umap(adata, neighbors_key='uncorrected_neighbors')


sc.pl.umap(adata, color='sample', title='Uncorrected UMAP', palette=list(matplotlib.colors.CSS4_COLORS.values()), 
           save="_atlas_uncorrected_sample.png")
sc.pl.umap(adata, color='study', title='Uncorrected UMAP', save="_atlas_uncorrected_study.png")


# Save uncorrected UMAP
adata.obsm[f"X_uncorrected_umap"] = adata.obsm['X_umap']
umap = adata.obsm[f"X_umap"]
umap = ad.AnnData(umap)
umap.write_h5ad(umap_embedding)

# Run TSNE
print("t-SNE")
sc.tl.tsne(adata, n_pcs=30)
sc.pl.tsne(adata, color='sample', title='Uncorrected t-SNE', palette=list(matplotlib.colors.CSS4_COLORS.values()), 
           save="_atlas_uncorrected_sample.png")
sc.pl.tsne(adata, color='study', title='Uncorrected t-SNE', save="_atlas_uncorrected_study.png")

# Save TSNE
adata.obsm[f"X_uncorrected_tsne"] = adata.obsm['X_tsne']
tsne = adata.obsm[f"X_tsne"]
tsne = ad.AnnData(tsne)
tsne.write_h5ad(tsne_embedding)


#################################################################################################################################
########################################################## Clustering ###########################################################
#################################################################################################################################

sc.tl.leiden(adata, neighbors_key="uncorrected_neighbors", key_added="uncorrected_leiden")
sc.pl.umap(adata, color="uncorrected_leiden", save="_atlas_uncorrected_leiden.png")
sc.pl.tsne(adata, color="uncorrected_leiden", save="_atlas_uncorrected_leiden.png")


#################################################################################################################################
################################################# Uncorrected Analysis Complete #################################################
#################################################################################################################################

print("Saving Uncorrected AnnData")
adata.write_h5ad(uncorrected_file)
