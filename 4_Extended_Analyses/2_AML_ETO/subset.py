#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
###################################################### AML ATLAS: SUBSET ########################################################
###################################### Subset for Samples with RUNX1-RUNX1T1 Translocation ######################################
#################################################################################################################################

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

if not os.path.exists("outs"):
    os.makedirs("outs")

adata = sc.read_h5ad('../annotation/outs/scvi_annotated.h5ad')

print(adata.n_obs, adata.n_vars)


#################################################################################################################################
############################################# Subset for RUNX1-RUNX1T1 Samples Only #############################################
#################################################################################################################################

print(adata.obs["translocations"].value_counts())
adata = adata[adata.obs["translocations"]=="RUNX1-RUNX1T1"]
print(adata.n_obs, adata.n_vars)

sc.pl.umap(adata, color="main_celltype", save="_RUNX1_RUNX1T1_original_main_celltype.png")


#################################################################################################################################
############################################## Re-Compute UMAP from scVI Embedding ##############################################
#################################################################################################################################

# Re-Compute UMAP
sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_scVI", key_added="scvi_neighbors")
sc.tl.umap(adata, neighbors_key="scvi_neighbors")

sc.pl.umap(adata, color="main_celltype", save="_RUNX1_RUNX1T1_main_celltype.png")
sc.pl.umap(adata, color="main_age_group", save="_RUNX1_RUNX1T1_main_age.png")
sc.pl.umap(adata, color="fine_age_group", save="_RUNX1_RUNX1T1_fine_age.png")


# Save Subset
adata.write_h5ad("outs/RUNX1-RUNX1T1.h5ad")


#################################################################################################################################
######################################################## Gene Filtering #########################################################
#################################################################################################################################

# Select Genes
min_cell_thresh = max(20,round(0.0001*len(adata.obs)))
print(str(min_cell_thresh))
sc.pp.filter_genes(adata, min_cells=min_cell_thresh)
print(adata)
genes = pd.DataFrame(adata.var.index)
genes.to_csv("outs/genes_used.csv")


# Write output files
print(adata.n_obs, adata.n_vars)
adata.write_h5ad("outs/filtered_adata.h5ad")

row_attrs = { 
    "Gene": np.array(adata.var.index) ,
}
col_attrs = { 
    "CellID":  np.array(adata.obs.index) ,
    "nGene": np.array( np.sum(adata.layers["counts"].transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.layers["counts"].transpose() , axis=0)).flatten() ,
}
lp.create("outs/gene_filtered_raw.loom", 
          adata.layers["counts"].transpose(), row_attrs, col_attrs )


#################################################################################################################################
############################################## Downsample to 4000 cells per sample ##############################################
#################################################################################################################################

print(adata.obs["sample"].value_counts())

# To optimise the runtime of pyscenic, downsample to 4000 cells per sample
target_cells = 4000

print("Downsampling to " + str(target_cells) + " per group")

adatas = [adata[adata.obs["sample"].isin([group])] for group in adata.obs["sample"].cat.categories]

for dat in adatas:
    if dat.n_obs > target_cells:
         sc.pp.subsample(dat, n_obs=target_cells)

adata_downsampled = adatas[0].concatenate(*adatas[1:])

print(adata_downsampled.n_obs, adata_downsampled.n_vars)
print(adata_downsampled.obs["sample"].value_counts())


# Write output files
print(adata_downsampled.n_obs, adata_downsampled.n_vars)
adata_downsampled.write_h5ad("outs/downsampled.h5ad")

row_attrs = {
    "Gene": np.array(adata_downsampled.var.index) ,
}
col_attrs = {
    "CellID":  np.array(adata_downsampled.obs.index) ,
    "nGene": np.array( np.sum(adata_downsampled.layers["counts"].transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata_downsampled.layers["counts"].transpose() , axis=0)).flatten() ,
}
lp.create("outs/downsampled_raw.loom",
          adata_downsampled.layers["counts"].transpose(), row_attrs, col_attrs )