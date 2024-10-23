#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
################################################## AML ATLAS BATCH CORRECTION ###################################################
####################################### Batch Correction Benchmarking: Combine Embeddings #######################################
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

if not os.path.exists("outs"):
    os.makedirs("outs")


#################################################################################################################################
################################################ Combine AnnData and Embeddings #################################################
#################################################################################################################################

# Load Anndata
adata = sc.read_h5ad("../outs/filtered_atlasData_uncorrected.h5ad")

# Harmonise input celltypes
adata.obs["celltype"] = adata.obs["celltype"].replace({"Perivascular cell":"Stromal Cells","Megakaryocyte":"Megakaryocytes",
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

# Normal celltypes as label key
adata.obs["normal_celltype"] = "celltype"
subset = adata.obs[adata.obs["clinical"]!="normal"]
s = subset.index
adata.obs.loc[s, 'normal_celltype'] = "Unknown"


# Load Embeddings

#Harmony
harmony = sc.read_h5ad("../harmony/outs/harmony_embedding.h5ad")

#scVI Study
scvi_study = sc.read_h5ad("../scvi/outs/scvi_study_embedding.h5ad")

#scVI Sample
scvi_sample = sc.read_h5ad("../scvi/outs/scvi_sample_embedding.h5ad")

#scANVI Study
scanvi_study = sc.read_h5ad("../scvi/outs/scanvi_study_embedding.h5ad")

#scANVI Sample
scanvi_sample = sc.read_h5ad("../scvi/outs/scanvi_sample_embedding.h5ad")

print("Embeddings imported")


# Add Embeddings to Anndata
adata.obsm["X_uncorrected"] = adata.obsm["X_pca"]
adata.obsm["X_harmony"] = harmony.X
adata.obsm["X_scvi_study"] = scvi_study.X
adata.obsm["X_scvi_sample"] = scvi_sample.X
adata.obsm["X_scanvi_study"] = scanvi_study.X
adata.obsm["X_scanvi_sample"] = scanvi_sample.X

# Save Anndata
adata.write_h5ad("outs/combined_embeddings.h5ad")

print("AnnData Saved")
