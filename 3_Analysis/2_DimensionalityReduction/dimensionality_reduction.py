#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
############################################## AML ATLAS DIMENSIONALITY REDUCTION ###############################################
###################################################### scVI Corrected Atlas #####################################################
#################################################################################################################################


import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import anndata as ad

#configure scanpy settings
sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc._settings.ScanpyConfig.n_jobs:10
sc.settings.set_figure_params(dpi=80, facecolor='white')

# create output directories
if not os.path.exists("outs"):
    os.makedirs("outs")

# Set output files
results_file = 'outs/scvi_dimsred.h5ad'
pca_embedding = 'outs/uncorrected_pca_embedding.h5ad'
umap_embedding = 'outs/uncorrected_umap_embedding.h5ad'
tsne_embedding = 'outs/uncorrected_tsne_embedding.h5ad'

# Read in full anndata object
adata = sc.read_h5ad('../core_atlasData.h5ad')

# Add embeddings and hvgs
embedding = sc.read_h5ad('../batch-correction/outs/scvi_sample_corrected.h5ad')
adata.obsm['X_scVI'] = embedding.obsm['X_scVI']
adata.obsm['X_pca'] = embedding.obsm['X_pca']
adata.obsm['X_uncorrected_tsne'] = embedding.obsm['X_uncorrected_tsne']
adata.obsm['X_uncorrected_umap'] = embedding.obsm['X_uncorrected_umap']
adata.var["highly_variable"] = False
hvgs = embedding.var.index
print(len(hvgs))
adata.var.loc[hvgs, 'highly_variable'] = True

print(adata.n_obs, adata.n_vars)



#################################################################################################################################
####################################################### Cell Cycle Scoring ######################################################
#################################################################################################################################


#open cell cycle list (obtained from scanpy github) and split into two lists
cell_cycle_genes = [x.strip() for x in open('/data/stemcell/jwhittle/ref/regev_lab_cell_cycle_genes.txt')]
print(len(cell_cycle_genes))

# Split into 2 lists
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]

#print number of genes in each list
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
print(len(cell_cycle_genes))


# Counts need to be normalised/logarithmised and scaled to run this function - be sure to save normalised/logarithmised counts in the raw slot for future analysis

#calculate cell cycle scores
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

sc.pl.violin(adata, ['S_score', 'G2M_score'],
             jitter=0.4, rotation=45, save="_cell_cycle_scores.png")

#save S/G2M difference in case needed
adata.obs["S_G2M_difference"] = adata.obs["S_score"] - adata.obs["G2M_score"]


# save scaled counts - put normalised counts back
sc.pp.scale(adata)
ZST = adata.X
emb = ad.AnnData(ZST)
emb.write_h5ad("outs/zst_scaled_counts.h5ad")


#################################################################################################################################
######################################################## Visualisatiom QC #######################################################
#################################################################################################################################


sc.pp.neighbors(adata, use_rep="X_scVI", key_added="scvi_neighbors")
sc.tl.umap(adata, neighbors_key="scvi_neighbors")
adata.obsm["X_scVI_umap"] = adata.obsm["X_umap"]

#check for potential cell cycle bias
sc.pl.embedding(adata, basis="X_scVI", color='phase', save="_cellcycle_scvi.png")
sc.pl.embedding(adata, basis="X_pca", color='phase', save="_cellcycle_Xpca.png")

#UMAPs by sample
sc.pl.embedding(adata, basis="X_scVI_umap", color='sample', palette=list(matplotlib.colors.CSS4_COLORS.values()), 
           title="scVI Corrected UMAP", save="_scvi_sample.png")
sc.pl.embedding(adata, basis="X_uncorrected_umap", color='sample', palette=list(matplotlib.colors.CSS4_COLORS.values()), 
           title="Uncorrected UMAP", save="_uncorrected_sample.png")

#UMAPs by study
sc.pl.embedding(adata, basis="X_scVI_umap", color='study', 
           title="scVI Corrected UMAP", save="_scvi_study.png")
sc.pl.embedding(adata, basis="X_uncorrected_umap", color='study',  
           title="Uncorrected UMAP", save="_uncorrected_study.png")

#UMAPs by celltype
sc.pl.embedding(adata, basis="X_scVI_umap", color=("celltype"),
           title="scVI Corrected UMAP", save="_scvi_celltype.png")
sc.pl.embedding(adata, basis="X_uncorrected_umap", color=("celltype"),
           title="Uncorrected UMAP", save="_uncorrected_celltype.png")

#UMAP QC
sc.pl.embedding(adata, basis="X_scVI_umap", color=["phase", "pct_counts_mt", "pct_counts_ribo",
                                                     "total_counts"],save="_scvi_quality_control.png")
sc.pl.embedding(adata, basis="X_uncorrected_umap", color=["phase", "pct_counts_mt", "pct_counts_ribo", 
                                                          "total_counts"],save="_uncorrected_quality_control.png")



#################################################################################################################################
#################################################### Visualise Full Dataset #####################################################
#################################################################################################################################

#plot predefined 'housekeeping' genes of interest on umap
sc.pl.embedding(adata, basis = "X_uncorrected_umap", color=['IL7R', #CD4 T cells
                         'CD79A', 'MS4A1', # B cells
                         'CD8A', # CD8 T cells
                         'LYZ', 'CD14', #CD14 monos
                         'S100A8', # monos/neuts
                         'GNLY', 'KLRB1', # NK cells
                         'FCGR3A', #FCGR3A monos
                         'FCER1A', 'CST3', #DCs
                         'PPBP', 'PF4', #Mks
                         'HBB', 'HBA2', #Erys
                         'CD34', 'RUNX1', 'MEIS1','SIX1'],#other
          save="_uncorrected_markers.png")

sc.pl.embedding(adata, basis="X_scVI_umap", color=['IL7R', #CD4 T cells
                         'CD79A', 'MS4A1', # B cells
                         'CD8A', # CD8 T cells
                         'LYZ', 'CD14', #CD14 monos
                         'S100A8', # monos/neuts
                         'GNLY', 'KLRB1', # NK cells
                         'FCGR3A', #FCGR3A monos
                         'FCER1A', 'CST3', #DCs
                         'PPBP', 'PF4', #Mks
                         'HBB', 'HBA2', #Erys
                         'CD34', 'RUNX1', 'MEIS1','SIX1'],#other
          save="_scvi_markers.png")


#UMAPs by marker gene
#IL7R
sc.pl.embedding(adata, basis="X_scVI_umap", color='IL7R', save="_scvi_IL7R.png")
sc.pl.embedding(adata, basis="X_uncorrected_umap", color='IL7R',  save="_uncorrected_IL7R.png")

#CD79A
sc.pl.embedding(adata, basis="X_scVI_umap", color='CD79A', save="_scvi_CD79A.png")
sc.pl.embedding(adata, basis="X_uncorrected_umap", color='CD79A', save="_uncorrected_CD79A.png")

#CD14
sc.pl.embedding(adata, basis="X_scVI_umap", color='CD14', save="_scvi_CD14.png")
sc.pl.embedding(adata, basis="X_uncorrected_umap", color='CD14', save="_uncorrected_CD14.png")

#HBB
sc.pl.embedding(adata, basis="X_scVI_umap", color='HBB', save="_scvi_HBB.png")
sc.pl.embedding(adata, basis="X_uncorrected_umap", color='HBB', save="_uncorrected_HBB.png")

#CD34
sc.pl.embedding(adata, basis="X_scVI_umap", color='CD34', save="_scvi_CD34.png")
sc.pl.embedding(adata, basis="X_uncorrected_umap", color='CD34', save="_uncorrected_CD34.png")


#UMAPs by clinical details
sc.pl.embedding(adata, basis="X_scVI_umap", color='clinical', 
           title="scVI Corrected UMAP", save="_scvi_clinical.png")
sc.pl.embedding(adata, basis="X_uncorrected_umap", color='clinical',  
           title="Uncorrected UMAP", save="_uncorrected_clinical.png")



#################################################################################################################################
################################################## Visualise Key Data Subsets ###################################################
#################################################################################################################################

#split normal
normal = np.array([s in ["normal"] for s in adata.obs.clinical])

adata_aml = adata[~normal].copy()
adata_norm = adata[normal].copy()


# Visualise normal samples only
sc.pl.embedding(adata_norm, basis="X_uncorrected_umap", color=['IL7R', #CD4 T cells
                         'CD79A', 'MS4A1', # B cells
                         'CD8A', # CD8 T cells
                         'LYZ', 'CD14', #CD14 monos
                         'S100A8', # monos/neuts
                         'GNLY', 'KLRB1', # NK cells
                         'FCGR3A', #FCGR3A monos
                         'FCER1A', 'CST3', #DCs
                         'PPBP', 'PF4', #Mks
                         'HBB', 'HBA2', #Erys
                         'CD34', 'RUNX1', 'MEIS1','SIX1'],#other
          save="_uncorrected_norm_markers.png")

sc.pl.embedding(adata_norm, basis="X_scVI_umap", color=['IL7R', #CD4 T cells
                         'CD79A', 'MS4A1', # B cells
                         'CD8A', # CD8 T cells
                         'LYZ', 'CD14', #CD14 monos
                         'S100A8', # monos/neuts
                         'GNLY', 'KLRB1', # NK cells
                         'FCGR3A', #FCGR3A monos
                         'FCER1A', 'CST3', #DCs
                         'PPBP', 'PF4', #Mks
                         'HBB', 'HBA2', #Erys
                         'CD34', 'RUNX1', 'MEIS1','SIX1'],#other
          save="_scvi_norm_markers.png")



# Visualise AML samples only
#plot predefined 'housekeeping' genes of interest on umap
sc.pl.embedding(adata_aml, basis="X_scVI_umap", color=['IL7R', #CD4 T cells
                         'CD79A', 'MS4A1', # B cells
                         'CD8A', # CD8 T cells
                         'LYZ', 'CD14', #CD14 monos
                         'S100A8', # monos/neuts
                         'GNLY', 'KLRB1', # NK cells
                         'FCGR3A', #FCGR3A monos
                         'FCER1A', 'CST3', #DCs
                         'PPBP', 'PF4', #Mks
                         'HBB', 'HBA2', #Erys
                         'CD34', 'RUNX1', 'MEIS1','SIX1'],#other
          save="_scvi_aml_markers.png")

sc.pl.embedding(adata_aml, basis="X_uncorrected_umap", color=['IL7R', #CD4 T cells
                         'CD79A', 'MS4A1', # B cells
                         'CD8A', # CD8 T cells
                         'LYZ', 'CD14', #CD14 monos
                         'S100A8', # monos/neuts
                         'GNLY', 'KLRB1', # NK cells
                         'FCGR3A', #FCGR3A monos
                         'FCER1A', 'CST3', #DCs
                         'PPBP', 'PF4', #Mks
                         'HBB', 'HBA2', #Erys
                         'CD34', 'RUNX1', 'MEIS1','SIX1'],#other
          save="_uncorrected_aml_markers.png")


#UMAPs by clinical
sc.pl.embedding(adata_aml, basis="X_scVI_umap", color='clinical', 
           title="scVI Corrected UMAP", save="_scvi_aml_clinical.png")
sc.pl.embedding(adata_aml, basis="X_uncorrected_umap", color='clinical',  
           title="Uncorrected UMAP", save="_uncorrected_aml_clinical.png")


del(adata_norm)
del(adata_aml)



# Visualise Celltypes

# Tidy Up Annotations
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

# Get Annotated Cells Only
adata_subset = adata[adata.obs["main_original_celltype"]!="Unknown"]

# Plot on UMAP
sc.pl.embedding(adata_subset, basis="X_scVI_umap", color=("main_original_celltype"),
           title="scVI Corrected UMAP", save="_scvi_subset_celltype.png")
sc.pl.embedding(adata_subset, basis="X_uncorrected_umap", color=("main_original_celltype"),
           title="Uncorrected UMAP", save="_uncorrected_subset_celltype.png")

# Plot on New UMAP
sc.pp.neighbors(adata_subset, use_rep="X_scVI", key_added="scvi_neighbors")
sc.tl.umap(adata_subset, neighbors_key="scvi_neighbors")
adata_subset.obsm["X_scVI_umap"] = adata_subset.obsm["X_umap"]

sc.pl.embedding(adata_subset, basis="X_scVI_umap", color=("main_original_celltype"),
           title="scVI Corrected UMAP", save="_scvi_subset_celltype2.png")


#################################################################################################################################
############################################### Dimensinoality Reduction Complete ###############################################
#################################################################################################################################

adata.write_h5ad(results_file)




