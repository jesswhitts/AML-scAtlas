#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
################################################ AML ATLAS SAMPLE SEX ANNOTATION ################################################
#################################################### Infer Sample XIST/ChrY #####################################################
#################################################################################################################################

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
sc.set_figure_params(scanpy=True, fontsize=10, dpi=180, dpi_save=500, frameon=False, format='png')
plt.rcParams.update({'font.weight': 'bold'})


#################################################################################################################################
######################################################## Prepare Dataset ########################################################
#################################################################################################################################

adata = sc.read_h5ad('../annotation/scvi_annotated_updated.h5ad')
adata.X = adata.layers['counts']

# Using biomart to identify relevant genes
annot = sc.queries.biomart_annotations(
        "hsapiens",
        ["ensembl_gene_id", "external_gene_name", "start_position", "end_position", "chromosome_name"],
    ).set_index("external_gene_name")


#################################################################################################################################
###################################################### Annotate XIST/ChrY #######################################################
#################################################################################################################################

chrY_genes = adata.var_names.intersection(annot.index[annot.chromosome_name == "Y"])
adata.obs['percent_chrY'] = np.sum(
    adata[:, chrY_genes].X, axis=1) / np.sum(adata.X, axis=1) * 100
adata.obs["XIST-counts"] = adata.X[:,adata.var_names.str.match('XIST')].toarray()

sc.pl.scatter(adata, x='XIST-counts', y='percent_chrY', color="cluster", save='full_XIST_chrY_scatter.png')

sc.pl.violin(adata, "XIST-counts", jitter=0.4, groupby = 'cluster', rotation= 45, save='full_XIST_violin.png')

sc.pl.violin(adata,  "percent_chrY", jitter=0.4, groupby = 'cluster', rotation= 45, save='full_chrY_violin.png')

# Select Relevant Samples
samples = ['AML11D', 'AML12D', 'AML14D', 'AML20D']

subset = np.array([s in samples for s in adata.obs.sample])
subset = adata[subset].copy()

sc.pl.scatter(subset, x='XIST-counts', y='percent_chrY', color="sample", save='selected_XIST_chrY_scatter.png')

sc.pl.violin(subset, "XIST-counts", jitter=0.4, groupby = 'sample', rotation= 45, save='selected_XIST_violin.png')

sc.pl.violin(subset,  "percent_chrY", jitter=0.4, groupby = 'sample', rotation= 45, save='selected_chrY_violin.png')
