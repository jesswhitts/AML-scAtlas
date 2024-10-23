#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
################################################## AML ATLAS CLUSTER MARKERS ####################################################
################################################ Uncorrected: Wilcoxon Rank Sum  ################################################
#################################################################################################################################

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import warnings

# Configure scanpy settings
sc.settings.verbosity = 0   ### suppress output of scanpy
warnings.filterwarnings('ignore') # suppress Warnings
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Set results output file
markers_file = 'outs/scvi_markers.h5ad'

# Read in anndata object
adata = sc.read_h5ad('../../clustering/outs/scvi_clustered.h5ad')
adata.uns['log1p']["base"] = None #temporary workaround to a known bug git #2181


#################################################################################################################################
################################################## Run Wilcoxon Rank Sum Test ###################################################
#################################################################################################################################

sc.tl.rank_genes_groups(adata, groupby='cluster', method='wilcoxon', key_added='wilcoxon_rank')

# Saving in case subsequent steps fail
adata.write_h5ad(markers_file)


#################################################################################################################################
######################################################### Save Results ##########################################################
#################################################################################################################################

gene_rank = sc.get.rank_genes_groups_df(adata, group=None, key='wilcoxon_rank')
gene_rank.to_csv("outs/gene_rank.csv")

pd.DataFrame(adata.uns['wilcoxon_rank']['names']).head(25)

result = adata.uns['wilcoxon_rank']
groups = result['names'].dtype.names
markers = pd.DataFrame(
          {group + '_' + key[:1]: result[key][group]
          for group in groups for key in ['names', 'pvals_adj']}).head(25)

markers.to_csv("outs/cluster_markers.csv")

adata.write_h5ad(markers_file)
