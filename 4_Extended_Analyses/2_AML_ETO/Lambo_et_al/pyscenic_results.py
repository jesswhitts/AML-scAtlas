#!/usr/bin/env python
# coding: utf-8


#################################################################################################################################
################################################## PYSCENIC: REGULON VALIDATION #################################################
#################################################################################################################################

import os
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial import distance
from scipy.cluster import hierarchy
import scipy.stats as stats
from adjustText import adjust_text #version 0.8 - newer ones have error
import gseapy as gp
from gseapy import GSEA
import pathlib

if not os.path.exists("plots"):
    os.makedirs("plots")

sc.settings.figdir = "plots"
sc.set_figure_params(scanpy=True, fontsize=10, dpi=180, dpi_save=500, frameon=False, format='png')
plt.rcParams.update({'font.weight': 'bold'})

lambo_auc = '../AML-ETO/pyscenic/gene_filtered/validation/lambo/outs/auc_mtx_validation_data.csv'
lambo_anno = '../AML-ETO/pyscenic/gene_filtered/validation/lambo/outs/cellAnnot_validation_data.csv'

def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    #for idx, (name, color) in enumerate(zip(names, colors)):
    #    ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f

save_name = 'aml_eto'
plots = 'plots/aml_eto/validation/'
if not os.path.exists("plots/aml_eto/validation"):
    os.makedirs("plots/aml_eto/validation")
sc.settings.figdir = "plots/aml_eto/validation"

## Lambo et al Data
auc_mtx = pd.read_csv(lambo_auc, index_col=0)
cellAnnot = pd.read_csv(lambo_anno, index_col=0)

# Select AML samples only
cellAnnot = cellAnnot[cellAnnot['sample'].str.contains('AML')==True]
# View Cell Types
print(cellAnnot.main_celltype.value_counts())

# Update Age Groups
cellAnnot["age_group"] = cellAnnot["age_group"].replace(to_replace="0-9", value="Pediatric")
cellAnnot["age_group"] = cellAnnot["age_group"].replace(to_replace="10-19", value="Adolescent")
cellAnnot["age_group"] = cellAnnot["age_group"].replace(to_replace="20-29", value="Adult")
cellAnnot["age_group"] = cellAnnot["age_group"].replace(to_replace="30-49", value="Adult")
cellAnnot["age_group"] = cellAnnot["age_group"].replace(to_replace="50-69", value="Adult")
cellAnnot["age_group"] = cellAnnot["age_group"].replace(to_replace="80+", value="Adult")
cellAnnot["age_group"] = cellAnnot["age_group"].replace(to_replace="unknown", value="Unknown")

# HSPCs only
progenitor_cells = ['HSC','Progenitor']
subset = np.array([s in progenitor_cells for s in cellAnnot.main_celltype])
cellAnnot = cellAnnot[subset]

cells = cellAnnot.index
auc_mtx = auc_mtx.T[cells].T

cm = {
    'HSC':'#fd8d3c',
    'Progenitor':'#843c39',
}
cm = dict(sorted(cm.items()))
    
lambo_age = {
    'Adolescent':'#029e73ff',
    'Pediatric':'#0173b2ff'
}

# Create Z Score matrix to enable comparison between regulons
auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
for col in list(auc_mtx.columns):
    auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)
    
# Create heatmap
samp = sorted(list(set(cellAnnot['sample'])))
colors = sns.color_palette('tab20',n_colors=len(samp) )
colorsd = dict( zip( samp, colors ))
colormap = [ colorsd[x] for x in cellAnnot['sample'] ]

sns.set()
sns.set(font_scale=0.8)
fig = palplot( colors, samp, size=1.0)
plt.savefig(plots + "lambo_heatmap_HSPC_legend_1_" + save_name + ".png", dpi=600, bbox_inches = "tight")

ct = sorted(list(set(cellAnnot['clinical_subtype'])))
#ct_colours = []
#for key in cm:
#    ct_colours.append(cm[key])
colors = sns.color_palette('deep',n_colors=len(ct) )
colorsd = dict( zip( ct, colors ))
colormap2 = [ colorsd[x] for x in cellAnnot['clinical_subtype'] ]

sns.set()
sns.set(font_scale=0.8)
fig = palplot( colors, ct, size=1.0)
plt.savefig(plots + "lambo_heatmap_HSPC_legend_2_" + save_name + ".png", dpi=600, bbox_inches = "tight")

cats = sorted(list(set(cellAnnot['age_group'])))
age_colours = []
for key in lambo_age:
    age_colours.append(lambo_age[key])
colors = sns.color_palette(age_colours,n_colors=len(cats) )
colorsd = dict( zip( cats, colors ))
colormap3 = [ colorsd[x] for x in cellAnnot['age_group'] ]

sns.set()
sns.set(font_scale=0.8)
fig = palplot( colors, cats, size=1.0)
plt.savefig(plots + "lambo_heatmap_HSPC_legend_3_" + save_name + ".png", dpi=600, bbox_inches = "tight")

sns.set(font_scale=1.4)
plt.rcParams.update({'font.weight': 'bold'})
g = sns.clustermap(auc_mtx_Z[key_regulons], annot=False,  square=False,  linecolor='gray',
    yticklabels=False, xticklabels=True, vmin=-5, vmax=5, row_colors=[colormap,colormap3],
    cmap='seismic', figsize=(18,16) )
g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')
plt.savefig(plots + "lambo_heatmap_HSPC" + save_name + ".png", dpi=600, bbox_inches = "tight")