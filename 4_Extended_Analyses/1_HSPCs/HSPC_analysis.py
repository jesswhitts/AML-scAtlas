#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
################################################### AML ATLAS: HSPC Analysis ####################################################
#################################################################################################################################

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation

sc.set_figure_params(scanpy=True, fontsize=10, dpi=180, dpi_save=500, frameon=False, format='png')
plot_dir = 'figures/'
plt.rcParams.update({'font.weight': 'bold'})
if not os.path.exists("outs"):
    os.makedirs("outs")
if not os.path.exists("figures"):
    os.makedirs("figures")

save_name = "HSPC"
adata = sc.read_h5ad('../annotation/scvi_annotated.h5ad')
adata.obs['pruned.labels'] = adata.obs['vgalen_pruned_labels']


#################################################################################################################################
############################################## Check CD34 Expression Across Clusters ############################################
#################################################################################################################################

print(adata.n_obs, adata.n_vars)

# Calculate Expression Per Cluster
obs_df = sc.get.obs_df(adata, keys=['cluster','CD34'], use_raw=True)
grouped = obs_df.groupby("cluster", observed=True)
mean, var = grouped.mean(), grouped.var()

celltypes = sc.get.obs_df(adata, keys=['cluster','main_celltype'], use_raw=True) #layer='normalised_corrected_counts')
celltypes = celltypes.drop_duplicates(subset='cluster')
celltypes = celltypes.set_index('cluster')

df = mean.join(celltypes)
df.to_csv("outs/cluster_CD34_celltype.csv")

# Select Clusters Above Threshold
M = df['CD34']
nmads = 1.5
CD34_high = (M < np.median(M) - nmads * median_abs_deviation(M)) | (np.median(M) + nmads * median_abs_deviation(M) < M)
select = df[CD34_high==True]

clusters = list(select.index)
subset = np.array([s in clusters for s in adata.obs.cluster])
adata = adata[subset].copy()

# Visualise Selected Clusters
sc.pl.umap(adata, color="main_celltype", save="_" + save_name + "_original_main_celltype.png")

print(adata.n_obs, adata.n_vars)


#################################################################################################################################
############################################## Re-Compute UMAP from scVI Embedding ##############################################
#################################################################################################################################

# Re-Compute UMAP
sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_scVI", key_added="scvi_neighbors")
sc.tl.umap(adata, neighbors_key="scvi_neighbors")

sc.pl.umap(adata, color="main_celltype", save="_" + save_name + "_main_celltype.png")
sc.pl.umap(adata, color="main_age_group", save="_" + save_name + "_main_age.png")
sc.pl.umap(adata, color="fine_age_group", save="_" + save_name + "_fine_age.png")
sc.pl.umap(adata, color="pruned.labels", save="_" + save_name + "_hspc_celltype.png")

# Save Subset
adata.write_h5ad("outs/" + save_name + ".h5ad")


#################################################################################################################################
################################################### Calculate LSC17 and LSC6 ####################################################
#################################################################################################################################

# Some LSC signature scores from literature - LSC17 for adults, LSC6 is the LSC17 adapted for paediatric AML
LSC17_genes = ['DNMT3B', 'ZBTB46', 'NYNRIN', 'ARHGAP22', 'LAPTM4B', 'MMRN1', 'DPYSL3', 'KIAA0125', 'CDK6','CPXM1', 
               'SOCS2', 'SMIM24', 'EMP1', 'NGFRAP1', 'CD34', 'AKR1C3', 'GPR56']
def LSC17(x):
    return x[0]*0.0874+x[1]*-0.0347+x[2]*0.00865+x[3]*-0.0138+x[4]*0.00582+x[5]*0.0258+x[6]\
           *0.0284+x[7]*0.0196+x[8]*-0.0704+x[9]*-0.0258+x[10]*0.0271+x[11]*-0.0226+x[12]\
           *0.0146+x[13]*0.0465+x[14]*0.0338+x[15]*-0.0402+x[16]*0.0501

LSC6_genes = ["DNMT3B","CD34","ADGRG1","SOCS2","SPINK2","FAM30A"]
def LSC6(x):
    return x[0]*0.0171+x[1]*0.109+x[2]*0.141+x[3]*0.0516+x[4]*0.054+x[5]*0.189

# LSC17
LSC17_adata = adata[:, LSC17_genes]
LSC17_genes = LSC17_adata.var.index
LSC17_cells = LSC17_adata.obs.index
LSC17_counts = pd.DataFrame(LSC17_adata.X, columns=LSC17_genes, index=LSC17_cells)

# LSC6
LSC6_adata = adata[:, LSC6_genes]
LSC6_genes = LSC6_adata.var.index
LSC6_cells = LSC6_adata.obs.index
LSC6_counts = pd.DataFrame(LSC6_adata.X, columns=LSC6_genes, index=LSC6_cells)

adata.obs['LSC17'] = "NA"
adata.obs['LSC6'] = "NA"
cells = adata.obs.index

for c in cells:
    df = LSC17_counts.loc[c]
    score = LSC17(df)
    adata.obs.loc[c, 'LSC17'] = score
    df = LSC6_counts.loc[c]
    score = LSC6(df)
    adata.obs.loc[c, 'LSC6'] = score
    
adata.obs['LSC17'] = adata.obs['LSC17'].astype(float)
adata.obs['LSC6'] = adata.obs['LSC6'].astype(float)

sc.pl.umap(adata, color=['LSC17', 'LSC6'], save='LSC_scores_umaps.png')
plt.rc("figure", figsize=(3,4))
sc.pl.violin(adata, keys=['LSC17', 'LSC6'], rotation=90, save='LSC_scores_violin.png')

plt.rc("figure", figsize=(16,4))
sc.pl.violin(adata,
             keys='LSC17',
             groupby="subset_cluster",
             inner="box",
             save='LSC17_cluster.png'
             )

plt.rc("figure", figsize=(16,4))
sc.pl.violin(adata,
             keys=['LSC6'],
             groupby="subset_cluster",
             inner="box",
             save='LSC6_cluster.png'
             )


#################################################################################################################################
#################################################### Correlate with Annotations #################################################
#################################################################################################################################

cm = {
    "T":"#1f77b4",
    "NK":"#4c31ad",
    "Plasma":"#0818A8",
    "CTL":"#052b6b",
    "Prog":"#01153e",
    "Prog-like":"#aec7e8",
    "HSC":"#fd8d3c",
    "HSC-like":'#fbdd7e',
    "LSPC-Primed":"#A95C68",
    "LSPC-Quiescent":"#f6b092",
    "LSPC-Cycle":"#d62728",
    "GMP":"#e7ba52",
    "GMP-like":"#843c39",
    "Mono":"#8ca252",
    "ProMono":"#054907",
    'ProMono-like':'#5F8575',
    "Mono-like":"#98df8a",
    "B":"#9467bd",
    "ProB":"#c5b0d5",
    "earlyEry":"#ff9896",
    "lateEry":"#e377c2",
    "cDC":"#31a354",
    "cDC-like":"#20B2AA",
    "pDC":'#79eb28',
    'Other':'#dcdcdc'
}

sc.set_figure_params(scanpy=True, fontsize=10, dpi=180, dpi_save=500, frameon=False, format='png')
plt.rcParams.update({'font.weight': 'bold'})
sc.pl.umap(adata, color='pruned.labels', palette=cm, save='hspc_singler.png')

plt.rc("figure", figsize=(8,8))
ax = sc.pl.scatter(adata, x='LSC17', y='LSC6', color='pruned.labels', size=10, show=False)
ax.grid(None)
plt.savefig(plot_dir + 'scatter_hspc_lsc_scores.png', dpi=500)


# Set LSC thresholds
x = adata.obs['LSC6']
t6 = np.median(x) + 2 * x.mad()

plt.hist(x, bins=200)
plt.axvline(t6, color='red')
plt.savefig(plot_dir + 'LSC6_threshold.png')

x = adata.obs['LSC17']
t17 = np.median(x) + 2 * x.mad()

plt.hist(x, bins=200)
plt.axvline(t17, color='red')
plt.savefig(plot_dir + 'LSC17_threshold.png')

celltypes = list(set(adata.obs['pruned.labels']))
celltypes = [x for x in celltypes if str(x) != 'nan']
for c in celltypes:
    adata.obs['highlight_celltype'] = "Other"
    ad = adata[adata.obs['pruned.labels']==c]
    sc.pl.umap(ad, color='pruned.labels', palette=cm, save=c + '_singler_thresholds.png')
    ax = sc.pl.scatter(ad, x='LSC17', y='LSC6', color='pruned.labels', size=10, show=False)
    ax.grid(None)
    plt.axvline(t17, color='red')
    plt.axhline(t6, color='red')
    plt.savefig(plot_dir + c + '_lsc_scores_thresholds.png', dpi=500)
    

#################################################################################################################################
############################################################ Select LSCs ########################################################
#################################################################################################################################

# Select relevant celltypes only
relevant_celltypes = ['Prog-like', 'Prog', 'HSC', 'HSC-like', 'LSPC-Quiescent', 'LSPC-Cycle',
                      'LSPC-Primed']
s = np.array([s in relevant_celltypes for s in subset.obs['pruned.labels']])
subset = subset[s].copy()


# Make sure there are no healthy donor cells mixed in
subset = subset[subset.obs['ELN_Classification']!='Healthy_Donor']

print(subset.obs['pruned.labels'].value_counts())


# Annotate LSCs
lsc_cells = list(subset.obs.index)
adata.obs['LSC'] = 'No'
adata.obs.loc[lsc_cells, 'LSC'] = "Yes"
print(adata.obs['LSC'].value_counts())


adata.write_h5ad('outs/HSPC_annotated.h5ad')