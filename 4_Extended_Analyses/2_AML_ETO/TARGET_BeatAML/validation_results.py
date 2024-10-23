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

tcga_adata = '../AML-ETO/pyscenic/gene_filtered/validation/TCGA/outs/AML_ETO_TARGET_BEATAML.h5ad'
tcga_auc = '../AML-ETO/pyscenic/gene_filtered/validation/TCGA/outs/auc_mtx_TARGET_BAML.csv'
tcga_anno = '../AML-ETO/pyscenic/gene_filtered/validation/TCGA/outs/cellAnnot_TARGET_BAML.csv'
counts_out = "../AML-ETO/pyscenic/gene_filtered/validation/TCGA/outs/TCGA_signature_groups_counts.csv"
groups_out = "../AML-ETO/pyscenic/gene_filtered/validation/TCGA/outs/TCGA_signature_groups_allocation.csv"
deseq = "../AML-ETO/pyscenic/gene_filtered/validation/TCGA/DESeq2_signature_Prenatal_vs_Postnatal.csv"
vds_cts = "../AML-ETO/pyscenic/gene_filtered/validation/TCGA/vst_filtered_signature_groups_counts.csv"
db1 = "/data/stemcell/jwhittle/ref/gene_sets/c2.all.v2023.1.Hs.symbols.gmt"
db2 = "/data/stemcell/jwhittle/ref/gene_sets/AML_drug_sensitivity.gmt"
target_outcomes = '../AML-ETO/pyscenic/gene_filtered/validation/TCGA/data/TARGET_outcomes.csv'
beataml_outcomes = '../AML-ETO/pyscenic/gene_filtered/validation/TCGA/data/BEATAML_outcomes.csv'


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


#################################################################################################################################
########################################################## Import Data ##########################################################
#################################################################################################################################

adata = sc.read_h5ad(tcga_adata)

# Define age groups
adata.obs["age_group"]="unknown"
subset = adata[adata.obs["age"]!="unknown"]
subset.obs["age"] = subset.obs["age"].astype('float')
group = subset[subset.obs["age"]>19]
sub = group.obs.index
adata.obs.loc[sub, 'age_group'] = "Adult"
group = subset[subset.obs["age"]<20]
group = group[group.obs["age"]>9]
sub = group.obs.index
adata.obs.loc[sub, 'age_group'] = "Adolescent"
group = subset[subset.obs["age"]<10]
sub = group.obs.index
adata.obs.loc[sub, 'age_group'] = "Pediatric"

# TCGA Pyscenic
cellAnnot = adata.obs
auc_mtx = pd.read_csv(tcga_auc, index_col=0)

# Find Columns with All 0 Values to Filter Out
remove = list(auc_mtx.loc[:, ~auc_mtx.any()].columns)
print(remove)
auc_mtx = auc_mtx.drop(remove, axis=1)


#################################################################################################################################
################################################# Visualise Signature Regulons ##################################################
#################################################################################################################################

# Set colourmap
age_cm = {
    'Adult':'#de8f05ff',
    'Pediatric':'#0173b2ff',
    'Adolescent':'#029e73ff'
}
age_cm = dict(sorted(age_cm.items()))

prenatal = ['TRIM28_(+)','CTCF_(+)','RAD21_(+)','SOX4_(+)','TAL1_(+)','MYB_(+)','FOXN3_(+)',
            'JUND_(+)','BCLAF1_(+)','ZBTB7A_(+)','IKZF1_(+)','MAZ_(+)','REST_(+)','YY1_(+)',
            'CUX1_(+)','KDM5A_(+)']
postnatal = ['YBX1_(+)','ENO1_(+)','HDAC2_(+)','GATA1_(+)','POLE3_(+)','TFDP1_(+)','MYBL2_(+)',
             'E2F4_(+)','KLF1_(+)','IRF1_(+)','STAT1_(+)','IRF7_(+)','MAFF_(+)','ATF4_(+)','TAGLN2_(+)',
             'SPI1_(+)','KLF2_(+)']
key_regulons = prenatal + postnatal

# Create Z Score matrix to enable comparison between regulons
auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
for col in list(auc_mtx.columns):
    auc_mtx_Z[ col ] = ( auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)

# Create heatmap
cats = sorted(list(set(cellAnnot['age_group'])))
age_colours = []
for key in age_cm:
    age_colours.append(age_cm[key])
colors = sns.color_palette(age_colours,n_colors=len(cats) )
colorsd = dict( zip( cats, colors ))
colormap2 = [ colorsd[x] for x in cellAnnot['age_group'] ]

sns.set()
sns.set(font_scale=0.8)
fig = palplot( colors, cats, size=1.0)
plt.savefig(plots + "TCGA_heatmap_legend_" + save_name + ".png", dpi=600, bbox_inches = "tight")

sns.set(font_scale=1.4)
plt.rcParams.update({'font.weight': 'bold'})
g = sns.clustermap(auc_mtx_Z[key_regulons], annot=False,  square=False, linecolor='gray',
    yticklabels=False, xticklabels=True, vmin=-5, vmax=5, row_colors=[colormap2],
    cmap='seismic', figsize=(18,16) )
g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')
plt.savefig(plots + "TCGA_heatmap_" + save_name + ".png", dpi=600, bbox_inches = "tight")


#################################################################################################################################
################################################### Select Signature Clusters ###################################################
#################################################################################################################################

# Select from clustermap
data = auc_mtx_Z[key_regulons]
reordered_rows = data.index[g.dendrogram_row.reordered_ind]
reordered_cols = data.columns[g.dendrogram_col.reordered_ind]
clustered_data = data.loc[reordered_rows,reordered_cols]

# Select Sample Clusters for Each Signature
# Plot Sample Cluster Dendogram
Z = g.dendrogram_row.linkage

plt.style.use('default')
plt.rcParams.update({'font.weight': 'bold'})
plt.figure(figsize=(25, 6))
#plt.title('Hierarchical Clustering', fontsize=20)
plt.xlabel('')
plt.ylabel('Euclidean Distance', fontsize=15, fontweight='bold')
plt.tight_layout()
plt.grid(False)

dendrogram(
    Z,
    labels=np.array(data.index),
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=6.,  # font size for the x axis labels
)

plt.savefig(plots + "heatmap_dendrogram_" + save_name + ".png", dpi=600, bbox_inches = "tight")

# Adjust Y Axis to match dendrogram
sns.set(font_scale=1.2)
g = sns.clustermap(auc_mtx_Z[key_regulons], annot=False,  square=False, linecolor='gray',
    yticklabels=False, xticklabels=True, vmin=-5, vmax=5, row_colors=[colormap2],
    cmap='seismic', figsize=(18,16) )

g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')
g.ax_heatmap.invert_yaxis()
g.ax_row_dendrogram.invert_yaxis()
g.ax_row_colors.invert_yaxis()
plt.savefig(plots + "TCGA_heatmap_" + save_name + ".png", dpi=600, bbox_inches = "tight")

postnatal_samps = ['TARGET-20-PANKEF-09A-01R','TARGET-20-PATDHA-09A-01R',
                   'TARGET-20-PANZXI-09A-01R','TARGET-21-PASLZE-09A-01R',
                   'TARGET-21-PASLZE-41A-02R','TARGET-20-PARUUB-09A-01R',
                   'TARGET-20-PASLDL-09A-02R','TARGET-20-PASLDL-09A-01R',
                   'BA3413R','TARGET-20-PAPVKG-09A-01R','BA2172R',
                   'TARGET-20-PASHYZ-04A-01R','TARGET-20-PANLIZ-04A-02R',
                   'TARGET-20-PASJEJ-09A-02R','TARGET-20-PARCZL-04A-02R',
                   'TARGET-20-PAKRUP-09A-02R','TARGET-20-PAPWHS-09A-03R',
                   'TARGET-20-PANDER-09A-02R','TARGET-20-PALGKX-09A-03R',
                   'TARGET-20-PARGVC-03A-03R','TARGET-20-PARBFZ-03A-02R',
                   'TARGET-20-PANJGR-09A-03R','TARGET-20-PANLIZ-09A-03R',
                   'TARGET-20-PARURW-09A-02R','TARGET-20-PANDER-09A-01R',
                   'TARGET-20-PARUUB-09A-02R','TARGET-20-PASHYZ-09A-03R']

prenatal_samps = ['TARGET-20-PASFIT-09A-01R','TARGET-20-PARENB-09A-01R',
                  'TARGET-20-PARENB-09A-03R','TARGET-20-PABHET-03A-02R',
                  'TARGET-20-PARLMY-09A-01R','TARGET-20-PARCZL-09A-02R',
                  'TARGET-20-PARZUU-04A-01R','BA3335R','TARGET-20-PANJTK-09A-01R',
                  'TARGET-20-PARIZR-03A-02R','TARGET-20-PASWPT-09A-01R',
                  'TARGET-20-PASXYG-04A-01R','TARGET-20-PABHET-09A-01R',
                  'TARGET-20-PATAST-04A-01R','TARGET-20-PARCZL-04A-01R',
                  'TARGET-20-PASHYZ-09A-01R','TARGET-20-PASADG-04A-01R',
                  'TARGET-20-PASPKE-09A-01R','TARGET-20-PARCVS-09A-01R',
                  'TARGET-20-PARZUU-09A-01R','TARGET-20-PARCZL-09A-01R',
                  'TARGET-20-PASADG-09A-01R','TARGET-20-PASXYG-09A-01R',
                  'TARGET-20-PARKUY-09A-01R','BA2449R','TARGET-20-PARDDA-09A-02R',
                  'TARGET-20-PARDDA-09A-01R','TARGET-20-PASKDS-03A-01R','BA3353R','BA2789R']

# Select Samples for Comparison
cellAnnot['signature'] = 'Intermediate'
cellAnnot.loc[postnatal_samps, 'signature'] = "Postnatal"
cellAnnot.loc[prenatal_samps, 'signature'] = "Prenatal"

cellAnnot['signature'].value_counts()

# Visualise Heatmap
g1 = sorted(list(set(cellAnnot['signature'])))
sig_cols = {
    'Postnatal':'#ff7c00',
    'Prenatal':'#001c7f',
    'Intermediate':'#D3D3D3'
}
sig_cols = dict(sorted(sig_cols.items()))
g1_colours = []
for key in sig_cols:
    g1_colours.append(sig_cols[key])
colors = sns.color_palette(g1_colours,n_colors=len(g1) )
colorsd = dict( zip( g1, colors ))
colormap1 = [ colorsd[x] for x in cellAnnot['signature'] ]

sns.set()
sns.set(font_scale=0.8)
fig = palplot( colors, g1, size=1.0)
plt.title("Signature", fontsize=12)
plt.savefig(plots + "heatmap_TCGA_signature_legend_1_" + save_name + ".png", dpi=600, bbox_inches = "tight")

age_colours = []
for key in age_cm:
    age_colours.append(age_cm[key])
colors = sns.color_palette(age_colours,n_colors=len(cats) )
colorsd = dict( zip( cats, colors ))
colormap3 = [ colorsd[x] for x in cellAnnot['age_group'] ]

sns.set()
sns.set(font_scale=0.8)
fig = palplot( colors, cats, size=1.0)
plt.title("Age Group", fontsize=12)
plt.savefig(plots + "heatmap_TCGA_signature_legend_3_" + save_name + ".png", dpi=600, bbox_inches = "tight")

sns.set(font_scale=1.4)
plt.rcParams.update({'font.weight': 'bold'})
g = sns.clustermap(auc_mtx_Z[key_regulons], annot=False,  square=False,  linecolor='gray',
    yticklabels=False, xticklabels=True, vmin=-5, vmax=5, row_colors=[colormap1,colormap3],
    cmap='seismic', figsize=(18,16) )
g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')
g.ax_heatmap.invert_yaxis()
g.ax_row_dendrogram.invert_yaxis()
plt.savefig(plots + "TCGA_heatmap_signature_" + save_name + ".png", dpi=600, bbox_inches = "tight")


#################################################################################################################################
################################################# Differential Gene Expression ##################################################
#################################################################################################################################

print(adata.obs['signature'].value_counts())

# Export Matrix to Run DESeq2
select = ['Prenatal','Postnatal']
subset = np.array([s in select for s in adata.obs.signature])
sub_adata = adata[subset].copy()

counts = pd.DataFrame(data=sub_adata.X, index=sub_adata.obs_names, columns=sub_adata.var_names).T
counts.to_csv(counts_out)

metadata = pd.DataFrame(sub_adata.obs['signature'])
metadata.to_csv(groups_out)

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
plt.rcParams.update({'font.weight': 'bold'})
sc.tl.rank_genes_groups(adata, 'signature', groups=['Prenatal'], reference='Postnatal', 
                        method='wilcoxon', key_added = "wilcoxon")
plt.rc("figure", figsize=(16,4)) 
sc.pl.rank_genes_groups(adata, n_genes=100, sharey=False, key="wilcoxon")

select = ['Prenatal','Postnatal']
subset = np.array([s in select for s in adata.obs.signature])
adata = adata[subset].copy()

tfs = []
for r in key_regulons:
    tf = r.replace("_(+)", "")
    if tf in adata.var_names:
        tfs.append(tf)

# Set Group Colours
adata.uns['signature_colors'] = np.array(['#dd8452','#4c72b0'])

# Plot Heatmap
sc.pl.heatmap(adata, tfs, groupby='signature', swap_axes=True, cmap='Spectral', save='TCGA_TF_expression.png')

dedf = sc.get.rank_genes_groups_df(adata, group='Prenatal', key="wilcoxon")
dedf = dedf.set_index('names')
marker_de = dedf.loc[tfs]
marker_de = marker_de[marker_de['pvals_adj']<0.01]

marker_de.sort_values(by='logfoldchanges')


##### DESeq2 - Run in R - Read Results Here
ds = pd.read_csv(deseq, index_col=0)
vds = pd.read_csv(vds_cts, index_col=0)

plt.style.use('default')
plt.rcParams.update({'font.weight': 'bold'})
plt.figure(figsize=(7, 5))
plt.scatter(x=ds['log2FoldChange'],y=ds['padj'].apply(lambda x:-np.log10(x)),s=1,label="Not significant",color='#C0C0C0')

# highlight down- or up- regulated genes
down = ds[(ds['log2FoldChange']<=-0.5)&(ds['padj']<=0.01)]
up = ds[(ds['log2FoldChange']>=0.5)&(ds['padj']<=0.01)]

plt.scatter(x=down['log2FoldChange'],y=down['padj'].apply(lambda x:-np.log10(x)),s=3,label="Down-regulated",color="#A7C7E7")
plt.scatter(x=up['log2FoldChange'],y=up['padj'].apply(lambda x:-np.log10(x)),s=3,label="Up-regulated",color="#FAA0A0")

# plot gene names for selected genes
genes_up = list( set(tfs) & set(up.index) )
print(genes_up)
texts=[]
for i,r in up.loc[genes_up].iterrows():
    texts.append(plt.text(x=r['log2FoldChange'],y=-np.log10(r['padj']),
                          fontsize=7.5, s=i, fontweight='bold'))
    plt.scatter(x=r['log2FoldChange'],y=-np.log10(r['padj']),s=4,color="black")
adjust_text(texts,arrowprops=dict(arrowstyle="-", color='black', lw=1))

genes_down = list( set(tfs) & set(down.index) )
print(genes_down)
texts=[]
for i,r in down.loc[genes_down].iterrows():
    texts.append(plt.text(x=r['log2FoldChange'],y=-np.log10(r['padj']),
                          fontsize=7.5, s=i, fontweight='bold'))
    plt.scatter(x=r['log2FoldChange'],y=-np.log10(r['padj']),s=4,color="black")
adjust_text(texts,arrowprops=dict(arrowstyle="-", color='black', lw=1))

plt.xlabel("log2FoldChange")
plt.ylabel("-logPAdj")
plt.axvline(-0.5,color="grey",linestyle="--",lw=0.8)
plt.axvline(0.5,color="grey",linestyle="--", lw=0.8)
plt.axhline(0.5,color="grey",linestyle="--", lw=0.8)
plt.legend(title="",
           fontsize=6,
           loc="center left",
           bbox_to_anchor=(1, 0, 0.5, 1),
           frameon=False)
plt.savefig(plots + "volcano_prenatal_vs_postnatal.png", format='png', bbox_inches='tight')

# Review Results
print(up.loc[genes_up].sort_values('log2FoldChange', ascending=False))
print(down.loc[genes_down].sort_values('log2FoldChange', ascending=True))
down = down.sort_values('log2FoldChange', ascending=True)
up = up.sort_values('log2FoldChange', ascending=False)


# Heatmap
genes = list(up.head(40).index) + list(down.head(40).index)
df = vds.T

genes = genes_up + genes_down

g1 = sorted(list(set(metadata['signature'])))
sig_cols = {
    'Postnatal':'#ff7c00',
    'Prenatal':'#001c7f',
}
sig_cols = dict(sorted(sig_cols.items()))
g1_colours = []
for key in sig_cols:
    g1_colours.append(sig_cols[key])
colors = sns.color_palette(g1_colours,n_colors=len(g1) )
colorsd = dict( zip( g1, colors ))
colormap = [ colorsd[x] for x in metadata['signature'] ]

sns.set(font_scale=1.4)
plt.rcParams.update({'font.weight': 'bold'})

g = sns.clustermap(df[genes], annot=False,  square=False,  linecolor='gray',
    yticklabels=False, xticklabels=True, vmin=4, vmax=15, row_colors=[colormap],
    cmap='Spectral_r', figsize=(12,16) )
g.cax.set_visible(True)
g.ax_heatmap.set_ylabel('')
g.ax_heatmap.set_xlabel('')
plt.savefig(plots + "heatmap_prenatal_vs_postnatal.png", format='png', bbox_inches='tight')


#################################################################################################################################
######################################################### GSEA on DEGs ##########################################################
#################################################################################################################################

bg = list(set(adata.var.index))

# Prepare gene rank
down = ds[(ds['log2FoldChange']<=-0.5)&(ds['padj']<=0.01)]
up = ds[(ds['log2FoldChange']>=0.5)&(ds['padj']<=0.01)]
gene_rank = pd.concat([up,down])

gene_rank = gene_rank.sort_values('log2FoldChange', ascending=False)
gene_rank = gene_rank[['log2FoldChange']]


# Select pathways to plot
db = db1
res = gp.prerank(rnk=gene_rank, gene_sets=db, background=bg)
filtered_results = res.res2d[res.res2d['FDR q-val']<0.05]
path = pathlib.Path(db)
path = path.stem

terms = ['GAL_LEUKEMIC_STEM_CELL_DN', 'SHEN_SMARCA2_TARGETS_UP',
         'AFFAR_YY1_TARGETS_DN','RAMALHO_STEMNESS_DN']
plt.style.use('default')
plt.rcParams.update({'font.weight': 'bold'})
cols_1 = ['#e377c2','#8c564b','#bcbd22','#17becf']
plt.rcParams["axes.prop_cycle"] = plt.cycler(color=cols_1)
plt.figure(figsize=(4, 3))
axs = res.plot(terms, show_ranking=False, legend_kws={'loc': (1.05, 0)}, )
axs.savefig(plots + "db1_prenatal_vs_postnatal_selected.png", format='png', bbox_inches='tight')


# Select pathways to plot
db = db2
res = gp.prerank(rnk=gene_rank, gene_sets=db, background=bg)
filtered_results = res.res2d[res.res2d['FDR q-val']<0.05]
path = pathlib.Path(db)
path = path.stem

terms = ['Venetoclax_Resistance_Combined_Zhang_2020', 'Cytarabine_Response_Xu_2019', 
         'Azacitidine_Response_Unnikrishnan_2017','Daunorubicin_Resistance_Williams_2019']
plt.style.use('default')
plt.rcParams.update({'font.weight': 'bold'})
cols_1 = ['#d95f02', '#7570b3', '#e7298a', '#66a61e']
plt.rcParams["axes.prop_cycle"] = plt.cycler(color=cols_1)
plt.figure(figsize=(4, 3))
axs = res.plot(terms, show_ranking=False, legend_kws={'loc': (1.05, 0)}, )
axs.savefig(plots + "db2_prenatal_vs_postnatal_selected.png", format='png', bbox_inches='tight')