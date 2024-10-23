#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
###################################################### AML ATLAS ANNOTATION #####################################################
################################################### Curate Atlas Annotations ####################################################
#################################################################################################################################

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns
import scanpy as sc

if not os.path.exists("figures"):
    os.makedirs("figures")

sc.set_figure_params(scanpy=True, fontsize=10, dpi=180, dpi_save=500, frameon=False, format='png')


#################################################################################################################################
########################################################## Read AnnData #########################################################
#################################################################################################################################

adata = sc.read_h5ad("annotate_celltypes/outs/scvi_markers.h5ad")
print(adata)

anno = pd.read_csv("annotate_celltypes/outs/full_anno.csv", index_col=0)
anno = anno.drop_duplicates(keep=False)

anno = anno[["predicted_labels", "majority_voting", "sr_main_labels", "sr_clust_main_labels",
             "sctype_labels"]]

adata.obs = adata.obs.join(anno)

# Add Van Galen Annotations
vgalen = pd.read_csv("annotate_celltypes/vgalen_annotation/outs/SingleR_custom_predictions_full.csv", index_col=0)
vgalen = vgalen[['labels','pruned.labels']]
vgalen = vgalen.rename(columns={'labels':'vgalen_labels','pruned.labels':'vgalen_pruned_labels'})
adata.obs = adata.obs.join(vgalen)

zeng = pd.read_csv("annotate_celltypes/vgalen_annotation/zeng_annotations/outs/SingleR_custom_predictions_full.csv", index_col=0)
zeng = zeng[['labels','pruned.labels']]
zeng = zeng.rename(columns={'labels':'zeng_labels','pruned.labels':'zeng_pruned_labels'})
adata.obs = adata.obs.join(zeng)

# Visualise Cluster
tmp = adata
clusters = list(tmp.obs["cluster"].value_counts().index)
for clust in clusters:
    tmp.obs['cluster_highlight'] = "FALSE"
    sub = tmp[tmp.obs['cluster'] == str(clust)]
    c = sub.obs.index
    adata.obs.loc[c, 'cluster_highlight'] = "TRUE"
    sc.pl.umap(adata, color="cluster_highlight", title="Cluster " + str(clust), legend_loc=None,
          save="_cluster_" + str(clust) + "_highlight.png")


#################################################################################################################################
#################################### Calculate Percentage of Healthy Donor Cells per Cluster ####################################
#################################################################################################################################

# Calculate Percentage of Healthy Donor Cells
adata.obs["pct_healthy_donor"] = adata.obs["cluster"].astype(str)
adata.obs["malignant_class"] = adata.obs["cluster"].astype(str)

clusters = list(set(adata.obs["cluster"]))

for c in clusters:
    subset = adata[adata.obs['cluster']==c]
    df = pd.DataFrame(subset.obs['clinical'].value_counts())
    if any(ele in 'normal' for ele in list(df.index)) == True:
        count = int(df.loc['normal'][0])
        total = int(len(subset.obs))
        pct = round(count/total*100)
    else:
        pct = 0
    cells = subset.obs.index
    adata.obs.loc[cells, 'pct_healthy_donor'] = pct
    if pct < 5:
        adata.obs.loc[cells, 'malignant_class'] = 'Leukaemic'
    elif pct > 40:
        adata.obs.loc[cells, 'malignant_class'] = 'Non-Leukaemic'
    else:
        adata.obs.loc[cells, 'malignant_class'] = 'Intermediate'

# Plot pct healthy donor
df = adata.obs[['cluster','pct_healthy_donor']]
df['cluster'] = df['cluster'].astype(int)
df = df.sort_values(by='cluster', ascending=True)
df['cluster'] = df['cluster'].astype(str)
df = df.set_index('cluster')

ax = df[['pct_healthy_donor']].plot(kind='bar', rot=0,
                                    width=0.8, edgecolor = "black", figsize=(12,4))
plt.title("Percentage of Healthy Donor Cells Per Cluster", fontsize=12)
plt.savefig("figures/hsc_hsc_like_pct.png", dpi=500, format='png', bbox_inches='tight', pad_inches=0.5)

#################################################################################################################################
############################################### Plot Key Marker Genes per Cluster ###############################################
#################################################################################################################################

# Plot some marker gene expression per cluster
markers = ['CD34','ELANE','LYZ','IL7R','CD8A','CD79A','MS4A1','FCGR3A',
           'CD14','GATA1','GATA2','KLF1','PF4','HBB','ITGAX','HLA-DRA',
           'CALD1','COL1A2','COL3A1']

for m in markers:
    plt.rc("figure", figsize=(16,4))
    sc.pl.violin(adata, m, groupby='cluster', inner='box',
                 save= m + '_expr_cluster.png')

sc.set_figure_params(scanpy=True, fontsize=10, dpi=180, dpi_save=500, frameon=False, format='png')


#################################################################################################################################
################################## Review Automated Annotation Outputs, GSEA, and Marker Genes ##################################
#################################################################################################################################

# Add Obs Categories
adata.obs["main_celltype"] = adata.obs["cluster"]


# Cluster 0
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="0", value="T")

# Cluster 1
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="1", value="HSPC")

# Cluster 2
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="2", value="CD14+ Mono")

# Cluster 3
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="3", value="CMP")#??

# Cluster 4
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="4", value="ProMono")

# Cluster 5
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="5", value="T")

# Cluster 6
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="6", value="HSPC")

# Cluster 7
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="7", value="HSPC")

# Cluster 8
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="8", value="B")

# Cluster 9
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="9", value="NK")

# Cluster 10
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="10", value="HSPC")

# Cluster 11
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="11", value="CD14+ Mono")

# Cluster 12
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="12", value="Erythroid")

# Cluster 13
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="13", value="MEP")

# Cluster 14
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="14", value="cDC")

# Cluster 15
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="15", value="HSPC")#MEP

# Cluster 16
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="16", value="Erythroid")

# Cluster 17
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="17", value="HSPC")

# Cluster 18
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="18", value="CD16+ Mono")

# Cluster 19
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="19", value="GMP")

# Cluster 20
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="20", value="HSPC")

# Cluster 21
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="21", value="HSPC")#MEP?

# Cluster 22
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="22", value="HSPC")

# Cluster 23
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="23", value="CD14+ Mono")

# Cluster 24
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="24", value="CMP")# GMP?

# Cluster 25
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="25", value="HSPC")

# Cluster 26
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="26", value="CD14+ Mono")

# Cluster 27
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="27", value="CD14+ Mono")

# Cluster 28
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="28", value="CD14+ Mono")

# Cluster 29
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="29", value="cDC")

# Cluster 30
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="30", value="CD14+ Mono")

# Cluster 31
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="31", value="MEP")

# Cluster 32
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="32", value="T")

# Cluster 33
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="33", value="Erythroid")

# Cluster 34
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="34", value="ProB")

# Cluster 35
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="35", value="HSPC")

# Cluster 36
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="36", value="pDC")

# Cluster 37
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="37", value="ProB")

# Cluster 38
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="38", value="HSPC")

# Cluster 39
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="39", value="CD14+ Mono")

# Cluster 40
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="40", value="T")

# Cluster 41
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="41", value="Plasma")

# Cluster 42
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="42", value="Erythroid")

# Cluster 43
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="43", value="HSPC")

# Cluster 44
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="44", value="HSPC")

# Cluster 45
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="45", value="HSPC")

# Cluster 46
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="46", value="HSPC")

# Cluster 47
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="47", value="CMP")

# Cluster 48
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="48", value="GMP")

# Cluster 49
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="49", value="HSPC")

# Cluster 50
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="50", value="GMP")

# Cluster 51
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="51", value="HSPC")

# Cluster 52
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="52", value="HSPC")

# Cluster 53
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="53", value="CD14+ Mono")

# Cluster 54
adata.obs["main_celltype"] = adata.obs["main_celltype"].replace(to_replace="54", value="Stromal")


cm = {
  "T":"#1f77b4",
  "HSPC":"#fd8d3c",
  "CMP":"#6E260E",
  "GMP":"#FFC000",
  "ProMono":"#C4B454",
  "CD14+ Mono":"#008000",
  "CD16+ Mono":"#054907",
  "B":"#9467bd",
  "ProB":"#c5b0d5",
  "cDC":"#B4C424",
  "pDC":"#20B2AA",
  "Plasma":"#7DF9FF",
  "MEP":"#d62728",
  "NK":"#aec7e8",
  "Erythroid":"#e377c2",
  "Stromal":"#98FB98"
}

sc.pl.umap(adata, color="main_celltype", palette=cm)
sc.pl.umap(adata, color="main_celltype", save="_main_celltype.png")


#################################################################################################################################
######################################################## Tidy Up Metadata #######################################################
#################################################################################################################################

adata.obs["pct_healthy_donor"] = adata.obs["pct_healthy_donor"].astype(str)

# Rename clinical
adata.obs["clinical"] = adata.obs["clinical"].replace(to_replace="normal", value="healthy_donor")
adata.obs["clinical"] = adata.obs["clinical"].replace(to_replace="M0", value="AML-M0")
adata.obs["clinical"] = adata.obs["clinical"].replace(to_replace="M1", value="AML-M1")
adata.obs["clinical"] = adata.obs["clinical"].replace(to_replace="M2", value="AML-M2")
adata.obs["clinical"] = adata.obs["clinical"].replace(to_replace="M3", value="AML-M3")
adata.obs["clinical"] = adata.obs["clinical"].replace(to_replace="M4", value="AML-M4")
adata.obs["clinical"] = adata.obs["clinical"].replace(to_replace="M5", value="AML-M5")
adata.obs["clinical"] = adata.obs["clinical"].replace(to_replace="M6", value="AML-M6")
adata.obs["clinical"] = adata.obs["clinical"].replace(to_replace="mixed_phenotype", value="AML-mixed_phenotype")
adata.obs["clinical"] = adata.obs["clinical"].replace(to_replace="unknown", value="AML-unknown")

adata.obs = adata.obs.rename(columns={"clinical":"clinical_subtype"})

# Rename mutational subgroup
adata.obs["mutational_subgroup"] = adata.obs["mutational_subgroup"].replace(to_replace="NA", value="Healthy_Donor")
adata.obs["mutational_subgroup"] = adata.obs["mutational_subgroup"].replace(to_replace="MLL-AF6", value="MLL")
adata.obs["mutational_subgroup"] = adata.obs["mutational_subgroup"].replace(to_replace="MLL-AF9", value="MLL")
adata.obs["mutational_subgroup"] = adata.obs["mutational_subgroup"].replace(to_replace="MLL-AF10", value="MLL")
adata.obs["mutational_subgroup"] = adata.obs["mutational_subgroup"].replace(to_replace="MLL-ELL", value="MLL")
adata.obs["mutational_subgroup"] = adata.obs["mutational_subgroup"].replace(to_replace="MLL-PDS5A", value="MLL")
adata.obs["mutational_subgroup"] = adata.obs["mutational_subgroup"].replace(to_replace="RUNX1,BCOR", value="RUNX1")
adata.obs["mutational_subgroup"] = adata.obs["mutational_subgroup"].replace(to_replace="RUNX1-RUNXT1", value="RUNX1-RUNX1T1")
adata.obs["mutational_subgroup"] = adata.obs["mutational_subgroup"].replace(to_replace="TP55", value="TP53")
adata.obs["mutational_subgroup"] = adata.obs["mutational_subgroup"].replace(to_replace="none_detected", value="Unknown")
adata.obs["mutational_subgroup"] = adata.obs["mutational_subgroup"].replace(to_replace="unknown", value="Unknown")
adata.obs["mutational_subgroup"] = adata.obs["mutational_subgroup"].replace(to_replace="none_detected", value="Unknown")

sc.pl.umap(adata, color="mutational_subgroup", save="_mutational_subgroup.png")

# Fix Typos
adata.obs["mutations"] = adata.obs["mutations"].replace(to_replace="FLT3,MPL1", value="FLT3-ITD,MPL")

adata.obs["mutations"] = adata.obs["mutations"].replace(to_replace="FLT3,CEBPA,MYO15A", value="FLT3-ITD,CEBPA,MYO15A")

subset = adata[adata.obs['sample'].str.contains('AML870')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'translocations'] = "MLL-AF9"

subset = adata[adata.obs['sample'].str.contains('AMLL3266')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'translocations'] = "MLL-AF10"

adata.obs["translocations"] = adata.obs["translocations"].replace(to_replace="RUNX1-RUNXT1", value="RUNX1-RUNX1T1")


#################################################################################################################################
######################################################## Add New Metadata #######################################################
#################################################################################################################################

# Add Age Group Annotations
adata.obs["main_age_group"]="unknown"

subset = adata[adata.obs["age"]!="unknown"]
subset.obs["age"] = subset.obs["age"].astype('float')

group = subset[subset.obs["age"]>17]
sub = group.obs.index
adata.obs.loc[sub, 'main_age_group'] = "adult"

group = subset[subset.obs["age"]<18]
sub = group.obs.index
adata.obs.loc[sub, 'main_age_group'] = "paediatric"

adata.obs["fine_age_group"]="unknown"

group = subset[subset.obs["age"]<10]
sub = group.obs.index
adata.obs.loc[sub, 'fine_age_group'] = "0-9"

group = subset[subset.obs["age"]>9]
group = group[group.obs["age"]<20]
sub = group.obs.index
adata.obs.loc[sub, 'fine_age_group'] = "10-19"

group = subset[subset.obs["age"]>19]
group = group[group.obs["age"]<30]
sub = group.obs.index
adata.obs.loc[sub, 'fine_age_group'] = "20-29"

group = subset[subset.obs["age"]>29]
group = group[group.obs["age"]<50]
sub = group.obs.index
adata.obs.loc[sub, 'fine_age_group'] = "30-49"

group = subset[subset.obs["age"]>49]
group = group[group.obs["age"]<70]
sub = group.obs.index
adata.obs.loc[sub, 'fine_age_group'] = "50-69"

group = subset[subset.obs["age"]>69]
sub = group.obs.index
adata.obs.loc[sub, 'fine_age_group'] = "Over 70"

sc.pl.umap(adata, color="main_age_group", save="_main_age_group.png")
sc.pl.umap(adata, color="fine_age_group", save="_fine_age_group.png")


# WHO Class
adata.obs["ELN_Classification"] = "Other"
adata.obs["ELN_Risk_Group"] = "Unclassified"


# Healthy Donors
subset = adata[adata.obs["mutational_subgroup"]=="Healthy_Donor"]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "Healthy_Donor"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Healthy_Donor"

subset = adata[adata.obs["clinical_subtype"]=="healthy_donor"]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "Healthy_Donor"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Healthy_Donor"


# Favourable
subset = adata[adata.obs['mutations'].str.contains('NPM1')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "NPM1"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Favourable"

subset = adata[adata.obs['mutations'].str.contains('CEBPA')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "CEBPA"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Favourable"

subset = adata[adata.obs['mutational_subgroup'].str.contains('inv\(16\)')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "inv(16)"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Favourable"
subset = adata[adata.obs['cytogenetics'].str.contains('inv\(16\)')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "inv(16)"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Favourable"

subset = adata[adata.obs['translocations'].str.contains('CBFB-MYH11')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "CBFB-MYH11"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Favourable"

subset = adata[adata.obs['mutational_subgroup'].str.contains('RUNX1-RUNX1T1')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "RUNX1-RUNX1T1"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Favourable"
subset = adata[adata.obs['translocations'].str.contains('RUNX1-RUNX1T1')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "RUNX1-RUNX1T1"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Favourable"

# Intermediate
subset = adata[adata.obs["ELN_Classification"]=="NPM1"]
subset = subset[subset.obs['mutations'].str.contains('FLT3-ITD')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "NPM1 with FLT3-ITD"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Intermediate"

subset = adata[adata.obs["ELN_Classification"]!="NPM1"]
subset = subset[subset.obs["ELN_Classification"]!="NPM1 with FLT3-ITD"]
subset = subset[subset.obs['mutations'].str.contains('FLT3-ITD')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "FLT3-ITD"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Intermediate"


# Adverse
#should not be used as an adverse prognostic marker if they co-occur with favorable-risk AML subtypes (ELN 2022)
subset = adata[adata.obs["ELN_Risk_Group"]!="Favourable"]
subset = subset[subset.obs['mutations'].str.contains('ASXL1') | subset.obs['mutations'].str.contains('BCOR') |
                subset.obs['mutations'].str.contains('EZH2') | subset.obs['mutations'].str.contains('RUNX1') |
                subset.obs['mutations'].str.contains('SF3B1') | subset.obs['mutations'].str.contains('SRSF2') |
                subset.obs['mutations'].str.contains('STAG2') | subset.obs['mutations'].str.contains('U2AF1') |
                subset.obs['mutations'].str.contains('ZRSR2')]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "MDS Related Genes"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Adverse"

subset = adata[adata.obs['cytogenetics'].str.contains('-7')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "-7"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Adverse"

subset = adata[adata.obs['cytogenetics'].str.contains('-5')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "-5"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Adverse"

subset = adata[adata.obs["mutational_subgroup"]=="CK"]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "CK"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Adverse"

subset = adata[adata.obs['mutational_subgroup'].str.contains('MLL')==True]
subset = subset[subset.obs["translocations"]!='MLL-AF9']
subset = subset[subset.obs["translocations"]!='none_detected']
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "MLLr-Other"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Adverse"

subset = adata[adata.obs['mutational_subgroup'].str.contains('BCR-ABL')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "BCR-ABL"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Adverse"

subset = adata[adata.obs['mutational_subgroup'].str.contains('KAT6A-CREBBP')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "KAT6A-CREBBP"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Adverse"

subset = adata[adata.obs['mutational_subgroup'].str.contains('RUNX1-EVI1')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "MECOM"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Adverse"

subset = adata[adata.obs['mutations'].str.contains('TP53')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "TP53"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Adverse"


# Intermediate (overrules other mutations)
subset = adata[adata.obs['translocations'].str.contains('MLL-AF9')==True]
cells = subset.obs.index
adata.obs.loc[cells, 'ELN_Classification'] = "MLL-AF9"
adata.obs.loc[cells, 'ELN_Risk_Group'] = "Intermediate"


print(adata.obs['ELN_Classification'].value_counts())
print(adata.obs['ELN_Risk_Group'].value_counts())


# Annotate XIST/ChrY
# Using biomart to identify relevant genes
adata.X = adata.layers['counts'] # Use raw counts for this
annot = sc.queries.biomart_annotations(
        "hsapiens",
        ["ensembl_gene_id", "external_gene_name", "start_position", "end_position", "chromosome_name"],
    ).set_index("external_gene_name")

chrY_genes = adata.var_names.intersection(annot.index[annot.chromosome_name == "Y"])
adata.obs['percent_chrY'] = np.sum(
    adata[:, chrY_genes].X, axis=1) / np.sum(adata.X, axis=1) * 100
adata.obs["XIST-counts"] = adata.X[:,adata.var_names.str.match('XIST')].toarray()

sc.pl.scatter(adata, x='XIST-counts', y='percent_chrY', color="cluster", save='full_XIST_chrY_scatter.png')
sc.pl.violin(adata, "XIST-counts", jitter=0.4, groupby = 'cluster', rotation= 45, save='full_XIST_violin.png')
sc.pl.violin(adata,  "percent_chrY", jitter=0.4, groupby = 'cluster', rotation= 45, save='full_chrY_violin.png')


# # Save Annotated AnnData
adata.X = adata.layers['normalised_counts']
del(adata.layers['normalised_counts'])
adata.write_h5ad("scvi_annotated.h5ad")
adata.obs.to_csv("annotated_obs.csv")
