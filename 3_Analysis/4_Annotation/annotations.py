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

adata = sc.read_h5ad("outs/scvi_markers.h5ad")
print(adata)

anno = pd.read_csv("outs/full_anno.csv", index_col=0)
anno = anno.drop_duplicates(keep=False)

anno = anno[["predicted_labels", "majority_voting", "sr_main_labels", "sr_clust_main_labels",
             "sctype_labels"]]

adata.obs = adata.obs.join(anno)
del(anno)

# Add Van Galen Annotations
vgalen = pd.read_csv("vgalen_annotation/outs/SingleR_custom_predictions_full.csv", index_col=0)
vgalen = vgalen[['labels','pruned.labels']]
vgalen = vgalen.rename(columns={'labels':'vgalen_labels','pruned.labels':'vgalen_pruned_labels'})
adata.obs = adata.obs.join(vgalen)

zeng = pd.read_csv("vgalen_annotation/zeng_annotations/outs/SingleR_custom_predictions_full.csv", index_col=0)
zeng = zeng[['labels','pruned.labels']]
zeng = zeng.rename(columns={'labels':'zeng_labels','pruned.labels':'zeng_pruned_labels'})
adata.obs = adata.obs.join(zeng)

# Visualise Cluster
clusters = adata.obs["cluster"].unique()

for clust in clusters:
    sc.pl.umap(
        adata,
        color="cluster",
        groups=[clust],
        title=f"Cluster {clust}",
        legend_loc=None,
        save=f"_cluster_{clust}_highlight.png"
    )

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

adata.obs["pct_healthy_donor"] = adata.obs["pct_healthy_donor"].astype(str)

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

mapping_dict = {
    0:'CD14+ Mono',
    1:'T',
    2:'HSPC',
    3:'ProMono',
    4:'CMP',
    5:'T',
    6:'T',
    7:'HSPC',
    8:'B',
    9:'HSPC',
    10:'T',
    11:'HSPC',
    12:'MEP',
    13:'CMP',
    14:'Erythroid',
    15:'cDC',
    16:'NK',
    17:'HSPC',
    18:'GMP',
    19:'HSPC',
    20:'Erythroid',
    21:'HSPC',
    22:'HSPC',
    23:'HSPC',
    24:'GMP',
    25:'CD14+ Mono',
    26:'GMP',
    27:'CD14+ Mono',
    28:'CD16+ Mono',
    29:'GMP',
    30:'CD16+ Mono',
    31:'HSPC',
    32:'ProB',
    33:'MEP',
    34:'HSPC',
    35:'GMP',
    36:'Erythroid',
    37:'cDC',
    38:'T',
    39:'pDC',
    40:'HSPC',
    41:'CD14+ Mono',
    42:'ProB',
    43:'HSPC',
    44:'HSPC',
    45:'Plasma',
    46:'GMP',
    47:'GMP',
    48:'ProB',
    49:'Erythroid',
    50:'CD14+ Mono',
    51:'HSPC',
    52:'T',
    53:'GMP',
    54:'HSPC',
    55:'CD14+ Mono'
}

adata.obs['main_celltype'] = adata.obs['cluster'].astype(int)
adata.obs['main_celltype'] = adata.obs['main_celltype'].replace(mapping_dict)
adata.obs['main_celltype'].value_counts()

cm = {
  "T":"#4682B4",
  "HSPC":"#fd8d3c",
  "CMP":"#6E260E",
  "GMP":"#FFC000",
  "ProMono":"#097969",
  "CD14+ Mono":"#8ca252",
  "CD16+ Mono":"#054907",
  "B":"#9467bd",
  "ProB":"#c5b0d5",
  "cDC":"#4CBB17",
  "pDC":"#40E0D0",
  "Plasma":"#191970",
  "MEP":"#d62728",
  "NK":"#aec7e8",
  "Erythroid":"#e377c2",
  "Granulocytes":"#9C7C38"
}

sc.pl.umap(adata, color="main_celltype", palette=cm)
sc.pl.umap(adata, color="main_celltype", save="_main_celltype.png")

adata.write_h5ad('scvi_celltype_annotations.h5ad')

#################################################################################################################################
######################################################## Tidy Up Metadata #######################################################
#################################################################################################################################

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

# Update Caron et al ages - found in paper correction table
mask = adata.obs['sample'] == 'PBMMC1'
adata.obs.loc[mask, "age"] = '1'

mask = adata.obs['sample'] == 'PBMMC2'
adata.obs.loc[mask, "age"] = '4'

mask = adata.obs['sample'] == 'PBMMC3'
adata.obs.loc[mask, "age"] = '2'


#################################################################################################################################
######################################################## Add New Metadata #######################################################
#################################################################################################################################

# Age Groups

# Create a temporary numeric age
age_num = pd.to_numeric(adata.obs["age"], errors="coerce")
known_age = age_num.notna()

# Main Age Group
adata.obs["main_age_group"] = "unknown"

adata.obs.loc[known_age & (age_num >= 18), "main_age_group"] = "adult"
adata.obs.loc[known_age & (age_num < 18), "main_age_group"] = "paediatric"

# Fine Age Group
adata.obs["fine_age_group"] = "unknown"

adata.obs.loc[known_age & (age_num < 10), "fine_age_group"] = "0-9"
adata.obs.loc[known_age & (age_num >= 10) & (age_num < 20), "fine_age_group"] = "10-19"
adata.obs.loc[known_age & (age_num >= 20) & (age_num < 30), "fine_age_group"] = "20-29"
adata.obs.loc[known_age & (age_num >= 30) & (age_num < 50), "fine_age_group"] = "30-49"
adata.obs.loc[known_age & (age_num >= 50) & (age_num < 70), "fine_age_group"] = "50-69"
adata.obs.loc[known_age & (age_num >= 70), "fine_age_group"] = "Over 70"

# Memory Saving
adata.obs["main_age_group"] = adata.obs["main_age_group"].astype("category")
adata.obs["fine_age_group"] = adata.obs["fine_age_group"].astype("category")

sc.pl.umap(adata, color="main_age_group", save="_main_age_group.png")
sc.pl.umap(adata, color="fine_age_group", save="_fine_age_group.png")



# WHO Class
adata.obs["ELN_Classification"] = "Other"
adata.obs["ELN_Risk_Group"] = "Unclassified"
obs = adata.obs


# Healthy Donors
healthy = (
    (obs["mutational_subgroup"] == "Healthy_Donor") |
    (obs["clinical_subtype"] == "healthy_donor")
)

obs.loc[healthy, ["ELN_Classification", "ELN_Risk_Group"]] = "Healthy_Donor"


# Favourable
# NPM1
mask = obs["mutations"].str.contains("NPM1", na=False)
obs.loc[mask, "ELN_Classification"] = "NPM1"
obs.loc[mask, "ELN_Risk_Group"] = "Favourable"

# CEBPA
mask = obs["mutations"].str.contains("CEBPA", na=False)
obs.loc[mask, "ELN_Classification"] = "CEBPA"
obs.loc[mask, "ELN_Risk_Group"] = "Favourable"

# inv(16)
mask = (
    obs["mutational_subgroup"].str.contains(r"inv\(16\)", na=False) |
    obs["cytogenetics"].str.contains(r"inv\(16\)", na=False)
)
obs.loc[mask, "ELN_Classification"] = "inv(16)"
obs.loc[mask, "ELN_Risk_Group"] = "Favourable"

# CBFB-MYH11
mask = obs["translocations"].str.contains("CBFB-MYH11", na=False)
obs.loc[mask, "ELN_Classification"] = "CBFB-MYH11"
obs.loc[mask, "ELN_Risk_Group"] = "Favourable"

# RUNX1-RUNX1T1
mask = (
    obs["mutational_subgroup"].str.contains("RUNX1-RUNX1T1", na=False) |
    obs["translocations"].str.contains("RUNX1-RUNX1T1", na=False)
)
obs.loc[mask, "ELN_Classification"] = "RUNX1-RUNX1T1"
obs.loc[mask, "ELN_Risk_Group"] = "Favourable"


# Intermediate
# NPM1 with FLT3-ITD
mask = (
    (obs["ELN_Classification"] == "NPM1") &
    obs["mutations"].str.contains("FLT3-ITD", na=False)
)

obs.loc[mask, "ELN_Classification"] = "NPM1 with FLT3-ITD"
obs.loc[mask, "ELN_Risk_Group"] = "Intermediate"

# FLT3-ITD (non-NPM1)
mask = (
    ~obs["ELN_Classification"].isin(["NPM1", "NPM1 with FLT3-ITD"]) &
    obs["mutations"].str.contains("FLT3-ITD", na=False)
)

obs.loc[mask, "ELN_Classification"] = "FLT3-ITD"
obs.loc[mask, "ELN_Risk_Group"] = "Intermediate"


# Adverse
# MDS-related gene mutations
# should not be used as an adverse prognostic marker if they co-occur with favorable-risk AML subtypes (ELN 2022)
mds_genes = (
    obs["mutations"].str.contains(
        r"ASXL1|BCOR|EZH2|RUNX1|SF3B1|SRSF2|STAG2|U2AF1|ZRSR2",
        na=False))

mask = (obs["ELN_Risk_Group"] != "Favourable") & mds_genes

obs.loc[mask, "ELN_Classification"] = "MDS Related Genes"
obs.loc[mask, "ELN_Risk_Group"] = "Adverse"

# -7/-5
mask = obs["cytogenetics"].str.contains(r"-7", na=False)
obs.loc[mask, ["ELN_Classification", "ELN_Risk_Group"]] = ["-7", "Adverse"]

mask = obs["cytogenetics"].str.contains(r"-5", na=False)
obs.loc[mask, ["ELN_Classification", "ELN_Risk_Group"]] = ["-5", "Adverse"]

# Complex karyotype
mask = obs["mutational_subgroup"] == "CK"
obs.loc[mask, ["ELN_Classification", "ELN_Risk_Group"]] = ["CK", "Adverse"]

# MLL rearrangements (excluding favourable)
mask = (
    obs["mutational_subgroup"].str.contains("MLL", na=False) &
    ~obs["translocations"].isin(["MLL-AF9", "none_detected"])
)
obs.loc[mask, ["ELN_Classification", "ELN_Risk_Group"]] = ["MLLr-Other", "Adverse"]

# Specific adverse translocations / subgroups
mask = obs["mutational_subgroup"].str.contains("BCR-ABL", na=False)
obs.loc[mask, ["ELN_Classification", "ELN_Risk_Group"]] = ["BCR-ABL", "Adverse"]

mask = obs["mutational_subgroup"].str.contains("KAT6A-CREBBP", na=False)
obs.loc[mask, ["ELN_Classification", "ELN_Risk_Group"]] = ["KAT6A-CREBBP", "Adverse"]

mask = obs["mutational_subgroup"].str.contains("RUNX1-EVI1", na=False)
obs.loc[mask, ["ELN_Classification", "ELN_Risk_Group"]] = ["MECOM", "Adverse"]

# TP53
mask = obs["mutations"].str.contains("TP53", na=False)
obs.loc[mask, ["ELN_Classification", "ELN_Risk_Group"]] = ["TP53", "Adverse"]

# Intermediate (overrules other mutations)
mask = obs["translocations"].str.contains("MLL-AF9", na=False)
obs.loc[mask, ["ELN_Classification", "ELN_Risk_Group"]] = ["MLL-AF9", "Intermediate"]

# Memory Saving
adata.obs["ELN_Classification"] = adata.obs["ELN_Classification"].astype("category")
adata.obs["ELN_Risk_Group"] = adata.obs["ELN_Risk_Group"].astype("category")

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
adata.write_h5ad("scvi_annotated.h5ad.gz", compression='gzip')
adata.obs.to_csv("annotated_obs.csv")
