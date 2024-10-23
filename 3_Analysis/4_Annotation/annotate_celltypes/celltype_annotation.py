#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
################################################ AML ATLAS CELL TYPE ANNOTATION #################################################
################################################ Automated Cell Type Annotations ################################################
#################################################################################################################################

import os
import sys
import numpy as np
import gseapy as gp
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import statistics
import celltypist
from celltypist import models

#configure R - note ensure the correct R lib is used! "R home found: /Library/Frameworks/R.framework/Resources"
import rpy2
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
from rpy2.robjects.packages import importr
# import R's "base" package
base = importr('base')
# import R's "utils" package
utils = importr('utils')

#configure scanpy settings
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

get_ipython().run_cell_magic('R', '', 'suppressMessages(library(SingleCellExperiment))\nsuppressMessages(library(celldex))\nsuppressMessages(library(SingleR))\nsuppressMessages(library(HGNChelper))\nsuppressMessages(library(openxlsx))\nsuppressMessages(library(ggraph))\nsuppressMessages(library(igraph))\nsuppressMessages(library(tidyverse))\nsuppressMessages(library(data.tree))\nsuppressMessages(library(data.table))')

if not os.path.exists("outs"):
    os.makedirs("outs")

save_list = "full"
results_file = "outs/" + str(save_list) + "_annotation.h5ad"


#################################################################################################################################
####################################################### Load AnnData Object #####################################################
#################################################################################################################################

print("Reading anndata file...")

# Read in filtered anndata object - using uncorrected counts
adata = sc.read_h5ad("../../clustering/outs/scvi_clustered.h5ad")

adata = adata[:, adata.var["highly_variable"]]

print(adata)

#################################################################################################################################
################################################# Run Celltypist and Create SCE #################################################
#################################################################################################################################

print("Running Celltypist...")

#models.download_models(force_update = True)
model = models.Model.load(model = 'Immune_All_Low.pkl')
predictions = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = True)
adata = predictions.to_adata()

sc.pl.umap(adata, color = ['predicted_labels', 'majority_voting', "conf_score"], 
           save="_celltypist_" + str(save_list) + ".png")


print("Creating SCE input tables...")

#feature/metadata tables
var = pd.DataFrame(adata.var)
obs = pd.DataFrame(adata.obs[["sample","cluster"]])
obs["sample"] = obs["sample"].astype('string')
#need gene and cell names for expression matrix
genes = var.index
cells = obs.index
#make expression matrix
counts = pd.DataFrame(adata.X, columns=genes, index = cells)
#create umap matrix
UMAP = np.asmatrix(adata.obsm["X_umap"])

print("Transposing count matrix... This can take a while.")
get_ipython().run_cell_magic('R', '-i counts', '\n#transpose for compatibility with sce\ncounts <- t(counts)')


print("Creating SCE Object...")


get_ipython().run_cell_magic('R', '-i obs,UMAP', '\nsce <- SingleCellExperiment(\n  list(logcounts=as.matrix(counts)),\n  colData = obs, reducedDims = list(UMAP = UMAP))\n\nsce')

get_ipython().run_cell_magic('R', '', '\n saveRDS(sce, "outs/sce_hvg.rds")')

#################################################################################################################################
################################################# Run SingleR with Haem Dataset #################################################
#################################################################################################################################

print("Running SingleR...")

get_ipython().run_cell_magic('R', '', '#import dataset to use as labelling reference\nsuppressMessages(haem.data <- celldex::NovershternHematopoieticData(ensembl = FALSE)) ')

get_ipython().run_cell_magic('R', '', '#run SingleR \n\n#main labels\n#at individual cell level\nsingler <- SingleR(test = sce, ref = haem.data, labels =   haem.data@colData@listData$label.fine)\n\n#at cluster level\nclust_singler <- SingleR(test = sce, ref = haem.data, clusters=sce@colData@listData$cluster, \n                         labels=haem.data@colData@listData$label.fine)')

get_ipython().run_cell_magic('R', '', '#main labels\n#extract the assigned scores for each cell\nsingler_scores <- singler@listData$scores %>% as_tibble() %>%\n  dplyr::mutate(assigned_score = NA)\n\nfor ( i in seq_len(nrow(singler_scores)) ) {\n  singler_scores$assigned_score[i] <- singler_scores[[singler@listData$labels[i]]][i]}\n\n#assign labels and scores to sce object\nsce@colData$sr_main_labels <- singler@listData$labels\nsce@colData$sr_main_score <- singler_scores$assigned_score')

get_ipython().run_cell_magic('R', '', '\n#main labels\n#assign labels and scores to each cluster\nsce@colData$sr_clust_main_labels <- ""\nfor(j in unique(clust_singler@rownames)){\n  cl_type = clust_singler[clust_singler@rownames==j,]; \n  sce@colData$sr_clust_main_labels[sce@colData@listData$cluster == j] <- as.character(cl_type$labels[1])}\n\nsce@colData$sr_clust_main_score <- ""\nfor(j in unique(clust_singler@rownames)){\n  cl_type = clust_singler[clust_singler@rownames==j,]; \n  sce@colData$sr_clust_main_score[sce@colData@listData$cluster == j] <- max(cl_type@listData[["scores"]])}')

get_ipython().run_cell_magic('R', '-i save_list', '#visualise\n#main by cell\nplot <- plotScoreHeatmap(\n    singler,\n    show.labels = TRUE)\n\nsaveName <- paste0("figures/heatmap_singleR_byCell_", save_list, ".png")\n\nggsave(saveName, plot ,height = 11, width = 14)')

get_ipython().run_cell_magic('R', '', '#visualise\n#main by clust\nplot <- plotScoreHeatmap(\n    clust_singler,\n    show.labels = TRUE)\n\nsaveName <- paste0("figures/heatmap_singleR_byClust_", save_list, ".png")\n\nggsave(saveName, plot ,height = 11, width = 14)')



#################################################################################################################################
########################################################## Run scType ###########################################################
#################################################################################################################################

print("Running scType...")

get_ipython().run_cell_magic('R', '', '#import functions from github\nsource("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); \nsource("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")')

get_ipython().run_cell_magic('R', '', '# get cell-type-specific gene sets from our in-built database (DB)\ndb_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";\ntissue <- "Immune system" \n# e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus \n\n# prepare gene sets\ngs_list <- gene_sets_prepare(db_, tissue)')

get_ipython().run_cell_magic('R', '', '#extract count matrix\ndata <- sce@assays@data@listData[["logcounts"]]\n\n#input into sctype function\nes.max <- sctype_score(scRNAseqData = data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)')

get_ipython().run_cell_magic('R', '', '# merge by cluster\nclust_results <- do.call("rbind", lapply(unique(sce@colData@listData$cluster), function(cl){\n    es.max.cl = sort(rowSums(es.max[ ,rownames(sce@colData[sce@colData@listData$cluster==cl, ])]), decreasing = !0)\n    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, \n                    ncells = sum(sce@colData@listData$cluster==cl)), 10)}))\nsctype_scores <- clust_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  \n\n# set low-confident (low ScType score) clusters to "unknown"\nsctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"\nprint(sctype_scores, n=40)')

get_ipython().run_cell_magic('R', '', '#assign labels to sce\nsce@colData$sctype_labels <- ""\nfor(j in unique(sctype_scores$cluster)){\n  cl_type = sctype_scores[sctype_scores$cluster==j,]; \n  sce@colData$sctype_labels[sce@colData@listData$cluster == j] <- as.character(cl_type$type[1])}\n\n#also add sctype scores\nsce@colData$sctype_score <- ""\nfor(j in unique(sctype_scores$cluster)){\n  cl_type = sctype_scores[sctype_scores$cluster==j,]; \n  sce@colData$sctype_score[sce@colData@listData$cluster == j] <- as.numeric(cl_type$scores[1])}')

get_ipython().run_cell_magic('R', '', '#prepare edges\nclust_results=clust_results[order(clust_results$cluster),]; edges = clust_results; edges$type = \npaste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = \nedges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL')

get_ipython().run_cell_magic('R', '', '#prepare nodes\nnodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); \nnodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); \nnodes_lvl2 = c(); \nccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", \n          "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")\n\nfor (i in 1:length(unique(clust_results$cluster))){\n  dt_tmp = clust_results[clust_results$cluster == unique(clust_results$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, \n           data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, \n                      realname = dt_tmp$type))\n}\n\nnodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;\nfiles_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, \n                          files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)\nnodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", \n                                                                          "Colour", "ord", "shortName", "realname")]\n\nmygraph <- graph_from_data_frame(edges, vertices=nodes)')

get_ipython().run_cell_magic('R', '', '\n# Make the graph\ngggr<- ggraph(mygraph, layout = \'circlepack\', weight=I(ncells)) + \n  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + \n  geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +\n  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), \n  fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  \n  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")\n\nsaveName <- paste0("figures/bubblePlot_scType_", save_list, ".png")\nggsave(saveName, gggr ,height = 11, width = 14)')

get_ipython().run_cell_magic('R', '', '\n saveRDS(sce, "outs/sce_hvg_annotated.rds")')

print("Annotations Finished!")



#################################################################################################################################
######################################################## Collate Results ########################################################
#################################################################################################################################

get_ipython().run_cell_magic('R', '-o obs', '\nobs <- as.data.frame(colData(sce))')

get_ipython().run_cell_magic('R', '', '\nprint(sessionInfo())')

# Save results
cols = ['sr_main_labels','sr_main_score', 'sr_clust_main_labels', 'sr_clust_main_score',
        'sctype_labels', 'sctype_score']

obs[cols].to_csv("outs/R_obs.csv")

cols2 = ['predicted_labels', 'majority_voting', 'conf_score']

adata.obs[cols2].to_csv("outs/celltypist_obs.csv")


# Add to AnnData
adata.obs = pd.concat([adata.obs, obs[cols]], axis=1)
adata.obs.to_csv("adata_obs.csv")

# Make scores numeric for plotting
adata.obs['sr_main_score'] = pd.to_numeric(adata.obs['sr_main_score'])
adata.obs['sr_clust_main_score'] = pd.to_numeric(adata.obs['sr_clust_main_score'])
adata.obs['sctype_score'] = pd.to_numeric(adata.obs['sctype_score'])

cols = ["cluster", 'predicted_labels', 'majority_voting', 'conf_score', 'sr_main_labels',
       'sr_main_score', 'sr_clust_main_labels', 'sr_clust_main_score',
       'sctype_labels', 'sctype_score']

anno = adata.obs[cols]
anno.to_csv("outs/" + str(save_list) + "_anno.csv")

# Plot celltype labels
sc.pl.umap(adata, color=["cluster", "study", "majority_voting", "conf_score", "sr_main_labels", "sr_main_score", 
                         "sr_clust_main_labels", "sr_clust_main_score", "sctype_labels", "sctype_score"],
           legend_fontsize = "6", save="_celltypeAnno_" + str(save_list) + ".png")


print("FINISHED")


adata.write_h5ad("outs/annotated_hvg.h5ad")
