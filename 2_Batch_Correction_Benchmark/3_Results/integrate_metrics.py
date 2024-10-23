#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
################################################## AML ATLAS BATCH CORRECTION ###################################################
############################################## Integrate Results from scIB Metrics ##############################################
#################################################################################################################################


import pandas as pd
import matplotlib.pyplot as plt
import re
import dataframe_image as dfi
import sys
import os

groups = ["all-data","hvg2000","hvg4000",'hvg6000','hvg8000','hvg10000']

full_results = pd.DataFrame(columns=["Embedding","Silhouette label","cLISI","Silhouette batch","iLISI","KBET",
                                     "Graph connectivity","PCR comparison","Batch correction","Bio conservation","Total"])

for g in groups:

    result = pd.read_csv("../" + g + "/scib_metrics/outs/sorted_sample_full_metrics.csv", index_col=0)

    value = len(result.index)
    values = list(range(0,value,1))

    for v in values:
        embedding1 = result.index[v]
        embedding2 = re.sub(r'^.*?X_', '', embedding1)
        embedding2 = g + "_" + embedding2
        result = result.rename(index={embedding1:embedding2})

    full_results = pd.concat([full_results, result])

full_results.to_csv("full_results.csv")

# Scale for plotting
metrics = full_results
metrics_scaled = (metrics - metrics.min()) / (metrics.max() - metrics.min())
metrics_scaled = metrics_scaled.sort_values(by=["Total"], ascending=False)
metrics_styled = metrics_scaled.style.background_gradient(cmap="Blues")
dfi.export(metrics_styled,"scaled_metrics_table.png", dpi=300, table_conversion="matplotlib")
metrics_scaled.style.background_gradient(cmap="Blues")

fig, ax = plt.subplots()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
metrics_scaled.plot.scatter(
    x="Batch correction",
    y="Bio conservation",
    c=range(len(metrics_scaled)),
    ax=ax,
)

for k, v in metrics_scaled[["Batch correction", "Bio conservation"]].iterrows():
    ax.annotate(
        k,
        v,
        xytext=(6, -3),
        textcoords="offset points",
        family="sans-serif",
        fontsize=8,
    )
plt.savefig("metrics_scatter_plot.png")
