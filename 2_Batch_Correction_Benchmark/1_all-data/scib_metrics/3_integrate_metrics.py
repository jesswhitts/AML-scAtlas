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

batch = sys.argv[1]
celltypes = sys.argv[2] #either 'normal' or 'full'
hvg = sys.argv[3]

# Load concatenated results
df = pd.read_csv("outs/" + batch + "_" + celltypes + "_metrics.csv", index_col=0)
df = df.drop_duplicates(keep=False)
df = df.T.add_prefix(hvg + "_").T
df = df.astype('float')
df.to_csv("outs/sorted_" + batch + "_" + celltypes + "_metrics.csv")

# Simplify rownames
value = len(df.index)
values = list(range(0,value,1))
for v in values:
    embedding1 = df.index[v]
    embedding2 = re.sub(r'^.*?X_', '', embedding1)
    df = df.rename(index={embedding1:embedding2})

# Scale for plotting
metrics = df
metrics_scaled = (metrics - metrics.min()) / (metrics.max() - metrics.min())
metrics_scaled = metrics_scaled.sort_values(by=["Total"], ascending=False)
metrics_styled = metrics_scaled.style.background_gradient(cmap="Blues")
dfi.export(metrics_styled,"outs/scaled_" + batch + "_" + celltypes + "_metrics_table.png", dpi=300, table_conversion="matplotlib")
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
plt.savefig('outs/' + batch + "_" + celltypes + "_metrics_scatter_plot.png")
