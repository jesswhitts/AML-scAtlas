#!/usr/bin/env python
# coding: utf-8


#################################################################################################################################
################################################## PYSCENIC: INTERPRET RESULTS ##################################################
######################################### Quick Prelim SCENIC Analysis - For Validation #########################################
#################################################################################################################################

import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import loompy as lp
import json
import base64
import zlib
import warnings

from pyscenic.export import add_scenic_metadata
from pyscenic.cli.utils import load_signatures
from pyscenic.utils import load_motifs
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
from adjustText import adjust_text #version 0.8 - newer ones have error

import operator as op
from cytoolz import compose
from IPython.display import HTML, display
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

warnings.filterwarnings('ignore')
sc.settings.verbosity = 0 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.set_figure_params(frameon=False, dpi=600, fontsize=10, dpi_save=600)

# Set Analysis Input Variables
analysis = "age_group"
auc_loom = 'outs/AML_ETO_TARGET_BEATAML_AUCell.loom'
input_adata = 'outs/AML_ETO_TARGET_BEATAML.h5ad'
save_name = "TARGET_BAML"
save_folder = "outs/"


#################################################################################################################################
####################################################### Define Functions ########################################################
#################################################################################################################################

def dfToNamedMatrix(df):
    arr_ip = [tuple(i) for i in df.values]
    dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
    arr = np.array(arr_ip, dtype=dtyp)
    return arr

def derive_regulons(motifs, db_names=['hg38__refseq-r80__10kb_up_and_down_tss.genes_vs_motifs.rankings',
                                      'hg38__refseq-r80__500bp_up_and_100bp_down_tss.genes_vs_motifs.rankings']):
    motifs.columns = motifs.columns.droplevel(0)

    def contains(*elems):
        def f(context):
            return any(elem in context for elem in elems)
        return f
    

    # For the creation of regulons we only keep the 10-species databases and the activating modules. We also remove the
    # enriched motifs for the modules that were created using the method 'weight>50.0%' (because these modules are not part
    # of the default settings of modules_from_adjacencies anymore.
    motifs = motifs[
        np.fromiter(map(compose(op.not_, contains('weight>50.0%')), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains(*db_names), motifs.Context), dtype=bool) & \
        np.fromiter(map(contains('activating'), motifs.Context), dtype=bool)]

    # We build regulons only using enriched motifs with a NES of 3.0 or higher; we take only directly annotated TFs or TF annotated
    # for an orthologous gene into account; and we only keep regulons with at least 10 genes.
    regulons = list(filter
                    (lambda r: len(r) >= 10, df2regulons(
                        motifs[(motifs['NES'] >= 3.0) & ((motifs['Annotation'] == 'gene is directly annotated')
                                                         | (motifs['Annotation'].str.startswith('gene is orthologous to')
                                                            & motifs['Annotation'].str.endswith('which is directly annotated for motif')))
                                                                     ])))
    
    # Rename regulons, i.e. remove suffix.
    return list(map(lambda r: r.rename(r.transcription_factor), regulons))

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)
    
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
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f


#################################################################################################################################
################################################### Add Metadata to Loom File ###################################################
#################################################################################################################################

lf = lp.connect( auc_loom, mode='r+', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
print("Regulons Detected in Full File: " + str(len(auc_mtx.T)))
regulons = lf.ra.Regulons

auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(')
regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )

# Regulon thresholds
rt = meta['regulonThresholds']
for i,x in enumerate(rt):
    tmp = x.get('regulon').replace("(","_(")
    x.update( {'regulon': tmp} )
    
# Update AnnData
adata = sc.read_h5ad(input_adata)
adata.obs["sample"] = adata.obs["patient"]

# Set Age Groups
adata.obs["main_age_group"]="unknown"
adata.obs["fine_age_group"]="unknown"
adata.obs["age_at_index"] = adata.obs["age_at_index"].astype('float')

# Main age groups
group = adata[adata.obs["age_at_index"]>18]
sub = group.obs.index
adata.obs.loc[sub, 'main_age_group'] = "adult"

group = adata[adata.obs["age_at_index"]<18]
sub = group.obs.index
adata.obs.loc[sub, 'main_age_group'] = "paediatric"

# Fine age groups
group = adata[adata.obs["age_at_index"]<10]
sub = group.obs.index
adata.obs.loc[sub, 'fine_age_group'] = "0-9"

group = adata[adata.obs["age_at_index"]>9]
group = group[group.obs["age_at_index"]<20]
sub = group.obs.index
adata.obs.loc[sub, 'fine_age_group'] = "10-19"

group = adata[adata.obs["age_at_index"]>19]
group = group[group.obs["age_at_index"]<30]
sub = group.obs.index
adata.obs.loc[sub, 'fine_age_group'] = "20-29"

group = adata[adata.obs["age_at_index"]>29]
group = group[group.obs["age_at_index"]<40]
sub = group.obs.index
adata.obs.loc[sub, 'fine_age_group'] = "30-39"

group = adata[adata.obs["age_at_index"]>39]
group = group[group.obs["age_at_index"]<50]
sub = group.obs.index
adata.obs.loc[sub, 'fine_age_group'] = "40-49"

group = adata[adata.obs["age_at_index"]>49]
group = group[group.obs["age_at_index"]<60]
sub = group.obs.index
adata.obs.loc[sub, 'fine_age_group'] = "50-59"

group = adata[adata.obs["age_at_index"]>59]
group = group[group.obs["age_at_index"]<70]
sub = group.obs.index
adata.obs.loc[sub, 'fine_age_group'] = "60-69"

group = adata[adata.obs["age_at_index"]>69]
group = group[group.obs["age_at_index"]<80]
sub = group.obs.index
adata.obs.loc[sub, 'fine_age_group'] = "70-79"

group = adata[adata.obs["age_at_index"]>79]
sub = group.obs.index
adata.obs.loc[sub, 'fine_age_group'] = "80+"

# Create new age groups
adata.obs["age_group"] = adata.obs["fine_age_group"]
adata.obs["age_group"] = adata.obs["age_group"].replace(to_replace="30-39", value="30-49")
adata.obs["age_group"] = adata.obs["age_group"].replace(to_replace="40-49", value="30-49")
adata.obs["age_group"] = adata.obs["age_group"].replace(to_replace="50-59", value="50-69")
adata.obs["age_group"] = adata.obs["age_group"].replace(to_replace="60-69", value="50-69")

# Add Relevant Metadata
metaJson = {}

metaJson["annotations"] = [
    {
        "name": "sample",
        "values": list(set(adata.obs['sample'].values))
    },
    {
        "name": "gender",
        "values": list(set(adata.obs['gender'].values))
    },
    {
        "name": "age_at_index",
        "values": list(set(adata.obs['age_at_index'].values))
    },
    {
     	"name": "age_group",
        "values": list(set(adata.obs['age_group'].values))
    },
    {
     	"name": "main_age_group",
        "values": list(set(adata.obs['main_age_group'].values))
    },
    {
     	"name": "fine_age_group",
        "values": list(set(adata.obs['fine_age_group'].values))
    }
]

# SCENIC regulon thresholds:
metaJson["regulonThresholds"] = rt

col_attrs = {
    "CellID": np.array(adata.obs.index),
    "sample": np.array(adata.obs['sample'].values),
    "main_age_group": np.array(adata.obs['main_age_group'].values),
    "fine_age_group": np.array(adata.obs['fine_age_group'].values),
    "gender": np.array(adata.obs['gender'].values),
    "age_at_index": np.array(adata.obs['age_at_index'].values),
    "age_group": np.array(adata.obs['age_group'].values),
    "RegulonsAUC": dfToNamedMatrix(auc_mtx),
}

row_attrs = {
    "Gene": lf.ra.Gene,
    "Regulons": regulons,
}

attrs = {
    "title": "sampleTitle",
    "MetaData": json.dumps(metaJson, cls=NpEncoder),
    "Genome": 'hg38',
    "SCopeTreeL1": "",
    "SCopeTreeL2": "",
    "SCopeTreeL3": ""
}

# Compress the metadata field:
attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson, cls=NpEncoder).encode('ascii'))).decode('ascii')


# Create new loom 
lp.create(
    filename = save_folder + "scenic_out_" + save_name + ".loom",
    layers=lf[:,:],
    row_attrs=row_attrs, 
    col_attrs=col_attrs, 
    file_attrs=attrs
)

lf.close()



#################################################################################################################################
################################################## Extract Data from Loom File ##################################################
#################################################################################################################################

# Extract Data from Loom File
lf = lp.connect(save_folder + "scenic_out_" + save_name + ".loom", mode='r', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
auc_mtx.to_csv(save_folder + "auc_mtx_" + save_name + ".csv")

regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
    regulons[i] =  list(r[r==1].index.values)
    
cellAnnot = pd.concat(
    [
        pd.DataFrame( lf.ca.sample, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.gender, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.age_at_index, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.main_age_group, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.fine_age_group, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.age_group, index=lf.ca.CellID ),
    ],
    axis=1
)
cellAnnot.columns = [
 'sample',
 'gender',
 'age_at_index',
 'age_group',
 'main_age_group',
 'fine_age_group']

cellAnnot.to_csv(save_folder + "cellAnnot_" + str(save_name) + ".csv")

lf.close()
