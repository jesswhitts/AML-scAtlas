#!/usr/bin/env python
# coding: utf-8


#################################################################################################################################
################################################## PYSCENIC: INTERPRET RESULTS ##################################################
############################################# Interpret and Visualise SCENIC Results ############################################
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

from pyscenic.plotting import plot_binarization
from pyscenic.export import add_scenic_metadata
from pyscenic.cli.utils import load_signatures
from pyscenic.utils import load_motifs
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
from adjustText import adjust_text #version 0.8 - newer ones have error
from pyscenic.binarization import binarize

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
auc_loom = 'outs/downsampled_raw_AUCell.loom'
save_name = "scAtlas_aml_eto"
save_folder = "outs/"

if not os.path.exists("outs"):
    os.makedirs("outs")


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


# Rename cell barcodes to remove suffix (added
for row in range(0,len(auc_mtx.index),1):
    cellid = auc_mtx.index[row]
    sep = '-'
    count = cellid.count(sep)
    res = cellid.split(sep, count)
    res = res[:-1]
    new_cellid = sep.join(res)
    auc_mtx = auc_mtx.rename(index = {cellid:new_cellid})

auc_mtx.to_csv(save_folder + save_name + "_auc_mtx.csv")

lf.close()