#!/usr/bin/env python
# coding: utf-8


import scanpy as sc
import anndata as ad
import pandas as pd
import loompy as lp
import numpy as np

counts = 'data/AML_ETO_sample_counts.csv'
tgt = 'data/TARGET_AML_ETO_obs.csv'
baml = 'data/BEATAML_obs.csv'
lm = 'outs/AML_ETO_TARGET_BEATAML.loom'
ad = 'outs/AML_ETO_TARGET_BEATAML.h5ad'


df = pd.read_csv(counts, index_col=0).T
target = pd.read_csv(tgt, index_col=0)
beataml = pd.read_csv(baml, index_col=0)

target['age'] = target['age_at_index']
beataml['age'] = beataml['age_at_diagnosis_years']
beataml['barcode'] = beataml['sample']

obs = pd.concat([target,beataml], axis=0)
obs = obs.set_index('barcode')


adata = sc.AnnData(df)
adata.obs = obs
adata.write_h5ad(ad)


df = df.T

row_attrs = { 
    "Gene": np.array(df.index) ,
}

col_attrs = { 
    "CellID":  np.array(df.columns) ,
}

lp.create(lm, 
          df.to_numpy(), row_attrs, col_attrs )

