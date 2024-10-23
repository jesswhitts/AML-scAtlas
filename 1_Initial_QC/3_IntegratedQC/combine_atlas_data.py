#!/usr/bin/env python
# coding: utf-8

#################################################################################################################################
################################################### AML ATLAS PRE-PROCESSING ####################################################
######### Concatenate All Atlas scRNA-seq Datasets, Add Original Authors Cell Types, and Perform Cell/Gene QC Filtering #########
#################################################################################################################################

import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from anndata import AnnData
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

#configure scanpy settings
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

# create output directories
if not os.path.exists("outs"):
    os.makedirs("outs")

#set results output file
results_file = 'outs/filtered_atlasData_original_celltypes.h5ad.gz'


#################################################################################################################################
################################################### Concatenate AnnData Files ###################################################
#################################################################################################################################

# AML Samples

print("Petti")
#petti - cell filtering performed
aml508084 = sc.read_h5ad("../1-petti/aml508084/outs/AML508084_cell_filtered.h5ad")
aml548327 = sc.read_h5ad("../1-petti/aml548327/outs/AML548327_cell_filtered.h5ad")
aml721214 = sc.read_h5ad("../1-petti/aml721214/outs/AML721214_cell_filtered.h5ad")
aml782328 = sc.read_h5ad("../1-petti/aml782328/outs/AML782328_cell_filtered.h5ad")
aml809653 = sc.read_h5ad("../1-petti/aml809653/outs/AML809653_cell_filtered.h5ad")

data = [aml508084,aml548327,aml721214,aml782328,aml809653]

adata = ad.concat(data, merge="same")
print(adata)

print("Zheng")
#zheng - cell filtering performed
aml027 = sc.read_h5ad("../2-zheng/aml027/outs/AML027_cell_filtered.h5ad")
aml035 = sc.read_h5ad("../2-zheng/aml035/outs/AML035_cell_filtered.h5ad")

data = [aml027,aml035]
adata = ad.concat(data, merge="same")
print(adata)

print("Van Galen")
#van galen - cells already filtered
aml1012 = sc.read_h5ad("../3-van-galen/aml1012/outs/AML1012_unfiltered.h5ad")
aml210A = sc.read_h5ad("../3-van-galen/aml210A/outs/AML210A_unfiltered.h5ad")
aml314 = sc.read_h5ad("../3-van-galen/aml314/outs/AML314_unfiltered.h5ad")
aml328 = sc.read_h5ad("../3-van-galen/aml328/outs/AML328_unfiltered.h5ad")
aml329 = sc.read_h5ad("../3-van-galen/aml329/outs/AML329_unfiltered.h5ad")
aml371 = sc.read_h5ad("../3-van-galen/aml371/outs/AML371_unfiltered.h5ad")
aml419A = sc.read_h5ad("../3-van-galen/aml419A/outs/AML419A_unfiltered.h5ad")
aml420B = sc.read_h5ad("../3-van-galen/aml420B/outs/AML420B_unfiltered.h5ad")

aml475 = sc.read_h5ad("../3-van-galen/aml475/outs/AML475_unfiltered.h5ad")
aml556 = sc.read_h5ad("../3-van-galen/aml556/outs/AML556_unfiltered.h5ad")
aml707B = sc.read_h5ad("../3-van-galen/aml707B/outs/AML707B_unfiltered.h5ad")
aml722B = sc.read_h5ad("../3-van-galen/aml722B/outs/AML722B_unfiltered.h5ad")
aml870 = sc.read_h5ad("../3-van-galen/aml870/outs/AML870_unfiltered.h5ad")
aml916 = sc.read_h5ad("../3-van-galen/aml916/outs/AML916_unfiltered.h5ad")
aml921A = sc.read_h5ad("../3-van-galen/aml921A/outs/AML921A_unfiltered.h5ad")
aml997 = sc.read_h5ad("../3-van-galen/aml997/outs/AML997_unfiltered.h5ad")

data = [aml1012,aml210A,aml314,aml328,aml329,aml371,aml419A,aml420B,aml475,aml556,aml707B,aml722B,aml870,aml916,aml921A,aml997]
adata = ad.concat(data, merge="same")
print(adata)

print("Stetson")
#stetson - cells already filtered
aml681 = sc.read_h5ad("../4-stetson/aml681/outs/AML681_unfiltered.h5ad")
aml682 = sc.read_h5ad("../4-stetson/aml682/outs/AML682_unfiltered.h5ad")
aml683 = sc.read_h5ad("../4-stetson/aml683/outs/AML683_unfiltered.h5ad")
aml684 = sc.read_h5ad("../4-stetson/aml684/outs/AML684_unfiltered.h5ad")
aml685 = sc.read_h5ad("../4-stetson/aml685/outs/AML685_unfiltered.h5ad")

data = [aml681,aml682,aml683,aml684,aml685]
adata = ad.concat(data, merge="same")
print(adata)

print("Wu")
# #wu - cell filtering performed
# #samples 1-20
amlP01 = sc.read_h5ad("../5-wu/amlP01/outs/AMLP01_cell_filtered.h5ad")
amlP02 = sc.read_h5ad("../5-wu/amlP02/outs/AMLP02_cell_filtered.h5ad")
amlP03 = sc.read_h5ad("../5-wu/amlP03/outs/AMLP03_cell_filtered.h5ad")
amlP04 = sc.read_h5ad("../5-wu/amlP04/outs/AMLP04_cell_filtered.h5ad")
amlP05 = sc.read_h5ad("../5-wu/amlP05/outs/AMLP05_cell_filtered.h5ad")

amlP06 = sc.read_h5ad("../5-wu/amlP06/outs/AMLP06_cell_filtered.h5ad")
amlP07 = sc.read_h5ad("../5-wu/amlP07/outs/AMLP07_cell_filtered.h5ad")
amlP08 = sc.read_h5ad("../5-wu/amlP08/outs/AMLP08_cell_filtered.h5ad")
amlP09 = sc.read_h5ad("../5-wu/amlP09/outs/AMLP09_cell_filtered.h5ad")
amlP10 = sc.read_h5ad("../5-wu/amlP10/outs/AMLP10_cell_filtered.h5ad")

amlP11 = sc.read_h5ad("../5-wu/amlP11/outs/AMLP11_cell_filtered.h5ad")
amlP12 = sc.read_h5ad("../5-wu/amlP12/outs/AMLP12_cell_filtered.h5ad")
amlP13 = sc.read_h5ad("../5-wu/amlP13/outs/AMLP13_cell_filtered.h5ad")
amlP14 = sc.read_h5ad("../5-wu/amlP14/outs/AMLP14_cell_filtered.h5ad")
amlP15 = sc.read_h5ad("../5-wu/amlP15/outs/AMLP15_cell_filtered.h5ad")

amlP16 = sc.read_h5ad("../5-wu/amlP16/outs/AMLP16_cell_filtered.h5ad")
amlP17 = sc.read_h5ad("../5-wu/amlP17/outs/AMLP17_cell_filtered.h5ad")
amlP18 = sc.read_h5ad("../5-wu/amlP18/outs/AMLP18_cell_filtered.h5ad")
amlP19 = sc.read_h5ad("../5-wu/amlP19/outs/AMLP19_cell_filtered.h5ad")
amlP20 = sc.read_h5ad("../5-wu/amlP20/outs/AMLP20_cell_filtered.h5ad")

# #wu
# #samples 21-42
amlP21 = sc.read_h5ad("../5-wu/amlP21/outs/AMLP21_cell_filtered.h5ad")
amlP22 = sc.read_h5ad("../5-wu/amlP22/outs/AMLP22_cell_filtered.h5ad")
amlP23 = sc.read_h5ad("../5-wu/amlP23/outs/AMLP23_cell_filtered.h5ad")
amlP24 = sc.read_h5ad("../5-wu/amlP24/outs/AMLP24_cell_filtered.h5ad")
amlP25 = sc.read_h5ad("../5-wu/amlP25/outs/AMLP25_cell_filtered.h5ad")

amlP26 = sc.read_h5ad("../5-wu/amlP26/outs/AMLP26_cell_filtered.h5ad")
amlP27 = sc.read_h5ad("../5-wu/amlP27/outs/AMLP27_cell_filtered.h5ad")
amlP28 = sc.read_h5ad("../5-wu/amlP28/outs/AMLP28_cell_filtered.h5ad")
amlP29 = sc.read_h5ad("../5-wu/amlP29/outs/AMLP29_cell_filtered.h5ad")
amlP30 = sc.read_h5ad("../5-wu/amlP30/outs/AMLP30_cell_filtered.h5ad")

amlP31 = sc.read_h5ad("../5-wu/amlP31/outs/AMLP31_cell_filtered.h5ad")
amlP32 = sc.read_h5ad("../5-wu/amlP32/outs/AMLP32_cell_filtered.h5ad")
amlP33 = sc.read_h5ad("../5-wu/amlP33/outs/AMLP33_cell_filtered.h5ad")
amlP34 = sc.read_h5ad("../5-wu/amlP34/outs/AMLP34_cell_filtered.h5ad")
amlP35 = sc.read_h5ad("../5-wu/amlP35/outs/AMLP35_cell_filtered.h5ad")

amlP36 = sc.read_h5ad("../5-wu/amlP36/outs/AMLP36_cell_filtered.h5ad")
amlP37 = sc.read_h5ad("../5-wu/amlP37/outs/AMLP37_cell_filtered.h5ad")
amlP38 = sc.read_h5ad("../5-wu/amlP38/outs/AMLP38_cell_filtered.h5ad")
amlP39 = sc.read_h5ad("../5-wu/amlP39/outs/AMLP39_cell_filtered.h5ad")
amlP40 = sc.read_h5ad("../5-wu/amlP40/outs/AMLP40_cell_filtered.h5ad")

amlP41 = sc.read_h5ad("../5-wu/amlP41/outs/AMLP41_cell_filtered.h5ad")
amlP42 = sc.read_h5ad("../5-wu/amlP42/outs/AMLP42_cell_filtered.h5ad")

data = [amlP01,amlP02,amlP03,amlP04,amlP05,amlP06,amlP07,amlP08,amlP09,amlP10,
        amlP11,amlP12,amlP13,amlP14,amlP15,amlP16,amlP17,amlP18,amlP19,amlP20,amlP21,
        amlP22,amlP23,amlP24,amlP25,amlP26,amlP27,amlP28,amlP29,amlP30,amlP31,amlP32,amlP33,
        amlP34,amlP35,amlP36,amlP37,amlP38,amlP39,amlP40,amlP41,amlP42]
adata = ad.concat(data, merge="same")
print(adata)

print("Zhai")
#zhai - cells already filtered
amlS220 = sc.read_h5ad("../6-zhai/amlS220/outs/AMLS220_unfiltered.h5ad")
amlS2275 = sc.read_h5ad("../6-zhai/amlS2275/outs/AMLS2275_unfiltered.h5ad")
amlS232 = sc.read_h5ad("../6-zhai/amlS232/outs/AMLS232_unfiltered.h5ad")
amlS292 = sc.read_h5ad("../6-zhai/amlS292/outs/AMLS292_unfiltered.h5ad")
amlS3432 = sc.read_h5ad("../6-zhai/amlS3432/outs/AMLS3432_unfiltered.h5ad")
amlS914 = sc.read_h5ad("../6-zhai/amlS914/outs/AMLS914_unfiltered.h5ad")

data = [amlS220,amlS2275,amlS232,amlS292,amlS3432,amlS914]
adata = ad.concat(data, merge="same")
print(adata)

print("Velten")
#velten - cells already filtered
amlV01 = sc.read_h5ad("../7-velten/amlV01/outs/AMLV01_unfiltered.h5ad")
amlV02 = sc.read_h5ad("../7-velten/amlV02/outs/AMLV02_unfiltered.h5ad")
amlV03 = sc.read_h5ad("../7-velten/amlV03/outs/AMLV03_unfiltered.h5ad")
amlV04 = sc.read_h5ad("../7-velten/amlV04/outs/AMLV04_unfiltered.h5ad")

data = [amlV01,amlV02,amlV03,amlV04]
adata = ad.concat(data, merge="same")
print(adata)

print("Jiang")
#jiang - cells already filtered
aml013 = sc.read_h5ad("../8-jiang/aml013/outs/AML013_unfiltered.h5ad")
aml016 = sc.read_h5ad("../8-jiang/aml016/outs/AML016_unfiltered.h5ad")
aml049 = sc.read_h5ad("../8-jiang/aml049/outs/AML049_unfiltered.h5ad")
aml060 = sc.read_h5ad("../8-jiang/aml060/outs/AML060_unfiltered.h5ad")
aml068 = sc.read_h5ad("../8-jiang/aml068/outs/AML068_unfiltered.h5ad")
aml070 = sc.read_h5ad("../8-jiang/aml070/outs/AML070_unfiltered.h5ad")
aml072 = sc.read_h5ad("../8-jiang/aml072/outs/AML072_unfiltered.h5ad")
aml076 = sc.read_h5ad("../8-jiang/aml076/outs/AML076_unfiltered.h5ad")
aml101 = sc.read_h5ad("../8-jiang/aml101/outs/AML101_unfiltered.h5ad")

data = [aml013,aml016,aml049,aml060,aml068,aml070,aml072,aml076,aml101]
adata = ad.concat(data, merge="same")
print(adata)

print("Johnston")
#johnston - cell filtering performed
aml7809 = sc.read_h5ad("../9-johnston/aml7809/outs/AML7809_cell_filtered.h5ad")
aml7810 = sc.read_h5ad("../9-johnston/aml7810/outs/AML7810_cell_filtered.h5ad")

data = [aml7809,aml7810]
adata = ad.concat(data, merge="same")
print(adata)

print("Beneyto-Calabuig")
#Beneyto-Calabuig - cells already filtered for A samples, filtering performed for B samples
amlA01 = sc.read_h5ad("../10-beneyto-calabuig/amlA01/outs/AMLA01_unfiltered.h5ad")
amlA02 = sc.read_h5ad("../10-beneyto-calabuig/amlA02/outs/AMLA02_unfiltered.h5ad")
amlA03 = sc.read_h5ad("../10-beneyto-calabuig/amlA03/outs/AMLA03_unfiltered.h5ad")
amlA04 = sc.read_h5ad("../10-beneyto-calabuig/amlA04/outs/AMLA04_unfiltered.h5ad")

amlA05 = sc.read_h5ad("../10-beneyto-calabuig/amlA05/outs/AMLA05_unfiltered.h5ad")
amlA06 = sc.read_h5ad("../10-beneyto-calabuig/amlA06/outs/AMLA06_unfiltered.h5ad")
amlA07 = sc.read_h5ad("../10-beneyto-calabuig/amlA07/outs/AMLA07_unfiltered.h5ad")
amlA08 = sc.read_h5ad("../10-beneyto-calabuig/amlA08/outs/AMLA08_unfiltered.h5ad")
amlA09 = sc.read_h5ad("../10-beneyto-calabuig/amlA09/outs/AMLA09_unfiltered.h5ad")

amlA10 = sc.read_h5ad("../10-beneyto-calabuig/amlA10/outs/AMLA10_unfiltered.h5ad")
amlA11 = sc.read_h5ad("../10-beneyto-calabuig/amlA11/outs/AMLA11_unfiltered.h5ad")
amlA12 = sc.read_h5ad("../10-beneyto-calabuig/amlA12/outs/AMLA12_unfiltered.h5ad")
amlA13 = sc.read_h5ad("../10-beneyto-calabuig/amlA13/outs/AMLA13_unfiltered.h5ad")
amlA14 = sc.read_h5ad("../10-beneyto-calabuig/amlA14/outs/AMLA14_unfiltered.h5ad")
amlA15 = sc.read_h5ad("../10-beneyto-calabuig/amlA15/outs/AMLA15_unfiltered.h5ad")

amlB01 = sc.read_h5ad("../10-beneyto-calabuig/amlB01/outs/AMLB01_cell_filtered.h5ad")
amlB02 = sc.read_h5ad("../10-beneyto-calabuig/amlB02/outs/AMLB02_cell_filtered.h5ad")
amlB03 = sc.read_h5ad("../10-beneyto-calabuig/amlB03/outs/AMLB03_cell_filtered.h5ad")
amlB04 = sc.read_h5ad("../10-beneyto-calabuig/amlB04/outs/AMLB04_cell_filtered.h5ad")

data = [amlA01,amlA02,amlA03,amlA04,amlA05,amlA06,amlA07,amlA08,amlA09,amlA10,amlA11,amlA12,amlA13,amlA14,amlA15,amlB01,amlB02,amlB03,amlB04]
adata = ad.concat(data, merge="same")
print(adata)

print("Li")
#li - cell filtering performed
amlPT03 = sc.read_h5ad("../11-li/amlPT03/outs/AMLPT03_cell_filtered.h5ad")
amlPT08 = sc.read_h5ad("../11-li/amlPT08/outs/AMLPT08_cell_filtered.h5ad")
amlPT09 = sc.read_h5ad("../11-li/amlPT09/outs/AMLPT09_cell_filtered.h5ad")
amlPT10 = sc.read_h5ad("../11-li/amlPT10/outs/AMLPT10_cell_filtered.h5ad")

data = [amlPT03,amlPT08,amlPT09,amlPT10]
adata = ad.concat(data, merge="same")
print(adata)

print("Fiskus")
#fiskus - cell filtering performed
aml8799 = sc.read_h5ad("../12-fiskus/aml8799/outs/AML8799_cell_filtered.h5ad")

data = [aml8799]
adata = ad.concat(data, merge="same")
print(adata)

print("Zhang")
#zhang - cell filtering performed
amlP105 = sc.read_h5ad("../13-zhang/amlP105/outs/AMLP105_cell_filtered.h5ad")
amlP106 = sc.read_h5ad("../13-zhang/amlP106/outs/AMLP106_cell_filtered.h5ad")
amlP108 = sc.read_h5ad("../13-zhang/amlP108/outs/AMLP108_cell_filtered.h5ad")
amlP114 = sc.read_h5ad("../13-zhang/amlP114/outs/AMLP114_cell_filtered.h5ad")
amlP115 = sc.read_h5ad("../13-zhang/amlP115/outs/AMLP115_cell_filtered.h5ad")
amlP116 = sc.read_h5ad("../13-zhang/amlP116/outs/AMLP116_cell_filtered.h5ad")
amlP117 = sc.read_h5ad("../13-zhang/amlP117/outs/AMLP117_cell_filtered.h5ad")
amlP118 = sc.read_h5ad("../13-zhang/amlP118/outs/AMLP118_cell_filtered.h5ad")
amlP119 = sc.read_h5ad("../13-zhang/amlP119/outs/AMLP119_cell_filtered.h5ad")
amlP120 = sc.read_h5ad("../13-zhang/amlP120/outs/AMLP120_cell_filtered.h5ad")
amlP122 = sc.read_h5ad("../13-zhang/amlP122/outs/AMLP122_cell_filtered.h5ad")
amlP123 = sc.read_h5ad("../13-zhang/amlP123/outs/AMLP123_cell_filtered.h5ad")
amlP124 = sc.read_h5ad("../13-zhang/amlP124/outs/AMLP124_cell_filtered.h5ad")

data = [amlP105,amlP106,amlP108,amlP114,amlP115,amlP116,amlP117,amlP118,amlP119,amlP120,amlP122,amlP123,amlP124]
adata = ad.concat(data, merge="same")
print(adata)

print("Mumme")
#mumme - cell filtering performed
aml01D = sc.read_h5ad("../14-mumme/aml01D/outs/AML01D_cell_filtered.h5ad")
aml02D = sc.read_h5ad("../14-mumme/aml02D/outs/AML02D_cell_filtered.h5ad")
aml03D = sc.read_h5ad("../14-mumme/aml03D/outs/AML03D_cell_filtered.h5ad")
aml04D = sc.read_h5ad("../14-mumme/aml04D/outs/AML04D_cell_filtered.h5ad")
aml05D = sc.read_h5ad("../14-mumme/aml05D/outs/AML05D_cell_filtered.h5ad")
aml06D = sc.read_h5ad("../14-mumme/aml06D/outs/AML06D_cell_filtered.h5ad")
aml07D = sc.read_h5ad("../14-mumme/aml07D/outs/AML07D_cell_filtered.h5ad")
aml08D = sc.read_h5ad("../14-mumme/aml08D/outs/AML08D_cell_filtered.h5ad")
aml09D = sc.read_h5ad("../14-mumme/aml09D/outs/AML09D_cell_filtered.h5ad")
aml10D = sc.read_h5ad("../14-mumme/aml10D/outs/AML10D_cell_filtered.h5ad")

aml11D = sc.read_h5ad("../14-mumme/aml11D/outs/AML11D_cell_filtered.h5ad")
aml12D = sc.read_h5ad("../14-mumme/aml12D/outs/AML12D_cell_filtered.h5ad")
aml13D = sc.read_h5ad("../14-mumme/aml13D/outs/AML13D_cell_filtered.h5ad")
aml14D = sc.read_h5ad("../14-mumme/aml14D/outs/AML14D_cell_filtered.h5ad")
aml16D = sc.read_h5ad("../14-mumme/aml16D/outs/AML16D_cell_filtered.h5ad")
aml17D = sc.read_h5ad("../14-mumme/aml17D/outs/AML17D_cell_filtered.h5ad")
aml18D = sc.read_h5ad("../14-mumme/aml18D/outs/AML18D_cell_filtered.h5ad")
aml19D = sc.read_h5ad("../14-mumme/aml19D/outs/AML19D_cell_filtered.h5ad")
aml20D = sc.read_h5ad("../14-mumme/aml20D/outs/AML20D_cell_filtered.h5ad")

data = [aml01D,aml02D,aml03D,aml04D,aml05D,aml06D,aml07D,aml08D,aml09D,aml10D,aml11D,aml12D,aml13D,aml14D,aml16D,aml17D,aml18D,aml19D,aml20D]
adata = ad.concat(data, merge="same")
print(adata)

print("Pei")
# pei
aml43363 = sc.read_h5ad("../15-pei/aml43363/outs/AML43363_cell_filtered.h5ad")

data = [aml43363]
adata = ad.concat(data, merge="same")
print(adata)

print("Naldini")
# naldini
amlM03 = sc.read_h5ad("../16-naldini/amlM03/outs/AMLM03_cell_filtered.h5ad")
amlM04 = sc.read_h5ad("../16-naldini/amlM04/outs/AMLM04_cell_filtered.h5ad")
amlM07 = sc.read_h5ad("../16-naldini/amlM07/outs/AMLM07_cell_filtered.h5ad")
amlM08 = sc.read_h5ad("../16-naldini/amlM08/outs/AMLM08_cell_filtered.h5ad")

amlM21 = sc.read_h5ad("../16-naldini/amlM21/outs/AMLM21_cell_filtered.h5ad")
amlM22 = sc.read_h5ad("../16-naldini/amlM22/outs/AMLM22_cell_filtered.h5ad")
amlM25 = sc.read_h5ad("../16-naldini/amlM25/outs/AMLM25_cell_filtered.h5ad")
amlM28 = sc.read_h5ad("../16-naldini/amlM28/outs/AMLM28_cell_filtered.h5ad")

amlM32 = sc.read_h5ad("../16-naldini/amlM32/outs/AMLM32_cell_filtered.h5ad")
amlM33 = sc.read_h5ad("../16-naldini/amlM33/outs/AMLM33_cell_filtered.h5ad")
amlM36 = sc.read_h5ad("../16-naldini/amlM36/outs/AMLM36_cell_filtered.h5ad")
amlM40 = sc.read_h5ad("../16-naldini/amlM40/outs/AMLM40_cell_filtered.h5ad")

amlM41 = sc.read_h5ad("../16-naldini/amlM41/outs/AMLM41_cell_filtered.h5ad")
amlM48 = sc.read_h5ad("../16-naldini/amlM48/outs/AMLM48_cell_filtered.h5ad")
amlM84 = sc.read_h5ad("../16-naldini/amlM84/outs/AMLM84_cell_filtered.h5ad")
amlM88 = sc.read_h5ad("../16-naldini/amlM88/outs/AMLM88_cell_filtered.h5ad")

data = [amlM03,amlM04,amlM07,amlM08,amlM21,amlM22,amlM25,amlM28,amlM32,amlM33,amlM36,amlM40,amlM41,amlM48,amlM84,amlM88]
adata = ad.concat(data, merge="same")
print(adata)

print("Lasry")
#lasry - some additional cell filtering performed
AMLL001 = sc.read_h5ad("../17-lasry/AML001/outs/AML001_cell_filtered.h5ad")
AMLL0024 = sc.read_h5ad("../17-lasry/AML0024/outs/AML0024_cell_filtered.h5ad")
AMLL0048 = sc.read_h5ad("../17-lasry/AML0048/outs/AML0048_cell_filtered.h5ad")
AMLL0102 = sc.read_h5ad("../17-lasry/AML0102/outs/AML0102_cell_filtered.h5ad")

AMLL012_etc = sc.read_h5ad("../17-lasry/AML012_AML005_AML055_AML048/outs/AML012_AML005_AML055_AML048_cell_filtered.h5ad")
AMLL026_AMLL051 = sc.read_h5ad("../17-lasry/AML026_AML051/outs/AML026_AML051_cell_filtered.h5ad")
AMLL038_etc = sc.read_h5ad("../17-lasry/AML038_AML008_AML043_AML028_AML056/outs/AML038_AML008_AML043_AML028_AML056_cell_filtered.h5ad")
AMLL0612_etc = sc.read_h5ad("../17-lasry/AML0612_AML3762_AML0160_AML0310_AML3133/outs/AML0612_AML3762_AML0160_AML0310_AML3133_cell_filtered.h5ad")

AMLL0693 = sc.read_h5ad("../17-lasry/AML0693/outs/AML0693_cell_filtered.h5ad")
AMLL073_etc = sc.read_h5ad("../17-lasry/AML073_AML006_AML025_AML003/outs/AML073_AML006_AML025_AML003_cell_filtered.h5ad")
AMLL1371_AMLL4340 = sc.read_h5ad("../17-lasry/AML1371_AML4340/outs/AML1371_AML4340_cell_filtered.h5ad")
AMLL2123 = sc.read_h5ad("../17-lasry/AML2123/outs/AML2123_cell_filtered.h5ad")
AMLL2451_Control0004 = sc.read_h5ad("../17-lasry/AML2451_Control0004/outs/AML2451_Control0004_cell_filtered.h5ad")

AMLL3266 = sc.read_h5ad("../17-lasry/AML3266/outs/AML3266_cell_filtered.h5ad")
AMLL3730 = sc.read_h5ad("../17-lasry/AML3730/outs/AML3730_cell_filtered.h5ad")
AMLL3948 = sc.read_h5ad("../17-lasry/AML3948/outs/AML3948_cell_filtered.h5ad")
AMLL4897 = sc.read_h5ad("../17-lasry/AML4897/outs/AML4897_cell_filtered.h5ad")

Control0005_AMLL009 = sc.read_h5ad("../17-lasry/Control0005_AML009/outs/Control0005_AML009_cell_filtered.h5ad")
Control0058_PAXMIJ_PAUMTZ = sc.read_h5ad("../17-lasry/Control0058_PAXMIJ_PAUMTZ/outs/Control0058_PAXMIJ_PAUMTZ_cell_filtered.h5ad")
Control0082_etc = sc.read_h5ad("../17-lasry/Control0082_AML052_AML022/outs/Control0082_AML052_AML022_cell_filtered.h5ad")
Control5_etc = sc.read_h5ad("../17-lasry/Control5_AML2910_AML3050_AML0361/outs/Control5_AML2910_AML3050_AML0361_cell_filtered.h5ad")
PAWWEE = sc.read_h5ad("../17-lasry/PAWWEE/outs/PAWWEE_cell_filtered.h5ad")

data = [AMLL001,AMLL0024,AMLL0048,AMLL0102,AMLL012_etc,AMLL026_AMLL051,AMLL038_etc,AMLL0612_etc,AMLL0693,AMLL073_etc,AMLL1371_AMLL4340,AMLL2123,AMLL2451_Control0004,AMLL3266,AMLL3730,AMLL3948,AMLL4897,Control0005_AMLL009,Control0058_PAXMIJ_PAUMTZ,Control0082_etc,Control5_etc,PAWWEE]
adata = ad.concat(data, merge="same")
print(adata)


# Healthy Donor Samples

print("Oetjen")
#oetjen - cell filtering performed
bm6161 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6161/outs/BM6161_cell_filtered.h5ad")
bm6162 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6162/outs/BM6162_cell_filtered.h5ad")
bm6163 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6163/outs/BM6163_cell_filtered.h5ad")
bm6164 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6164/outs/BM6164_cell_filtered.h5ad")
bm6165 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6165/outs/BM6165_cell_filtered.h5ad")
bm6166 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6166/outs/BM6166_cell_filtered.h5ad")

bm6167 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6167/outs/BM6167_cell_filtered.h5ad")
bm6168 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6168/outs/BM6168_cell_filtered.h5ad")
bm6169 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6169/outs/BM6169_cell_filtered.h5ad")
bm6170 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6170/outs/BM6170_cell_filtered.h5ad")
bm6171 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6171/outs/BM6171_cell_filtered.h5ad")
bm6172 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6172/outs/BM6172_cell_filtered.h5ad")

bm6173 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6173/outs/BM6173_cell_filtered.h5ad")
bm6174 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6174/outs/BM6174_cell_filtered.h5ad")
bm6175 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6175/outs/BM6175_cell_filtered.h5ad")
bm6176 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6176/outs/BM6176_cell_filtered.h5ad")
bm6177 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6177/outs/BM6177_cell_filtered.h5ad")
bm6178 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6178/outs/BM6178_cell_filtered.h5ad")

bm6179 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6179/outs/BM6179_cell_filtered.h5ad")
bm6180 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6180/outs/BM6180_cell_filtered.h5ad")
bm6181 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6181/outs/BM6181_cell_filtered.h5ad")
bm6182 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6182/outs/BM6182_cell_filtered.h5ad")
bm6183 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6183/outs/BM6183_cell_filtered.h5ad")
bm6184 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6184/outs/BM6184_cell_filtered.h5ad")
bm6185 = sc.read_h5ad("../healthy_donors/1-oetjen/bm6185/outs/BM6185_cell_filtered.h5ad")

data = [bm6161,bm6162,bm6163,bm6164,bm6165,bm6166,bm6167,bm6168,bm6169,bm6170,bm6171,bm6172,bm6173,bm6174,
        bm6175,bm6176,bm6177,bm6178,bm6179,bm6180,bm6181,bm6182,bm6183,bm6184,bm6185]
adata = ad.concat(data, merge="same")
print(adata)

print("Wu")
#wu - cell filtering performed
n01 = sc.read_h5ad("../healthy_donors/2-wu/N01/outs/N01_cell_filtered.h5ad")
n02 = sc.read_h5ad("../healthy_donors/2-wu/N02/outs/N02_cell_filtered.h5ad")

data = [n01,n02]
adata = ad.concat(data, merge="same")
print(adata)

print("Van Galen")
#van galen - cell filtering NOT performed
bm1 = sc.read_h5ad("../healthy_donors/3-van-galen/BM1/outs/BM1_unfiltered.h5ad")
bm2 = sc.read_h5ad("../healthy_donors/3-van-galen/BM2/outs/BM2_unfiltered.h5ad")
bm3 = sc.read_h5ad("../healthy_donors/3-van-galen/BM3/outs/BM3_cell_filtered.h5ad")
bm4 = sc.read_h5ad("../healthy_donors/3-van-galen/BM4/outs/BM4_cell_filtered.h5ad")
bm5a = sc.read_h5ad("../healthy_donors/3-van-galen/BM5a/outs/BM5a_unfiltered.h5ad")
bm5b = sc.read_h5ad("../healthy_donors/3-van-galen/BM5b/outs/BM5b_unfiltered.h5ad")

data = [bm1,bm2,bm3,bm4,bm5a,bm5b]
adata = ad.concat(data, merge="same")
print(adata)

print("Setty")
#setty HCA data - cell filtering NOT performed
setty = sc.read_h5ad("../healthy_donors/4-setty/outs/setty_unfiltered.h5ad")
print(setty)

print("Beneyto-Calabuig")
# ben-c - cell filtering not needed
bmA00 = sc.read_h5ad("../healthy_donors/5-ben-c/bmA00/outs/BMA00_unfiltered.h5ad")
print(bmA00)

print("Human Cell Atlas (Immune)")
# hca - cell filtering performed
hca = sc.read_h5ad("../healthy_donors/6-hca/outs/hca_cell_filtered.h5ad")
print(hca)

print("Zhang")
# zhang - cell filtering performed
N004 = sc.read_h5ad("../healthy_donors/7-zhang/N004/outs/N004_cell_filtered.h5ad")
print(N004)

print("Caron")
# caron - cell filtering performed
PBMMC1 = sc.read_h5ad("../healthy_donors/8-caron/PBMMC1/outs/PBMMC1_cell_filtered.h5ad")
PBMMC2 = sc.read_h5ad("../healthy_donors/8-caron/PBMMC2/outs/PBMMC2_cell_filtered.h5ad")
PBMMC3 = sc.read_h5ad("../healthy_donors/8-caron/PBMMC3/outs/PBMMC3_cell_filtered.h5ad")

data = [PBMMC1,PBMMC2,PBMMC3]
adata = ad.concat(data, merge="same")
print(adata)


# Merge All Samples

data = [aml508084,aml548327,aml721214,aml782328,aml809653,
        aml027,aml035,
        aml1012,aml210A,aml314,aml328,aml329,aml371,aml419A,aml420B,aml475,aml556,aml707B,aml722B,aml870,aml916,aml921A,aml997,
        aml681,aml682,aml683,aml684,aml685,
        amlP01,amlP02,amlP03,amlP04,amlP05,amlP06,amlP07,amlP08,amlP09,amlP10,
        amlP11,amlP12,amlP13,amlP14,amlP15,amlP16,amlP17,amlP18,amlP19,amlP20,amlP21,
        amlP22,amlP23,amlP24,amlP25,amlP26,amlP27,amlP28,amlP29,amlP30,amlP31,amlP32,amlP33,
        amlP34,amlP35,amlP36,amlP37,amlP38,amlP39,amlP40,amlP41,amlP42,
        amlS220,amlS2275,amlS232,amlS292,amlS3432,amlS914,
        amlV01,amlV02,amlV03,amlV04,
        aml013,aml016,aml049,aml060,aml068,aml070,aml072,aml076,aml101,
        aml7809,aml7810,
        amlA01,amlA02,amlA03,amlA04,amlA05,amlA06,amlA07,amlA08,amlA09,amlA10,amlA11,amlA12,
        amlA13,amlA14,amlA15,amlB01,amlB02,amlB03,amlB04,
        amlPT03,amlPT08,amlPT09,amlPT10,
        aml8799,
        amlP105,amlP106,amlP108,amlP114,amlP115,amlP116,amlP117,amlP118,amlP119,
        amlP120,amlP122,amlP123,amlP124,
        aml01D,aml02D,aml03D,aml04D,aml05D,aml06D,aml07D,aml08D,aml09D,aml10D,
        aml11D,aml12D,aml13D,aml14D,aml16D,aml17D,aml18D,aml19D,aml20D,
        aml43363,
        amlM03,amlM04,amlM07,amlM08,amlM21,amlM22,amlM25,amlM28,amlM32,amlM33,amlM36,amlM40,
        amlM41,amlM48,amlM84,amlM88,
        AMLL001,AMLL0024,AMLL0048,AMLL0102,AMLL012_etc,AMLL026_AMLL051,AMLL038_etc,AMLL0612_etc,
        AMLL0693,AMLL073_etc,AMLL1371_AMLL4340,AMLL2123,AMLL2451_Control0004,AMLL3266,AMLL3730,
        AMLL3948,AMLL4897,Control0005_AMLL009,Control0058_PAXMIJ_PAUMTZ,Control0082_etc,Control5_etc,PAWWEE,
        bm6161,bm6162,bm6163,bm6164,bm6165,bm6166,bm6167,bm6168,bm6169,bm6170,
        bm6171,bm6172,bm6173,bm6174,bm6175,bm6176,bm6177,bm6178,bm6179,bm6180,
        bm6181,bm6182,bm6183,bm6184,bm6185,
        n01,n02,
        bm1,bm2,bm3,bm4,bm5a,bm5b,
        setty,
        bmA00,
        hca,
        N004,
        PBMMC1,PBMMC2,PBMMC3]


#concetenate files
adata = ad.concat(data, join="outer", merge="same")
print(adata)


del(aml508084,aml548327,aml721214,aml782328,aml809653,
        aml027,aml035,
        aml1012,aml210A,aml314,aml328,aml329,aml371,aml419A,aml420B,aml475,aml556,aml707B,aml722B,aml870,aml916,aml921A,aml997,
        aml681,aml682,aml683,aml684,aml685,
        amlP01,amlP02,amlP03,amlP04,amlP05,amlP06,amlP07,amlP08,amlP09,amlP10,
        amlP11,amlP12,amlP13,amlP14,amlP15,amlP16,amlP17,amlP18,amlP19,amlP20,amlP21,
        amlP22,amlP23,amlP24,amlP25,amlP26,amlP27,amlP28,amlP29,amlP30,amlP31,amlP32,amlP33,
        amlP34,amlP35,amlP36,amlP37,amlP38,amlP39,amlP40,amlP41,amlP42,
        amlS220,amlS2275,amlS232,amlS292,amlS3432,amlS914,
        amlV01,amlV02,amlV03,amlV04,
        aml013,aml016,aml049,aml060,aml068,aml070,aml072,aml076,aml101,
        aml7809,aml7810,
        amlA01,amlA02,amlA03,amlA04,amlA05,amlA06,amlA07,amlA08,amlA09,amlA10,amlA11,amlA12,
        amlA13,amlA14,amlA15,amlB01,amlB02,amlB03,amlB04,
        amlPT03,amlPT08,amlPT09,amlPT10,
        aml8799,
        amlP105,amlP106,amlP108,amlP114,amlP115,amlP116,amlP117,amlP118,amlP119,
        amlP120,amlP122,amlP123,amlP124,
        aml01D,aml02D,aml03D,aml04D,aml05D,aml06D,aml07D,aml08D,aml09D,aml10D,
        aml11D,aml12D,aml13D,aml14D,aml16D,aml17D,aml18D,aml19D,aml20D,
        aml43363,
        amlM03,amlM04,amlM07,amlM08,amlM21,amlM22,amlM25,amlM28,amlM32,amlM33,amlM36,amlM40,
        amlM41,amlM48,amlM84,amlM88,
        AMLL001,AMLL0024,AMLL0048,AMLL0102,AMLL012_etc,AMLL026_AMLL051,AMLL038_etc,AMLL0612_etc,
        AMLL0693,AMLL073_etc,AMLL1371_AMLL4340,AMLL2123,AMLL2451_Control0004,AMLL3266,AMLL3730,
        AMLL3948,AMLL4897,Control0005_AMLL009,Control0058_PAXMIJ_PAUMTZ,Control0082_etc,Control5_etc,PAWWEE,
        bm6161,bm6162,bm6163,bm6164,bm6165,bm6166,bm6167,bm6168,bm6169,bm6170,
        bm6171,bm6172,bm6173,bm6174,bm6175,bm6176,bm6177,bm6178,bm6179,bm6180,
        bm6181,bm6182,bm6183,bm6184,bm6185,
        n01,n02,
        bm1,bm2,bm3,bm4,bm5a,bm5b,
        setty,
        bmA00,
        hca,
        N004,
        PBMMC1,PBMMC2,PBMMC3)


print("Data Combined")
print(adata.n_obs, adata.n_vars)
adata.obs = adata.obs.drop(columns=['doublet_score','predicted_doublets','doublet_info','van-galen-celltype',
                                   'samples', 'Broad_cell_identity', 'Cell_type_identity', 'run', "outlier", "mt_outlier"])

adata.write_h5ad(results_file, compression="gzip")

#################################################################################################################################
################################################### Add Cell Type Annotations ###################################################
#################################################################################################################################


print("Add Celltype Annotations")
# Include as original labels and as harmonised labels across datasets

celltypes = pd.read_csv('./original_celltypes.csv', index_col=0, header=None)
original_celltypes = celltypes

celltypes = celltypes.replace({"B":"B Cells", "cDC":"Dendritic Cells", "cDC-like":"Dendritic Cells", 
                                         "CTL":"T Cells", "earlyEry":"Erythroid Progenitor Cells", "GMP-like":"GMP", 
                                         "HSC-like":"HSC", "lateEry":"Erythroid Cells", "Mono":"Monocytes", 
                                         "Mono-like":"Monocytes", "NK":"NK Cells","pDC":"Plasmacytoid Dendritic Cells",
                                         "Plasma":"Plasma Cells", "ProB":"B Cell Precursors", "Prog":"HSC/MPPs", 
                                         "Prog-like":"HSC/MPPs", "ProMono":"Monocyte Progenitor Cells", 
                                         "ProMono-like":"Monocyte Progenitor Cells", "T":"T Cells"}) 
#update van galen labels

celltypes = celltypes.replace({"BASO":"Granulocytes", "B-CELL":"B Cells", "DEND(L)":"Dendritic Cells", 
                               "DEND(M)":"Dendritic Cells", "EOS":"Granulocytes", "ERY":"Erythroid Cells", 
                               "ERY(CD34+)":"Erythroid Progenitor Cells", "GRAN": "Granulocytes", 
                               "MEGA":"Megakaryocytes", "MONO":"Monocytes", "NK":"NK Cells", "NKT":"NK T Cells", 
                               "PRE-B-CELL":"B Cell Precursors", "T-CELL":"T Cells"
}) 
#update petti labels

celltypes = celltypes.replace({"Myeloid-AM":"Unknown", "Bcell-B cells":"B Cells", "Myeloid-CD34+CD117bri":"CD34+ Blasts", 
                               "Myeloid-CD34+CD117bri-G2M":"CD34+ Blasts", "Myeloid-CD34+CD117bri-S":"CD34+ Blasts", 
                               "Myeloid-CD34+CD117dim":"CD34+ Blasts", "Myeloid-DC":"Dendritic Cells", 
                               "IMERY-Immature Ery":"Erythroid Progenitor Cells", 
                               "Myeloid-Monocyte Precursors":"Monocyte Progenitor Cells",
                               "Tcell-T cells":"T Cells", "Myeloid-Monocytes":"Monocytes"

}) 
#update jiang labels

celltypes = celltypes.replace({"Aberrant erythroid":"Erythroid Cells", "CD11c+ memory B cells":"B Cells", 
                               "CD4+ cytotoxic T cells":"CD4+ T Cells", "CD4+ memory T cells":"CD4+ T Cells", 
                               "CD4+ naive T cells":"CD4+ T Cells", "CD56brightCD16- NK cells":"NK Cells", 
                               "CD56dimCD16+ NK cells":"NK Cells", "CD69+PD-1+ memory CD4+ T cells":"CD4+ T Cells", 
                               "CD8+CD103+ tissue resident memory T cells":"CD8+ T Cells", 
                               "CD8+ central memory T cells":"CD8+ T Cells", 
                               "CD8+ effector memory T cells":"CD8+ T Cells", "CD8+ naive T cells":"CD8+ T Cells", 
                               "Classical Monocytes":"Monocytes", "Class switched memory B cells":"B Cells", 
                               "Conventional dendritic cell 1":"Dendritic Cells", 
                               "Conventional dendritic cell 2":"Dendritic Cells", 
                               "Early erythroid progenitor":"Erythroid Progenitor Cells", 
                               "Early promyelocytes":"Promyelocytes", 
                               "Eosinophil-basophil-mast cell progenitors":"EBMP", 
                               "Erythro-myeloid progenitor":"CMP", "GammaDelta T cells":"T Cells", 
                               "HSCs & MPPs":"HSC/MPPs","Immature B cells":"B Cell Precursors", 
                               "Late erythroid progenitor":"Erythroid Cells", "Late promyelocytes":"Promyelocytes", 
                               "Lymphomyeloid prog":"Lympho-Myeloid Progenitor Cells", 
                               "Mature naive B cells":"B Cells", "NK cell progenitors":"NK Cells", 
                               "Non-classical monocytes":"Monocytes", "Nonswitched memory B cells":"B Cells", 
                               "Plasmacytoid dendritic cell progenitors":"Plasmacytoid Dendritic Cells", 
                               "Pre-B cells":"B Cell Precursors", "Pre-pro-B cells":"B Cell Precursors", 
                               "Pro-B cells":"B Cell Precursors", "Small pre-B cell":"B Cell Precursors"


}) 
#update BC labels

celltypes = celltypes.replace({"Bcellsandprogenitors":"B Cells", "CD34+Blasts":"CD34+ Blasts", 
                               "CD34+BlastsandHSPCs":"CD34+ Blasts", "CD34-Blasts(Calprotectin-AZU1+)":"CD34- Blasts",
                               "CD34-Blasts(Calprotectin+AZU1-)":"CD34- Blasts", 
                               "CD34-Blasts(Calprotectin+AZU1+)":"CD34- Blasts", 
                               "CD34-Blasts(Intermediate)":"CD34- Blasts", "CD34-Blasts(Unclear)":"Unknown", 
                               "CD34+HBZ+Blasts":"CD34+ Blasts", "CentralmemoryT-cells":"T Cells", 
                               "CytotoxicT-cells":"T Cells", "EffectormemoryT-cells":"T Cells", 
                               "Erythroidprecursors":"Erythroid Progenitor Cells","MitoticHSPCs(G2/M)":"HSC/MPPs", 
                               "Neutrophilprecursors":"Granulocytes", "NKandTcells":"Unknown", 
                               "OtherT/Nkcells":"Unknown"



}) 
#update velten labels

celltypes = celltypes.replace({"Dendtiric progenitor cell":"Dendritic Progenitor Cell", 
                               "MK progenitor":"Megakaryocyte Progenitor Cells"
}) 
#update setty labels

celltypes = celltypes.replace({"B cell T cell doublet":"Lymphocyte", 
                               "CD14+ monocyte type 1":"Monocytes", "CD14+ monocyte type 2":"Monocytes", 
                               "CD16+ monocyte":"Monocytes", "CD4+ naive T cell":"CD4+ T Cells", 
                               "conventional dendritic cell":"Dendritic Cells", 
                               "cytotoxic T cell type 1":"T Cells", "cytotoxic T cell type 2":"T Cells", 
                               "erythroid cell type 1":"Erythroid Cells", "erythroid cell type 2":"Erythroid Cells",
                               "hematopoietic stem cell":"HSC", "megakaryocyte":"Megakaryocytes", 
                               "memory B cell":"B Cells", "mesenchymal stem cell":"Stromal Cells", 
                               "naive B cell":"B Cells", "naive CD8+ T cell":"CD8+ T Cells", 
                               "natural killer cell":"NK Cells", "plasma cell":"Plasma Cells", 
                               "plasmacytoid dendritic cell":"Plasmacytoid Dendritic Cells", 
                               "precursor B cell":"B Cell Precursors", "pro-B cell":"B Cell Precursors",
                               "T-helper cell":"T Cells"
}) 
#update hca labels

celltypes = celltypes.replace({"Mesenchymal cells_1":"Stromal Cells", 
                               "NK T cells":"NK T Cells",
                               "OtherT/NKcells":"Unknown",
                               "Monocyte-like blasts":"Unknown",
                               "Erythroid progenitor":"Erythroid Progenitor Cells",
                               "NKcells":"NK Cells",
                               "Plasmacytoid dendritic cells":"Plasmacytoid Dendritic Cells",
                               "Megakaryocyte progenitors":"Megakaryocyte Progenitor Cells",
                               "Erythro-myeloid progenitors":"Erythro-Myeloid Progenitors",
                               "Plasma cells":"Plasma Cells"
                               
}) 
# update zhang labels

celltypes = celltypes.replace({"proB-like":"B Cell Precursors", "erythrocyte-like":"Erythroid Cells",
                               "MK":"Megakaryocytes", "macrophage":"Macrophage", "stromal":"Stromal Cells",
                               "proB":"B Cell Precursors", "CLP-like":"CLP", "E/B/M":"EBMP",
                               "pDC":"Plasmacytoid Dendritic Cells","pDC-like":"Plasmacytoid Dendritic Cells",
                               "cDC":"Dendritic Cells","MEP-like":"MEP","neutrophil-like":"Granulocytes",
                               "CTL":"T Cells","monocyte-like":"Monocytes","naiveT":"T Cells","LMPP-like":"LMPP",
                               "HSC-like":"HSC","GMP-like":"GMP", "plasmaB":"Plasma Cells","erythrocyte":"Erythroid Cells",
                               "monocyte":"Monocytes","E/B/M-like":"EBMP", "B":"B Cells", "cDC-like":"Dendritic Cells",
})
#update missed labels/fix capitalisation


celltypes = celltypes.replace({"Perivascular cell":"Stromal Cells","Megakaryocyte":"Megakaryocytes",
                               "Plasma cell":"Plasma Cells","Plasmablast":"Plasma Cells",
                               "cDC1":"Conventional Dendritic Cells","cDC2":"Conventional Dendritic Cells",
                               "DC precursor":"Dendritic Progenitor Cell","Ery":"Erythroid Cells",
                               "Dendritic Cells":"Conventional Dendritic Cells","MAIT":"T Cells",
                               "CD8+ T":"CD8+ T Cells","CD4+ T":"CD4+ T Cells","gd T":"T Cells",
                               "Pre-B":"B Cell Precursors","Pro-B":"B Cell Precursors",
                               "Myelocytes":"Granulocytes","Granulocyte":"Granulocytes",
                               "Promyelocytes":"Granulocytes","HLA-II+ monocyte":"Monocytes",
                               "HSC":"HSC/MPPs","MPP":"HSC/MPPs","CD14+ monocyte":"Monocytes"
}) 


celltypes = celltypes.replace({"CD8+ IFN+":"T Cells", "CD4+ Activated":"T Cells",
                               "CD8+ NK-like":"T Cells", "Treg":"T Cells",
                               "CD8+ TRM":"T Cells", "CD8+ Naive":"T Cells", 
                               "CD8+ GZMK+":"T Cells", "CD4+ TCM":"T Cells",
                               "CD4+ Naive":"T Cells", "CD8+ Cytotoxic":"T Cells",
                               "CD8+ effector T cells":"T Cells", "CD16+ monocytes":"Monocytes",
                               "CD10+ B cells":"B Cells", "Natural killer cells":"NK Cells",
                               "Late erythroid progenitors":"Erythroid Cells",
                               "HSPCs":"HSC/MPPs", "Late erythrocytes":"Erythroid Cells",
                               "Early erythroid progenitors":"Erythroid Cells",
                               "Late erythroid progenitors":"Erythroid Cells",
                               "Early erythrocytes":"Erythroid Cells", "CD20+ B cells":"B Cells",
                               "Erythroid Progenitor Cells":"Erythroid Cells",
                               "CD14+ monocytes":"Monocytes", "CD8+ effector T cells":"T Cells",
                               "CD4+ T Cells":"T Cells", "CD8+ T Cells":"T Cells",
                               "Monocyte progenitors":"Monocyte Progenitor Cells",
                               "CD11c+":"Unknown"
}) 

print(celltypes[1].value_counts())

# Merge to AnnData Obs
annotations = pd.merge(original_celltypes, celltypes, how='left', left_index=True, right_index=True)

annotations = annotations.rename(columns={"1_x": "author_original_celltype", "1_y": "main_original_celltype"})

adata.obs = pd.merge(adata.obs, annotations, how="left", left_index=True, right_index=True)


adata.obs[["author_original_celltype","main_original_celltype"]] = adata.obs[["author_original_celltype","main_original_celltype"]].fillna("Unknown")
adata.obs["celltype"] = adata.obs["main_original_celltype"]


#################################################################################################################################
###################################################### Run PreFiltering QC ######################################################
#################################################################################################################################


print("Pre-Filtering QC")

#annotate mitochondrial/ribosomal/hb genes
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL"))
adata.var['hb'] = adata.var_names.str.contains(("^HB[^(P)]"))


#run qc function
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], 
                           inplace=True, percent_top=[20], log1p=True)


#find proportion of mito genes
mito_genes = adata.var_names.str.startswith('MT-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mt2'] = np.sum(
    adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)


# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1)


#plot qc
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo', 'pct_counts_hb'],
             jitter=0.4, groupby = 'sample', rotation= 45, save="_prefilt_initial_qc.png")
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save="_prefilt_mt_by_counts.png")
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color="pct_counts_mt", 
              save="_prefilt_nfeatures_by_counts.png")


#plot qc
fig = sns.displot(adata.obs["total_counts"], bins=100, kde=False).fig
fig.savefig("figures/hist_prefilt_total_counts.png") 

fig = sns.displot(adata.obs["n_genes_by_counts"], bins=100, kde=False).fig
fig.savefig("figures/hist_prefilt_ngenes_counts.png") 

#improve visibility for lower counts
sub = adata[adata.obs["total_counts"]<5000]
fig = sns.displot(sub.obs["total_counts"], bins=100, kde=False).fig
fig.savefig("figures/hist_prefilt_total_counts_lowcounts.png")

fig = sns.displot(sub.obs["n_genes_by_counts"], bins=100, kde=False).fig
fig.savefig("figures/hist_prefilt_ngenes_counts_lowcounts.png")

# Gene QC
sc.pl.highest_expr_genes(adata, n_top=20, save="_highest_expr_genes.png")



#################################################################################################################################
###################################################### Cell/Gene Filtering ######################################################
#################################################################################################################################

# Cell Filtering
print("Remove cells with fewer than 1000 counts and expressing fewer than 300 genes")
sc.pp.filter_cells(adata, min_counts=1000)
sc.pp.filter_cells(adata, min_genes=300)

print(adata.n_obs, adata.n_vars)

print("Remove cells with over 10% mitochondrial gene counts")
adata = adata[adata.obs['pct_counts_mt'] < 10, :]

print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")


# Gene Filtering
print("Removing Genes Detected in Fewer than 50 Cells")
sc.pp.filter_genes(adata, min_cells=50)
print(adata.n_obs, adata.n_vars)

print("Removing Artefactual Genes")
malat1 = adata.var_names.str.startswith('MALAT1')

remove = np.array(malat1)
keep = np.invert(remove)
adata = adata[:,keep]
print(adata.n_obs, adata.n_vars)



#################################################################################################################################
###################################################### Run PostFiltering QC #####################################################
#################################################################################################################################


print("Post-Filtering QC")

#plot qc
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt','pct_counts_ribo', 'pct_counts_hb'],
             jitter=0.4, groupby = 'sample', rotation= 45, save="_postfilt_initial_qc.png")
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save="_postfilt_mt_by_counts.png")
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color="pct_counts_mt", 
              save="_postfilt_nfeatures_by_counts.png")


#plot qc
fig = sns.displot(adata.obs["total_counts"], bins=100, kde=False).fig
fig.savefig("figures/hist_postfilt_total_counts.png") 

fig = sns.displot(adata.obs["n_genes_by_counts"], bins=100, kde=False).fig
fig.savefig("figures/hist_postfilt_ngenes_counts.png") 

#improve visibility for lower counts
sub = adata[adata.obs["total_counts"]<5000]
fig = sns.displot(sub.obs["total_counts"], bins=100, kde=False).fig
fig.savefig("figures/hist_postfilt_total_counts_lowcounts.png")

fig = sns.displot(sub.obs["n_genes_by_counts"], bins=100, kde=False).fig
fig.savefig("figures/hist_postfilt_ngenes_counts_lowcounts.png")

sc.pl.highest_expr_genes(adata, n_top=20, save="_filtered_highest_expr_genes.png")



#################################################################################################################################
##################################################### PreProcessing Complete ####################################################
#################################################################################################################################

print("Filtered Dataset")
print(adata.n_obs, adata.n_vars)


adata.write_h5ad(results_file, compression="gzip")
