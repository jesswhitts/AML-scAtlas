<!-- PROJECT LOGO -->
<br />

  <h3 align="center">Single-Cell Atlas of Acute Myeloid Leukemia (AML scAtlas)</h3>

  <p align="center">
    Analysis for the creation of AML scAtlas.
    <br />
    <a href="https://cellxgene.bmh.manchester.ac.uk/AML/"><strong>Explore the Data »</strong></a>
    <br />
    <br />
    <a href="#code">View Code</a>
    ·
    <a href="https://bioarchivelink.md">Read Preprint</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#description">Description</a>
    <li><a href="#code">Code</a>
      <ul>
        <li><a href="#initial-qc">Initial Quality Control</a></li>
        <li><a href="#batch-correction">Batch Correction Benchmark</a></li>
        <li><a href="#aml-scatlas-analysis">AML scAtlas Analysis</a></li>
        <li><a href="#further-analysis">Further Analysis</a></li>
      </ul>
    </li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## Description

[![Data Visualisation][data-image]](https://cellxgene.bmh.manchester.ac.uk/AML)

Using large scale data integration of publicly available single-cell data, we created a single-cell transcriptomic atlas of acute myeloid leukaemia (AML). The data is hosted for easy gene expression exploration, and is downloadable as an AnnData object. 

[data-image]: https://github.com/jesswhitts/AML-scAtlas/blob/main/images/celltype_umap.png
<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CODE -->
## Code

All code is included in this repository, divided into the main analysis steps:

#### Initial Quality Control

QC was completed for each scRNA-seq sample individually and on the combined dataset, using a [standardised workflow](https://github.com/jesswhitts/AML-scAtlas/tree/main/1_Initial_QC). 

1. Sample QC - Example notebook shows preprocessing steps performed on one sample
2. Batch Assessment - Batch effects were quantified within and between studies
3. Integrated QC - Combining all data and QC steps
4. Uncorrected Analysis - Dimensionality reduction on dataset without batch correction

#### Batch Correction

To handle batch effects, we performed [benchmarking](https://github.com/jesswhitts/AML-scAtlas/tree/main/2_Batch_Correction_Benchmark) of 3 batch correction methods scalable to single-cell atlas integration tasks - Harmony, scVI, and scANVI. 

1. all-data - batch correction methods implemented using all genes
2. hvg2000 - batch correction methods implemented using 2000 highly variable genes. The same process was also performed using 4000, 6000, 8000 and 10000 highly variable genes as part of the benchmarking.
3. Results - scripts/notebooks used to combine and visualise the benchmarking results.

#### AML scAtlas Analysis

Scripts used to carry out the [main analysis](https://github.com/jesswhitts/AML-scAtlas/tree/main/3_Analysis) steps on the complete AML scAtlas, split into relevant sub-folders. 

1. Batch Correction
2. Dimensionality reduction
3. Clustering
4. Annotation

#### Further Analysis

1. Code used to [identify leukemic stem cell (LSC) populations](https://github.com/jesswhitts/AML-scAtlas/tree/main/4_Extended_Analyses/1_HSPCs) within AML scAtlas.
2. Analysis performed on the t(8;21) data of the AML scAtlas to [identify age-associated GRN in t(8;21) AML](https://github.com/jesswhitts/AML-scAtlas/tree/main/4_Extended_Analyses/2_AML_ETO), and subsequent [validation with the TARGET/BeatAML cohorts](https://github.com/jesswhitts/AML-scAtlas/tree/main/4_Extended_Analyses/2_AML_ETO/TARGET_BeatAML), and [Lambo et al data](https://github.com/jesswhitts/AML-scAtlas/tree/main/4_Extended_Analyses/2_AML_ETO/Lambo_et_al)

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!--Contact -->
## Contact

* Jessica Whittle (author) - jessica.whittle@cruk.manchester.ac.uk
* Mudassar Iqbal (corresponding author) - mudassar.iqbal@manchester.ac.uk
* Georges Lacaud (corresponding author) - georges.lacaud@cruk.manchester.ac.uk
* Syed Murtuza Baker (corresponding author) - syed.murtuzabaker@manchester.ac.uk
<p align="right">(<a href="#readme-top">back to top</a>)</p>

