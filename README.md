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
    <a href="https://github.com/othneildrew/Best-README-Template/issues/new?labels=enhancement&template=feature-request---.md">Read Preprint</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-project">About The Project</a>
    <li><a href="#code">Code</a>
      <ul>
        <li><a href="#initial-qc">Initial Quality Control</a></li>
        <li><a href="#batch-correction">Batch Correction Benchmark</a></li>
        <li><a href="#aml-scatlas-analysis">AML scAtlas Analysis</a></li>
        <li><a href="#further-analysis">Further Analysis</a></li>
      </ul>
    </li>
    <li><a href="#reference-mapping">AML scAtlas as a Single-Cell Reference</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

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

[Scripts](https://github.com/jesswhitts/AML-scAtlas/tree/main/3_Analysis) used to carry out the main analysis steps on the complete AML scAtlas, split into relevant sub-folders. 

1. Batch Correction
2. Dimensionality reduction
3. Clustering
4. Annotation - automated annotation tool scripts and code used to identify marker genes, used for manual curation

#### Further Analysis

Code used to identify leukemic stem cell (LSC) populations within AML scAtlas. t(8;21)

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## AML scAtlas as a Single-Cell Reference

AML scAtlas can be used as a single-cell reference for new single-cell datasets. An example notebook for this is included.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!--Contact -->
## Contact

* Jessica Whittle (author) - jessica.whittle@cruk.manchester.ac.uk
* Mudassar Iqbal (co-corresponding author) - 
* Georges Lacaud (co-corresponding author) - 
<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!--Acknowledgements -->
## Acknowledgements

AML scAtlas can be used as a single-cell reference for new single-cell datasets. An example notebook for this is included.

<p align="right">(<a href="#readme-top">back to top</a>)</p>
