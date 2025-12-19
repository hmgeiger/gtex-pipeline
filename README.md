# Analysis pipelines for the [GTEx Consortium](https://www.gtexportal.org) and [TOPMed](https://www.nhlbi.nih.gov/science/trans-omics-precision-medicine-topmed-program)

This repository contains analysis pipelines for:

* RNA-seq alignment, quantification, and quality control
* eQTL mapping and annotation
* Allele-specific expression quantification
* Generation of the collapsed annotation used for gene-level expression quantification

Pipeline components are available as docker images, and execution scripts are provided in WDL (to add links to NYGC WDL files). All Python scripts are written in Python 3.

A detailed description of the RNA-seq pipeline settings used for TOPMed is provided in [TOPMed_RNAseq_pipeline.md](https://github.com/hmgeiger/gtex-pipeline/blob/patch-1/TOPMed_RNAseq_pipeline.md).
