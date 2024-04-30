# BrainFANS

## Project Overview

Increased understanding about the functional complexity of the genome has led to growing recognition about the role of non-sequence-based epigenetic variation in health and disease. Current genomic analyses of the human brain, however, are limited by the use of “bulk” tissue, comprising a heterogeneous mix of different neural cell-types. Because epigenetic processes play a critical role in determining cell type specific patterns of gene regulation it is critically important to consider cellular composition in regulatory genomics studies of human post-mortem tissue. Furthermore, the valuable nature of human post-mortem tissue means it is important to use methods that maximize the amount of genomic data generated on each sample.

We have developed a [method using fluorescence-activated nuclei sorting (FANS)](https://dx.doi.org/10.17504/protocols.io.bmh2k38e) to isolate and profile nuclei from different human brain cell-types. We demonstrate that this can be used to robustly purify populations of neuronal, oligodendrocytes and other (glial) origin nuclei from adult and fetal post-mortem frozen brain, with each tissue sample yielding purified populations of nuclei amenable to simultaneous analysis of i) DNA modification analysis (5mC and 5hmC), ii) histone modification analysis, iii) open chromatin analysis, and iv) RNA sequencing. Our overarching aim is to generate an integrated resource of genomic, regulatory, epigenomic and transcriptomic landscapes across both developing and adult human brains to explore molecular mechanisms involved in schizophrenia and autism.

This repository contains analysis pipelines for the pre-processing and analysis of the data generated as part of this project. They are developed for use with a HPC system with the SLURM job scheduler. 

## Publications

The following publications are based on analyses performed using scripts from this repository. Links to the relevant folders can be found below. 

[Hannon, E., Dempster, E.L., Davies, J.P. et al. Quantifying the proportion of different cell types in the human cortex using DNA methylation profiles. BMC Biol 22, 17 (2024). https://doi.org/10.1186/s12915-024-01827-y]( https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-024-01827-y#citeas)

* Pipeline version: [v1.0.0](https://github.com/ejh243/BrainFANS/releases/tag/v1.0.0) 
* [QC Scripts](https://github.com/ejh243/BrainFANS/tree/master/array/DNAm/preprocessing)
* [Analysis](https://github.com/ejh243/BrainFANS/tree/master/array/DNAm/analysis/neuralCellComposition)

