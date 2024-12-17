# Neural cell-specific EWAS of schizophrenia   

## Overview

The scripts in this folder relate to the analyses performed in XXX. 

## Requirements

Scripts rely on the following packages being installed

* CETYGO
* ggplot2
* ggpubr
* qqman
* tidyr
* dplyr
* IlluminaHumanMethylationEPICanno.ilm10b4.hg19


## Script Orientation

Order below loosely corresponds to order analyses presented in manuscript. Scripts are designed to be submitted from the command line with paths to files/ folder added as arguments. Note that the command line arguments need to be specified in order listed in the table below. 

For example to run the **summariseDataset.r** script you might execute

` Rscript summariseDataset.r <path to folder> `

| Filename | Description | Required Arguments | 
| --- | ----------- | ----------- |
| summariseDataset.r | Create summary tables and plots of demographic variables | <ol><li> path to RDS file with normalised dataset, containing  phenotype matrix (`QCmetrics`) |
| testCellComposition.r | Calculate celluar composition for each sample and test against case control status | <ol><li> path to RDS file with normalised dataset, containing  normalised beta matrix (`celltypeNormbeta`) and phenotype matrix (`QCmetrics`) |



## Data Availablity

Raw and processed DNA methylation data are available from GEO under accession number GSE279509. 

## Citation

## Contact

<E.J.Hannon@exeter.ac.uk>
