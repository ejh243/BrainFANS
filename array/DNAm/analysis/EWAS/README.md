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
* lme4
* lmerTest
* doParallel
* devtools
* cdegUtilities
* pwr

## Script Orientation

Order below loosely corresponds to order analyses presented in manuscript. Scripts are designed to be submitted from the command line with paths to files/ folder added as arguments. Note that the command line arguments need to be specified in order listed in the table below. 

For example to run the **summariseDataset.r** script you might execute

` Rscript summariseDataset.r <path to folder> `

| Filename | Description | Required Arguments | 
| --- | ----------- | ----------- |
| summariseDataset.r | Create summary tables and plots of demographic variables | <ol><li> path to RDS file with normalised dataset, containing  phenotype matrix (`QCmetrics`) </ol></li> |
| testCellComposition.r | path to project data folder where normalised dataset, containing  normalised beta matrix (`celltypeNormbeta`) and phenotype matrix (`QCmetrics`) is within a folder called 3_normalised </ol></li> |
| testAgeAcceleration.r | Calculate epigenetic age for each sample and test against chronological age and age acceleration residuals against case control status | <ol><li> path to project data folder where normalised dataset, containing  normalised beta matrix (`celltypeNormbeta`) and phenotype matrix (`QCmetrics`) is within a folder called 3_normalised </li><li> path to folder containing cs files with clock coefficients </ol></li> |
| lmWithinCT.r | Performs case control EWAS for a single cell type | <ol><li> path to project data folder where normalised dataset, containing  normalised beta matrix (`celltypeNormbeta`) and phenotype matrix (`QCmetrics`) is within a folder called 3_normalised </li><li> Cell Type to be analysed </ol></li>  |
| mlm.r | Performs case control EWAS across all cell types | <ol><li> path to project data folder where normalised dataset, containing  normalised beta matrix (`celltypeNormbeta`) and phenotype matrix (`QCmetrics`) is within a folder called 3_normalised </ol></li>  |
| summmariseWithinCTDMPs.r | Loads cell specific EWAS for all cell types and produces summary tables and plots | <ol><li> path to project data folder where normalised dataset, containing  normalised beta matrix (`celltypeNormbeta`) and phenotype matrix (`QCmetrics`) is within a folder called 3_normalised </li><li> Path to folder where EPIC annotation files are stored </ol></li>  |
| calcPowerinOtherCT.r | Performs power analysis for typical DMPs | <ol><li> path to project data folder where normalised dataset, containing  normalised beta matrix (`celltypeNormbeta`) and phenotype matrix (`QCmetrics`) is within a folder called 3_normalised </li><li> Path to folder where EPIC annotation files are stored </ol></li>  |


## Data Availability

Raw and processed DNA methylation data are available from GEO under accession number GSE279509. 

## Citation

## Contact

<E.J.Hannon@exeter.ac.uk>
