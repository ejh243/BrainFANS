# Training and testing reference panels for quantifying cellular heterogeneity from DNA methylation profiles from bulk tissues  

## Overview

The scripts in this folder relate to the analyses performed in [Hannon, E., Dempster, E.L., Davies, J.P. et al. Quantifying the proportion of different cell types in the human cortex using DNA methylation profiles. BMC Biol 22, 17 (2024). https://doi.org/10.1186/s12915-024-01827-y](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-024-01827-y). 

## Requirements

Scripts rely on the following packages being installed

* paletteer
* ggplot2
* ggpubr
* Rmisc
* tidyverse
* reshape
* CETYGO
* tidyr
* corrplot

Scripts rely on the following array annotation files being available in the same folder

AllProbeIlluminaAnno.Rdata
CrossHybridisingProbesPriceORWeksberg.csv
SNPsinProbesAnno.csv



## Script Orientation

Order below loosely corresponds to order analyses presented in manuscript. Scripts are designed to be submitted from the command line with paths to files/ folder added as arguments. Note that the command line arguments need to be specified in order listed in the table below. 

For example to run the **summariseDataset.r** script you might execute

` Rscript summariseDataset.r normData.rdata plots`

| Filename | Description | Required Arguments | 
| --- | ----------- | ----------- |
| summariseDataset.r | Create summary tables of demographics and perform PCA analysis | <ol><li> path to RDS file with normalised dataset, containing betas matrix object (`norm.all`) and phenotype matrix (`pheno.all`) </li><li> path to folder to save plots </li></ol> |
| fitModelsSimulations.r | Trains and tests prediction of specified combination of neural cell types against reconstructed bulk profiles |  <ol><li> path to RDS file with normalised dataset, containing betas matrix object (`norm.all`) and phenotype matrix (`pheno.all`) </li><li>  path to folder with array annotation files </li><li> path to csv file with reference panels specified </li><li>  index of which reference panel to train (row number of previous file) needs to be rerun for each reference panel </li></ol> |
| summariseSimulations.r | Aggregate and visualise the simulation results comparing different reference panels  | <ol><li> path to RDS objects with output from `fitModelsSimulations.r` </ol></li>|
| fitModelsAll.r | Trains a series of models to predict different combinations of neural cell types  | <ol><li> path to RDS file with normalised dataset, containing betas matrix object (`norm.all`) and phenotype matrix (`pheno.all`) </li><li>  path to folder with array annotation files </li><li> path to csv file with reference panels specified </li></ol>|
| profileCellCompBulkBrainSamples.r | Estimates cellular composition of bulk brain DNAm  profiles using pretrained models and tests against biological factors  | <ol><li> path to RDS file with normalised dataset, containing betas matrix object (`norm.all`) and phenotype matrix (`pheno.all`) </li><li>  path pretrained models, output of `fitModelsAll.r` </li><li> path to folder to save plots </li></ol> |
| testCellCompPathology.r | Tests best estimates of cellular composition against AD neuropathology measured by Braak stage | <ol><li> path to RDS file with estimated cellular composition </li><li> path to folder to save plots </li></ol> |

## Data Availablity

Raw and processed DNA methylation data are available from GEO under accession number GSE279509. 

## Citation

If you use these scripts please cite the following manuscript: [Hannon, E., Dempster, E.L., Davies, J.P. et al. Quantifying the proportion of different cell types in the human cortex using DNA methylation profiles. BMC Biol 22, 17 (2024). https://doi.org/10.1186/s12915-024-01827-y](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-024-01827-y). 

## Contact

<E.J.Hannon@exeter.ac.uk>
