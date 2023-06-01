# Training and testing reference panels for quantifying cellular heterogeneity from DNA methylation profiles from bulk tissues  

## Overview

The scripts in this folder relate to the analyses performed in XXX. 

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

| Filename | Description | Command line arguments | 
| --- | ----------- | ----------- |
| summariseDataset.r | Create summary tables of demographics and perform PCA analysis | <ul><li> path to RDS file with normalised dataset, containing betas matrix object (`norm.all`) and phenotype matrix (`pheno.all`) </li><li> path to folder to save plots </li></ul> |
| fitModelsSimulations.r | Trains and tests prediction of specified combination of neural cell types against reconstructed bulk profiles |  <ul><li> path to RDS file with normalised dataset, containing betas matrix object (`norm.all`) and phenotype matrix (`pheno.all`) </li><li>  path to folder with array annotation files </li><li> path to csv file with reference panels specified </li><li>  index of which reference panel to train (row number of previous file) needs to be rerun for each reference panel </li></ul> |
| summariseSimulations.r | Aggregate and visualise the simulation results comparing different reference panels  | <ul><li> path to RDS objects with output from `fitModelsSimulations.r` </ul></li>|
| fitModelsAll.r | Trains a series of models to predict different combinations of neural cell types  | <ul><li> path to RDS file with normalised dataset, containing betas matrix object (`norm.all`) and phenotype matrix (`pheno.all`) </li><li>  path to folder with array annotation files </li><li> path to csv file with reference panels specified </li></ul>|
| profileCellCompBulkBrainSamples.r | Estimates cellular composition of bulk brain DNAm  profiles using pretrained models and tests against biological factors  | <ul><li> path to RDS file with normalised dataset, containing betas matrix object (`norm.all`) and phenotype matrix (`pheno.all`) </li><li>  path pretrained models, output of `fitModelsAll.r` </li><li> path to folder to save plots </li></ul> |
| testCellCompPathology.r | Test estimated cell composition against AD neuropathology | <ul><li> path to RDS file with estimated cellular composition </li><li> path to folder to save plots </li></ul> |

## Data Availablity

## Citation

## Contact

<E.J.Hannon@exeter.ac.uk>
