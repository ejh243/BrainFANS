# ATAC-seq analysis pipeline 

The Assay for Transposase-Accessible Chromatin followed by sequencing (ATAC-seq) experiment provides genome-wide profiles of chromatin accessibility. This pipeline has been developed to provide end-to-end quality control, processing and analysis of paired-end ATAC-seq datasets. This pipeline has been developed so that can be run end-to-end, starting from raw FASTQ files all the way to peak calling and perform differential accessibility analysis. 

## Pipeline's structure

This scripts necessary for running this pipeline should be found in 3 different subfolders:
- jobSubmission: includes main scripts of pipeline, one for each step. These are the scripts meant to be directly run by the user.
- subScripts: includes subscripts used from the main scripts of the pipeline. These performs some of the substeps that can be specified from the main scripts.
- Rscripts: includes the R and Rmarkdown scripts used from the main scripts of the pipeline.

## Configuration files

In order to use the ATAC-seq pipeline, two main configuration files need to be set up. For examples of these config files, go to ATACSeq/config
- config.txt: This is the main configuration file that will source the pipeline.
  - This file should be found in the main directory of the project/dataset to be analysed.
  - The main variables to be changed in this configuration file are:
    - `PROJECT`: name of the project's folder.
    - `MAIN_DIR`: full path to the project's folder.
    - `REFERENCES_DIR`: full path to directory with all references used throughout the pipeline.
    - Specify the modules versions to be loaded and the full paths to the softwares used in the pipeline (e.g. picard).
  - An alphabetically ordered array of cell types to which samples belong to: `CELLTYPES`. 
- config.r: This is the configuration file for running R scripts.
  - The `dir` variable should be changed to be the full path to the project's directory. This should match `MAIN_DIR` variable in the config.txt file.
  - Other parameters or threshold should be changed as required.
  - An alphabetically ordered list of cell types to which samples belong to needs to be specified as *cellTypes*.
- Additional source files of R packages that are installed manually will also be located in this folder.
- A temporary directory can be either specified in the system's _.bashrc_ file or in the *config.txt* as an additional variable `TMPDIR`. Please make sure you have full permissions for the directory specified.

## Requisites:

- The full path to the project folder must be specified as the first argument on the command line for all scripts.
- Requires config files (file paths should be specified in config.txt):
  - config.txt: with variables for bash scripts.
  - config.r with parameters for R scripts
  - environment.yml : file with packages and their version to be installed by conda.
  - packagesPip.txt : list of packages and their version to be installed by pip in the conda environment. 
- Requires a samples.txt file in the META_DATA (0_metadata) folder with the names of the samples that will be used to run the pipeline.
- Requires a sampleSheet.csv file in the META_DATA (0_metadata) folder with information about the samples. More information below.
- Samples need to be in the RAWDATADIR folder in the project folder (1_raw)
- An alphabetically ordered array of cell types to which samples belong to needs to be specified in the *config.txt* file as `CELLTYPES` in order to perform peak calling by cell group (STEP 7). This also required in the *config.r* file.
- A Reference folder is needed with the following resources:
  - Reference genome in .fa format (hg38): [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
  - Index for Bowtie2 alignment: [Bowtie index](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer)
  - ENCODE blacklist regions: [hg38](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz)
  - Pseudoautosomal regions (PAR)
  - 1000 Genomes Project [data](https://www.ensembl.org/Homo_sapiens/Info/Index)
  - X chromosome gene list.

## Software 

Several software are needed to run the pipeline. Some of these can be used from a conda or pip environment. The set up script (0_setUp.sh) will create a conda environment where some of these packages will be installed either by conda or pip. As explained in *Requisites*, the packagesPip.txt and environment.yml files are needed for the conda environment to be set up. This can be done using either Anaconda or Miniconda, which local source directory needs to be specified in the *config.txt* file (`CONDA`). In order to download and install refer to: [Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)

Important software that will be installed in this environment are: MACS3 (3.0.2), SAMstats, SAMtools, BEDTools, Bowtie2, MultiQC, R, Python3, Pandoc, GATK, Picard. If you already have a conda environment, please specify its name or path in the *config.txt* file. 
  - Python version in this conda environment should be >= 3.12, as this is required by MACS3. For further details about MACS3 requirements: [MACS3 documentation](https://macs3-project.github.io/MACS/docs/INSTALL.html)
  - R version should be >= 4.2. Some packages are no longer supported in the regular CRAN repository and need to be installed from source. To do this, download the source file of the package [ptest R package](https://cran.r-project.org/src/contrib/Archive/ptest/), keep this in the *config* folder and this will be installed from source at the *0_setUp.sh* script.
Other software that needs to be available locally are (and path to these need to be specified in the *config.txt* file):
  - [VerifyBamID](https://github.com/Griffan/VerifyBamID)
  - [Phantompeakqualtools](https://github.com/kundajelab/phantompeakqualtools)
Important software that will be installed in this environment are: MACS3 (3.0.2) and samstats. If you already have a conda environment, please specify its name or path in the *config.txt* file. Please check the python version in this conda environment is <= 3.12, as this is required by MACS3. For further details about MACS3 requirements: [MACS3 documentation](https://macs3-project.github.io/MACS/docs/INSTALL.html)

Other software needed are: BEDTools, Bowtie2, Picard and R. Note R also needs to be installed in the conda environment, as there are some packages that are not available for later versions of R. 

## Sample sheet 

A sample sheet with the metadata about samples is essential to run the ATAC-seq pipeline. This should be located in the META_DATA (0_metadata) folder and should be named as sampleSheet.csv. 
Samples must be in the same order as in the samples.txt file. An example of information about a sample is shown below. 

| sequencingBatch | sampleID | sampleCode | cohort | fraction | individualID | age | gender | clinical | vcfID |
| --------------- | -------- | ---------- |------- | -------- | ------------ | --- | ------ | -------- | ----- |
| 11031	| 11031_EX161_NEUN_S1	| EX161_NEUN | BRISTOL_A	| NEUN | 859 | 85	| F |	5	| 201023670019_R08C02_EX161 |

- `sequencingBatch`: sequencing batch number
- `sampleID`: name of sample in samples.txt and name of fastq file corresponding to the sample. This is usually formed by the sequencingBatchNumber_sampleID_fraction_SNumber.
- `sampleCode`: ID of sample.
- `cohort`: brain bank or cohort.
- `fraction`: cell type
- `individualID`: ID of the individual to which the sample belongs to. All samples that belong to the same individual must have the same ID.
- `age`
- `gender`: F or M.
- `clinical`: braak stage
- `vcfID`: ID of imputed array data to perform genotype comparison and confirm sample's identity.

## STEPS

Parameters in [] are optional. If no step is specified, all steps in the script will be run. For more information about the inputs needed for each step and the outputs generated, please refer to each script.

### 0. Set-up

`sbatch --array=0 ./jobSubmission/0_setUp.sh (project directory) `

This script checks for required files, packages or libraries that are needed later in the pipeline.

##### -scripts executed-
- ./Rscripts/installLibraries.r : checks for already installed R libraries and install if not found.

##### -parameters-
- `--array`: should be 0 for this script as general analysis is run rather than individual samples.
- `(project-name)` project's directory.

##### -requires-
- packagesPip.txt: list of packages to be installed by pip in the conda environment.
- environment.yml: list of conda packages to be installed in conda.
- config.r: R config file with data paths and thresholds used in analysis.

#### -outputs-
- creates data folder if not found: 2_trimmed, 3_aligned, 4_calledPeaks
- creates conda environment with required software for later subprocesses of the pipeline.
- installs R libraries.


### 1. Pre-analysis (QC and alignment)

 `sbatch --array=<number of batch jobs> ./jobSubmission/1_batchRunPreAnalysis.sh (project directory) [STEPS]`

Performs pre-analysis of ATAC-seq data, including pre-alignment quality control, alignment and post-alignment quality control of samples. 

##### -scripts executed-
- ./subScripts/fastqc.sh : FastQC for pre-alignment quality control.
- ./subScripts/fastp.sh : FastP for trimming samples.
- ./subScripts/alignment.sh : Alingment of samples to reference genome using Bowtie2.
- ./subScripts/calcENCODEQCMetrics.sh : ENCODE QC metrics are calculated on aligned samples.
  
##### -parameters-
- `--array`: number of batch jobs, each number matches a sample.
- `(project-name)` project's directory.
-optional-
- `[STEPS]` Option to combine steps with desired steps included as single string, i.e. `FASTQC,TRIM`. Default if left blank is to run all steps. Options:
  - `FASTQC`: Perform only quality control. Run only fastqc.sh on samples.
  - `TRIM`: Perform only trimming on samples. Run only fastp.sh on samples.
  - `ALIGN`: Align samples. Run only alignment.sh on samples.
  - `ENCODE`: Calculation of ENCODE QC metrics.

### 2. Post-alignment processing

 `sbatch --array=<number of batch jobs/10> ./jobSubmission/2_batchCalcQCMetrics.sh  (project directory)`

This script uses the ATACseqQC R package to generate the fragment distribution and calculate some summary statistics. It splits the samples into groups of 10 to run in parallel
If the number of samples is less than 10, set batch number to 0.

##### -scripts executed-  
- ./Rscripts/fragmentDistribution.r (R config file) (array-number) : fragment distribution of samples specified by array. 

##### -parameters-
- `--array`: number of batch jobs, should be the number of samples divided by 10. e.g. If there are 50 samples, batch number should be 0-5, producing 5 batches of 10 samples each.
- `(R config file)` directory of R config file


### 3. Peak calling

  `sbatch --array=<number of batch jobs> ./jobSubmission/3_batchRunPeakCalling.sh (project directory) [STEPS]`

This script performs the core part of the ATACseq pipeline: calling peaks. In this stage, peak calling is performed at sample level using the Paired-end mode of MACS3. After this, FRIP statistics are calculated and collated in a single file.

##### -scripts executed-
- ./subScripts/shiftAlignedReads.sh : takes a filtered bam file converts to a tagalign file, calculate CC scores and shifts reads ready for peak calling
    * This step is not needed if peak calling is performed using PE mode or with HMMRATAC
- ./subScripts/samplePeaks.sh : runs MAC version 3 peak calling with BAM files paired end reads
    * it then filters the peaks to exclude those that overlap with blacklisted regions, sex chromosomes 
    * peaks are sorted by chromosome
- ./subScripts/collateCalcFrip.sh : calculates fraction of reads in peaks, number of reads and peaks for peak calling at sample level
    * a single --array job number should be specified.

##### -parameters-
- `--array`: number of batch jobs, each number matches a sample.
- `(project-name)` project's directory.
-optional-
- `[STEPS]` Option to combine steps with desired steps included as single string, i.e. `FASTQC,TRIM`. Default if left blank is to run all steps. Options:
  - `SHIFT`: Perform only shifting of reads.
  - `PEAKS`: Perform only peak calling on samples. 
  - `FRIP`: Calculate peak calling results for all samples.
  
### 4. Stage 1 QC metrics summary

   `sbatch ./jobSubmission/4_collateStage1QCMetrics.sh  (project directory) [STEPS]`

This scripts uses MultiQC to collate the output of fastqc and bowtie2 alginment, as well as peak calling results. 

##### -scripts executed-
- ./subScripts/progressReport.sh : identifies how many samples have been successful at each stage of the processing pipeline and for each fastq file, how far through the process it has progressed.
- ./subScripts/countMTcollateFS.sh : collates counts of the number of reads aligned to MT chromosome and collates flagstat summary of- aligned sorted reads.
- ./Rscripts/collateDataQualityStats.Rmd : generates Rmarkdown report summarising stage 1 qc metrics: raw reads metrics, trimming, alignment and peak calling
    * outputs a list of samples that pass Stage 1 of QC: passStage1SampleList.txt

##### -parameters-
- `--array` : a single job number should be specified.
- `(project-name)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. MULIQC,COLLATE. Default if left blank is to run all of them.
  - `MULTIQC`: Collates fastqc and alignment statistics in a single report.
  - `COLLATE`: Collate results from all samples to show progress of each sample through the pipeline in order to avoid missing steps.
  - `SUMMARY`: Collates stage 1 QC and peak calling results in a single Rmarkdown report.

### 5. Sex chromosomes 

  `sbatch --array=<number of batch jobs> ./jobSubmission/5_batchFormatSexChrs.sh (project directory) [STEPS]`
  
This script creates bam files for each sample containing only reads aligned to sex chromosomes, perform peak calling on these and output results.

##### -scripts executed-
- ./subScripts/subsetSexChrs.sh : subsets reads of only X and Y chromosomes for input sample in bam file. 
- ./subScripts/sexChrPeaks.sh : performs peak calling on the sex chromomes, filter and read counts using MACS3 in Single-end mode.
- ./Rscripts/collateSexChecks.r : collates peak calling results on sex chromosomes and uses it to check the assigned sex of the sample.

##### -parameters-
- `--array` : number of batch jobs, each number matches a sample.
- `(project-directory)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. MULIQC,COLLATE. Default if left blank is to run all of them.
  - `SPLIT`: Generate only samples with X and Y chromosomes reads.
  - `PEAKS`: Perform only peak calling on new sex chromosomes reads samples.
  - `CHECK`: Collate peak calling on sex chromosomes results and check assigned sex of each sample.

### 6. Genotype check and Stage 2 QC metrics summary

  `sbatch --array=<number of batch jobs> ./jobSubmission/6_batchRunGenotypeCheck.sh (project directory) [STEPS]`

This script will detect possible DNA contamination in order to ensure high quality sequence reads. If any contamination is detected, possible swaps will be suggested. It will also collate results from stage 2 of QC.

##### -scripts executed-
- ./subScripts/compareBamWithGenotypes.sh  : prepares bam file for comparison of assigned genotype of each sample and runs verifyBamID 
  * requires a file in METADATA folder called matchedVCFIDs.txt which lists the samples with their matched vcfID.
- ./Rscripts/collateSampleChecks.Rmd : collates results from previous steps (sex and genotype check).
- ./subScripts/searchBestGenoMatch.sh : Outputs a summary of stats from previous step and finds any sample that might be contaminated.
  * Contaminated samples will go to a created file potentialSwitches.txt. If this file is not empty, searchBestGenoMatch.sh is executed and an alternative Genotype search is done for the sample.
- ./subScripts/progressReportS2.sh : identifies what samples have been through stage 2 of QC in order to avoid missing processing samples

##### -parameters-
- `--array` : number of batch jobs, each number matches a sample. Starting from 2
- `(project-directory)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. GENCHECK,COMPARE. Default if left blank is to run all of them.
  - `COMPARE`: Compare only samples genotype with matched VCF file genotype.
  - `COLLATE`: Collate results from all samples to show progress of each sample through stage 2 of QC pipeline in order to avoid missing steps or samples.
  - `GENCHECK`: Sex and gencheck results are collated in a report.
  - `SWITCH`: Indentify only mismatched genotype samples and output potential switches and seach for alternative genotype.


### 7. Group peak calling

 `sbatch ./jobSubmission/7_batchPeakCallingByGroup.sh (project directory) [STEPS] [GROUP]`

This script groups samples per cell type, performs group level peak calling and gets read counts in peaks called. Samples used are those that have passed previous quality control stages 1 and 2. Cell fractions to which samples belong to need to be specified in the *config.txt* file.

##### -scripts executed-
- ./Rscripts/samplesForGroupAnalysis.r: creates a file with samples of cell group input.
- ./Rscripts/fragDistributionForPeaks.r: outputs the fragment distribution of samples that are used for group peak calling, in order to check their ATAC-seq quality.
- ./subScripts/groupPeaks.sh : performs peak calling on samples that belong to the cell group chosen by the user.
- ./Rscripts/countsInPeaks.r: : gets read counts in peaks called on each cell group.

##### -parameters-
- `(project-directory)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. FRAGSIZE,PEAKS. Default if left blank is to run all of them.
  - `FRAGSIZE`: Get fragment size distribution of samples chosen for peak calling.
  - `PEAKS`: Performs peak calling on samples chosen for peak calling using MACS3 in paired-end mode.
  - `COUNTS`: Gets read counts in peaks called at cell group level.
- `[GROUP]` Select cell-type of group of samples to perform peak calling on. This is only needed for the `PEAKS` step.

### 8. Stage 3 QC metrics summary

 `sbatch ./jobSubmission/8_collateStage3QCMatrics.sh (project directory) [STEPS] [peakSet]`

This script normalises counts in peaks and uses results to check cell-type identity of samples, as this might be of poor quality or mislabelled.

##### -scripts executed-
- ./Rscripts/normCounts.r: outputs results of Variance Partition Analysis (VPA) on raw counts, normalising counts in peaks and results of VPA in these.
- ./Rscripts/collateCellTypeCheck.Rmd: outputs a summary report of VPA on raw and normalised counts, the normalisation results and identifies samples that fail cell-type check.

##### -parameters-
- `(project-directory)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. FRAGSIZE,PEAKS. Default if left blank is to run all of them.
  - `NORM`: performs Variance Partition Analysis (VPA) on raw counts, normalises counts in peaks and performs VPA again in these.
  - `CTCHECK`: Collates results from before in a Rmarkdown report and identifies samples that fail cell-type identity check.
- `[peakSet]` : as counts in peaks have been produced for only promoters or all peaks, user can choose whether this stage is performed using either of those peak sets.