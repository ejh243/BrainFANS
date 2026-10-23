# PIPPA
Pipeline for Integrative Peak calling in Population-scale ATAC-seq (**PIPPA**) provides end-to-end quality control and processing of ATAC-seq datasets and curated of reproducible peak sets. The pipeline has been developed to work from raw FASTQ files all the way to curated peak sets that can be used for downstream analyses.

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
  - Cutoff values for peak calling should be specified by the user. If not specified, default will be used (0.05). 
    - `macs_sample_thres` : Sample-level peak calling q-value cutoff
    - `macs_group_thres` : Cell fraction-level peak calling q-value cutoff
    - `macs_chr_thres` : Peak calling on sex chromosomes q-value cutoff
- config.r: This is the configuration file for running R scripts.
  - The `mainDir` variable should be changed to be the full path to the project's directory. This should match `MAIN_DIR` variable in the config.txt file.
  - Other parameters or threshold should be changed as required.
  - An alphabetically ordered list of cell types to which samples belong to needs to be specified as *cellTypes*.
- Additional source files of R packages that are installed manually will also be located in this folder.
- A temporary directory can be either specified in the system's _.bashrc_ file or in the *config.txt* as an additional variable `TMPDIR`. Please make sure you have full permissions for the directory specified.
- For each main bash script, make sure the SBATCH details are filled in/removed if not apply.

## Requisites:

- The full path to the project folder must be specified as the first argument on the command line for all scripts.
- Requires config files (file paths should be specified in config.txt):
  - config.txt: with variables for bash scripts.
  - config.r with parameters for R scripts
  - condaEnv.yml : file with packages and their version to be installed by conda as an environment.
- Requires a samples.txt file in the META_DIR (0_metadata) folder with the names of the samples that will be used to run the pipeline. 
- Requires a sampleSheet.csv file in the META_DIR (0_metadata) folder with information about the samples. Name given/path should be specified as `sampleSheet` in the config.r file. More information below.
- Samples need to be in the RAWDATADIR folder in the project folder (1_raw)
- An alphabetically ordered array of cell types to which samples belong to needs to be specified in the *config.txt* file as `CELLTYPES` in order to perform peak calling by cell group (STEP 8). This also required in the *config.r* file. 
- A Reference folder is needed with the following resources:
  - Reference genome in .fa format (hg38): [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/)
  - Index for Bowtie2 alignment: [Bowtie index](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer)
  - ENCODE blacklist regions: [hg38](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz)
  - Pseudoautosomal regions (PAR) (hg38) [UCSC](https://hgdownload.soe.ucsc.edu/gbdb/hg38/html/description.html)
  - X chromosome gene list (FIRRE and XIST) : chrX	73820656	73852723	XIST	10	- , chrX  131688779	131830604	FIRRE	10	-
  - RefSeq file: [NCBI Reference Sequence Database](https://www.ncbi.nlm.nih.gov/refseq/)

## Software 

Several software are needed to run the pipeline. Some of these can be used from a conda. The set up script (0_setUpPipeline.sh) will create a conda environment where some of these packages will be installed by conda. As explained in *Requisites*, the condaEnv.yml file can be used to set up the conda environment. This can be done using either Anaconda or Miniconda, which local source directory needs to be specified in the *config.txt* file (`CONDA`). In order to download and install refer to: [Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://docs.anaconda.com/miniconda/miniconda-install/)

Important software that will be installed in this environment are: MACS3 (3.0.2), SAMstats, SAMtools, BEDTools, Bowtie2, MultiQC, R, VerifyBamID, Python3, Pandoc, GATK, Picard. If you already have a conda environment, please specify its name or path in the *config.txt* file. 
  - Python version in this conda environment should be >= 3.12, as this is required by MACS3. For further details about MACS3 requirements: [MACS3 documentation](https://macs3-project.github.io/MACS/docs/INSTALL.html)
  - R version should be >= 4.3. Most packages are installed through the conda environment. A couple of packages would be installed when running the *0_setUpPipeline.sh* script. Some packages are no longer supported in the regular CRAN repository and need to be installed from source. To do this, download the source file of the package [ptest R package](https://cran.r-project.org/src/contrib/Archive/ptest/), keep this in the *config* folder and this will be installed from source when running the *0_setUpPipeline.sh* script.

## Sample sheet 

A sample sheet with the metadata about samples is essential to run the ATAC-seq pipeline. This should be located in the META_DIR (0_metadata) folder. The name of this file should be specified in the config files as `SAMPLE_SHEET` and `sampleSheet` in the *config.txt* and *config.r* files, respectively.
Samples must be in the same order as in the samples.txt file. An example of information about a sample is shown below. 

| sequencingBatch | sampleID | sampleCode | cohort | fraction | individualID | age | gender | clinical | vcfID |
| --------------- | -------- | ---------- |------- | -------- | ------------ | --- | ------ | -------- | ----- |
| X	| X_A61_NEUN_S1	| A61_NEUN | Y	| NEUN | 859 | 85	| F |	5	| XXXX_YYYYY_A61 |

- `sequencingBatch`: sequencing batch number
- `sampleID`: name of sample in samples.txt and name of fastq file corresponding to the sample. This is usually formed by the sequencingBatchNumber_sampleID_fraction_SNumber.
- `sampleCode`: ID of sample.
- `cohort`: brain bank or cohort.
- `fraction`: cell type
- `individualID`: ID of the individual to which the sample belongs to. All samples that belong to the same individual must have the same ID.
- `age`
- `gender`: F or M.
- `clinical`: clinical group of samples (optional)
- `vcfID`: ID of imputed array data to perform genotype comparison and confirm sample's identity (optional)

## NOTES

This pipeline has been originally designed to be executed with cell type sorted ATAC-seq datasets. However, with the right modifications, it can also be used on bulk tissue ATAC-seq datasets. Instead of a list of cell types in the `CELLTYPES` variable, use "BULK" for example. This has not been tested yet.

## STEPS

Parameters in [] are optional. If no step is specified, all steps in the script will be run. For more information about the inputs needed for each step and the outputs generated, please refer to each script.

### 0. Set-up (optional)

`sbatch ./jobSubmission/0_setUpPipeline.sh (project directory) `

This script checks for required files, packages or libraries that are needed in the pipeline.

##### -scripts executed-
- ./Rscripts/installLibraries.r : checks for already installed R libraries and install if not found.

##### -parameters-
- `(project-name)` project's directory.

##### -requires-
- condaEnv.yml: list of conda packages to be installed in conda.
- config.r: R config file with data paths and thresholds used in analysis.

### 1. Pre-analysis (stage 1 QC) 

 `sbatch --array=<number of batch jobs> ./jobSubmission/1_batchRunPreAnalysis.sh (project directory) [STEPS]`

Performs pre-analysis of ATAC-seq data, including pre-alignment quality control, alignment and post-alignment quality control of samples. 

##### -scripts executed-
- ./subScripts/fastQC.sh (array-number): FastQC for pre-alignment quality control.
- ./subScripts/fastp.sh (array-number): FastP for trimming samples.
- ./subScripts/alignment.sh (array-number): Alingment of samples to reference genome using Bowtie2.
- ./subScripts/calcENCODEQCMetrics.sh (array-number): ENCODE QC metrics are calculated on aligned samples.
  
##### -parameters-
- `--array`: number of batch jobs, each number matches a sample, starting from 0.
- `(project-name)` project's directory.
-optional-
- `[STEPS]` Option to combine steps with desired steps included as single string, i.e. `FASTQC,TRIM`. Default if left blank is to run all steps. Options:
  - `FASTQC`: Perform only quality control. Run only fastQC.sh on sample.
  - `TRIM`: Perform only trimming on samples. Run only fastp.sh on sample.
  - `ALIGN`: Align samples. Run only alignment.sh on sample.
  - `ENCODE`: Calculation of ENCODE QC metrics. Run only calcENCODEQCMetrics.sh on sample.

##### -requires-
- file with adapters to be trimmed at `STEP 1.2 TRIM`. Should be in fasta format .fa and specified in *config.txt* file as `ADAPTERS_FILE`. 

### 2. ATAC-seq QC metrics (stage 1 QC)

 `sbatch --array=<number of batch jobs> ./jobSubmission/2_batchCalcQCMetrics.sh  (project directory) [STEPS] (GROUP_SAMPS)`

Performs analysis of ATAC-seq-specific characteristics, including fragment size distribution and TSS enrichment scores.

##### -scripts executed-  
- ./Rscripts/fragmentDistribution.r (R config file) (array-number) : Calculation of fragment distribution of samples specified by array. 
- ./subScripts/prepareBigWig.sh (array-number) : Coversion of sample bam file to bigWig format and shift.
- ./subScripts/calcTssEnrich.sh (array-number) : Calculation of TSS enrichment score of sample.
- ./subScripts/collateTssScores.sh (array-number) : Collation of TSS enrichment score of all samples.

##### -parameters-
- `--array`: number of batch jobs, each number matches a sample, starting from 0. Only required for BW and SCORE.
- `(project-name)` project's directory.
- `[STEPS]` Option to combine steps with desired steps included as single string, i.e. `FASTQC,TRIM`. Default if left blank is to run all steps. Options:
  - `FRAG`: generate the fragment distribution and calculate some summary statistics. Run fragmentDistribution.r on group of 10 samples.
  - `BW`: Transform bam file to bigwig format and shift. Run only prepareBigWig.sh on sample.
  - `SCORE`: Calculate TSS enrichment score. Run only calcTssEnrich.sh on sample.
  - `COLLATE`: Collate TSS enrichment scores of all samples. Run only collateTSSscores.
- `GROUP_SAMPS` : for step `FRAG` , should be the number of samples divided by 10. e.g. If there are 50 samples, batch number should be 0-5, producing 5 batches of 10 samples each.

##### -requires-
- RefSeq file for TSS scores calculation. Should be specified as `REFTSS` in the *config.txt* file. 

### 3. Samples peak calling (stage 1 QC)

  `sbatch --array=<number of batch jobs> ./jobSubmission/3_batchRunPeakCalling.sh (project directory) [STEPS]`

This script performs the core part of the ATACseq pipeline: calling peaks. In this stage, peak calling is performed at sample level using the paired-end mode of MACS3. After this, FRIP statistics are calculated and collated in a single csv file.

##### -scripts executed-

- ./subScripts/samplePeaks.sh (array-number): MACS version 3 peak calling with BAM files paired-end reads
    * it then filters the peaks to exclude those that overlap with blacklisted regions and sex chromosomes.
    * peaks are sorted by chromosome.
- ./subScripts/collateCalcFrip.sh (array-number): Calculation of fraction of reads in peaks, number of reads and peaks for peak calling at sample level
    * a single --array job number should be specified (any number).

##### -parameters-
- `--array`: number of batch jobs, each number matches a sample. No needed for `FRIP`.
- `(project-name)` project's directory.
-optional-
- `[STEPS]` Option to combine steps with desired steps included as single string, i.e. `FASTQC,TRIM`. Default if left blank is to run all steps. Options:
  - `PEAKS`: Perform only peak calling on samples. Run only samplePeaks.sh on sample.
  - `FRIP`: Calculate peak calling results for all samples. Run only collateCalcFrip.sh on all samples found.
  
### 4. Stage 1 QC metrics summary

   `sbatch ./jobSubmission/4_collateStage1QCMetrics.sh (project directory) [STEPS]`

This scripts collates results of Stage 1 of QC pipeline to confirm data quality. It uses MultiQC to collate the output of fastQC, bowtie2 alignment, ATAC-seq QC metrics and sample peak calling results. 

##### -scripts executed-
- ./subScripts/progressReport.sh : Identification of samples that have been successful at each stage of the processing pipeline and for each fastq file, how far through the process it has progressed.
- ./subScripts/countMTcollateFS.sh : Counting of the number of reads aligned to MT chromosome and collates flagstat summary of- aligned sorted reads.
- ./Rscripts/collateDataQualityStats.Rmd (R config file): Output Rmarkdown report summarising stage 1 QC metrics: raw reads metrics, trimming, alignment and peak calling.

##### -parameters-
- `--array` : no need to be specified.
- `(project-name)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. MULIQC,COLLATE. Default if left blank is to run all of them.
  - `MULTIQC`: Collates fastQC, trimming and alignment statistics in a single report. Run multiQC to collate fastQC, trimming and alignment results.
  - `COLLATE`: Collate results from all samples to show progress of each sample through the pipeline in order to avoid missing steps. Run only progressReport.sh and countMTcollateFS.sh.
  - `SUMMARY`: Collates stage 1 QC and peak calling results in a single Rmarkdown report. Run only collateDataQualityStats.Rmd.

### 5. Sex check (stage 2 QC)

  `sbatch --array=<number of batch jobs> ./jobSubmission/5_batchPeaksSexChrs.sh (project directory) [STEPS]`
  
This script checks the sex of samples as part of their identity check by creating bam files for each sample containing only reads aligned to sex chromosomes and performing peak calling on these.

##### -scripts executed-
- ./subScripts/subsetSexChrs.sh (array-number): Subsetting reads of only X and Y chromosomes in bam file for input sample . 
- ./subScripts/sexChrPeaks.sh : Peak calling on the sex chromomes, filtering and read counts using MACS3 in paired-end mode.

##### -parameters-
- `--array` : number of batch jobs, each number matches a sample. Only for step 5.1 `SPLIT`.
- `(project-directory)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. MULIQC,COLLATE. Default if left blank is to run all of them.
  - `SPLIT`: Generate only samples with X and Y chromosomes reads. Run only subsetSexChrs.sh on sample.
  - `PEAKS`: Perform only peak calling on sex chromosomes reads. Run only sexChrPeaks.sh.

##### -requires-
- List of pseudoautosomal regions specified in *config.txt* as `PAR`.
- X chromosome gene list (FIRRE and XIST) specified in *config.txt* as `XCHRBED`.
  
### 6. Genotype check (stage 2 QC)

  `sbatch --array=<number of batch jobs> ./jobSubmission/6_batchRunGenotypeCheck.sh (project directory) [STEPS]`

This script checks the genotype of the samples as part of their identity check by detecting possible DNA contamination in order to ensure high quality sequence reads. If any contamination is detected, possible swaps will be suggested. 

##### -scripts executed-
- ./subScripts/compareBamWithGenotypes.sh (array-number): Preparation of bam file for comparison of assigned genotype of each sample and runs verifyBamID 
  * requires a file in METADATA folder called matchedVCFIDs.txt which lists the samples with their matched vcfID.
- ./subScripts/searchBestGenoMatch.sh : Output a summary of stats from previous step and finds any sample that might be contaminated.
  * Contaminated samples will go to a created file potentialSwitches.txt. If this file is not empty, searchBestGenoMatch.sh is executed and an alternative Genotype search is done for the sample.

##### -parameters-
- `--array` : number of batch jobs, each number matches a sample. 
- `(project-directory)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. GENCHECK,COMPARE. Default if left blank is to run all of them.
  - `COMPARE`: Compare only samples genotype with matched VCF file genotype. Run only compareBamWithGenotypes.sh on sample.
  - `SWITCH`: Indentify only mismatched genotype samples and output potential switches and seach for alternative genotype. Run only searchBestGenoMatch.sh.

### 7. Summary stage 2 QC

   `sbatch ./jobSubmission/7_collateStage2QCMetrics.sh (project directory) [STEPS]`

This scripts collates results of Stage 2 of QC pipeline to confirm samples identity. Sex and genotype check results are collated.

##### -scripts executed-

- ./subScripts/progressReportS2.sh : Identification of samples that have been through stage 2 of QC in order to avoid missing processing samples
- ./Rscripts/collateSampleChecks.Rmd (R config file): Output results from sex and genotype check in a Rmarkdown report.

##### -parameters-
- `--array` : no need to be specified.
- `(project-name)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. MULIQC,COLLATE. Default if left blank is to run all of them.
  - `COLLATE`: Collate results from all samples to show progress of each sample through stage 2 of QC pipeline in order to avoid missing steps or samples.
  - `SUMMARY`: Collates stage 2 QC results in a single Rmarkdown report. Run only collateSampleChecks.Rmd.

### 8. Group peak calling (Stage 3 QC)

 `sbatch ./jobSubmission/8_batchRunGroupPeakCalling.sh (project directory) [STEPS] [GROUP]`

This script performs stage 3 of QC to check cell type of samples by grouping samples per cell type, performing group level peak calling and gets read counts in peaks called. Samples used are those that have passed previous QC stages 1 and 2. Cell fractions to which samples belong to need to be specified in the *config.txt* file.

##### -scripts executed-
- ./Rscripts/samplesForGroupAnalysis.r (R config file) (GROUP): Creation of a file with a list of samples of cell group input.
- ./subScripts/groupPeaks.sh (GROUP): Peak calling on merged bam file for cell type group or on samples that belong to the cell group.
- ./Rscripts/countsInPeaksGroup.r (R config file): Obtain read counts in peaks called on each cell group and performs differential analysis to find topmost differential peaks.
- ./Rscripts/collateCellTypeCheck.Rmd (R config file): Output results of cell-groups peak calling and counts in peaks to perform cell type check.

##### -parameters-
- `(project-directory)` project's directory.
- `(GROUP)` Select cell type of group of samples to perform peak calling on. This is needed for the `PEAKS` and `MERGE` steps.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. FRAGSIZE,PEAKS. Default if left blank is to run all of them.
  - `MERGE`: Collates samples bam files and creates a sorted, merged bam file for cell type group.
  - `PEAKS`: Performs peak calling on bam file or samples for chosen cell type using MACS3 in paired-end mode. Run only groupPeaks.sh.
  - `COUNTS`: Gets read counts in peaks called at cell group level and analysis for topmost differential peaks. Run only countsInPeaksGroup.r.
  - `CTCHECK` : Collates results from before in a Rmarkdown report and identifies samples that fail cell type identity.

### 9. Peak set curation 

 `sbatch ./jobSubmission/9_batchPeakSetCuration.sh (project directory) [STEPS] [GROUP] [THRES]`

This script performs the analysis and curation of peak sets using high-quality samples in order to produce a curated peak set that can be used for downstream analyses. Cell fractions to which samples belong to need to be specified in the *config.txt* file.

##### -scripts executed-

- ./Rscripts/samplesForDownAnalysis.r (R config file) (GROUP): Creation of a file with a list of samples of cell group input.
- ./subScripts/peaksForCuration.sh (GROUP): Peak calling on merged bam file for cell type group or on samples that belong to the cell group.
- ./Rscripts/cutoffAnalysis.r (R config file) (GROUP): Perform segmented analysis to decide optimal peak calling threshold.
- ./subScripts/peaksStrThres.sh (GROUP) (THRES): Peak calling on merged bam file for cell type group or on samples that belong to the cell group with a chosen stringent significance threshold.

##### -parameters-
- `(project-directory)` project's directory.
- `(GROUP)` Select cell type of group of samples to perform peak calling on. This is needed for the `PEAKS`, `PEAKSTR`, and `MERGE` steps.
- `(THRES)` Select significance threshold for peak calling in step 9.5 `STRPEAKS` as suggested by results in `CUTOFF` step 9.4.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. FRAGSIZE,PEAKS. Default if left blank is to run all of them.
  - `MERGE`: Collates samples bam files and creates a sorted, merged bam file for cell type group.
  - `BW`: Converts merged bam file into bigwig file for genome visualisation.
  - `PEAKS`: Performs peak calling on bam file or samples for chosen cell type using MACS3 in paired-end mode. Run only peaksForCuration.sh.
  - `CUTOFF`: Performs segmented analysis to decide an optimal peak calling significance threshold. Run only cutoffAnalysis.r.
  - `PEAKSTR` : Performs peak calling on bam file or samples for chosen cell type using MACS3 in paired-end mode with a stringent significance threshold. Run only peaksStrThres.sh.

### 10. Downstream analysis preparation

 `sbatch ./jobSubmission/9_batchDownstreamAnalyses.sh (project directory) [STEPS] [GROUP] [THRES]`

This scripts prepares inputs for downstream analyses, especially counts in peaks matrices. Cell fractions to which samples belong to need to be specified in the *config.txt* file. 

##### -scripts executed-
- ./Rscripts/curatedPeaksChecks.Rmd (R config file): Output results of curated cell-groups peak calling.
- ./Rscripts/countsInCleanPeaksCT.r (R config file): Obtain read counts in peaks called on each cell group and performs differential analysis to find topmost differential peaks from peak sets used for curation (5e-2). This would be with only samples that passed all QC stages. Note STEP 9.3 should have run before this and created the peak sets.
- ./Rscripts/countsInCuratedPeaks.r (R config file) (GROUP) (THRES): Obtain read counts in peaks called in specified cell type and significance threshold. Note STEP 9.5 PEAKSTR should have run before this and created the peak set.

##### -parameters-
- `(project-directory)` project's directory.
- `(GROUP)` Select cell type of group of samples to perform peak calling on. This is needed for the `COUNTS` step.
- `(THRES)` Select significance threshold of peaks for counts in peaks step 10.3 `COUNTS`. 
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. FRAGSIZE,PEAKS. Default if left blank is to run all of them.
  - `SUMMARY`: Collates curated peak sets generated across cell types in a Rmarkdown report with  properties of peaks: number, genomic annotation, width, etc. Run curatedPeaksChecks.Rmd only.
  - `CTCLEAN`: Gets read counts in peaks called at cell group level and analysis for topmost differential peaks. Run only countsInCleanPeaksCT.r.
  - `COUNTS`: Gets read counts in peaks called in curated peak set in chosen cell type. Run only countsInCuratedPeaks.r.

# Citation

If you use PIPPA, please cite: 
