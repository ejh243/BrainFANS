# ATAC-seq analysis pipeline 

This readme explains the content of the scripts for running the ATAC analysis pipeline

REQUISITES:
- Scripts submitted from the scripts/sequencing folder.
- The name of the project folder must be specified as the first argument on the command line.
- Requires config files (file paths should be specified in the first file config.txt):
  - config.txt with variables for bash scripts
  - config.r with parameters for R scripts
  - packagesPip.txt : list of packages and their version to be installed in the created pip virtual environment.
  - environment.yml : file with packages and their version to be installed by conda.
- Requires a samples.txt file in the METADATA folder with the names of the samples that will be used to run the pipeline
- Samples need to be in the RAW folder in the project folder (1_raw)
- Project directory needs to be specified in the command line (project directory)

Parameters in [] are optional.


### 0. Set-up

 `sbatch --array=<number of batch jobs> ATACSeq/jobSubmission/0_setUp.sh (project directory)`
 
This script checks for required files, packages or libraries that are needed later in the pipeline.

##### -scripts executed-
- ATACSeq/preprocessing/intallLibraries.sh : checks for already installed R libraries and install if not found.

##### -parameters-
- `--array`: should be 0 for this script as general analysis is run rather than individual samples.
- `(project-name)` project's directory.

##### -requires-
- packagesPip.txt, environment.yml and config.r file paths to be specified in the config.txt file.

### 1. Pre-analysis (QC and alignment)

 `sbatch --array=<number of batch jobs> ATACSeq/jobSubmission/1_batchRunPreAnalysis.sh (project directory) [STEPS]`

Performs pre-analysis of ATAC-seq data, including pre-alignment quality control, alignment and post-alignment quality control of samples. 

##### -scripts executed-
- preScripts/fastqc.sh : FastQC for pre-alignment quality control.
- preScripts/fastp.sh : FastP for trimming samples.
- ATACSeq/preprocessing/alignment.sh : Alingment of samples to reference genome using Bowtie2.
- ATACSeq/preprocessing/calcENCODEQCMetrics.sh : ENCODE QC metrics are calculated on aligned samples.
  
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

 `sbatch --array=<number of batch jobs/10> ATACSeq/jobSubmission/2_batchCalcQCMetrics.sh  (project directory)`

This script uses the ATACseqQC R package to generate the fragment distribution and calculate some summary statistics. It splits the samples into groups of 10 to run in parallel
If the number of samples is less than 10, set batch number to 0.

##### -scripts executed-  
- ATACSeq/preprocessing/fragmentDistribution.r (R config file) (array-number) : fragment distribution of samples specified by array. 

##### -parameters-
- `--array`: number of batch jobs, should be the number of samples divided by 10. e.g. If there are 50 samples, batch number should be 0-5, producing 5 batches of 10 samples each.
- `(R config file)` directory of R config file


### 3. Peak calling

  `sbatch --array=<number of batch jobs> ATACSeq/jobSubmission/3_batchRunPeakCalling.sh (project directory) [STEPS]`

This script performs the core part of the ATACseq pipeline: calling peaks. In this stage, peak calling is performed at sample level using the Paired-end mode of MACS3. After this, FRIP statistics are calculated and collated in a single file.

##### -scripts executed-
- ATACSeq/preprocessing/shiftAlignedReads.sh : takes a filtered bam file converts to a tagalign file, calculate CC scores and shifts reads ready for peak calling
    * This step is not needed if peak calling is performed using PE mode or with HMMRATAC
- ATACSeq/preprocessing/samplePeaks.sh : runs MAC version 3 peak calling with BAM files paired end reads
		* it then filters the peaks to exclude those that overlap with blacklisted regions, sex chromosomes 
    * peaks are sorted by chromosome
- ATACSeq/preprocessing/collateCalcFrip.sh : calculates fraction of reads in peaks, number of reads and peaks for peak calling at sample level

##### -parameters-
- `--array`: number of batch jobs, each number matches a sample.
- `(project-name)` project's directory.
-optional-
- `[STEPS]` Option to combine steps with desired steps included as single string, i.e. `FASTQC,TRIM`. Default if left blank is to run all steps. Options:
  - `SHIFT`: Perform only shifting of reads.
  - `PEAKS`: Perform only peak calling on samples. 
  - `FRIP`: Calculate peak calling results for all samples.
  
### 4. Stage 1 QC metrics summary

   `sbatch --array=0 ATACSeq/jobSubmission/4_collateStage1QCMetrics.sh  (project directory) [STEPS]`

This scripts uses MultiQC to collate the output of fastqc and bowtie2 alginment, as well as peak calling results. 

##### -scripts executed-
- ATACSeq/preprocessing/progressReport.sh : identifies how many samples have been successful at each stage of the processing pipeline and for each fastq file, how far through the process it has progressed.
- ATACSeq/preprocessing/countMTReads.sh : collates counts of the number of reads aligned to MT chromosome
- ATACSeq/preprocessing/collateFlagStatOutput.sh : collates flagstat summary of aligned sorted reads.
- ATACSeq/preprocessing/collateDataQualityStats.Rmd : generates Rmarkdown report summarising stage 1 qc metrics: raw reads metrics, trimming, alignment and peak calling

##### -parameters-
- `--array` : should be 0 for this script as general analysis is run rather than individual samples.
- `(project-name)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. MULIQC,COLLATE. Default if left blank is to run all of them.
  - `MULTIQC`: Collates fastqc and alignment statistics in a single report.
  - `COLLATE`: Collate results from all samples to show progress of each sample through the pipeline in order to avoid missing steps.
  - `SUMMARY`: Collates stage 1 QC and peak calling results in a single Rmarkdown report.
  - `BATCH`: If different batches of samples are found, results related to these are compared.


### 5. Sex chromosomes 

  `sbatch --array=<number of batch jobs> ATACSeq/jobSubmission/5_batchFormatSexChrs.sh (project directory) [STEPS]`
  
This script creates bam files for each sample containing only reads aligned to sex chromosomes, perform peak calling on these and output results.

##### -scripts executed-
- ATACSeq/preprocessing/subsetSexChrs.sh : subsets reads of only X and Y chromosomes for input sample in bam file. 
	* requires a file in METADATA folder called passStage1SampleList.txt which lists the samples to be included for sex check
- ATACSeq/preprocessing/sexChrPeaks.sh : performs peak calling on the sex chromomes, filter and read counts using MACS3 in Single-end mode.
- ATACSeq/preprocessing/collateSexChecks.r : collates peak calling results on sex chromosomes and uses it to check the assigned sex of the sample.

##### -parameters-
- `--array` : number of batch jobs, each number matches a sample.
- `(project-directory)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. MULIQC,COLLATE. Default if left blank is to run all of them.
  - `SPLIT`: Generate only samples with X and Y chromosomes reads.
  - `PEAKS`: Perform only peak calling on new sex chromosomes reads samples.
  - `CHECK`: Collate peak calling on sex chromosomes results and check assigned sex of each sample.

### 6. Genotype check

  `sbatch --array=<number of batch jobs> ATACSeq/jobSubmission/6_batchRunGenotypeConcordance.sh (project directory) [STEPS]`

This script will detect possible DNA contamination in order to ensure high quality sequence reads. If any contamination is detected, possible swaps will be suggested. 

##### -scripts executed-
- ATACSeq/preprocessing/compareBamWithGenotypes.sh  : prepares bam file for comparison of assigned genotype of each sample and runs verifyBamID 
  * requires a file in METADATA folder called matchedVCFIDs.txt which lists the samples with their matched vcfID.
- ATACSeq/preprocessing/collateSampleChecks.Rmd : collates results from previous steps (sex and genotype check).
- ATACSeq/preprocessing/searchBestGenoMatch.sh : Outputs a summary of stats from previous step and finds any sample that might be contaminated.
  * Contaminated samples will go to a created file potentialSwitches.txt. If this file is not empty, searchBestGenoMatch.sh is executed and an alternative Genotype search is done for the sample.
  
##### -parameters-
- `--array` : number of batch jobs, each number matches a sample. Starting from 2
- `(project-directory)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. GENCHECK,COMPARE. Default if left blank is to run all of them.
  - `COMPARE`: Compare only samples genotype with matched VCF file genotype.
  - `GENCHECK`: Sex and gencheck results are collated in a report.
  - `SWITCH`: Indentify only mismatched genotype samples and output potential switches and seach for alternative genotype.


#### Still to be developed ####


### 7. Group peak calling

 7.0) `sbatch --array=<number of batch jobs> ATACSeq/jobSubmission/7_helperTools.sh (project directory) [GROUPS] [STEPS]`
 
This script includes several MACS functionalities/tools and samtools to perform different functionalities, such as produce samples with a subset of reads or create pseudo-replicates. 


 7.1) `sbatch --array=<number of batch jobs> ATACSeq/jobSubmission/7_1_batchPeakCallingByGroup.sh (project directory) [GROUPS] [STEPS]`

This script repeats peak calling but on groups of samples, which can be grouped by cell fraction or if they have passed the first stage of quality control.  

##### -scripts executed-
- /general/processing/makeGroupAnalysisFile.r : creates a txt file (samplesForGroupAnalysis.txt) with samples classified by fraction/tissue
    * tissue can be specified, default is "prefrontal cortex|PFC"
- ATACSeq/preprocessing/groupPeaks.sh : performs peak calling on the subsets of samples specified 
- ATACSeq/preprocessing/calcFripGroup.sh : calculates fraction of reads in peaks in subsets of samples specified. 

##### -parameters-
- `--array` : cell fraction corresponding number.
- `(project-directory)` project's directory.
-optional-
- `[STEPS]` They may be combined, with desired steps included as single string, i.e. GENCHECK,COMPARE. Default if left blank is to run all of them.

  7.2) `sbatch --array=<number of batch jobs> ATACSeq/jobSubmission/7_2_batchPeakCallingExtra.sh (project directory) [GROUPS] [STEPS]`

### 8. Advanced analysis

  8.1) `sbatch --array=<number of batch jobs> ./sequencing/ATACSeq/jobSubmission/8_1_idrAnalysis.sh <projectName> [GROUPS] [STEPS]`
  
  This script performs IDR analysis.
  
  8.2) `sbatch --array=<number of batch jobs> ./sequencing/ATACSeq/jobSubmission/8_2_diffAnalysis.sh <projectName> [GROUPS] [STEPS]`
  
  This script performs counts in peaks, various ways of normalisation of peaks and differential accessibility analysis.
  

