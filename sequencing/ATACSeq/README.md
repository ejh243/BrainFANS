
# ATAC-seq analysis pipeline 

This readme explains the content of the scripts for running the ATAC analysis pipeline

REQUISITES:
- Scripts submitted from the scripts/sequencing folder.
- Users require a username folder in the ATACSeq/logFiles directory.
- The name of the project folder must be specified as the first argument on the command line.
- Requires two config files in the ATACSeq/config directory:
  - config.txt with variables for bash scripts
  - config.r with parameters for R scripts
- Requires a samples.txt file in the METADATA folder with the names of the samples that will be used to run the pipeline
- Samples need to be in the RAW folder in the project folder (1_raw)

Parameters in [] are optional.
### 1. Pre-analysis


#### 1.1 `sbatch --array="number of batch jobs" ATACSeq/jobSubmission/1_batchRunAlignment.sh (project name) [STEPS]`

Performs pre-analysis of ATAC-seq data, including pre-alignment quality control, alignment and post-alignment quality control of samples. 

##### -scripts executed-
- preScripts/fastqc.sh : FastQC for pre-alignment quality control.
- preScripts/fastp.sh : FastP for trimming samples.
- ATACSeq/preprocessing/alignment.sh : Alingment of samples to reference genome using Bowtie2.
- ATACSeq/preprocessing/calcENCODEQCMetrics.sh : ENCODE QC metrics are calculated on aligned samples.
  
##### -parameters-
- `--array`: number of batch jobs, each number matches a sample.
- `(project-name)` the name of the project folder in the main data type data directory
-optional-
- `[STEPS]` Option to combine steps with desired steps included as single string, i.e. `FASTQC,TRIM`. Default if left blank is to run all steps. Options:
  - `FASTQC`: Perform only quality control. Run only fastqc.sh on samples.
  - `TRIM`: Perform only trimming on samples. Run only fastp.sh on samples.
  - `ALIGN`: Align samples. Run only alignment.sh on samples.
  - `ENCODE`: Calculation of ENCODE QC metrics.

### 2. Post-alignment processing

#### 2.1 `sbatch --array=<number of batch jobs/10> ATACSeq/jobSubmission/2_batchCalcQCMetrics.sh  (project-name)`

This script uses the ATACseqQC R package to generate the fragment distribution and calculate some summary statistics. It splits the samples into groups of 10 to run in parallel, so batch number should be the number of samples divided by 10.
e.g. If there are 50 samples, batch number should be 0-5, producing 5 batches of 10 samples each.
If the number of samples is less than 10, set batch number to 0.
* executes ATACSeq/preprocessing/3_fragmentDistribution.r (aligned-dir) (array-number)
 
### 3. Peak calling by sample

#### 3.1 `sbatch --array=<number of batch jobs> ATACSeq/jobSubmission/3_batchRunPeakCalling.sh (project-name) [STEPS]`

This script searchs within aligned data folders for all bam files
Then submits a batch script to process one sample
* executes ATACSeq/preprocessing/4_shiftAlignedReadsPE.sh 
	* takes a filtered bam file converts to a tagalign file, calculate CC scores and shifts reads ready for peak calling
* executes ATACSeq/preprocessing/5_samplePeaks.sh
	* which runs MAC peak calling with shifted tagAlign files and with BAM files and paired end reads
	* it then filters the peaks to exclude those that overlap with blacklisted regions
* executes ATACSeq/preprocessing/6_calcFrip.sh 
	* calculates fraction of reads in peaks for each set of peaks
	
##### -optional- 

[STEPS] option to either SHIFT, PEAKS or FRIP. The first one shift reads, the second one performs peak calling and the third calculates fraction of reads in peaks. 

#### 3.2 `sbatch ATACSeq/jobSubmission/4_collateStage1QCMetrics.sh  (project-name) [STEPS]`

This scripts uses MultiQC to collate the output of fastqc and bowtie2 alginment. It also runs the following utility scripts.
* executes ATACSeq/preprocessing/7_progressReport.sh 
	* identifies how many samples have been successful at each stage of the processing pipeline and for each fastq file, how far through the process it has progressed.

* executes ATACSeq/preprocessing/9_countMTReads.sh
	* collates counts of the number of reads aligned to MT chromosome

* executes ATACSeq/preprocessing/9_collateFlagStatOutput.sh 
	* collates flagstat summary of aligned sorted reads.

* executes ATACSeq/preprocessing/collateDataQualityStats.Rmd
	* generates Rmarkdown report summarising stage 1 qc metrics

  !FILTER option not available 
	#* executes ATACSeq/preprocessing/10.2_filterOnS1SumStats.r [min. no reads] [min. Alignment rate] [min. no filtered reads]
		#* outputs a list file of samples that have passed the set thresholds. If optional parameters not specified, defaults are 10, 80 and 20.

##### -optional- 

[STEPS] option of either MULTIQC, COLLATE, SUMMARY. They may be combined, with desired steps included as single string, i.e. MULIQC,COLLATE. Default if left blank is to run all of them.

#### 3.3 `sbatch --array=<number of batch jobs> ATACSeq/jobSubmission/5_batchFormatSexChrs.sh <project ID>`

	* executes ./ATACSeq/preprocessing/subsetSexChrs.sh for all samples that pass stage 1 filtering
		* subsets tagalign files to X and Y chromosome reads
	* requires a file in METADATA folder called passStage1SampleList.txt which lists the samples to be included for sex check
 
#### 3.4 `sbatch ATACSeq/jobSubmission/6_sexCheck.sh  <project ID>`
	* executes ATACSeq/preprocessing/sexChrPeaks.sh which performs peak calling on the sex chromomes, filter and read counts
 
### 4. Genotype concordance

#### 4.1 `7.1 sbatch --array=<number of batch jobs> ATACSeq/jobSubmission/7_batchRunGenotypeConcordance.sh <project ID> [OPTIONS]`
  * executes /ATACSeq/preprocessing/compareBamWithGenotypes.sh which prepares bam file for comparison with verifyBamID
    * creates an recalibration file for base quality scores to be adjusted.
  * requires a file in METADATA folder called matchedVCFIDs.txt which lists the samples with their matched vcfID.
  * --array: starts from 2
  
  -optional-
  
  [OPTIONS] option of GENCHECK and COMPARE. COMPARE performs comparison of bam file with genotype. GENCHEK to collate results from the previous steps.
  
  * executes /ATACSeq/preprocessing/collateSampleChecks.Rmd which collates the results from previous sex check and Genotype check
    
#### 4.2 `7.2 sbatch --array=<number of batch jobs> ATACSeq/jobSubmission/7_2_batchRunGenotypeSearch.sh <project ID> `
  * Outputs a summary of stats from previous step and finds any sample that might be contaminated.
  * Contaminated samples will go to a created file potentialSwitches.txt. If this file is not empty, searchBestGenoMatch.sh is executed and an alternative Genotype search is done for the sample.

### 5. Peak calling by group

#### 5.1 `8. sbatch --array=<number of cell fractions> ATACSeq/jobSubmission/8_batchPeakCallingByGroup.sh <project ID> [GROUPS] [STEPS]`
* executes /general/processing/makeGroupAnalysisFile.r which creates a txt file (samplesForGroupAnalysis.txt) with samples classified by fraction/tissue
    * tissue can be specified, default is "prefrontal cortex|PFC"
    
* executes /ATACSeq/preprocessing/groupPeaks.sh which performs peak calling on the subsets of samples specified (<PEAK> option)
  
* executes /ATACSeq/preprocessing/calcFripGroup.sh which calculates fraction of reads in peaks in subsets of samples specified (<FRIP> option)
    
  -optional-
  
  [GROUPS] Peak calling by group can be done either by PASS (all samples that passed Stage 1 Quality control) or by FRACTION (samples grouped by their fraction).
  [STEPS] option to either PEAK or FRIP. The first performs peak calling, the second calculates fraction of reads in peaks between and within subsets.
  
### 6. Advanced analysis

#### Still to be developed

7. idr analysis against sample peaks

8. pca analysis of peaks 

#### Misc data preprocessing

mergeFastq.sh
Where a sample was procesed multiple times, the fastq files are merged. 

#### OUTPUT files

| Script | Filename | Description | Downstream processing |
| ------ | -------- | ----------- | --------------------- |
| fastqc.sh | * | FastQC report | Input for QC summary |
| fastp.sh | *_trimmed.f*q | Trimmed raw sequence file output from fastp input | Input for alignment |
| 1_alignment.sh | ${ALIGNEDDIR}/*_statsperchr.txt | Counts of number of reads mapped to each chromosome before any filtering | Input to CountMTReads.sh |
| 1_alignment.sh | ${ALIGNEDDIR}/*_sorted.bam | Aligned reads, sorted and indexed | Input to calcATACQCMetrics.r |
| 1_alignment.sh | ${ALIGNEDDIR}/*_sorted_chr1.bam | Aligned chr 1 reads, sorted and indexed | Input to calcENCODEQCMetricsPE.sh |
| 1_alignment.sh | 
${ALIGNEDDIR}/*_dupMetrics.txt | Output of Picard with duplication statistics | ? |
| 1_alignment.sh | ${ALIGNEDDIR}/*_depDup_q30.bam | Aligned reads with MT chr, low quality (q < 30),improperly paired reads, pcr optical duplicates, secondary alignments excluded, sorted and indexed  | Input to peak calling |
| 1_alignment.sh | ${ALIGNEDDIR}/*_postFilter_statsperchr.txt | Counts of number of reads mapped to each chromosome after filtering  |  |
 
