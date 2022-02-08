This readme explains the content of the scripts for running the WGBS analysis pipeline

Scripts should be submitted from the scripts/sequencing folder. 
The ${DIR} variable of the config file should be edited to be specific to the dataset.

#### Data pre-processing:

1. 0_preparatory.sh

2. sbatch --array=<number of batch jobs> WGBS/jobSubmissionScripts/batchRunAlignment.sh <config.txt> <STEPS>

	This script searchs within raw data folders for all fastq files
	Then submits a batch script to process one set of paired samples at a time

	<STEPS>=option of either FASTQC or ALIGN. The former runs only the script 1_qcRawData.sh, fastqc-ing and trimming files. The latter runs only the alignment. Default if left blank is to run the whole script. 

	This script:
		* executes preScripts/1_qcRawData.sh
			* finds paired fastq files
			* runs fastqc
			* runs fastp or trimgalore dependent on data type

		* executes WGBS/preprocessing/2_alignment.sh 

3. sbatch 


