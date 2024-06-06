This readme explains the content of the scripts for running the genotype quality control pipeline.
It assumes the idat files have been processed through Genome Studio and exported as plink format
For CNV calling it assumes that the log R and B allele ratio have been exported and are included in the final report 

It is requires a config file. There is a template for this in the config folder.

#### Data pre-processing:

1. sbatch runSNPQC.sh /path/to/configfile
Executes the following scripts
	1_QC.sh performs basic QC of samples and variants
	2_CheckEthnicity.sh takes QCd SNP data, merges with 1000 genomes, classifies samples to super populations
	3_CheckRelatedness.sh 
	4_formatForImputation.sh prepares QCd data for imputation via Michegan imputation server

#### Post processing of imputed genotypes:

2. sbatch liftoverImputation.sh /path/to/configfile
Executes for each chromosome in parallel (as a separate job)
	5_liftoverhg38.sh requires output from Michegan imputation server to liftover to hg 38


3. sbatch runSNPImputationQC.sh /path/to/configfile
Executes
	6_summarizeImputation.r creates plots to summarise quality metrics of imputation
	7_combineImputationOutput.sh filters variants to those iwth rsq > 0.3 and maf > 0.01 and collates QC'd, lifted over, imputed output into a single plink file and single vcf file
	
#### CNV calling

4. sbatch runCNVCalling.sh  /path/to/configfile
Executes
  8_pennCNV.sh takes final report and reformats for Penn CNV and runs CNV callings
  9_filterCNVCalls.sh performs both CNV and sample filtering of CNV calls and annotates with genes
  10_summarizeCNVCalls.rmd creates report of CNV calls