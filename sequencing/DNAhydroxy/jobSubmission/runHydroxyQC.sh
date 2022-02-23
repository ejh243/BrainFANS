#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --error=LogFiles/CEGX5hmCQC.err # error file
#SBATCH --output=LogFiles/CEGX5hmCQC.log # output file
#SBATCH --job-name=CEGX5hmCQC



## print start date and time
echo Job started on:
date -u

## load config parameters
source $1

module load MACS2
module load Miniconda2
source activate epic2
module load HTSeq

##  call peaks on sex chromosomes across all samples and do peak quanitification
./sexChrPeaks.sh

## calculate reads in gene body
./geneBodyCounts.sh 

## Run QC report
module load Pandoc/2.5
module load R/3.6.3-foss-2020a

Rscript calcQCMetrics.r config.r 

Rscript -e "rmarkdown::render('QCReport.rmd', params = list(config='config.r'), output_file='QCJuly2021.html')"


## print finish date and time
echo Job finished on:
date -u