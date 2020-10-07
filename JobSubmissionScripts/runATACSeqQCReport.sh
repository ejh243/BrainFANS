#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=150:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=e.j.hannon@exeter.ac.uk # email me at job completion
#SBATCH --error=LogFiles/ATAQQCReport.err # error file
#SBATCH --output=LogFiles/ATAQQCReport.log # output file
#SBATCH --job-name=ATAQQCReport

## needs to be executed from the scripts folder


## print start date and time
echo Job started on:
date -u


####### 
## NOTE: Do not store confidential information in this file use the config file
######

echo "Changing Folder to: "
echo $SLURM_SUBMIT_DIR

cd $SLURM_SUBMIT_DIR


## generate QC metrics
module purge
module load R/3.6.3-foss-2020a
Rscript ATACSeq/calcATACSeqQCMetrics.r ATACSeq/config.r

## produce QC report
module load Pandoc
Rscript -e "rmarkdown::render('ATACSeq/ATACSeqQCReport.rmd', params = list(configFile='ATACSeq/config.r'), output_file='../ATACSeq/RawData/QCOutput/QC.html')"