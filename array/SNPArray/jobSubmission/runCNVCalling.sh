#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=SNPArray/logFiles/CNVcalling.o
#SBATCH --error=SNPArray/logFiles/CNVcalling.e
#SBATCH --job-name=CNVcalling


## print start date and time
echo Job started on:
date -u


####### 
## NOTE: Do not store confidential information in this file use the config file
######

source $1

module load Perl/5.32.1-GCCcore-10.3.0
module load R

cd ${DATADIR}/scripts ## all scripts are written to be executed from the scripts folder

sh ${SCRIPTDIR}/preprocessing/8_pennCNV.sh
echo CNVs called on full sample

cd ${SCRIPTDIR}/
sh preprocessing/9_filterCNVCalls.sh

module load Pandoc
R -e "rmarkdown::render('summarizeCNVCalls.rmd', params = list(fileName = '../../SNPdata/CNV/PennCNVOutput/SCZ_GCModel_MergedFiltered_AnnoGencodev29.rawcnv', subset = '../../SNPdata/Merged/MRCSCZsamples.txt'), output_file = '../../SNPdata/CNV/SummarizeCNVCallsBackgroundAll_SCZSamples.html')"
R -e "rmarkdown::render('summarizeCNVCalls.rmd', params = list(fileName = '../../SNPdata/CNV/PennCNVOutput/SCZ_GCModel_MergedFiltered_AnnoGencodev29.rawcnv', subset = '../../SNPdata/Merged/SFARIFoetalSamples.txt'), output_file = '../../SNPdata/CNV/SummarizeCNVCallsBackgroundAll_SFARISamples.html')"
R -e "rmarkdown::render('summarizeCNVCalls.rmd', params = list(fileName = '../../SNPdata/CNV/PennCNVOutput/EURonly/SCZ_GCModel_MergedFiltered_AnnoGencodev29.rawcnv', subset = '../../SNPdata/Merged/MRCSCZsamples.txt'), output_file = '../../SNPdata/CNV/SummarizeCNVCallsBackgroundEUR_SCZSamples.html')"
R -e "rmarkdown::render('summarizeCNVCalls.rmd', params = list(fileName = '../../SNPdata/CNV/PennCNVOutput/EURonly/SCZ_GCModel_MergedFiltered_AnnoGencodev29.rawcnv', subset = '../../SNPdata/Merged/SFARIFoetalSamples.txt'), output_file = '../../SNPdata/CNV/SummarizeCNVCallsBackgroundEUR_SFARISamples.html')"

## print finish date and time
echo Job finished on:
date -u
