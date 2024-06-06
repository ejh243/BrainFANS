#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=CNVcalling_%j.log
#SBATCH --error=CNVcalling_%j.err
#SBATCH --job-name=CNVcalling


## print start date and time
echo Job started on:
date -u


source $1 || exit

module load Perl/5.32.1-GCCcore-10.3.0
module load $RVERS

cd ${SCRIPTDIR}/preprocessing 

sh 8_pennCNV.sh
echo CNVs called on full sample

sh 9_filterCNVCalls.sh

module load Pandoc
R -e "rmarkdown::render('10_summarizeCNVCalls.rmd', output_file = 'SummarizeCNVCalls.html')" --args ${CNVDIR}/PennCNVOutput/${FILEPREFIX}_GCModel_MergedFiltered_AnnoGencodev39.rawcnv NULL ${CNVLIST} $USER

mv SummarizeCNVCalls.html ${CNVDIR}PennCNVOutput/${FILEPREFIX}_SummarizeCNVCalls.html

## print finish date and time
echo Job finished on:
date -u

mkdir -p ${DATADIR}/logFiles

mv "${SLURM_SUBMIT_DIR}/CNVcalling_${SLURM_JOB_ID}.log" \
"${DATADIR}/logFiles/CNVcalling_${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/CNVcalling_${SLURM_JOB_ID}.err" \
"${DATADIR}/logFiles/CNVcalling_${SLURM_JOB_ID}.err"