#!/bin/bash
#SBATCH --export=ALL #export all enviroment variables to the batch job
#SBATCH -p mrcq #submit to the serial queue
#SBATCH --time=24:00:00 ##maximum wall time for the job
#SBATCH -A Research_Project-MRC190311 #research project to submit under
#SBATCH --nodes=1 #specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --error=QCDNAdata_%j.err # error file
#SBATCH --output=QCDNAdata_%j.log # output file
#SBATCH --job-name=QCDNAdata

#------------------------------------------------------

# 1. the first input is the full path to config file.

#-----------------------------------------------------

print_error_message() {
    IFS=$'\n'
    error_message=$*
    RED='[0;31m'
    NO_COLOUR='[0m'

cat 1>&2 << MESSAGE
${RED}
ERROR:
${error_message}
${NO_COLOUR}
MESSAGE
    exit 1
}

## print start date and time
echo Job started on:
date -u

echo Log files intially stored in: ${SLURM_SUBMIT_DIR}/QCDNAdata_${SLURM_JOB_ID}.log and ${SLURM_SUBMIT_DIR}/QCDNAdata_${SLURM_JOB_ID}.err

config_file=$1
source "${config_file}" || print_error_message \
    "The provided config file was not sourced correctly." \
    "Please check the path you gave ('$config_file') exists, exiting..." 

echo "Processing data located in :" ${DATADIR}

## load modules
echo "Loading R module :" $RVERS
module load Pandoc
module load $RVERS   # load specified R version

cd ${SCRIPTSDIR}/array/DNAm/preprocessing/

Rscript checkRconfigFile.r ${DATADIR}
config_malformed=$?
if [[ "${config_malformed}" -ne 0 ]]; then
    print_error_message \
        "Malformed config file has been identified." \
        "Please check the error logs for more information, exiting..."
fi

Rscript installLibraries.r ${DATADIR}
library_did_not_install=$?
if [[ "${library_did_not_install}" -ne 0 ]]; then
    print_error_message \
        "A required library did not install properly." \
        "Please check the error logs as to why this happened." \
        "If the problem is not easily fixed, consider opening an issue." \
        "https://github.com/ejh243/BrainFANS/issues/new/choose" \
        "Exiting..."
fi

Rscript checkColnamesSampleSheet.r ${DATADIR}
sample_sheet_malformed=$?
if [[ "${sample_sheet_malformed}" -ne 0 ]]; then
    print_error_message \
        "Malformed sample sheet has been identified." \
        "Please check the error logs for more information, exiting..."
fi

mkdir -p ${GDSDIR}/QCmetrics

Rscript loadDataGDS.r ${DATADIR}
gds_problem_identified=$?
if [[ "${gds_problem_identified}" -ne 0 ]]; then
    print_error_message \
        "A problem with the GDS data has been identified." \
        "Please check the error logs for more information, exiting..."
fi

chmod 755 ${DATADIR}/2_gds/raw.gds

Rscript calcQCMetrics.r ${DATADIR} ${REFDIR}

Rscript clusterCellTypes.r ${DATADIR} ${REFDIR}

Rscript -e "rmarkdown::render('QC.rmd', output_file='QC.html')" --args ${DATADIR} ${REFDIR} ${RCONFIG} $USER

## mv markdown report to correct location
mv QC.html ${GDSDIR}/QCmetrics/

mkdir -p ${DATADIR}/3_normalised

Rscript normalisation.r ${DATADIR} ${REFDIR}
chmod 755 ${DATADIR}/2_gds/rawNorm.gds

mkdir -p ${GDSDIR}/QCmetrics/CETYGO

Rscript CETYGOdeconvolution.r ${DATADIR}

## print finish date and time
echo Job finished on:
date -u

mkdir -p ${DATADIR}/logFiles

mv "${SLURM_SUBMIT_DIR}/QCDNAdata_${SLURM_JOB_ID}.log" \
"${DATADIR}/logFiles/QCDNAdata_${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/QCDNAdata_${SLURM_JOB_ID}.err" \
"${DATADIR}/logFiles/QCDNAdata_${SLURM_JOB_ID}.err"

echo Log files moved to: ${DATADIR}/logFiles/
