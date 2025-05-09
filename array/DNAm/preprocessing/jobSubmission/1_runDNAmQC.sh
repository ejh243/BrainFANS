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


if [[ -z "${DNAM_CONDA_ENVIRONMENT}" ]]; then
    echo "Conda environment does not exist, please run the setup script first."
    echo "The setup script can be found at ${SCRIPTSDIR}/array/DNAm/Setup/setup.sh"
    exit 1
fi
    

echo "Processing data located in :" ${DATADIR}

source "${CONDA_SHELL}"
conda activate "${DNAM_CONDA_ENVIRONMENT}"

cd ${SCRIPTSDIR}/array/DNAm/preprocessing/

Rscript installPackages.r

Rscript checkRconfigFile.r ${RCONFIG}
config_malformed=$?
if [[ "${config_malformed}" -ne 0 ]]; then
    print_error_message \
        "Malformed config file has been identified." \
        "Please check the error logs for more information, exiting..."
fi

Rscript checkColnamesSampleSheet.r ${RCONFIG} ${DATADIR}
sample_sheet_malformed=$?
if [[ "${sample_sheet_malformed}" -ne 0 ]]; then
    print_error_message \
        "Malformed sample sheet has been identified." \
        "Please check the error logs for more information, exiting..."
fi

mkdir -p ${GDSDIR}/QCmetrics

Rscript loadDataGDS.r ${DATADIR} ${RCONFIG}
gds_problem_identified=$?
if [[ "${gds_problem_identified}" -ne 0 ]]; then
    print_error_message \
        "A problem with the GDS data has been identified." \
        "Please check the error logs for more information, exiting..."
fi

chmod 755 ${DATADIR}/2_gds/raw.gds

Rscript calcQCMetrics.r ${DATADIR} ${REFDIR} ${RCONFIG}

Rscript clusterCellTypes.r ${DATADIR} ${REFDIR} ${RCONFIG}

most_recent_git_tag=$(git describe --tags --abbrev=0)
current_commit_hash=$(git rev-parse HEAD)
Rscript -e "rmarkdown::render('QC.rmd', output_file='QC.html')" \
    --args "${DATADIR}" "${REFDIR}" "${RCONFIG}" "${USER}" \
    "${most_recent_git_tag}" "${current_commit_hash}"

## mv markdown report to correct location
mv QC.html ${GDSDIR}/QCmetrics/

mkdir -p ${DATADIR}/3_normalised

Rscript normalisation.r ${DATADIR} ${REFDIR} ${RCONFIG}
chmod 755 ${DATADIR}/2_gds/rawNorm.gds

mkdir -p ${GDSDIR}/QCmetrics/CETYGO

Rscript CETYGOdeconvolution.r ${DATADIR} ${RCONFIG}

conda deactivate

## print finish date and time
echo Job finished on:
date -u

mkdir -p ${DATADIR}/logFiles

mv "${SLURM_SUBMIT_DIR}/QCDNAdata_${SLURM_JOB_ID}.log" \
"${DATADIR}/logFiles/QCDNAdata_${SLURM_JOB_ID}.log"
mv "${SLURM_SUBMIT_DIR}/QCDNAdata_${SLURM_JOB_ID}.err" \
"${DATADIR}/logFiles/QCDNAdata_${SLURM_JOB_ID}.err"

echo Log files moved to: ${DATADIR}/logFiles/
