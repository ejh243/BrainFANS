#!/bin/bash
## ===================================================================================================================##
##                               ATAC-seq pipeline STEP 1.2: Pre-analysis -- trimming                                 ##
## ===================================================================================================================##
## EXECUTION: sbatch ./subScripts/fastp.sh <sampleName> <R1> <R2>                                                     ||
## - execute from pipeline's subScripts directory                                                                     ||
##                                                                                                                    ||
## DESCRIPTION: This script performs trimming of raw reads and following quality control of the new trimmed reads.    ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - Variables in config file: RAWDATADIR, TRIM_DIR, ADAPTERS_FILE                                                    ||
## - Software: fastp, fastqc (in conda environment)                                                                   ||
## - Read 1 and Read 2 of the same sample in the same directory (RAWDATADIR)                                          ||
## - A fasta file (.fa) with all adapters or sequences to be trimmed. Path to file should be specified in config file ||
##   as "ADAPTERS_FILE"                                                                                               ||
## - The length of sequences should be specified in the config.txt file as "MAX_LEN_R1" and "MAX_LEN_R2" for the max  ||
##   length of reads and "MIN_LEN" for the minimum.                                                                   ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <sampleName> Name of sample specified in command line.                                                       ||
## $2 -> <R1> Read 1 name                                                                                             ||
## $3 -> <R2> Read 2 name                                                                                             ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           || 
##  *_fastp.html, *_trimmed.f, *_fastp.json                                                                           ||
##                                                                                                                    ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##


sampleName=$1

## Get read 1 and read 2 for specified sample name
f1=$(basename $2)
f2=$(basename $3)

cd ${RAWDATADIR}

## create output filenames
outf1=${f1/.f/_trimmed.f}
outf2=${f2/.f/_trimmed.f}

mkdir -p ${TRIM_DIR}/qc/

## ============ ##
##    TRIM      ##
## ============ ##

echo "TRIM on ${sampleName} is done with maximum length of ${MAX_LEN_R1} (R1) and ${MAX_LEN_R2} (R2)."
echo "File with adapters to trim: ${ADAPTERS_FILE}"

fastp --detect_adapter_for_pe \
    --length_required=${MIN_LEN} --thread=$(( (${SLURM_ARRAY_TASK_ID} % 16) + 1 )) \
    --in1=${f1} --in2=${f2} \
    --out1=${TRIM_DIR}/${outf1} --out2=${TRIM_DIR}/${outf2} \
    --html=${TRIM_DIR}/qc/${sampleName}_fastp.html \
    --json=${TRIM_DIR}/qc/${sampleName}_fastp.json \
    --adapter_fasta ${ADAPTERS_FILE} \
    --max_len1 ${MAX_LEN_R1} --max_len2 ${MAX_LEN_R2} \
    -g -c -x
    

if [[ ! -f ${TRIM_DIR}/${outf1} ]] && [[ ! -f ${TRIM_DIR}/${outf2} ]]
then 
   { echo "TRIM on ${sampleName} was not completed. Please check variables and inputs." ; exit 1; }
else
   echo "Trimming of ${sampleName} done."
fi

## ============ ##
##    QC        ##
## ============ ##

## QC is performed in trimmed reads
fastqc ${TRIM_DIR}/${outf1} ${TRIM_DIR}/${outf2} -t 8 -o ${TRIM_DIR}/qc/

