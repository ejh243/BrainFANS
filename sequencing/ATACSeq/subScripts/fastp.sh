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
## - Variables in config file: RAWDATADIR, TRIM_DIR                                                                   ||
## - Software: fastp, fastqc (in conda environment)                                                                   ||
## - Read 1 and Read 2 of the same sample in the same directory (RAWDATADIR)                                          ||
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

##Specific adapters for R1 and R2 need to be specified for correct trimming
adapterR1=${ADAP_R1}
adapterR2=${ADAP_R2}

## ============ ##
##    TRIM      ##
## ============ ##

fastp --detect_adapter_for_pe --adapter_sequence=$adapterR1 --adapter_sequence_r2=$adapterR2 --length_required=27 --thread=$(( (${SLURM_ARRAY_TASK_ID} % 16) + 1 )) --in1=${f1} --in2=${f2} --out1=${TRIM_DIR}/${outf1} --out2=${TRIM_DIR}/${outf2} --html=${TRIM_DIR}/qc/${sampleName}_fastp.html --json=${TRIM_DIR}/qc/${sampleName}_fastp.json

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