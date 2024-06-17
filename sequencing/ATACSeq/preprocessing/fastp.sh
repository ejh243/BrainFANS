#!/bin/bash
## ===================================================================================================================##
##                               ATAC-seq pipeline STEP 1.2: Pre-analysis -- trimming                                 ##
## ===================================================================================================================##
## EXECUTION: sbatch --array= ./sequencing/ATACSeq/preprocessing/fastp.sh <sampleName> <R1> <R2>                      ||
## - execute from scripts directory                                                                                   ||
##                                                                                                                    ||
## DESCRIPTION: This script performs trimming of raw reads and following quality control of the new trimmed reads.    ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - Variables in config file: RAWDATADIR, TRIMDIR                                                                    ||
## - Software: fastp, fastqc                                                                                          ||
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

echo
echo "Starting trimming on" ${sampleName} "at: "
date -u

## Get read 1 and read 2 for specified sample name
f1=$(basename $2)
f2=$(basename $3)

cd ${RAWDATADIR}

## create output filenames
outf1=${f1/.f/_trimmed.f}
outf2=${f2/.f/_trimmed.f}

mkdir -p ${TRIMDIR}/qc/

##Specific adapters for R1 and R2 need to be specified for correct trimming
adapterR1=CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
adapterR2=CTGTCTCTTATACACATCTGACGCTGCCGACGA 

## ============ ##
##    TRIM      ##
## ============ ##

fastp --detect_adapter_for_pe --adapter_sequence=$adapterR1 --adapter_sequence_r2=$adapterR2 --length_required=27 --thread=8 --in1=${f1} --in2=${f2} --out1=${TRIMDIR}/${outf1} --out2=${TRIMDIR}/${outf2} --html=${TRIMDIR}/qc/${sampleName}_fastp.html --json=${TRIMDIR}/qc/${sampleName}_fastp.json

## ============ ##
##    QC        ##
## ============ ##

fastqc ${TRIMDIR}/${outf1} ${TRIMDIR}/${outf2} -t 8 -o ${TRIMDIR}/qc/