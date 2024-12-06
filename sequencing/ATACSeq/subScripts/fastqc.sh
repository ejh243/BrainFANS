#!/bin/bash
## ===================================================================================================================##
##                    ATAC-seq pipeline STEP 1.1: Pre-analysis -- quality-control raw data                            ##
## ===================================================================================================================##
## EXECUTION: sbatch  ./subScripts/fastqc.sh <sampleName> <R1> <R2>                                                   ||
## - execute from pipeline's main directory                                                                           ||
##                                                                                                                    ||
## DESCRIPTION: This script calculates sequencing quality control qc metrics with fastqc for paired files for a       ||
## single sample.                                                                                                     ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - Variables in config file: RAWDATADIR, FASTQCDIR                                                                  ||
## - Software: fastqc (in conda environment)                                                                          ||
## - Paired-end files for sample: Read 1 and Read 2 in the same directory (RAWDATADIR)                                ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <sampleName> Name of sample specified in command line.                                                       ||
## $2 -> <R1> Read 1 name                                                                                             ||
## $3 -> <R2> Read 2 name                                                                                             ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           || 
##  sample_R1.fastq.gz, sample_R2.fastq.gz, sample_R1_fastqc.html, sample_R2_fastqc.html                              ||
##                                                                                                                    ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

sampleName=$1

f1=$(basename $2)
f2=$(basename $3)

## ========== ##
##    QC      ##
## ========== ##

fastqc ${RAWDATADIR}/${f1}  ${RAWDATADIR}/${f2} -t 8 -o ${FASTQCDIR}

ending1=".fastq.gz"
ending2=".fq.gz"
out1=$(echo "$f1" | sed -e "s/$ending1//g" -e "s/$ending2//g")
out2=$(echo "$f2" | sed -e "s/$ending1//g" -e "s/$ending2//g")

echo ${out1}

if [[ ! -f "${FASTQCDIR}/${out1}_fastqc.html" ]] && [[ ! -f "${FASTQCDIR}/${out2}_fastqc.html" ]]
then 
  { echo "FASTQC on ${sampleID} was not completed. Please check variables and inputs." ; exit 1; }
fi 

echo Job finished on:
date -u 