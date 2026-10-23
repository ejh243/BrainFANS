#!/bin/bash
## ===========================================================================================================##
##                                PIPPA STEP 5.1: Subset sex chromosomes                                      ##
## ===========================================================================================================##
## EXECUTION: sbatch ./subScripts/subsetSexChrs.sh <sampleName>                                               ||
##                                                                                                            ||
## DESCRIPTION: This script isolates reads from sex chromosomes to prepare for peak calling.                  ||
##                                                                                                            ||
## INPUTS:                                                                                                    || 
## $1 -> <sampleName> sample name to subset                                                                   ||
##                                                                                                            ||
## OUTPUTS:                                                                                                   ||
## ALIGNED_DIR/sexChr/<sampleName>.chrX.bam, ALIGNED_DIR/sexChr/<sampleName>.chrY.bam                         ||
##                                                                                                            ||
## REQUIRES:                                                                                                  ||
## - aligned filtered, no duplicated tagAlign reads: <sampleName>.filt.nodup.bam                              ||
## - Samtools in conda environment                                                                            ||
## - Variables in config file: ALIGNED_DIR                                                                    ||
##                                                                                                            ||
## ===========================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

sampleName=$1
echo "Split reads for sex chromosomes for sample: " ${sampleName}
cd ${ALIGNED_DIR}/

## ===================== ##
##  BAMPE sample split   ##
## ===================== ##

bamfile=${sampleName}.filt.nodup.bam

## extract X chr reads
samtools view -b $bamfile chrX > "${ALIGNED_DIR}/sexChr/${sampleName}.chrX.bam"

## extract Y chr reads
samtools view -b $bamfile chrY > "${ALIGNED_DIR}/sexChr/${sampleName}.chrY.bam"

echo "Split reads for sex chromosomes for sample: ${sampleName} done"

if [[ ! -f "${ALIGNED_DIR}/sexChr/${sampleName}.chrX.bam" ]] || [[ ! -f "${ALIGNED_DIR}/sexChr/${sampleName}.chrY.bam" ]]
then 
  { echo "Bam files with reads in sex chr created not created for sample ${sampleName}." ; exit 1; }
else
  echo "Bam files with reads in sex chr created for sample ${sampleName}."
  echo Job finished on:
  date -u 
fi 