#!/bin/bash
## ===================================================================================================================##
##                          ATAC-seq pipeline STEP 5.1: Subset sex chromosomes                                        ##
## ===================================================================================================================##
## EXECUTION: sh ./subScripts/subsetSexChrs.sh <sampleName>                                                           ||
## - execute from pipeline's main directory                                                                           ||
##                                                                                                                    ||
## DESCRIPTION: This script isolates reads from sex chromosomes to prepare for peak calling.                          ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <sampleName> sample name to subset                                                                           ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           ||
## ALIGNED_DIR/sexChr/<sampleName>.chrX.tn5.tagAlign.gz, ALIGNED_DIR/sexChr/<sampleName>.chrY.tn5.tagAlign.gz         ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - aligned filtered, no duplicated tagAlign reads: <sampleName>.tn5.tagAlign.gz                                     ||
## - Software: bedtools in conda environment                                                                          ||
## - Variables in config file: ALIGNED_DIR                                                                            ||
##                                                                                                                    ||
## ===================================================================================================================##


sampleName=$1
echo "Split reads for sex chromosomes for sample: " ${sampleName}
cd ${ALIGNED_DIR}/

tagalign=${sampleName}.tn5.tagAlign.gz

if [[ ! -f ${ALIGNED_DIR}/${tagalign} ]]
then
   { echo "Sample in tagAlign format not found. Please refer to STEP 3.1 to obtain this." ; exit 1 ;}
else
  echo "Sample in tagAlign format found: $tagalign. Reads in sex chromosomes will be subset."
fi

## ===================== ##
## TagAlign sample split ##
## ===================== ##

## extract X chr from tag align file
zcat $tagalign | grep "chrX"  | gzip -c > sexChr/${tagalign/tn5/chrX.tn5}

## extract Y chr from tag align file
zcat $tagalign | grep "chrY"  | gzip -c > sexChr/${tagalign/tn5/chrY.tn5}


if [[ ! -f ${ALIGNED_DIR}/sexChr/${sampleName}.chrY.tn5.tagAlign.gz ]] && [[ ! -f ${ALIGNED_DIR}/sexChr/${sampleName}.chrX.tn5.tagAlign.gz ]]
then
   { echo "Sample could not be subset. Please refer to STEP 3.1 to check tagAlign file of ${sampleName} exists." ; exit 1 ;}
else
  echo "Split reads for sex chromosomes for sample: ${sampleName} done"
fi


