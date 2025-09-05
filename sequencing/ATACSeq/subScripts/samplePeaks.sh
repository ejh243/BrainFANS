#!/bin/bash
## ===================================================================================================================##
##                             ATAC-seq pipeline STEP 3.2: Sample Peak Calling                                        ##
## ===================================================================================================================##
## EXECUTION: sh ./subScripts/samplePeaks.sh <sampleName>                                                             ||
## - execute from pipeline's main  directory                                                                          ||
##                                                                                                                    ||
## DESCRIPTION: This script calls peaks in an input sample using MACS3 in paired-end mode. Parameter choices guided by||
## https://github.com/taoliu/MACS/issues/331                                                                          ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <sampleName> sample name to call peaks on                                                                    ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           ||
## - PEAK_DIR/BAMPE/<sampleName>.narrowPeak.filt                                                                      ||
## - output will be peaks called in all chromosomes and not aligned to blacklist regions                              ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - aligned filtered, no duplicated bam file for MACS3 PE: <sampleName>.filt.nodup.bam                               ||
## - MACS3 installed in a conda environment, bedtools                                                                 ||
## - Variables in config file: ALIGNED_DIR, PEAK_DIR_SAMPLES, BLACKLIST,                                              ||
## - BLACKLIST: list of blacklist regions to exclude peaks called in these                                            ||
##                                                                                                                    ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

sampleName=$1

if [[ ! -s ${ALIGNED_DIR}/${sampleName}.filt.nodup.bam ]]
then 
  { echo "Aligned bam file for ${sampleName} not found. Peaks can't be called." ; exit 1; }
else
  echo "Aligned bam file for ${sampleName} found." 
fi
  
## ============================= ##
## MACS3 BAMPE mode peak calling ##
## ============================= ##

sampleBAM=${sampleName}.filt.nodup.bam
echo "Starting peak calling using MACS3 PE for ${sampleName} at:"
date -u
echo "Cutoff for broad peak calling is ${MACS_SAMPLE}"
macs3 callpeak -t ${ALIGNED_DIR}/${sampleBAM} --outdir ${PEAK_DIR_SAMPLES} -n ${sampleName} -f BAMPE -g 2.9e9 -q ${MACS_SAMPLE} --keep-dup all --nomodel

## exclude peaks aligned to blacklist regions, sorted by chromosome
bedtools intersect -v -a ${PEAK_DIR_SAMPLES}/${sampleName}_peaks.narrowPeak -b ${BLACKLIST} \
| sort -k1 \
| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' > ${PEAK_DIR_SAMPLES}/${sampleName}.narrowPeak.filt

rm ${PEAK_DIR_SAMPLES}/${sampleName}_peaks.narrowPeak
rm ${PEAK_DIR_SAMPLES}/${sampleName}_peaks.xls
rm ${PEAK_DIR_SAMPLES}/${sampleName}_peaks.gappedPeak  

if [[ ! -f ${PEAK_DIR_SAMPLES}/${sampleName}.narrowPeak.filt ]]
then
  { echo "Peak calling on ${sampleName} could not be completed. Please make sure STEP 1.3 ALIGN was properly run." ; exit 1; }
else
  echo "Peak calling on ${sampleName} done."
  echo Job ended on:
  date -u
fi

