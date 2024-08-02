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
## - PEAK_DIR/BAMPE/<sampleName>.broadPeak.filt                                                                       ||
## - output will be peaks called in autosomal chromosomes and not aligned to blacklist regions                        ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - aligned filtered, no duplicated bam file for MACS3 PE: <sampleName>.filt.nodup.bam                               ||
## - MACS3 installed in a pip environment, bedtools                                                                   ||
## - Variables in config file: ALIGNED_DIR, PEAK_DIR, BLACKLIST                                                       ||
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

echo "Calling peaks using MACS3 PE for"" ${sampleName}"
macs3 callpeak -t ${ALIGNED_DIR}/${sampleBAM} --outdir ${PEAK_DIR}/BAMPE -n ${sampleName} -f BAMPE -g 2.9e9 -q 5e-2 --keep-dup all --broad --broad-cutoff 5e-2

## exclude peaks aligned to blacklist regions, exclude peaks called in chr X and Y, sorted by chromosome
bedtools intersect -v -a ${PEAK_DIR}/BAMPE/${sampleName}_peaks.broadPeak -b ${BLACKLIST} \
| awk '!/^(chrY|chrX)/' \
| sort -k1 \
| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' > ${PEAK_DIR}/BAMPE/${sampleName}.broadPeak.filt

rm ${PEAK_DIR}/BAMPE/${sampleName}_peaks.broadPeak
rm ${PEAK_DIR}/BAMPE/${sampleName}_peaks.xls
rm ${PEAK_DIR}/BAMPE/${sampleName}_peaks.gappedPeak  

if [[ ! -f ${PEAK_DIR}/BAMPE/${sampleName}.broadPeak.filt ]]
then
  { echo "Peak calling on ${sampleName} could not be completed. Please make sure STEP 1.3 ALIGN was properly run." ; exit 1; }
else
  echo "Peak calling on ${sampleName} done."
  echo Job ended on:
  date -u
fi

