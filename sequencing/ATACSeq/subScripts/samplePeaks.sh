#!/bin/bash
## ===================================================================================================================##
##                                     PIPPA STEP 3.1: Sample Peak Calling                                            ##
## ===================================================================================================================##
## EXECUTION: sbatch ./subScripts/samplePeaks.sh <sampleName>                                                         ||
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
## - MACS3,bedtools installed in a conda environment,                                                                 ||
## - Variables in config file: ALIGNED_DIR, PEAK_DIR_SAMPLES, BLACKLIST, macs_sample_thres                            ||
## - BLACKLIST: list of blacklist regions to exclude peaks called in these                                            ||
##                                                                                                                    ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

sampleName=$1

INPUT_BAM="${ALIGNED_DIR}/${sampleName}.filt.nodup.bam"
OUTPUT_PEAKS="${PEAK_DIR_SAMPLES}/${sampleName}.narrowPeak.filt"

if [[ ! -s "${INPUT_BAM}" ]]
then 
  { echo "Aligned bam file for ${sampleName} not found. Peaks can't be called." ; exit 1; }
else
  echo "Aligned bam file for ${sampleName} found." 
fi
  
## ============================= ##
## MACS3 BAMPE mode peak calling ##
## ============================= ##

echo "Starting peak calling using MACS3 PE for sample ${sampleName} at:"
date -u
echo "Cutoff for narrow peak calling is ${macs_sample_thres}"

cd ${TMPDIR}
macs3 callpeak -t "${INPUT_BAM}" --outdir "${PEAK_DIR_SAMPLES}" \
  -n "${sampleName}" -f BAMPE -g 2.9e9 \
  -q "${macs_sample_thres}" --keep-dup all --nomodel

## exclude peaks aligned to blacklist regions, sorted by chromosome
bedtools intersect -v -a "${PEAK_DIR_SAMPLES}/${sampleName}_peaks.narrowPeak" -b "${BLACKLIST}" \
| sort -k1 \
| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' > "${OUTPUT_PEAKS}"

rm "${PEAK_DIR_SAMPLES}/${sampleName}_summits.bed"
rm "${PEAK_DIR_SAMPLES}/${sampleName}_peaks.narrowPeak"
rm "${PEAK_DIR_SAMPLES}/${sampleName}_peaks.xls"


if [[ ! -f "${OUTPUT_PEAKS}" ]]
then
  { echo "Peak calling on ${sampleName} could not be completed. Please make sure STEP 1.3 ALIGN was properly run." ; exit 1; }
else
  echo "Peak calling on ${sampleName} done."
  echo Job ended on:
  date -u
fi

