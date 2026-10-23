#!/bin/bash
## ===================================================================================================================##
##                                  PIPPA STEP 8.2: Cell group Peak Calling                                           ##
## ===================================================================================================================##
## EXECUTION: sh ./groupPeaks.sh <cell-group>                                                                         ||
##                                                                                                                    ||
## DESCRIPTION: This script calls peaks at cell group level, on all samples that belong to the same cell type using   ||
## MACS3 in paired-end mode.                                                                                          ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <cell-group> Cell group of samples to call peaks on                                                          ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           ||
## - ${PEAK_DIR_GROUPS}/${GROUP}/${PEAKS_NAME_GROUP}_peaks.narrowPeak : all peaks                                     ||
## - ${PEAK_DIR_GROUPS}/${GROUP}/${PEAKS_NAME_GROUP}.narrowPeak.filt : peaks in blacklisted regions filtered          ||
## - ${PEAK_DIR_GROUPS}/${GROUP}/${PEAKS_NAME_GROUP}.narrowPeak.noSC.filt : peaks in blacklisted regions and sex chrs ||
##   filtered                                                                                                         ||
## - ${PEAK_DIR_GROUPS}/${GROUP}/${PEAKS_NAME_GROUP}.narrowPeak.noY.filt : peaks in blacklisted regions and Y chr filt||
## - ${PEAK_DIR_GROUPS}/${GROUP}/${PEAKS_NAME_GROUP}.sorted.noSC.narrowPeak.filt : peaks in blacklisted regions and   ||
##   sex chromosomes filtered. Sorted by chromosome number                                                            ||
## - output will be peaks called in all samples that belong to the input cell type (that passed stage 1 and 2 of QC). ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - if available, merged bam file with samples that belong to cell-type input: ${GROUP}.sorted.merged.bam            ||
## - otherwise, aligned bam files. List of samples in samplesForGroupAnalysis_${GROUP}.txt                            ||
## - MACS3 and bedtools installed in a conda environment                                                              ||
## - Variables in config file: ALIGNED_DIR, PEAK_DIR_GROUPS, BLACKLIST,macs_group_thres                               ||
## - BLACKLIST: list of blacklist regions to exclude peaks called in these, specified in config.txt file              ||
##                                                                                                                    ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

PEAKS_NAME_GROUP="${GROUP}_${macs_group_thres}"
OUTPUT_DIR="${PEAK_DIR_GROUPS}/${GROUP}/"


echo "Starting peak calling using MACS3 PE on samples that belong to cell fraction ${GROUP} at: "
date -u

cd ${TMPDIR}

echo "Cutoff for narrow peak calling is $macs_group_thres"
echo "Peaks file name is ${PEAKS_NAME_GROUP}"

## ============================= ##
## MACS3 BAMPE mode peak calling ##
## ============================= ##

mergedBam="${ALIGNED_DIR}/merged/passS12/${GROUP}.sorted.merged.bam"

if [[ -f ${mergedBam} ]]
then
  echo "Merged bam file for $GROUP found. Using this for peak calling"
  
  macs3 callpeak -t ${mergedBam} \
  --outdir ${OUTPUT_DIR} -n ${PEAKS_NAME_GROUP} \
  -f BAMPE -g 2.9e9 -q $macs_group_thres \
  --keep-dup all --nomodel 
else
  echo "Merged bam file for $GROUP not found. Using aligned samples bam files instead."
  mapfile -t SAMPLES < "${META_DIR}/samplesForGroupAnalysis_${GROUP}.txt"

  ## Collate samples for peak calling
  bamFiles=()
  for sample in ${SAMPLES[@]}
  do
   bam_sample="${ALIGNED_DIR}/${sample}.filt.nodup.bam"
   bamFiles+=( $bam_sample )
  done
  
  macs3 callpeak -t ${bamFiles[@]} \
  --outdir ${OUTPUT_DIR} -n ${PEAKS_NAME_GROUP} \
  -f BAMPE -g 2.9e9 -q $macs_group_thres \
  --keep-dup all --nomodel 
fi

## exclude peaks aligned to blacklist regions
bedtools intersect -v -a "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}_peaks.narrowPeak" -b ${BLACKLIST} \
	| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' > "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}.narrowPeak.filt"
 
echo "Sorting peaks by chr"
sort -k1 "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}.narrowPeak.filt" > "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}.sorted.narrowPeak.filt"

rm "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}_peaks.narrowPeak"
rm "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}_peaks.xls" 

if [[ ! -f "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}.sorted.narrowPeak.filt" ]]
then
  { echo "Peak calling on ${GROUP} could not be completed. Please make sure STEP 8.0 or 8.1 were properly run." ; exit 1; }
else
  echo "Finished calling peaks using MACS3 PE on samples in group ${GROUP}"
  echo Job ended on:
  date -u
fi