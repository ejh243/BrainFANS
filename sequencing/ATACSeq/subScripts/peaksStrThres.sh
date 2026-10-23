#!/bin/bash
## ===================================================================================================================##
##                                PIPPA STEP 9.5: Cell group Peak Calling                                             ##
## ===================================================================================================================##
## EXECUTION: sh ./subScripts/peaskStrThres.sh <cell-group> <thres>                                                   ||
##                                                                                                                    ||
## DESCRIPTION: This script calls peaks at cell group level, on all samples that belong to the same cell type using   ||
## MACS3 in paired-end mode.                                                                                          ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <cell-group> Cell group of samples to call peaks on                                                          ||
## $2 -> <thres> Stringent threshold to use for peak calling                                                          ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           ||
## - ${CURATED_PEAKS}/${GROUP}/${PEAKS_NAME_GROUP}_peaks.narrowPeak : all peaks                                       ||
## - ${CURATED_PEAKS}/${GROUP}/${PEAKS_NAME_GROUP}.narrowPeak.filt : peaks in blacklisted regions filtered            ||
## - ${CURATED_PEAKS}/${GROUP}/${PEAKS_NAME_GROUP}.narrowPeak.noSC.filt : peaks in blacklisted regions and sex chrs   ||
##   filtered                                                                                                         ||
## - ${CURATED_PEAKS}/${GROUP}/${PEAKS_NAME_GROUP}.narrowPeak.noY.filt : peaks in blacklisted regions and Y chr filt  ||
## - ${CURATED_PEAKS}/${GROUP}/${PEAKS_NAME_GROUP}.sorted.noSC.narrowPeak.filt : peaks in blacklisted regions and sex ||
##   chromosomes filtered. Sorted by chromosome number                                                                ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - aligned filtered, no duplicated bam files for list of samples ${META_DIR}/downAnalysisSamples_${GROUP}.txt       ||
## - or merged bam file: /merged/passAllQC/${GROUP}.sorted.merged.bam                                                 ||
## - MACS3 and bedtools installed in a conda environment                                                              ||
## - Variables in config file: ALIGNED_DIR, PEAK_DIR, BLACKLIST,CURATED_PEAKS                                         ||
## - BLACKLIST: list of blacklist regions to exclude peaks called in these,specified in config.txt file               ||
##                                                                                                                    ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

PEAKS_NAME_GROUP="${GROUP}_${THRES}"
OUTPUT_DIR="${CURATED_PEAKS}/${GROUP}/"

echo "Starting peak calling using MACS3 PE on samples that belong to cell fraction ${GROUP} at: "
date -u

cd ${TMPDIR}

echo "Cutoff for narrow peak calling is ${THRES}"
echo "Output directory is: ${OUTPUT_DIR}"
echo "Peaks file name is ${PEAKS_NAME_GROUP}"

## ============================= ##
## MACS3 BAMPE mode peak calling ##
## ============================= ##

mergedBam="${ALIGNED_DIR}/merged/passAllQC/${GROUP}.sorted.merged.bam"

if [[ -f ${mergedBam} ]]
then
  echo "Merged bam file for $GROUP found. Using this for peak calling"
  
  macs3 callpeak -t ${mergedBam} --outdir ${OUTPUT_DIR} \
  -n ${PEAKS_NAME_GROUP} -f BAMPE -g 2.9e9 \
  -q ${THRES} --keep-dup all --nomodel 
else
  echo "Merged bam file for $GROUP not found. Using aligned samples bam files instead."
  mapfile -t SAMPLES < "${META_DIR}/downAnalysisSamples_${GROUP}.txt"

  ## Collate samples for peak calling
  bamFiles=()
  for sample in ${SAMPLES[@]}
  do
   bam_sample="${ALIGNED_DIR}/${sample}.filt.nodup.bam"
   bamFiles+=( $bam_sample )
  done
  
  macs3 callpeak -t ${bamFiles[@]} --outdir ${OUTPUT_DIR} \
  -n ${PEAKS_NAME_GROUP} -f BAMPE -g 2.9e9 \
  -q ${THRES} --keep-dup all --nomodel 
fi

## exclude peaks aligned to blacklist regions and peaks called in chr X and Y
bedtools intersect -v -a "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}_peaks.narrowPeak" -b ${BLACKLIST} \
	| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' > "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}.narrowPeak.filt"

## exclude peaks called in chr X and Y
awk '!/^(chrY|chrX)/' "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}.narrowPeak.filt" > "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}.narrowPeak.noSC.filt"

## exclude peaks called in chr Y
awk '!/^(chrY)/' "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}.narrowPeak.filt" > "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}.narrowPeak.noY.filt"

echo "Sorting peaks by chr"
sort -k1 "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}.narrowPeak.noSC.filt" > "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}.sorted.noSC.narrowPeak.filt"

rm "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}_peaks.xls" 
rm "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}_peaks.narrowPeak"

if [[ ! -f "${OUTPUT_DIR}/${PEAKS_NAME_GROUP}.sorted.noSC.narrowPeak.filt" ]]
then
  { echo "Peak calling on ${GROUP} could not be completed. Please make sure STEP 9.0 or 9.1 were properly run." ; exit 1; }
else
  echo "Finished calling peaks using MACS3 PE on samples in group ${GROUP}"
  echo Job ended on:
  date -u
fi