#!/bin/bash
## ===================================================================================================================##
##                          ATAC-seq pipeline STEP 7.1: Cell group Peak Calling                                       ##
## ===================================================================================================================##
## EXECUTION: sh ./subScripts/groupPeaks.sh <cell-group>                                                              ||
## - execute from pipeline's main  directory                                                                          ||
##                                                                                                                    ||
## DESCRIPTION: This script calls peaks at cell group level, on all samples that belong to the same cell type using   ||
## MACS3 in paired-end mode.                                                                                          ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <cell-group> Cell group of samples to call peaks on                                                          ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           ||
## - ${PEAK_DIR_GROUPS}/${GROUP}.sorted.narrowPeak.filt                                                               ||
## - output will be peaks called in all samples that belong to the input cell type (that passed stage 1 and 2 of QC.  ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - aligned filtered, no duplicated bam file for MACS3 PE: <sampleName>.filt.nodup.bam                               ||
## - MACS3 and bedtools installed in a conda environment                                                              ||
## - Variables in config file: ALIGNED_DIR, PEAK_DIR_GROUPS, BLACKLIST,MACS_GROUP                                     ||
## - BLACKLIST: list of blacklist regions to exclude peaks called in these                                            ||
##                                                                                                                    ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

mapfile -t SAMPLES < ${META_DIR}/samplesForGroupAnalysisOrdered_${GROUP}.txt

PEAKS_NAME_GROUP="${GROUP}_${MACS_GROUP}"

## Collate samples for peak calling
bamFiles=()
for sample in ${SAMPLES[@]}
do
  bam_sample=${ALIGNED_DIR}/${sample}.filt.nodup.bam
  bamFiles+=( $bam_sample )
done

echo ${#SAMPLES[@]}

## ============================= ##
## MACS3 BAMPE mode peak calling ##
## ============================= ##

echo "Starting peak calling using MACS3 PE on samples that belong to cell fraction ${GROUP} at: "
date -u

cd ${TMPDIR}

echo "Cutoff for narrow peak calling is $MACS_GROUP"
echo "Peaks file name is ${PEAKS_NAME_GROUP}"
macs3 callpeak -t  ${bamFiles[@]} --outdir ${PEAKS_NAME_GROUP} -n ${PEAKS_NAME_GROUP} -f BAMPE -g 2.9e9 -q $MACS_GROUP --keep-dup all --nomodel 

## exclude peaks aligned to blacklist regions and peaks called in chr X and Y
bedtools intersect -v -a ${PEAK_DIR_GROUPS}/${PEAKS_NAME_GROUP}_peaks.narrowPeak -b ${BLACKLIST} \
	| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
	| awk '!/^(chrY|chrX)/' > ${PEAK_DIR_GROUPS}/${PEAKS_NAME_GROUP}.narrowPeak.filt
 
rm ${PEAK_DIR_GROUPS}/${PEAKS_NAME_GROUP}_peaks.narrowPeak
rm ${PEAK_DIR_GROUPS}/${PEAKS_NAME_GROUP}_peaks.xls 

echo "Sorting peaks by chr for merged sample peaks"
sort -k1 ${PEAK_DIR_GROUPS}/${PEAKS_NAME_GROUP}.narrowPeak.filt > ${PEAK_DIR_GROUPS}/${PEAKS_NAME_GROUP}.sorted.narrowPeak.filt

if [[ ! -f ${PEAK_DIR_GROUPS}/${PEAKS_NAME_GROUP}.sorted.narrowPeak.filt ]]
then
  { echo "Peak calling on ${sampleName} could not be completed. Please make sure STEP 7.0 was properly run." ; exit 1; }
else
  echo "Finished calling peaks using MACS3 PE on samples in group ${GROUP}"
  echo Job ended on:
  date -u
fi
