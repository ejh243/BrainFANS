#!/bin/bash
## ===================================================================================================================##
##                          ATAC-seq pipeline STEP 7.2: Cell group Peak Calling                                       ##
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
## - ${PEAK_DIR}/MACS/BAMPE/group/${GROUP}.sorted.broadPeak.filt                                                      ||
## - output will be peaks called in all samples that belong to the input cell type (that passed stage 1 and 2 of QC.  ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - aligned filtered, no duplicated bam file for MACS3 PE: <sampleName>.filt.nodup.bam                               ||
## - MACS3 installed in a conda environment, bedtools                                                                 ||
## - Variables in config file: ALIGNED_DIR, PEAK_DIR, BLACKLIST                                                       ||
## - BLACKLIST: list of blacklist regions to exclude peaks called in these                                            ||
##                                                                                                                    ||
## ===================================================================================================================##

pip list
## ============ ##
##    SET UP    ##
## ============ ##

mapfile -t SAMPLES < ${META_DIR}/samplesForGroupAnalysisOrdered_${GROUP}.txt

## Collate samples for peak calling
bamFiles=()
for sample in ${SAMPLES[@]}
do
  bam_sample=${ALIGNED_DIR}/${sample}.filt.nodup.bam
  bamFiles+=( $bam_sample )
done

## ============================= ##
## MACS3 BAMPE mode peak calling ##
## ============================= ##

echo "Starting peak calling using MACS3 PE on samples that belong to cell fraction ${GROUP} at: "
date -u

cd ${TMPDIR}
macs3 callpeak -t  ${bamFiles[@]} --outdir ${PEAK_DIR}/MACS/BAMPE/group -n ${GROUP} -f BAMPE -g 2.9e9 -q 5e-2 --keep-dup all --nomodel --broad --broad-cutoff 5e-2

## exclude peaks aligned to blacklist regions and peaks called in chr X and Y
bedtools intersect -v -a ${PEAK_DIR}/MACS/BAMPE/group/${GROUP}_peaks.broadPeak -b ${BLACKLIST} \
	| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
	| awk '!/^(chrY|chrX)/' > ${PEAK_DIR}/MACS/BAMPE/group/${GROUP}.broadPeak.filt

echo "Sorting peaks by chr for merged sample peaks"
sort -k1 ${PEAK_DIR}/MACS/BAMPE/group/${GROUP}.broadPeak.filt > ${PEAK_DIR}/MACS/BAMPE/group/${GROUP}.sorted.broadPeak.filt

if [[ ! -f ${PEAK_DIR}/MACS/BAMPE/group/${GROUP}.sorted.broadPeak.filt ]]
then
  { echo "Peak calling on ${sampleName} could not be completed. Please make sure STEP 7.0 was properly run." ; exit 1; }
else
  echo "Finished calling peaks using MACS3 PE on samples in group ${GROUP}"
  echo Job ended on:
  date -u
fi
