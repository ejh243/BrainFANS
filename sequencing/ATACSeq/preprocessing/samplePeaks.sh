#!/bin/bash
## ===================================================================================================================##
##                             ATAC-seq pipeline STEP 3.2: Sample Peak Calling                                        ##
## ===================================================================================================================##
## EXECUTION: sh ./sequencing/ATACSeq/preprocessing/samplePeaks.sh <sampleName>                                       ||
## - execute from scripts directory                                                                                   ||
##                                                                                                                    ||
## DESCRIPTION: This script calls peaks in an input sample using MACS3 in paired-end mode. Parameter choices guided by||
## https://github.com/taoliu/MACS/issues/331                                                                          ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <sampleName> sample name to call peaks on                                                                    ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           ||
## BAMPE/.broadPeak.filt                                                                                              ||
## - output will be peaks called in autosomal chromosomes and not aligned to blacklist regions                        ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - aligned filtered, no duplicated bam file for MACS3 PE                                                            ||
## - MACS3 installed in a pip environment, bedtools                                                                   ||
## - ALIGNEDDIR, PEAKDIR, BLACKLIST                                                                                   ||
##                                                                                                                    ||
## ===================================================================================================================##

sampleName=$1

## ============================= ##
## MACS3 BAMPE mode peak calling ##
## ============================= ##


sampleBAM=${sampleName}.filt.nodup.bam
echo "Calling peaks using MACS3 BAMPE for"" ${sampleName}"
macs3 callpeak -t ${ALIGNEDDIR}/${sampleBAM} --outdir ${PEAKDIR}/MACS/BAMPE -n ${sampleName} -f BAMPE -g 2.9e9 -q 5e-2 --keep-dup all --broad --broad-cutoff 5e-2

## exclude peaks aligned to blacklist regions, exclude peaks called in chr X and Y, sorted by chromosome
bedtools intersect -v -a ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.broadPeak -b ${BLACKLIST} \
| awk '!/^(chrY|chrX)/' \
| sort -k1 \
| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' > ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt

rm ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.broadPeak
rm ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.xls
rm ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.gappedPeak  

echo "Peaks for ${sampleName} using MACS3 PE called"

