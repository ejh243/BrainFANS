## ===================================================================================================================##
##                        ATAC-seq pipeline STEP 3.3: Sample Peak Calling using HMMRATAC                              ##
## ===================================================================================================================##
## EXECUTION: sh ./sequencing/ATACSeq/preprocessing/sampleHMMRatac.sh <sampleName>                                    ||
## - execute from scripts directory                                                                                   ||
##                                                                                                                    ||
## DESCRIPTION: This script calls peaks in an input sample using MACS3 HMMRATAC                                       ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <sampleName> sample name to call peaks on                                                                    ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           ||
## HMMRATAC/.sorted.filteredPeaks.gappedPeak                                                                          ||
## - output will be peaks called in autosomal chromosomes and not aligned to blacklist regions                        ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - aligned filtered, no duplicated bam reads                                                                        ||
## - MACS3 installed in a pip environment, bedtools                                                                   ||
## - ALIGNEDDIR, PEAKDIR, BLACKLIST                                                                                   ||
##                                                                                                                    ||
## ===================================================================================================================##

sampleName=$1

## ================================ ##
## MACS3 HMMRATAC mode peak calling ##
## ================================ ##

echo "Calling peaks for sample ${sampleName} using MACS3 HMMRATAC"

macs3 hmmratac -i ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam -f BAMPE --outdir ${PEAKDIR}/HMMRATAC/ -n ${sampleName} --keep-duplicates

## exclude peaks aligned to blacklist regions, exclude peaks called in chr X and Y, sorted by chromosome
bedtools intersect -v -a ${PEAKDIR}/HMMRATAC/${sampleName}_accessible_regions.gappedPeak -b ${BLACKLIST} \
  | awk '!/^(chrY|chrX)/' \
  | sort -k1 > ${PEAKDIR}/HMMRATAC/${sampleName}.sorted.filteredPeaks.gappedPeak 

echo "Peaks for ${sampleName} using MACS3 HMMRATAC called"
