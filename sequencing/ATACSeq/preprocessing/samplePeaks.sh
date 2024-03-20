## ===================================================================================================================##
##                             ATAC-seq pipeline STEP 3.2: Sample Peak Calling                                        ##
## ===================================================================================================================##
## EXECUTION: sh ./sequencing/ATACSeq/preprocessing/samplePeaks.sh <sampleName>                                       ||
## - execute from scripts directory                                                                                   ||
##                                                                                                                    ||
## DESCRIPTION: This script calls peaks in an input sample using MACS3 in two different modes: single-end and         ||
## paired-end reads. Parameter choices guided by https://github.com/taoliu/MACS/issues/331                            ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <sampleName> sample name to call peaks on                                                                    ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           ||
## ShiftedTagAlign/.broadPeak.filt, BAMPE/.broadPeak.filt                                                             ||
## - output will be peaks called in autosomal chromosomes and not aligned to blacklist regions                        ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - aligned filtered, no duplicated bam file for MACS3 PE                                                            ||
## - shifted TagAlign files for MACS3 TA                                                                              ||
## - MACS3 installed in a pip environment, bedtools                                                                   ||
## - ALIGNEDDIR, PEAKDIR, BLACKLIST                                                                                   ||
##                                                                                                                    ||
## ===================================================================================================================##

cd ${ALIGNEDDIR}
sampleName=$1

## ========================== ##
## MACS3 TA mode peak calling ##
## ========================== ##

f_TA=${sampleName}.tn5.tagAlign.gz
echo "Calling peaks using MACS3 TA for"" ${sampleName}"
macs3 callpeak -t ${f_TA} --outdir ${PEAKDIR}/MACS/ShiftedTagAlign -n ${sampleName} -f BED -g 2.9e9 -q 5e-2 --keep-dup all --shift 100 --extsize 200 --nomodel --broad --broad-cutoff 5e-2

## exclude peaks aligned to blacklist regions, exclude peaks called in chr X and Y, and sort by chromosome
bedtools intersect -v -a ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}_peaks.broadPeak -b ${BLACKLIST} \
| awk '!/^(chrY|chrX)/' \
| sort -k1 \
| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' > ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}.broadPeak.filt
 
rm ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}_peaks.broadPeak
rm ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}_peaks.xls
rm ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}_peaks.gappedPeak  

echo "Peaks for ${sampleName} using MACS3 TA called"

## ============================= ##
## MACS3 BAMPE mode peak calling ##
## ============================= ##


f_PE=${sampleName}.filt.nodup.bam
echo "Calling peaks using MACS3 BAMPE for"" ${sampleName}"
macs3 callpeak -t ${f_PE} --outdir ${PEAKDIR}/MACS/BAMPE -n ${sampleName} -f BAMPE -g 2.9e9 -q 5e-2 --keep-dup all --broad --broad-cutoff 5e-2

## exclude peaks aligned to blacklist regions, exclude peaks called in chr X and Y, and sort by chromosome
bedtools intersect -v -a ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.broadPeak -b ${BLACKLIST} \
| awk '!/^(chrY|chrX)/' \
| sort -k1 \
| awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' > ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt

rm ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.broadPeak
rm ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.xls
rm ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.gappedPeak  

echo "Peaks for ${sampleName} using MACS3 PE called"

