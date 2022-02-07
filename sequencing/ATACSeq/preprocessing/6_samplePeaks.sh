## Performs sample level peak calling using MACS with 1. shifted tagAlign files 
## and 2. paired end BAM files
## parameter choices guided by this post: https://github.com/taoliu/MACS/issues/331


## EXECUTION
# sh ./ATACSeq/preprocessing/6_samplePeaks.sh <sampleName>
# where 
# <sampleName> is the file prefix
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR, BLACKLIST

## REQUIRES the following software
# bedtools, macs2

## INPUT
# shifted tag align file
#

## OUTPUT
# ${PEAKDIR}/MACS/ShiftedTagAlign/*.broadPeak (and other macs output)
# ${PEAKDIR}/MACS/BAMPE/*.broadPeak
# ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}.broadPeak.filt
# ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt



cd ${ALIGNEDDIR}
mkdir -p ${PEAKDIR}/MACS/ShiftedTagAlign
sampleName=$1
f=${sampleName}.tn5.tagAlign.gz

macs2 callpeak -t ${f} --outdir ${PEAKDIR}/MACS/ShiftedTagAlign -n ${sampleName} -f BED -g 2.9e9 -p 1e-3 --keep-dup all --shift 100 --extsize 200 --nomodel --broad


mkdir -p ${PEAKDIR}/MACS/BAMPE

f=${sampleName}.filt.nmsrt.nodup.bam
macs2 callpeak -t ${f} --outdir ${PEAKDIR}/MACS/BAMPE -n ${sampleName} -f BAMPE -g 2.9e9 -p 1e-3 --keep-dup all --nomodel --broad 

## exclude peaks aligned to blacklist regions
bedtools intersect -v -a ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}_peaks.broadPeak -b ${BLACKLIST} \
  | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
  | grep -P 'chr[\dXY]+[ \t]' > ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}.broadPeak.filt

rm ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}_peaks.broadPeak
  
bedtools intersect -v -a ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.broadPeak -b ${BLACKLIST} \
  | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
  | grep -P 'chr[\dXY]+[ \t]' > ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt
  
rm ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.broadPeak