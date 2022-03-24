## Performs sample level peak calling using MACS with 1. shifted tagAlign files 
## and 2. paired end BAM files
## parameter choices guided by this post: https://github.com/taoliu/MACS/issues/331


## EXECUTION
# sh ./ChIPSeq/preprocessing/3_samplePeaks.sh <sampleName>
# where 
# <sampleName> is the sample specific file prefix
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR, BLACKLIST

## REQUIRES the following software
# bedtools, macs2

## INPUT
# chip bam files
#

## OUTPUT
# ${PEAKDIR}/MACS/*.narrowPeak (and other macs output)
# ${PEAKDIR}/MACS/${sampleName}.narrowPeak.filt

cd ${ALIGNEDDIR}
mkdir -p ${PEAKDIR}/MACS
sampleName=$1

if [[ ! -z "$2" ]]
then 
  control=$2
  c=${control}_depDup_q30.bam
  echo "Control file:" ${c}
fi

f=${sampleName}_depDup_q30.bam
echo "Processing:" ${f} 

macs2 callpeak -t ${f} -c ${c} --outdir ${PEAKDIR}/MACS/ -n ${sampleName} -g 2.9e9  -B 2> ${PEAKDIR}/${sampleName}.macs2.log

## exclude peaks aligned to blacklist regions
bedtools intersect -v -a ${PEAKDIR}/MACS/${sampleName}_peaks.narrowPeak -b ${BLACKLIST} \
  | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
  | grep -P 'chr[\dXY]+[ \t]' > ${PEAKDIR}/MACS/${sampleName}.narrowPeak.filt
