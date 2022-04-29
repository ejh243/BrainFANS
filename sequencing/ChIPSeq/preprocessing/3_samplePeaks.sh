## Performs sample level peak calling using MACS with paired end BAM files
## parameter choices guided by this post: https://github.com/crazyhottommy/ChIP-seq-analysis/blob/master/part1.3_MACS2_peak_calling_details.md

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

## OUTPUT
# ${PEAKDIR}/*.broadPeak (and other macs output)
# ${PEAKDIR}/${sampleName}.broadPeak.filt
# ${PEAKDIR}/*.narrowPeak (and other macs output)
# ${PEAKDIR}/${sampleName}.narrowPeak.filt

## broad and narrow marks taken from https://www.encodeproject.org/chip-seq/histone/
broad=(H3F3A H3K27me3 H3K36me3 H3K4me1 H3K79me2 H3K79me3 H3K9me1 H3K9me2 H4K20me1)
narrow=(H2AFZ H3ac H3K27ac H3K4me2 H3K4me3 H3K9ac)

cd ${ALIGNEDDIR}
mkdir -p ${PEAKDIR}
mark=$1
sampleName=$2

if [[ ! -z "$3" ]]
then 
  control=$3
  c=${control}.filt.nodup.bam
  echo "Control file:" ${c}
fi

f=${sampleName}.filt.nodup.bam
echo "Processing:" ${f} 


if [[ ${broad[*]} =~  ${mark} ]]
then
  macs2 callpeak -t ${f} -c ${c} --outdir ${PEAKDIR} -n ${sampleName} -g 2.9e9 -q 5e-2 --broad --nomodel --extsize 147 --keep-dup all 2> ${PEAKDIR}/${sampleName}.macs2.log


  ## exclude peaks aligned to blacklist regions
  bedtools intersect -v -a ${PEAKDIR}/${sampleName}_peaks.broadPeak -b ${BLACKLIST} \
    | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
    | grep -P 'chr[\dXY]+[ \t]' > ${PEAKDIR}/${sampleName}.broadPeak.filt

  rm ${PEAKDIR}/${sampleName}_peaks.broadPeak
  rm ${PEAKDIR}/${sampleName}_peaks.xls
  rm ${PEAKDIR}/${sampleName}_peaks.gappedPeak   
else if [[ ${narrow[*]} =~  ${mark} ]]
then
  # calculate narrowPeak
  macs2 callpeak -t ${f} -c ${c} --outdir ${PEAKDIR} -n ${sampleName} -g 2.9e9 -q 1e-2 --nomodel --extsize 147 --keep-dup all 2> ${PEAKDIR}/${sampleName}.macs2.log

  ## exclude peaks aligned to blacklist regions
  bedtools intersect -v -a ${PEAKDIR}/${sampleName}_peaks.narrowPeak -b ${BLACKLIST} \
    | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
    | grep -P 'chr[\dXY]+[ \t]' > ${PEAKDIR}/${sampleName}.narrowPeak.filt

  rm ${PEAKDIR}/${sampleName}_peaks.narrowPeak
  rm ${PEAKDIR}/${sampleName}_peaks.xls
  rm ${PEAKDIR}/${sampleName}_peaks.gappedPeak 
else
  echo 'Target not found in either broad or narrow peak mark'
fi