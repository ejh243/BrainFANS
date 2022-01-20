## Written by Eilis
## Uses HMMRATAC 
## calls peaks per sample for QC purposes
## parameter choices guided by this post: https://github.com/taoliu/MACS/issues/331

cd ${ALIGNEDDIR}
mkdir -p ${PEAKDIR}/HMMRATAC

## create file paths from sample name
## requires sorted bam file with index from paired end data
sampleName=$1
BAM=${sampleName}_depDup_q30.bam # nb should be identical to ${sampleName}.filt.nmsrt.nodup.bam used for MACS
INDEX=${sampleName}_depDup_q30.bam.bai


java -jar ${HMMRATAC} -b ${BAM} -i ${INDEX} -g ${GENOMESIZE}  -o ${PEAKDIR}/HMMRATAC/${sampleName}

## exclude peaks aligned to blacklist regions
bedtools intersect -v -a ${PEAKDIR}/HMMRATAC/${sampleName}_peaks.gappedPeak -b ${BLACKLIST} \
  | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
  | grep -P 'chr[\dXY]+[ \t]'  | gzip -nc > ${PEAKDIR}/HMMRATAC/${sampleName}.gappedPeak.filt.gz
