## Written by Eilis
## Uses MACs in two methods
## 1. with shifted tagAlign files
## 2. using BAM files and paired end reads
## calls peaks per sample for QC purposes
## parameter choices guided by this post: https://github.com/taoliu/MACS/issues/331

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
  | grep -P 'chr[\dXY]+[ \t]'  | gzip -nc > ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}.broadPeak.filt.gz
  
bedtools intersect -v -a ${PEAKDIR}/MACS/BAMPE/${sampleName}_peaks.broadPeak -b ${BLACKLIST} \
  | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
  | grep -P 'chr[\dXY]+[ \t]'  | gzip -nc > ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt.gz
  
