## Written by Eilis
## Uses shifted tagAlign files
## calls peaks per sample
## parameter choices guided by this post: https://github.com/taoliu/MACS/issues/331

cd ${ALIGNEDDIR}
mkdir -p ${PEAKDIR}
TAGFILES=$(ls *.tn5.tagAlign.gz)

echo "Number of bam files found for peak calling:" "${#BAMFILES[@]}"

for f in ${TAGFILES};
do
   echo "Peak Calling" ${f}
  sampleName=${f%.tn5.tagAlign.gz}
  macs2 callpeak -t ${f} --outdir ${PEAKDIR} -n ${sampleName} -f BED -g 2.9e9  -B --keep-dup all --shift 100 --extsize 200 --nomodel --broad

  ## exclude peaks aligned to blacklist regions
  bedtools intersect -v -a ${PEAKDIR}/${sampleName}.narrowPeak -b ${BLACKLIST} \
  | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
  | grep -P 'chr[\dXY]+[ \t]'  | gzip -nc > ${PEAKDIR}/${sampleName}.narrowPeak.filt.gz

done

