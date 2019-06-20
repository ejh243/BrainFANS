## create sampleSheet for ChipQC


FQFILES=$(ls ${DATADIR}/H*/11_trimmed/*.fq.gz)

echo "SampleID,bamReads,Peaks" > ${DATADIR}/SampleSheet.csv

for f in ${FQFILES};
do
  basename=${f%.fq.gz}
  foldername=${basename//11_trimmed/alignedData}
  foldername=${foldername%/*}
  basename=${basename##*/}
  echo ${basename},${foldername}/${basename}.bam","${foldername}/MACS2/${basename}.bed >> ${DATADIR}/SampleSheet.csv
done
