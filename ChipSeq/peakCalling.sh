## Written by Eilis


BAMFILES=$(ls ${DATADIR}/*/alignedData/*_depDup_q30.bam)

echo "Number of bam files found for peak calling:" "${#distro[@]}"

for f in ${BAMFILES};
do
  echo "Peak Calling" ${f}
  basename=${f%.bam}
  foldername=${basename//alignedData/MACS2}
  foldername=${foldername%/*}
  basename=${basename##*/}
  
  echo ${DATADIR}/${f}

  mkdir -p ${foldername}

  macs2 callpeak -t ${f} --outdir ${foldername} -n ${basename} -g 2.9e9  -B

done


