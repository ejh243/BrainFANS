## Written by Eilis
## Takes advantage of paired end reads

BAMFILES=$(ls ${DATADIR}/alignedData/*_depDup_q30.bam)

echo "Number of bam files found for peak calling:" "${#BAMFILES[@]}"

for f in ${BAMFILES};
do
   echo "Peak Calling" ${f}
  basename=${f%.bam}
  foldername=${basename//alignedData/MACS2}
  foldername=${foldername%/*}
  basename=${basename##*/}

  mkdir -p ${foldername}
  macs2 callpeak -t ${f} --outdir ${foldername} -n ${basename} -f BAMPE -g 2.9e9  -B --keep-dup all
  
  ## exclude peaks aligned to blacklist regions
done

