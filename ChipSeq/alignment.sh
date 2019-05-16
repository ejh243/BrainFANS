## Written by Eilis
## Takes trimmed reads and aligns to genome with bowtie2

FQFILES=$(ls ${DATADIR}/*/11_trimmed/*.fq.gz)

echo "Number of .fq.gz files found for alignment:" "${#distro[@]}"

for f in ${FQFILES};
do
  echo "Aligning" ${f}
  basename=${f%.fq.gz}
  foldername=${basename//11_trimmed/alignedData}
  foldername=${foldername%/*}
  basename=${basename##*/}
  
  echo ${DATADIR}/${f}

  mkdir -p ${foldername}
  bowtie2 -x ${REFGENOME}/genome -U ${f} -S ${foldername}/${basename}.sam &> ${foldername}/${basename}.bowtie.log
 
done


