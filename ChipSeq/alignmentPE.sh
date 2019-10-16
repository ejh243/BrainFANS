## Written by Eilis	
## Takes trimmed reads and aligns to genome with bowtie2	
## converts to bam files	
## excludes duplicates 	
	
FQFILES=$(ls ${DATADIRPE}/M*/11_trimmed/*_r1.fq.gz)	
	
echo "Number of .fq.gz files found for alignment:"" ""${#FQFILES[@]}"""	
	
	
for f in ${FQFILES};	
do	
  echo "Aligning"" ${f}"	
  basename=${f%_r1.fq.gz}	
  foldername=${basename//11_trimmed/alignedData}
  foldername=${foldername%/*}
  basename=${basename##*/}

  if [ -f ${foldername}/${basename}_depDup_q30.bam ]	
  then
    echo "already aligned"
  else
    ## count uniqueness
	f2=${f//r1/r2}
	${BBMAP}/bbcountunique.sh in=${f} in2=${f2} out=${f//_r1.fq.gz/_hist.txt} interval=5000 overwrite=true cumulative=true count=t

  
    mkdir -p ${foldername}
	
	bowtie2 -x ${REFGENOME}/genome -1 ${f} -2 ${f2} -S ${foldername}/${basename}.sam &> ${foldername}/${basename}.bowtie.log
	## convert to sam files
	samtools view -bSo ${foldername}/${basename}.bam ${foldername}/${basename}.sam
	samtools sort ${foldername}/${basename}.bam > ${foldername}/${basename}_sorted.bam
	samtools index ${foldername}/${basename}_sorted.bam
	rm ${foldername}/${basename}.sam
	
	## remove duplicates
	java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${foldername}/${basename}_sorted.bam O=${foldername}/${basename}_depDuplicated.bam M=${foldername}/${basename}_dupMetrics.txt REMOVE_DUPLICATES=TRUE
	  
	## remove reads with q < 30 nb 
	samtools view -q 30 -h ${foldername}/${basename}_depDuplicated.bam > ${foldername}/${basename}_depDup_q30.bam
  fi	
done	
