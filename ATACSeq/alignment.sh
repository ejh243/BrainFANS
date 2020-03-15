## Written by Eilis	
## Takes trimmed reads and aligns to genome with bowtie2	
## converts to bam files	
## excludes duplicates 	

CHR=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)

FQFILES=$(ls ${DATADIR}/11_trimmed/*_r1.fq.gz)	
	
echo "Number of .fq.gz files found for alignment:"" ""${#FQFILES[@]}"""	
	
	
for f in ${FQFILES};	
do	
  echo "Aligning"" ${f}"	
  basename=${f%_r1.fq.gz}	
  foldername=${basename//11_trimmed/alignedData}
  foldername=${foldername%/*}
  basename=${basename##*/}

  #if [ -f ${foldername}/${basename}_depDup_q30.bam ]	
  #then
   # echo "already aligned"
  #else
    ## count uniqueness
	f2=${f//r1/r2}
	${BBMAP}/bbcountunique.sh in=${f} in2=${f2} out=${f//_r1.fq.gz/_hist.txt} interval=5000 overwrite=true cumulative=true count=t

    mkdir -p ${foldername}
	
	## alignment
	bowtie2  -p 10 -k ${multimap} -x ${REFGENOME}/genome -1 ${f} -2 ${f2} -S ${foldername}/${basename}.sam &> ${foldername}/${basename}.bowtie.log
	
	## convert to sam files, sort and index
	samtools view -bSo ${foldername}/${basename}.bam ${foldername}/${basename}.sam
	samtools sort ${foldername}/${basename}.bam > ${foldername}/${basename}_sorted.bam
	samtools index ${foldername}/${basename}_sorted.bam
  samtools idxstats ${foldername}/${basename}_sorted.bam > ${foldername}/${basename}_statsperchr.txt
	rm ${foldername}/${basename}.bam
	rm ${foldername}/${basename}.sam

	## remove reads from MT and other random chrs
	samtools view -b ${foldername}/${basename}_sorted.bam ${CHR[@]} > ${foldername}/${basename}_noMT.bam	
	
	## remove duplicates
	java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${foldername}/${basename}_noMT.bam O=${foldername}/${basename}_depDuplicated.bam M=${foldername}/${basename}_dupMetrics.txt REMOVE_DUPLICATES=TRUE
	  
	## remove reads with q < 30;only keep properly paired reads; exclude reads marked as pcr optical duplicate, or secondary alignment
	samtools view -f 0x2 -b -F 0x400 -F 0x100 -q 30 -h ${foldername}/${basename}_depDuplicated.bam > ${foldername}/${basename}_depDup_q30.bam
	rm ${foldername}/${basename}_depDuplicated.bam
	rm ${foldername}/${basename}_noMT.bam
	
	samtools index ${foldername}/${basename}_depDup_q30.bam
  samtools idxstats ${foldername}/${basename}_depDup_q30.bam > ${foldername}/${basename}_postFilter_statsperchr.txt
  
  ## shift reads prior to peak calling
  bedtools bamtobed -i ${foldername}/${basename}_depDup_q30.bam -bedpe | awk -v OFS="\t" '{if($9=="+"){print $1,$2+4,$6+4}else if($9=="-"){print $1,$2-5,$6-5}}' > ${foldername}/${basename}_shifted.bed

	
  #fi	
done	
