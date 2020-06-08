## Written by Eilis
## Takes trimmed reads and aligns to genome with bowtie2
## converts to bam files
## excludes duplicates 
## output: filtered bam file, tagAlign file

NTHREADS=8

cd ${DATADIR}
FQFILES=$(ls 11_trimmed/*q.gz)
mkdir -p ${ALIGNEDDIR}	

echo "Number of .fq.gz files found for alignment:" "${#FQFILES[@]}"

for f in ${FQFILES};
do
  echo "Aligning" ${f}
  fileName=$(basename ${f})
  sampleName=${fileName/_R*}
  if [ ! -f ${f//.*.gz/_hist.txt} ]
  then
	## count uniqueness 
	${BBMAP}/bbcountunique.sh in=${f} out=${f//.*.gz/_hist.txt} interval=5000 overwrite=true cumulative=true count=t k=31
  fi  
  if [ ! -f ${ALIGNEDDIR}/${sampleName}_depDup_q30.bam ]
  then
	  bowtie2 -x ${REFGENOME}/genome -p ${NTHREADS} -U ${f} -S ${ALIGNEDDIR}/${sampleName}.sam &> ${ALIGNEDDIR}/${sampleName}.bowtie.log

	  ## convert to sam files
	  samtools view -bSo ${ALIGNEDDIR}/${sampleName}.bam ${ALIGNEDDIR}/${sampleName}.sam
	  samtools sort ${ALIGNEDDIR}/${sampleName}.bam > ${ALIGNEDDIR}/${sampleName}_sorted.bam
	  ## index bam files for QC package
	  samtools index ${ALIGNEDDIR}/${sampleName}_sorted.bam
	  rm ${ALIGNEDDIR}/${sampleName}.sam
	  rm ${ALIGNEDDIR}/${sampleName}.bam

	  ## remove duplicates
	  java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${ALIGNEDDIR}/${sampleName}_sorted.bam O=${ALIGNEDDIR}/${sampleName}_depDuplicated.bam M=${ALIGNEDDIR}/${sampleName}_dupMetrics.txt REMOVE_DUPLICATES=TRUE
	  
	  ## remove reads with q < 30 nb 
	  samtools view -q 30 -h ${ALIGNEDDIR}/${sampleName}_depDuplicated.bam > ${ALIGNEDDIR}/${sampleName}_depDup_q30.bam
	  rm ${ALIGNEDDIR}/${sampleName}_depDuplicated.bam
  fi
  if [ ! -f ${ALIGNEDDIR}/${sampleName}.SE.tagAlign.gz ]		
  then	  
   	# Create tagAlign file nb this is what is used in ENCODE pipeline (although subject to additional filters)
    bedtools bamtobed -i ${ALIGNEDDIR}/${sampleName}_depDup_q30.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -nc > ${ALIGNEDDIR}/${sampleName}.SE.tagAlign.gz
  fi
done

