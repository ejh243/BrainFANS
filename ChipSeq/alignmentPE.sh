## Written by Eilis	
## Takes trimmed reads and aligns to genome with bowtie2	
## converts to bam files	
## excludes duplicates 	
## output: filtered bam file, tagAlign file (virtual single end), BEDPE file (with read pairs on each line) 

NTHREADS=8

cd ${DATADIRPE}
FQFILES=($(ls 11_trimmed/*[rR]1*q.gz))
mkdir -p ${ALIGNEDDIR}	

	
echo "Number of .fq.gz files found for alignment:"" ""${#FQFILES[@]}"""	

	
for f in ${FQFILES[@]};	
do	
  echo "Aligning"" ${f}"	
  fileName=$(basename ${f})
  sampleName=${fileName/_[rR]*}
  pairedFiles=($(ls 11_trimmed/${sampleName}*.gz))
  f1=${pairedFiles[0]}
  f2=${pairedFiles[1]}
	
  ## count uniqueness of fastq files
  if [ ! -f ${sampleName}_hist.txt ]		
  then
    ## count uniqueness
    ${BBMAP}/bbcountunique.sh in=${f1} in2=${f2} out=${sampleName}_hist.txt interval=5000 overwrite=true cumulative=true count=t
  fi	
  if [ ! -f ${ALIGNEDDIR}/${sampleName}_depDup_q30.bam ]		
  then
    bowtie2 -x ${REFGENOME}/genome -1 ${f1} -2 ${f2} -p ${NTHREADS} -S ${ALIGNEDDIR}/${sampleName}.sam &> ${ALIGNEDDIR}/${sampleName}.bowtie.log

    ## convert to sam file, sort and index
    samtools view -bSo ${ALIGNEDDIR}/${sampleName}.bam ${ALIGNEDDIR}/${sampleName}.sam
    samtools sort ${ALIGNEDDIR}/${sampleName}.bam > ${ALIGNEDDIR}/${sampleName}_sorted.bam
    samtools index ${ALIGNEDDIR}/${sampleName}_sorted.bam
    rm ${ALIGNEDDIR}/${sampleName}.sam
    rm ${ALIGNEDDIR}/${sampleName}.bam
    ## remove duplicates
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${ALIGNEDDIR}/${sampleName}_sorted.bam O=${ALIGNEDDIR}/${sampleName}_depDuplicated.bam M=${ALIGNEDDIR}/${sampleName}_dupMetrics.txt REMOVE_DUPLICATES=TRUE

    ## remove reads with q < 30 nb 
    samtools view -q 30 -h ${ALIGNEDDIR}/${sampleName}_depDuplicated.bam > ${ALIGNEDDIR}/${sampleName}_depDup_q30.bam
  fi
  #if [ ! -f ${ALIGNEDDIR}/${sampleName}.PE.tagAlign.gz ]		
  #then	
    # Create BEDPE file
    #samtools sort -n ${ALIGNEDDIR}/${sampleName}_depDup_q30.bam -o ${ALIGNEDDIR}/${sampleName}.filt.nmsrt
    #bedtools bamtobed -bedpe -mate1 -i ${ALIGNEDDIR}/${sampleName}.filt.nmsrt | gzip -nc > ${ALIGNEDDIR}/${sampleName}.bedpe.gz
    #zcat ${ALIGNEDDIR}/${sampleName}.bedpe.gz | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",$1,$2,$3,$9,$4,$5,$6,$10}' | gzip -nc > ${ALIGNEDDIR}/${sampleName}.PE.tagAlign.gz
  #fi
done	
