## Written by Eilis	
## Takes trimmed reads and aligns to genome with bowtie2	
## converts to bam files	
## excludes duplicates, mt reads, only keeps properly paired reads
## shifts read prior to peak calling

CHR=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)
cd ${DATADIRPE}

RAWDATADIR=($(find . -name '11_trimmed'))

for FOLDER in ${RAWDATADIR[@]}
do
	echo "Processing:"
	echo ${FOLDER}
	FQFILES=($(ls ${FOLDER}/*[rR]1*q.gz))	
	echo "Number of R1 .fq.gz files found for alignment:"" ""${#FQFILES[@]}"""	

	mkdir -p ${ALIGNEDDIR}	
		
	for f in ${FQFILES[@]};	
	do	
	  echo "Aligning"" ${f}"	
	  basename=$(basename $f)
	  sampleName=${basename/_trimmed_r1.fq.gz}
	  pairedFiles=($(ls ${FOLDER}/${sampleName}*.gz))
	  f1=${pairedFiles[0]}
	  f2=${pairedFiles[1]}

	  if [ ! -f ${sampleName}_hist.txt ]	
	  then
		## count uniqueness
		${BBMAP}/bbcountunique.sh in=${f1} in2=${f2} out=${sampleName}_hist.txt interval=5000 overwrite=true cumulative=true count=t
	  fi
	  
	  if [ ! -f ${ALIGNEDDIR}/${sampleName}_postFilter_statsperchr.txt ]	
	  then
		## alignment
		bowtie2  -p 10 -k ${multimap} -x ${REFGENOME}/genome -1 ${f1} -2 ${f2} -S ${ALIGNEDDIR}/${sampleName}.sam &> ${ALIGNEDDIR}/${sampleName}.bowtie.log
		
		## convert to sam files, sort and index
		samtools view -bSo ${ALIGNEDDIR}/${sampleName}.bam ${ALIGNEDDIR}/${sampleName}.sam
		samtools sort ${ALIGNEDDIR}/${sampleName}.bam > ${ALIGNEDDIR}/${sampleName}_sorted.bam
		samtools index ${ALIGNEDDIR}/${sampleName}_sorted.bam
		samtools idxstats ${ALIGNEDDIR}/${sampleName}_sorted.bam > ${ALIGNEDDIR}/${sampleName}_statsperchr.txt
		rm ${ALIGNEDDIR}/${sampleName}.bam
		rm ${ALIGNEDDIR}/${sampleName}.sam

		## remove reads from MT and other random chrs
		samtools view -b ${ALIGNEDDIR}/${sampleName}_sorted.bam ${CHR[@]} > ${ALIGNEDDIR}/${sampleName}_noMT.bam	
		
		## remove duplicates
		java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${ALIGNEDDIR}/${sampleName}_noMT.bam O=${ALIGNEDDIR}/${sampleName}_depDuplicated.bam M=${ALIGNEDDIR}/${sampleName}_dupMetrics.txt REMOVE_DUPLICATES=TRUE
		  
		## remove reads with q < 30;only keep properly paired reads; exclude reads marked as pcr optical duplicate, or secondary alignment
		samtools view -f 0x2 -b -F 0x400 -F 0x100 -q 30 -h ${ALIGNEDDIR}/${sampleName}_depDuplicated.bam > ${ALIGNEDDIR}/${sampleName}_depDup_q30.bam
		rm ${ALIGNEDDIR}/${sampleName}_depDuplicated.bam
		rm ${ALIGNEDDIR}/${sampleName}_noMT.bam
		
		samtools index ${ALIGNEDDIR}/${sampleName}_depDup_q30.bam
	  samtools idxstats ${ALIGNEDDIR}/${sampleName}_depDup_q30.bam > ${ALIGNEDDIR}/${sampleName}_postFilter_statsperchr.txt
	  fi	
	done
done
