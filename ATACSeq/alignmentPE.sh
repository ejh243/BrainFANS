## script to process a single sample
## Requires a raw (r1) fastq file provided on the command line 
## Takes trimmed reads and aligns to genome with bowtie2	
## converts to bam files	
## excludes duplicates, mt reads, only keeps properly paired reads
## shifts read prior to peak calling


CHR=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX)

f=$1 
cd ${DATADIRPE}

## extract sample names
FOLDER=$(dirname ${f})
f1=$(basename $f)
sampleName=${f1%_[rR]*}
## later samples have an additional _S[num] in the file name need to remove
sampleName=${sampleName%_S[0-9]*}
echo "Processing" ${sampleName}

## need trimmed files for alignment

## filename format changed sometime r1_trimmed sometimes trimmed_r1
if ls ${FOLDERTRIM}/${sampleName}*[rR]1*_trimmed.f* 1> /dev/null 2>&1;
then 
  f1=$(basename $(ls ${FOLDERTRIM}/${sampleName}*[rR]1*_trimmed.f*))
  f2=$(basename $(ls ${FOLDERTRIM}/${sampleName}*[rR]2*_trimmed.f*))
else
  f1=$(basename $(ls ${FOLDERTRIM}/${sampleName}*_trimmed_[rR]1*.f*))
  f2=$(basename $(ls ${FOLDERTRIM}/${sampleName}*_trimmed_[rR]2*.f*))
fi

echo "Found trimmed files:"
echo ${f1}
echo ${f2}

## checks if last file exists and is not empty
if [ ! -s ${ALIGNEDDIR}/${sampleName}_postFilter_statsperchr.txt ]	
	then
	echo "Running alignment for"" ${sampleName}"
	date -u	
	## alignment
	bowtie2  -p 10 -X 1000 -k ${multimap} -x ${REFGENOME}/genome -1 ${FOLDERTRIM}/${f1} -2 ${FOLDERTRIM}/${f2} -S ${ALIGNEDDIR}/${sampleName}.sam &> ${ALIGNEDDIR}/${sampleName}.bowtie.log

	## convert to sam files, sort and index
	echo "converting to sam file, sorting and indexing"
	samtools view -bSo ${ALIGNEDDIR}/${sampleName}.bam ${ALIGNEDDIR}/${sampleName}.sam
	samtools sort ${ALIGNEDDIR}/${sampleName}.bam > ${ALIGNEDDIR}/${sampleName}_sorted.bam
	samtools index ${ALIGNEDDIR}/${sampleName}_sorted.bam
	samtools idxstats ${ALIGNEDDIR}/${sampleName}_sorted.bam > ${ALIGNEDDIR}/${sampleName}_statsperchr.txt
	rm ${ALIGNEDDIR}/${sampleName}.bam
	rm ${ALIGNEDDIR}/${sampleName}.sam

	## extract chr1 for calculation of ENCODE metrics
	echo "extracting chr 1 for ENCODE metrics calculations"
	samtools view -b ${ALIGNEDDIR}/${sampleName}_sorted.bam chr1 > ${ALIGNEDDIR}/${sampleName}_sorted_chr1.bam
	
	## remove reads from MT and other random chrs
	echo "removing MT chr"
	samtools view -b ${ALIGNEDDIR}/${sampleName}_sorted.bam ${CHR[@]} > ${ALIGNEDDIR}/${sampleName}_noMT.bam	

	## remove duplicates
	echo "removing duplicates"
	java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${ALIGNEDDIR}/${sampleName}_noMT.bam O=${ALIGNEDDIR}/${sampleName}_depDuplicated.bam M=${ALIGNEDDIR}/${sampleName}_dupMetrics.txt REMOVE_DUPLICATES=TRUE
	  
	## remove reads with q < 30;only keep properly paired reads; exclude reads marked as pcr optical duplicate, or secondary alignment
	echo "filtering aligned reads"
	samtools view -f 0x2 -b -F 0x400 -F 0x100 -q 30 -h ${ALIGNEDDIR}/${sampleName}_depDuplicated.bam > ${ALIGNEDDIR}/${sampleName}_depDup_q30.bam
	rm ${ALIGNEDDIR}/${sampleName}_depDuplicated.bam
	rm ${ALIGNEDDIR}/${sampleName}_noMT.bam

	samtools index ${ALIGNEDDIR}/${sampleName}_depDup_q30.bam
	samtools idxstats ${ALIGNEDDIR}/${sampleName}_depDup_q30.bam > ${ALIGNEDDIR}/${sampleName}_postFilter_statsperchr.txt
fi

echo "Alignment and post filtering complete"
date -u

