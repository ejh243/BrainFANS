## Takes trimmed reads for single sample aligns with bowtie2 and filters out duplicates, mt reads, incorrectly paired reads

## EXECUTION
# sh ./ATACSeq/preprocessing/2_alignment.sh <fastq file>
# where 
# <fastq file> is the path to the "R1" fastq files which are expected to be compressed, and have either r1 or R1 in the filename
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# RAWDATADIR, TRIMDIR, ALIGNEDDIR, REFGENOME, multimap

## REQUIRES the following software
# bowtie2, samtools, picard

## INPUT
# 2 trimmed fastq files

## OUTPUT
# *_statsperchr.txt
# *_sorted.bam
# *_sorted_chr1.bam
# *_dupMetrics.txt
# *_depDup_q30.bam
# *_postFilter_statsperchr.txt

CHR=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

sampleName=$1

cd ${RAWDATADIR}

echo "Processing" ${sampleName}

## need trimmed files for alignment

## filename format changed sometime r1_trimmed sometimes trimmed_r1
if ls ${TRIMDIR}/${sampleName}*1*_trimmed.f* 1> /dev/null 2>&1;
then 
  f1=$(basename $(ls ${TRIMDIR}/${sampleName}*1*_trimmed.f*))
  f2=$(basename $(ls ${TRIMDIR}/${sampleName}*2*_trimmed.f*))
else
  f1=$(basename $(ls ${TRIMDIR}/${sampleName}*_trimmed_*1*.f*))
  f2=$(basename $(ls ${TRIMDIR}/${sampleName}*_trimmed_*2*.f*))
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
	bowtie2  -p 10 -X 2000 -k ${multimap} -x ${REFGENOME}/genome -1 ${TRIMDIR}/${f1} -2 ${TRIMDIR}/${f2} -S ${ALIGNEDDIR}/${sampleName}.sam &> ${ALIGNEDDIR}/${sampleName}.bowtie.log

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

	## remove reads with q < 30;only keep properly paired reads; exclude reads marked as pcr optical duplicate, or secondary alignment
	echo "filtering aligned reads"
	samtools view -f 0x2 -b -F 0x400 -F 0x100 -q 30 -h ${ALIGNEDDIR}/${sampleName}_noMT.bam > ${ALIGNEDDIR}/${sampleName}_q30.bam

	## remove duplicates
	echo "removing duplicates"
	java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${ALIGNEDDIR}/${sampleName}_q30.bam O=${ALIGNEDDIR}/${sampleName}_depDup_q30.bam M=${ALIGNEDDIR}/${sampleName}_dupMetrics.txt REMOVE_DUPLICATES=TRUE
	
 
	rm ${ALIGNEDDIR}/${sampleName}_q30.bam
	rm ${ALIGNEDDIR}/${sampleName}_noMT.bam

	samtools index ${ALIGNEDDIR}/${sampleName}_depDup_q30.bam
	samtools idxstats ${ALIGNEDDIR}/${sampleName}_depDup_q30.bam > ${ALIGNEDDIR}/${sampleName}_postFilter_statsperchr.txt
else
	{ echo "Aligned file found so not aligning"; exit 1;}

fi

echo "Alignment and post filtering complete"
date -u

