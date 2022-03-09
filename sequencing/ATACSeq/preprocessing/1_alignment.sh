## Takes trimmed reads for single sample aligns with bowtie2 and filters out duplicates, mt reads, incorrectly paired reads

## EXECUTION
# sh ./ATACSeq/preprocessing/1_alignment.sh <fastq file>
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

cd ${ALIGNEDDIR}

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
if [ ! -s ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam ]	
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

	## remove reads from MT and other random chrs
	echo "removing MT chr"
	samtools view -b ${ALIGNEDDIR}/${sampleName}_sorted.bam ${CHR[@]} > ${ALIGNEDDIR}/${sampleName}_noMT.bam	

	## get sum stats prior to filtering
	samtools stats ${ALIGNEDDIR}/${sampleName}_noMT.bam	> ${ALIGNEDDIR}/QCOutput/${sampleName}_noMT.stats	

	# remove reads with q < 30, unmapped, mate unmapped, secondary alignment, reads failing platform
	# only keep properly paired reads
	# Obtain name sorted BAM file
	echo "filtering aligned reads"
    samtools view -F 524 -f 2  -q 30 -u ${ALIGNEDDIR}/${sampleName}_noMT.bam | samtools sort -n /dev/stdin -o ${ALIGNEDDIR}/${sampleName}_q30.tmp.nmsrt.bam
    samtools view -h ${ALIGNEDDIR}/${sampleName}_q30.tmp.nmsrt.bam | $(which assign_multimappers.py) -k $multimap --paired-end | samtools fixmate -r /dev/stdin ${ALIGNEDDIR}/${sampleName}_q30.tmp.nmsrt.fixmate.bam
	
	# Remove orphan reads (pair was removed)
	# and read pairs mapping to different chromosomes
	# Obtain position sorted BAM
	samtools view -F 1804 -f 2 -u ${ALIGNEDDIR}/${sampleName}_q30.tmp.nmsrt.fixmate.bam | samtools sort /dev/stdin -o ${ALIGNEDDIR}/${sampleName}.filt.bam
	
	rm ${ALIGNEDDIR}/${sampleName}_q30.tmp.nmsrt.fixmate.bam
	rm ${ALIGNEDDIR}/${sampleName}_q30.tmp.nmsrt.bam
	rm ${ALIGNEDDIR}/${sampleName}_noMT.bam
	rm ${ALIGNEDDIR}/${sampleName}_sorted.bam
	
   # Mark duplicates
   java -Xmx4G -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=${ALIGNEDDIR}/${sampleName}.filt.bam OUTPUT=${ALIGNEDDIR}/${sampleName}.filt.dupmark.bam METRICS_FILE=${ALIGNEDDIR}/${sampleName}_dupMetrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
   
   rm ${ALIGNEDDIR}/${sampleName}.filt.bam

   # Remove duplicates
   samtools view -F 1804 -f 2 -b ${ALIGNEDDIR}/${sampleName}.filt.dupmark.bam > ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam
   
   # Index Final BAM file
   samtools index ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam
   
   ## get sum stats post to filtering
   samtools stats ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam	> ${ALIGNEDDIR}/QCOutput/${sampleName}filt.nodup.stats
else
	{ echo "Aligned file found so not aligning"; exit 1;}
fi

echo "Alignment and post filtering complete"
date -u

