## Takes trimmed reads for single sample aligns with bowtie2 and filters out duplicates, mt reads, incorrectly paired reads

## EXECUTION
# sh ./ChIPSeq/preprocessing/1_alignment.sh <fastq file>
# where 
# <fastq file> is the path to the "R1" fastq files which are expected to be compressed
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# RAWDATADIR, TRIMDIR, ALIGNEDDIR, REFGENOME, multimap

## REQUIRES the following software
# bowtie2, samtools, picard

## INPUT
# 2 trimmed fastq files

## OUTPUT
# *_statsperchr.txt
# *_dupMetrics.txt
# *_depDup_q30.bam
# *_postFilter_statsperchr.txt
# *.filt.nodup.bam
# *.filt.dupmark.bam
# *.filt.nodup.bam.bai 
# *_sorted.bam.bai


sampleName=$1
echo
echo "Starting alignment on" ${sampleName} "at: "
date -u

cd ${TRIMDIR}
f=($(ls ${sampleName}*trimmed*.f*))


if [ ! -f ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam ]		
then
  cd ${TRIMDIR}

  # Alignment
  bowtie2 -p 10 -x ${REFGENOME}/genome -1 ${f[0]} -2 ${f[1]}  -S ${ALIGNEDDIR}/${sampleName}.sam &> ${ALIGNEDDIR}/${sampleName}.bowtie.log

  ## Convert to sam file, sort and index
  samtools view -bSo ${ALIGNEDDIR}/${sampleName}.bam ${ALIGNEDDIR}/${sampleName}.sam
  samtools sort ${ALIGNEDDIR}/${sampleName}.bam > ${ALIGNEDDIR}/${sampleName}_sorted.bam

  samtools index ${ALIGNEDDIR}/${sampleName}_sorted.bam
  samtools idxstats ${ALIGNEDDIR}/${sampleName}_sorted.bam > ${ALIGNEDDIR}/${sampleName}_statsperchr.txt
  rm ${ALIGNEDDIR}/${sampleName}.bam
  rm ${ALIGNEDDIR}/${sampleName}.sam

  ## Get sum stats prior to filtering
  samtools stats ${ALIGNEDDIR}/${sampleName}_sorted.bam > ${ALIGNEDDIR}/QCOutput/${sampleName}_sorted.stats 

  # Remove reads with q < 30, unmapped, mate unmapped, secondary alignment, reads failing platform
  # only keep properly paired reads
  # Obtain name sorted bam
  echo "Filtering aligned reads"
  samtools view -q 30 -h ${ALIGNEDDIR}/${sampleName}_sorted.bam | samtools sort -n /dev/stdin -o ${ALIGNEDDIR}/${sampleName}_q30.tmp.nmsrt.bam
  samtools view -h ${ALIGNEDDIR}/${sampleName}_q30.tmp.nmsrt.bam | samtools fixmate -r /dev/stdin ${ALIGNEDDIR}/${sampleName}_q30.tmp.nmsrt.fixmate.bam
  
  # Remove orphan reads (pair was removed)
  # Obtain position sorted BAM
  samtools view -F 1804 -f 2 -u ${ALIGNEDDIR}/${sampleName}_q30.tmp.nmsrt.fixmate.bam | samtools sort /dev/stdin -o ${ALIGNEDDIR}/${sampleName}.filt.bam

  rm ${ALIGNEDDIR}/${sampleName}_q30.tmp.nmsrt.fixmate.bam
  rm ${ALIGNEDDIR}/${sampleName}_q30.tmp.nmsrt.bam
  rm ${ALIGNEDDIR}/${sampleName}_sorted.bam

  # Mark duplicates
  java -Xmx4G -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=${ALIGNEDDIR}/${sampleName}.filt.bam OUTPUT=${ALIGNEDDIR}/${sampleName}.filt.dupmark.bam METRICS_FILE=${ALIGNEDDIR}/${sampleName}_dupMetrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
  rm ${ALIGNEDDIR}/${sampleName}.filt.bam

  ## Remove duplicates
  samtools view -F 1804 -f 2 -b ${ALIGNEDDIR}/${sampleName}.filt.dupmark.bam > ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam

  # Index Final BAM file
  samtools index ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam
   
  ## get sum stats post to filtering
  samtools stats ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam  > ${ALIGNEDDIR}/QCOutput/${sampleName}.filt.nodup.stats
else
  { echo "Aligned file found so not aligning"; exit 1;}
fi


if [[ $? == 0 ]]
  then echo "Finished alignment"
fi