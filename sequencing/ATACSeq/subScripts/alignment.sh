#!/bin/bash
## ===================================================================================================================##
##                               ATAC-seq pipeline STEP 1.3: Pre-analysis -- alignment                                ##
## ===================================================================================================================##
## EXECUTION: sbatch ./subScripts/alignment.sh <sampleName>                                                           ||
## - execute from pipeline's main directory                                                                           ||
##                                                                                                                    ||
## DESCRIPTION: This scripts performs alignment of reads against the reference genome and subsequent processing of    ||  
## aligned reads, such as reads with quality < 30, unmapped, mate unmapped, secondary alignment, reads failing        || 
## platform, remove duplicates and reads mapped to a different chromosome                                             ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - Variables in config file: RAWDATADIR, TRIM_DIR, ALIGNED_DIR, REFGENOME, multimap                                 ||
## - Software: bowtie2, samtools, picard (in conda environment)                                                       ||
## - assign_multimappers.py file in ./subScripts folder.                                                              ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <sampleName> Name of sample specified in command line. Reads should be previously trimmed                    ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           || 
##  <sampleName>_filt.nodup.stats, <sampleName>_statsperchr.txt, <sampleName>_noMT.stats, <sampleName>_sorted.bam,    ||  
##  <sampleName>_sorted.bam.bai, <sampleName>_noMT.bam, <sampleName>_.filt.nodup.bam, <sampleName>_.filt.nodup.bam.bai|| 
##                                                                                                                    ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

## Sample to align
sampleName=$1
echo "Processing" ${sampleName}

## Check if trimmed files for alignment exist
cd ${TRIM_DIR}
f=($(ls ${sampleName}*trimmed*.f*))

if [[ ! -f ${TRIM_DIR}/${f[0]} ]] && [[ ! -f ${TRIM_DIR}/${f[1]} ]]
then 
  { echo "Trimmed reads for ${sampleName} not found. Please run STEP 1.2 TRIM first." ; exit 1; }
else
  echo "Found trimmed files:" ${f[0]} ${f[1]}
fi

## =============== ##
##    ALIGNMENT    ##
## =============== ##

CHR=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

echo "Running alignment for: ${sampleName}"
echo Job started on:
date -u	

bowtie2  -p 10 -X 2000 -k ${multimap} -x ${REFGENOME}/genome -1 ${f[0]} -2 ${f[1]} -S ${ALIGNED_DIR}/${sampleName}.sam &> ${ALIGNED_DIR}/${sampleName}.bowtie.log

## convert to sam files, sort and index
echo "converting to sam file, sorting and indexing"
samtools view -bSo ${ALIGNED_DIR}/${sampleName}.bam ${ALIGNED_DIR}/${sampleName}.sam
samtools sort ${ALIGNED_DIR}/${sampleName}.bam > ${ALIGNED_DIR}/${sampleName}_sorted.bam
samtools index ${ALIGNED_DIR}/${sampleName}_sorted.bam
samtools idxstats ${ALIGNED_DIR}/${sampleName}_sorted.bam > ${ALIGNED_DIR}/${sampleName}_statsperchr.txt
rm ${ALIGNED_DIR}/${sampleName}.bam
rm ${ALIGNED_DIR}/${sampleName}.sam

## remove reads from MT and other random chrs
echo "removing MT chr"
samtools view -b ${ALIGNED_DIR}/${sampleName}_sorted.bam ${CHR[@]} > ${ALIGNED_DIR}/${sampleName}_noMT.bam	

## get sum stats prior to filtering
samtools stats ${ALIGNED_DIR}/${sampleName}_noMT.bam	> ${ALIGNED_DIR}/QCOutput/${sampleName}_noMT.stats	

# remove reads with q < 30, unmapped, mate unmapped, secondary alignment, reads failing platform
# only keep properly paired reads
echo "filtering aligned reads"
samtools view -F 524 -f 2  -q 30 -u ${ALIGNED_DIR}/${sampleName}_noMT.bam | samtools sort -n /dev/stdin -o ${ALIGNED_DIR}/${sampleName}_q30.tmp.nmsrt.bam
samtools view -h ${ALIGNED_DIR}/${sampleName}_q30.tmp.nmsrt.bam | ${SUB_SCRIPTS_DIR}/assign_multimappers.py -k $multimap --paired-end | samtools fixmate -r /dev/stdin ${ALIGNED_DIR}/${sampleName}_q30.tmp.nmsrt.fixmate.bam

# Remove orphan reads (pair was removed)
# and read pairs mapping to different chromosomes
# Obtain position sorted BAM
samtools view -F 1804 -f 2 -u ${ALIGNED_DIR}/${sampleName}_q30.tmp.nmsrt.fixmate.bam | samtools sort /dev/stdin -o ${ALIGNED_DIR}/${sampleName}.filt.bam

rm ${ALIGNED_DIR}/${sampleName}_q30.tmp.nmsrt.fixmate.bam
rm ${ALIGNED_DIR}/${sampleName}_q30.tmp.nmsrt.bam
rm ${ALIGNED_DIR}/${sampleName}_noMT.bam

# Mark duplicates
samtools addreplacerg -r "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa" -o ${ALIGNED_DIR}/${sampleName}_filt.name.bam ${ALIGNED_DIR}/${sampleName}.filt.bam

picard MarkDuplicates -INPUT ${ALIGNED_DIR}/${sampleName}_filt.name.bam -OUTPUT ${ALIGNED_DIR}/${sampleName}.filt.dupmark.bam -METRICS_FILE ${ALIGNED_DIR}/${sampleName}_dupMetrics.txt \
-VALIDATION_STRINGENCY LENIENT -ASSUME_SORTED true -REMOVE_DUPLICATES false --TMP_DIR ${TMPDIR}

rm ${ALIGNED_DIR}/${sampleName}.filt.bam

# Remove duplicates
samtools view -F 1804 -f 2 -b ${ALIGNED_DIR}/${sampleName}.filt.dupmark.bam > ${ALIGNED_DIR}/${sampleName}.filt.nodup.bam

# Index Final BAM file
samtools index ${ALIGNED_DIR}/${sampleName}.filt.nodup.bam

## get sum stats post to filtering
samtools stats ${ALIGNED_DIR}/${sampleName}.filt.nodup.bam	> ${ALIGNED_DIR}/QCOutput/${sampleName}.filt.nodup.stats

if [[ ! -f ${ALIGNED_DIR}/${sampleName}.filt.nodup.bam ]] || [[ ! -f ${ALIGNED_DIR}/QCOutput/${sampleName}.filt.nodup.stats ]]
then 
  { echo "Alignment for ${sampleName} could not be completed." ; exit 1; }
else
  echo "Alignment for ${sampleName} is done."
  echo Job finished on:
  date -u 
fi 
