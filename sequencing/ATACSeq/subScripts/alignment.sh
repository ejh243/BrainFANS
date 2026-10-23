#!/bin/bash
## ===================================================================================================================##
##                                       PIPPA STEP 1.3: Pre-analysis -- alignment                                    ##
## ===================================================================================================================##
## EXECUTION: sbatch ./subScripts/alignment.sh <sampleName>                                                           ||
##                                                                                                                    ||
## DESCRIPTION: This scripts performs alignment of reads against the reference genome and subsequent processing of    ||  
## aligned reads, such as reads with quality < 30, unmapped, mate unmapped, secondary alignment, reads failing        || 
## platform, remove duplicates and reads mapped to a different chromosome                                             ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - Variables in config file: RAWDATADIR, TRIM_DIR, ALIGNED_DIR, REFGENOME, multimap, p_threads, max_frag_len, min_q ||
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
echo "Processing " ${sampleName}

## Check if trimmed files for alignment exist
cd ${TRIM_DIR}

R1=$(ls ${sampleName}*R1*trimmed*.f* | head -n 1)
R2=$(ls ${sampleName}*R2*trimmed*.f* | head -n 1)

if [[ ! -f ${TRIM_DIR}/$R1 ]] && [[ ! -f ${TRIM_DIR}/$R2 ]]
then 
 { echo "Trimmed reads for ${sampleName} not found. Please run STEP 1.2 TRIM first." ; exit 1; }
else
  echo "Found trimmed files:" ${f[0]} ${f[1]}
fi

if [[ ${#R1[@]} -ne 1 || ${#R2[@]} -ne 1 ]]; then
    echo "Error: missing or ambiguous FASTQ pairs"
    exit 1
fi

## =============== ##
##    ALIGNMENT    ##
## =============== ##

cat <<EOF
Running alignment for: ${sampleName}

Parameters specified in config file are: 
Number of threads: ${p_threads}
Maximum fragment length: ${max_frag_len}
Number of alignments reported for multimapping: ${multimap}
Job started on:
EOF
date -u	

cd ${TMPDIR}
echo "${TRIM_DIR}/$R1 ${TRIM_DIR}/$R2"
CHR=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
bowtie2 -p ${p_threads} -X ${max_frag_len} -k ${multimap} -x ${REFGENOME}/genome -1 "${TRIM_DIR}/$R1" -2 "${TRIM_DIR}/$R2" -S ${ALIGNED_DIR}/${sampleName}.sam &> ${ALIGNED_DIR}/${sampleName}.bowtie.log


## convert to sam files, sort and index
echo "Alignment completed. Converting to sam file, sorting and indexing"
samtools view -bSo "${ALIGNED_DIR}/${sampleName}.bam" "${ALIGNED_DIR}/${sampleName}.sam"
samtools sort "${ALIGNED_DIR}/${sampleName}.bam" > "${ALIGNED_DIR}/${sampleName}_sorted.bam"
samtools index "${ALIGNED_DIR}/${sampleName}_sorted.bam"
samtools idxstats "${ALIGNED_DIR}/${sampleName}_sorted.bam" > "${ALIGNED_DIR}/${sampleName}_statsperchr.txt"

## remove reads from MT and other random chrs
echo "Removing MT chr"
samtools view -b "${ALIGNED_DIR}/${sampleName}_sorted.bam" ${CHR[@]} > "${ALIGNED_DIR}/${sampleName}_noMT.bam"	
samtools idxstats "${ALIGNED_DIR}/${sampleName}_noMT.bam" > "${ALIGNED_DIR}/${sampleName}_statsperchr.txt"

## get sum stats prior to filtering
samtools stats "${ALIGNED_DIR}/${sampleName}_noMT.bam"	> "${ALIGNED_DIR}/QCOutput/${sampleName}_noMT.stats"

# remove reads with q < 30, unmapped, mate unmapped, secondary alignment, reads failing platform
# only keep properly paired reads
echo "Filtering aligned reads. Reads with quality below ${mapq_thres} will be excluded"
samtools view -F 524 -f 2 -q ${mapq_thres} -u "${ALIGNED_DIR}/${sampleName}_noMT.bam" | samtools sort -n /dev/stdin -o "${ALIGNED_DIR}/${sampleName}_q30.tmp.nmsrt.bam"

samtools view -h "${ALIGNED_DIR}/${sampleName}_q30.tmp.nmsrt.bam" \
  | "${SUB_SCRIPTS_DIR}/assign_multimappers.py" -k $multimap --paired-end \
  | samtools view -bS - \
  -o "${ALIGNED_DIR}/${sampleName}_q30.tmp.multimapped.bam"

samtools sort -n "${ALIGNED_DIR}/${sampleName}_q30.tmp.multimapped.bam" \
  -o "${ALIGNED_DIR}/${sampleName}_q30.tmp.multimapped.namesort.bam"

samtools fixmate -r  "${ALIGNED_DIR}/${sampleName}_q30.tmp.multimapped.namesort.bam" \
  "${ALIGNED_DIR}/${sampleName}_q30.tmp.nmsrt.fixmate.bam"

# Remove orphan reads (pair was removed)
# and read pairs mapping to different chromosomes
# Obtain position sorted BAM
samtools view -F 1804 -f 2 -u "${ALIGNED_DIR}/${sampleName}_q30.tmp.nmsrt.fixmate.bam" \
 | samtools sort /dev/stdin -o "${ALIGNED_DIR}/${sampleName}.filt.bam"

# Mark duplicates
samtools addreplacerg -r "@RG\tID:RG1\tSM:sampleName\tPL:Illumina\tLB:Library.fa" -o "${ALIGNED_DIR}/${sampleName}_filt.name.bam" \
 "${ALIGNED_DIR}/${sampleName}.filt.bam"

picard MarkDuplicates -I "${ALIGNED_DIR}/${sampleName}_filt.name.bam" \
  -O "${ALIGNED_DIR}/${sampleName}.filt.dupmark.bam" \
 --METRICS_FILE "${ALIGNED_DIR}/${sampleName}_dupMetrics.txt" \
  --VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true \
  --REMOVE_DUPLICATES false --TMP_DIR ${TMPDIR}

# Remove duplicates
samtools view -F 1804 -f 2 -b "${ALIGNED_DIR}/${sampleName}.filt.dupmark.bam" > "${ALIGNED_DIR}/${sampleName}.filt.nodup.bam"

# Index Final BAM file
samtools index "${ALIGNED_DIR}/${sampleName}.filt.nodup.bam"

## get sum stats post to filtering
samtools stats "${ALIGNED_DIR}/${sampleName}.filt.nodup.bam"	> "${ALIGNED_DIR}/QCOutput/${sampleName}.filt.nodup.stats"

## Outputs final reads metrics to output.log file
samtools flagstat "${ALIGNED_DIR}/${sampleName}.filt.nodup.bam"

## Remove intermediate files ##

rm "${ALIGNED_DIR}/${sampleName}.bam"
rm "${ALIGNED_DIR}/${sampleName}.sam"
rm "${ALIGNED_DIR}/${sampleName}_q30.tmp.nmsrt.fixmate.bam"
rm "${ALIGNED_DIR}/${sampleName}_q30.tmp.nmsrt.bam"
rm "${ALIGNED_DIR}/${sampleName}_q30.tmp.multimapped.namesort.bam"
rm "${ALIGNED_DIR}/${sampleName}_noMT.bam"
rm "${ALIGNED_DIR}/${sampleName}_q30.tmp.multimapped.bam"
rm "${ALIGNED_DIR}/${sampleName}.filt.bam"
rm "${ALIGNED_DIR}/${sampleName}_filt.name.bam"
rm "${ALIGNED_DIR}/${sampleName}_dupMetrics.txt"

if [[ ! -f "${ALIGNED_DIR}/${sampleName}.filt.nodup.bam" ]] || [[ ! -f "${ALIGNED_DIR}/QCOutput/${sampleName}.filt.nodup.stats" ]]
then 
  { echo "Alignment for ${sampleName} could not be completed." ; exit 1; }
else
  echo "Alignment for ${sampleName} is done."
  echo Job finished on:
  date -u 
fi 
