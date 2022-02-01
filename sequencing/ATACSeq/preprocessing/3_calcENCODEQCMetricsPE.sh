## calculates ENCODE library complexity metrics (taken from https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit#)

## EXECUTION
# sh ./ATACSeq/preprocessing/3_calcENCODEQCMetricsPE.sh <bam file>
# where 
# <bam file> is the path to a sorted bam file
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, multimap

## REQUIRES the following software
# samtools, samstats, picard, bedtools

## INPUT
# aligned sorted bam file

## OUTPUT
# *.pbc.qc
# ENCODEMetrics/*.bam
# *.flagstat.qc

RAW_BAM_FILE=$1 

echo "Calculating ENCODE QC metrics"
echo Job started on:
date -u

cd ${ALIGNEDDIR}

mkdir -p ENCODEMetrics

# =============================
# Remove  unmapped, mate unmapped
# not primary alignment, reads failing platform
# Only keep properly paired reads
# Obtain name sorted BAM file
# ==================

OFPREFIX=${RAW_BAM_FILE%_sorted_chr1.bam}	
FILT_BAM_PREFIX="${OFPREFIX}.filt"
FILT_BAM_FILE="ENCODEMetrics/${FILT_BAM_PREFIX}.bam"
TMP_FILT_BAM_PREFIX="${FILT_BAM_PREFIX}.tmp.nmsrt"
TMP_FILT_BAM_FILE="ENCODEMetrics/${TMP_FILT_BAM_PREFIX}.bam"
TMP_FILT_FIXMATE_BAM_FILE="ENCODEMetrics/${TMP_FILT_BAM_PREFIX}.fixmate.bam"
if [ ${multimap} -gt 0 ]
then 
  samtools view -F 524 -f 2 -u ${RAW_BAM_FILE} | samtools sort -n /dev/stdin -o ${TMP_FILT_BAM_FILE}
  samtools view -h ${TMP_FILT_BAM_FILE} | $(which assign_multimappers.py) -k $multimap --paired-end | samtools fixmate -r /dev/stdin ${TMP_FILT_FIXMATE_BAM_FILE}
else
  ## for multimap==0:
  MAPQ_THRESH=255
  samtools view -F 1804 -f 2 -q ${MAPQ_THRESH} -u ${RAW_BAM_FILE} | samtools sort -n /dev/stdin -o ${TMP_FILT_BAM_FILE} 
  samtools fixmate -r ${TMP_FILT_BAM_FILE} ${TMP_FILT_FIXMATE_BAM_FILE}
fi

echo "Step 1 complete on:"
date -u

# Remove orphan reads (pair was removed)
# and read pairs mapping to different chromosomes
# Obtain position sorted BAM
samtools view -F 1804 -f 2 -u ${TMP_FILT_FIXMATE_BAM_FILE} | samtools sort /dev/stdin -o ${FILT_BAM_FILE}

rm  ${TMP_FILT_FIXMATE_BAM_FILE}
rm ${TMP_FILT_BAM_FILE}

echo "Step 2 complete on:"
date -u

# =============
# Mark duplicates
# =============

TMP_FILT_BAM_FILE="ENCODEMetrics/${FILT_BAM_PREFIX}.dupmark.bam"
DUP_FILE_QC="ENCODEMetrics/${FILT_BAM_PREFIX}.dup.qc"

java -Xmx4G -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false

mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}

echo "Step 3 complete on:"
date -u

# ============================
# Remove duplicates
# Index final position sorted BAM
# Create final name sorted BAM
# ============================
FINAL_BAM_PREFIX="ENCODEMetrics/${OFPREFIX}.nodup"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${FINAL_BAM_FILE}.bai"
FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file

samtools view -F 1804 -f 2 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}

# Index Final BAM file
samtools index ${FINAL_BAM_FILE}

samtools sort -n --threads 10 ${FINAL_BAM_FILE} -O SAM  | SAMstats --sorted_sam_file -  --outf ${FINAL_BAM_FILE_MAPSTATS}

echo "Step 4 complete on:"
date -u

# =============================
# Compute library complexity
# =============================
# Sort by name
# convert to bedPE and obtain fragment coordinates
# sort by position and strand
# Obtain unique count statistics

PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"
# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair

samtools sort -n ${FILT_BAM_FILE} -o ${OFPREFIX}.srt.tmp.bam

bedtools bamtobed -bedpe -i ${OFPREFIX}.srt.tmp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{if (mt > 0) {printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2} else printf "%d\t%d\t%d\t%d\t%s\t%s\t%s\n",mt,m0,m1,m2,NA,NA,NA }'  > ${PBC_FILE_QC}

rm ${OFPREFIX}.srt.tmp.bam
rm ${FILT_BAM_FILE}

echo "QC metrics calculated:"
date -u