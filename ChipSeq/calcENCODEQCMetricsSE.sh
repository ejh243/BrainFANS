## taken from ENCODE pipeline https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#
## Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform,  reads mapping to chrM.
## Remove multi-mapped reads (i.e. those with MAPQ < 30, using -q in SAMtools)
## Remove PCR duplicates (using Picard’s MarkDuplicates or FixSeq)

## input: Raw BAM file ${RAW_BAM_FILE}

cd ${ALIGNEDDIR}
BAMFILES=($(ls *_sorted.bam))
MAPQ_THRESH=30
OUTFOLDER=encodeQCMetrics
mkdir -p ${OUTFOLDER}

for RAW_BAM_FILE in ${BAMFILES[@]};	
do
  OFPREFIX=${RAW_BAM_FILE%_sorted.bam}
  if [ ! -f ${OUTFOLDER}/${OFPREFIX}*pbc.qc ]
  then
    # =============================
    # Remove  unmapped, mate unmapped
    # not primary alignment, reads failing platform
    # Remove low MAPQ reads
    # ==================  
    FILT_BAM_PREFIX="${OFPREFIX}.filt.srt"
    FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
  
    samtools view -F 1804 -q ${MAPQ_THRESH} -b ${RAW_BAM_FILE} -o ${FILT_BAM_FILE}

    # ========================
    # Mark duplicates
    # ======================

    TMP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.dupmark.bam"
    DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc" # QC file

    java -Xmx4G -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false

    mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}

    # ============================
    # Remove duplicates
    # Index final position sorted BAM
    # ============================
    FINAL_BAM_PREFIX="${OFPREFIX}.filt.nodup.srt"
    FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
    FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai" # To be stored
    FINAL_BAM_FILE_MAPSTATS="${OUTFOLDER}/${FINAL_BAM_PREFIX}.flagstat.qc" # QC file

    samtools view -F 1804 -b ${FILT_BAM_FILE} -o ${FINAL_BAM_FILE}

    # Index Final BAM file
    samtools index ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX_FILE}
    samtools sort -n --threads 10 ${FINAL_BAM_FILE} -O SAM  | SAMstats --sorted_sam_file -  --outf ${FINAL_BAM_FILE_MAPSTATS}

    # =============================
    # Compute library complexity
    # =============================
    # sort by position and strand
    # Obtain unique count statistics

    PBC_FILE_QC="${OUTFOLDER}/${FINAL_BAM_PREFIX}.pbc.qc"

    # PBC File output
    # TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair

    bedtools bamtobed -i ${FILT_BAM_FILE} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}
    rm ${FILT_BAM_FILE}

  fi
 done
