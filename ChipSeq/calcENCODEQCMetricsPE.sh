## taken from ENCODE pipeline https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#
## Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform,  reads mapping to chrM.
## Retain properly paired reads -f 2
## Remove multi-mapped reads (i.e. those with MAPQ < 30, using -q in SAMtools)
## Remove PCR duplicates (using Picard’s MarkDuplicates or FixSeq)

## input: Raw BAM file ${RAW_BAM_FILE}

BAMFILES=($(ls ${ALIGNEDDIR}/*_sorted.bam))
MAPQ_THRESH=30
OUTFOLDER=${ALIGNEDDIR}/encodeQCMetrics
mkdir -p ${OUTFOLDER}

for RAW_BAM_FILE in ${BAMFILES[@]};	
do
  FOLDER=$(dirname ${RAW_BAM_FILE})
  BASENAME=$(basename ${RAW_BAM_FILE})  
  OFPREFIX=${BASENAME%_sorted.bam}
  if [ ! -f ${OFPREFIX}*pbc.qc ]
  then
	# =============================
	# Remove  unmapped, mate unmapped
	# not primary alignment, reads failing platform
	# Remove low MAPQ reads
	# Only keep properly paired reads
	# Obtain name sorted BAM file
	# ==================
	FILT_BAM_PREFIX="${OFPREFIX}.filt.srt"
	FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
	TMP_FILT_BAM_PREFIX="tmp.${FILT_BAM_PREFIX}.nmsrt"
	TMP_FILT_BAM_FILE="${TMP_FILT_BAM_PREFIX}.bam"
		
	samtools view -F 1804 -f 2 -q ${MAPQ_THRESH} -u ${RAW_BAM_FILE} | samtools sort -n /dev/stdin -o ${OUTFOLDER}/${TMP_FILT_BAM_FILE} -O bam        

	# Remove orphan reads (pair was removed)
	# and read pairs mapping to different chromosomes
	# Obtain position sorted BAM
	samtools fixmate -r ${OUTFOLDER}/${TMP_FILT_BAM_FILE} ${OUTFOLDER}/${OFPREFIX}.fixmate.tmp
	samtools view -F 1804 -f 2 -u ${OUTFOLDER}/${OFPREFIX}.fixmate.tmp | samtools sort /dev/stdin -o ${OUTFOLDER}/${FILT_BAM_FILE} -O bam 
	rm ${OUTFOLDER}/${OFPREFIX}.fixmate.tmp
	rm ${OUTFOLDER}/${TMP_FILT_BAM_FILE}
  
	# =============
	# Mark duplicates
	# =============
	TMP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.dupmark.bam"
	DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc"
	  
	java -Xmx4G -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=${OUTFOLDER}/${FILT_BAM_FILE} OUTPUT=${OUTFOLDER}/${TMP_FILT_BAM_FILE} METRICS_FILE=${OUTFOLDER}/${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
	  
	mv ${OUTFOLDER}/${TMP_FILT_BAM_FILE} ${OUTFOLDER}/${FILT_BAM_FILE}
	  
	# ============================
	# Remove duplicates
	# Index final position sorted BAM
	# Create final name sorted BAM
	# ============================
	FINAL_BAM_PREFIX="${OFPREFIX}.filt.srt.nodup"
	FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
	FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai"
	FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file
	FINAL_NMSRT_BAM_PREFIX="${OFPREFIX}.filt.nmsrt.nodup"
	FINAL_NMSRT_BAM_FILE="${FINAL_NMSRT_BAM_PREFIX}.bam" # To be stored

	samtools view -F 1804 -f 2 -b ${OUTFOLDER}/${FILT_BAM_FILE} > ${OUTFOLDER}/${FINAL_BAM_FILE}
	#samtools sort -n ${FINAL_BAM_FILE} > ${FINAL_NMSRT_BAM_FILE}
	#samtools index ${FINAL_NMSRT_BAM_FILE} 
	samtools sort -n --threads 10 ${OUTFOLDER}/${FINAL_BAM_FILE} -O SAM  | SAMstats --sorted_sam_file -  --outf ${OUTFOLDER}/${FINAL_BAM_FILE_MAPSTATS}
 
	# =============================
	# Compute library complexity
	# =============================
	# Sort by name
	# convert to bedPE and obtain fragment coordinates
	# sort by position and strand
	# Obtain unique count statistics
	PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"
	# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
	  
	samtools sort -n ${OUTFOLDER}/${FILT_BAM_FILE} -o ${OUTFOLDER}/${OFPREFIX}.srt.tmp.bam
	bedtools bamtobed -bedpe -i ${OUTFOLDER}/${OFPREFIX}.srt.tmp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${OUTFOLDER}/${PBC_FILE_QC}
	rm ${OUTFOLDER}/${OFPREFIX}.srt.tmp.bam
	rm ${OUTFOLDER}/${FILT_BAM_FILE}
  fi
done
