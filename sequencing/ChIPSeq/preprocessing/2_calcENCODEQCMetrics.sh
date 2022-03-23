## taken from ENCODE pipeline https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#
## Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform,  reads mapping to chrM.
## Retain properly paired reads -f 2
## Remove multi-mapped reads (i.e. those with MAPQ < 30, using -q in SAMtools)
## Remove PCR duplicates (using Picard’s MarkDuplicates or FixSeq)

## input: sampleName

sampleName=$1 

echo "Calculating ENCODE QC metrics on:" ${sampleName}
echo Job started on:
date -u

cd ${ALIGNEDDIR}

mkdir -p ENCODEMetrics

MAPQ_THRESH=30

f=$(find . -name ${sampleName}_sorted.bam)
 
if [ ! -f ${sampleName}*pbc.qc ]
then
	# =============================
	# Remove  unmapped, mate unmapped
	# not primary alignment, reads failing platform
	# Remove low MAPQ reads
	# Only keep properly paired reads
	# Obtain name sorted BAM file
	# ==================
	FILT_BAM_PREFIX="${sampleName}.filt.srt"
	FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
	TMP_FILT_BAM_PREFIX="tmp.${FILT_BAM_PREFIX}.nmsrt"
	TMP_FILT_BAM_FILE="${TMP_FILT_BAM_PREFIX}.bam"
		
	samtools view -F 1804 -f 2 -q ${MAPQ_THRESH} -u ${f} | samtools sort -n /dev/stdin -o ENCODEMetrics/${TMP_FILT_BAM_FILE} -O bam        

	# Remove orphan reads (pair was removed)
	# and read pairs mapping to different chromosomes
	# Obtain position sorted BAM
	samtools fixmate -r ENCODEMetrics/${TMP_FILT_BAM_FILE} ENCODEMetrics/${sampleName}.fixmate.tmp
	samtools view -F 1804 -f 2 -u ENCODEMetrics/${sampleName}.fixmate.tmp | samtools sort /dev/stdin -o ENCODEMetrics/${FILT_BAM_FILE} -O bam 
	rm ENCODEMetrics/${sampleName}.fixmate.tmp
	rm ENCODEMetrics/${TMP_FILT_BAM_FILE}
  
	# =============
	# Mark duplicates
	# =============
	TMP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.dupmark.bam"
	DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc"
	  
	java -Xmx4G -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=ENCODEMetrics/${FILT_BAM_FILE} OUTPUT=ENCODEMetrics/${TMP_FILT_BAM_FILE} METRICS_FILE=ENCODEMetrics/${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
	  
	mv ENCODEMetrics/${TMP_FILT_BAM_FILE} ENCODEMetrics/${FILT_BAM_FILE}
	  
	# ============================
	# Remove duplicates
	# Index final position sorted BAM
	# Create final name sorted BAM
	# ============================
	FINAL_BAM_PREFIX="${sampleName}.filt.srt.nodup"
	FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
	FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai"
	FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file
	FINAL_NMSRT_BAM_PREFIX="${sampleName}.filt.nmsrt.nodup"
	FINAL_NMSRT_BAM_FILE="${FINAL_NMSRT_BAM_PREFIX}.bam" # To be stored

	samtools view -F 1804 -f 2 -b ENCODEMetrics/${FILT_BAM_FILE} > ENCODEMetrics/${FINAL_BAM_FILE}
	samtools sort -n --threads 10 ENCODEMetrics/${FINAL_BAM_FILE} -O SAM  | SAMstats --sorted_sam_file -  --outf ENCODEMetrics/${FINAL_BAM_FILE_MAPSTATS}
 
	# =============================
	# Compute library complexity
	# =============================
	# Sort by name
	# convert to bedPE and obtain fragment coordinates
	# sort by position and strand
	# Obtain unique count statistics
	PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"

	# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
	  
	samtools sort -n ENCODEMetrics/${FILT_BAM_FILE} -o ENCODEMetrics/${sampleName}.srt.tmp.bam
	bedtools bamtobed -bedpe -i ENCODEMetrics/${sampleName}.srt.tmp.bam | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ENCODEMetrics/${PBC_FILE_QC}
	rm ENCODEMetrics/${sampleName}.srt.tmp.bam
	rm ENCODEMetrics/${FILT_BAM_FILE}

	if [[ $? == 0 ]]
		then echo 'Library complexity metrics calculated'
	fi
fi



