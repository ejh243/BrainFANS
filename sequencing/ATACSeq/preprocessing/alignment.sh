## ===================================================================================================================##
##                               ATAC-seq pipeline STEP 1.3: Pre-analysis -- alignment                                ##
## ===================================================================================================================##
## EXECUTION: sbatch --array= ./sequencing/ATACSeq/preprocessing/alignment.sh <sampleName>                            ||
## - execute from scripts directory                                                                                   ||
##                                                                                                                    ||
## DESCRIPTION: This scripts performs alignment of reads against the reference genome and subsequent processing of    ||  
## aligned reads, such as reads with quality < 30, unmapped, mate unmapped, secondary alignment, reads failing        || 
## platform, remove duplicates and reads mapped to a different chromosome                                             ||
##                                                                                                                    ||
## REQUIRES:                                                                                                          ||
## - File in ${METADIR}/samples.txt that lists sample names.                                                          ||
## - Variables in config file: RAWDATADIR, TRIMDIR, ALIGNEDDIR, REFGENOME, multimap                                   ||
## - Software: bowtie2, samtools, picard                                                                              ||
##                                                                                                                    ||
## INPUTS:                                                                                                            || 
## $1 -> <sampleName> Name of sample specified in command line. Reads should be previously trimmed                    ||
##                                                                                                                    ||
## OUTPUTS:                                                                                                           || 
##  *_filt.nodup.stats, *_statsperchr.txt, *_noMT.stats, *_sorted.bam, *_sorted.bam.bai, *_noMT.bam, *_.filt.nodup.bam||	
##  *_.filt.nodup.bam.bai                                                                                             || 
##                                                                                                                    ||
## ===================================================================================================================##

## ============ ##
##    SET UP    ##
## ============ ##

## Sample to align
sampleName=$1
echo "Processing" ${sampleName}

CHR=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

## Check if trimmed files for alignment exist
cd ${TRIMDIR}
f=($(ls ${sampleName}*trimmed*.f*))

echo "Found trimmed files:" ${f[0]} ${f[1]}

## =============== ##
##    ALIGNMENT    ##
## =============== ##

## checks if aligned file for input sample already exists
if [ ! -s ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam ]	
then
	echo "Running alignment for"" ${sampleName}"
	date -u	
	## alignment
	bowtie2  -p 10 -X 2000 -k ${multimap} -x ${REFGENOME}/genome -1 ${f[0]} -2 ${f[1]} -S ${ALIGNEDDIR}/${sampleName}.sam &> ${ALIGNEDDIR}/${sampleName}.bowtie.log

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
	
   # Mark duplicates
   java -Xmx4G -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=${ALIGNEDDIR}/${sampleName}.filt.bam OUTPUT=${ALIGNEDDIR}/${sampleName}.filt.dupmark.bam METRICS_FILE=${ALIGNEDDIR}/${sampleName}_dupMetrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
   
   rm ${ALIGNEDDIR}/${sampleName}.filt.bam

   # Remove duplicates
   samtools view -F 1804 -f 2 -b ${ALIGNEDDIR}/${sampleName}.filt.dupmark.bam > ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam
   
   # Index Final BAM file
   samtools index ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam
   
   ## get sum stats post to filtering
   samtools stats ${ALIGNEDDIR}/${sampleName}.filt.nodup.bam	> ${ALIGNEDDIR}/QCOutput/${sampleName}.filt.nodup.stats
else
	{ echo "Aligned file found so not aligning"; exit 1;}
fi
