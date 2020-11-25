
cd ${ALIGNEDDIR}

## check if required output directories exist, if not create
for d in sorted dup.marked dedup dup.marked.120bp dedup.120bp; do
if [ ! -d $d ]; then
	mkdir $d
fi
done

BAMFILES=($(ls *_aligned_reads.bam))
for f in ${BAMFILES[@]}
do
	sampleName=${f/_aligned_reads.bam}
		
	echo "Filtering unmapped fragments... ""$sampleName".bam
	date
	samtools view -bh -f 3 -F 4 -F 8 ${f} > sorted/"$sampleName".step1.bam

	echo "Sorting BAM... ""$sampleName".bam
	date
	java -jar $EBROOTPICARD/picard.jar SortSam \
	INPUT=sorted/"$sampleName".step1.bam OUTPUT=sorted/"$sampleName".bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT
	rm -rf sorted/"$sampleName".step1.bam

	echo "Marking duplicates... ""$sampleName".bam
	date
	java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
	INPUT=sorted/"$sampleName".bam OUTPUT=dup.marked/"$sampleName".bam VALIDATION_STRINGENCY=SILENT \
	METRICS_FILE=metrics."$sampleName".txt

	echo "Removing duplicates... ""$sampleName".bam
	date
	samtools view -bh -F 1024 dup.marked/"$sampleName".bam > dedup/"$sampleName".bam

	echo "Filtering to <120bp... ""$sampleName".bam
	date
	samtools view -h dup.marked/"$sampleName".bam |LC_ALL=C awk -f $CUTRUNPATH/filter_below.awk |samtools view -Sb - > dup.marked.120bp/"$sampleName".bam
	samtools view -h dedup/"$sampleName".bam |LC_ALL=C awk -f $CUTRUNPATH/filter_below.awk |samtools view -Sb - > dedup.120bp/"$sampleName".bam

	echo "Creating bam index files... ""$sampleName".bam
	date
	samtools index sorted/"$sampleName".bam
	samtools index dup.marked/"$sampleName".bam
	samtools index dedup/"$sampleName".bam
	samtools index dup.marked.120bp/"$sampleName".bam
	samtools index dedup.120bp/"$sampleName".bam
done
