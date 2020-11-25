## Written by Eilis
## script is based on cutruntools pipeline and utilties
## https://bitbucket.org/qzhudfci/cutruntools/src/master/
## this script performs the two trimmming steps and aligns using bowtie2


cd ${DATADIR}
FQFILES=($(ls *[rR]1*q.gz))
for f in ${FQFILES[@]}
do
	sampleName=${f/_[rR]*}
	if [ ! -f $ALIGNEDDIR/"$sampleName"_aligned_reads.bam ]		
	then	
		## Run trimmomatic
		echo "Trimming file $sampleName ..."
		java -jar ${TRIMMOMATICPATH}trimmomatic-0.39.jar PE -threads 1 -phred33 $DATADIR/"$sampleName"_r1.fq.gz $DATADIR/"$sampleName"_r2.fq.gz $TRIMDIR/"$sampleName"_1.paired.fastq.gz $TRIMDIR/"$sampleName"_1.unpaired.fastq.gz $TRIMDIR/"$sampleName"_2.paired.fastq.gz $TRIMDIR/"$sampleName"_2.unpaired.fastq.gz ILLUMINACLIP:$CUTRUNPATH/adapters/Truseq3.PE.fa:2:15:4:4:true LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:25

		echo "Second stage trimming $sampleName ..."
		# This tool further trims the reads by 6 nt to get rid of the problem of possible adapter run-through. 
		date
		$KSEQPATH/kseq_test $TRIMDIR/"$sampleName"_1.paired.fastq.gz $len $TRIMDIR2/"$sampleName"_1.paired.fastq.gz
		$KSEQPATH/kseq_test $TRIMDIR/"$sampleName"_2.paired.fastq.gz $len $TRIMDIR2/"$sampleName"_2.paired.fastq.gz

		echo "Aligning file $sampleName ..."
		date
		(bowtie2 -p 2 --dovetail --phred33 -x ${REFGENOME}/genome -1 $TRIMDIR2/"$sampleName"_1.paired.fastq.gz -2 $TRIMDIR2/"$sampleName"_2.paired.fastq.gz) 2> $ALIGNEDDIR/"$sampleName".bowtie2 | samtools view -bS - > $ALIGNEDDIR/"$sampleName"_aligned_reads.bam
	fi
done
