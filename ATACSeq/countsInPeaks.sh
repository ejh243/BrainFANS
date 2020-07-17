 


#mkdir -p ${PEAKDIR}/GFF
#mkdir -p ${PEAKDIR}/CountsInPeaks

cd ${PEAKDIR}
PEAKFILES=($(ls *.broadPeak))
for peakSet in ${PEAKFILES[@]};
do
	gffFile=${peakSet/broadPeak/gff}
	## create gff files
	awk -v OFS='\t' '{print $1,"MACS2","BroadPeak",$2,$3,$5,$6,".","ID="$4";"}' ${peakSet} > GFF/${gffFile}
	## count reads in peaks
	htseq-count -f bam -t BroadPeak -i ID --stranded=no ${BAMFILES[@]} GFF/${gffFile} > CountsInPeaks/${gffFile}
done

