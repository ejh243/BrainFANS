## Calls peaks on sex chromosomes and performs quantification of these peaks
## filters peaks to exclude pseudoautosomal regions


## EXECUTION
# sh ./ATACSeq/preprocessing/sexChrPeaks.sh
# where 

# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR, PEAKCOUNTDIR, PAR

## REQUIRES the following software
# macs2, epic2, featureCounts, bedtools

## INPUT
# indexed, bam files filtered to sex chromosomes

## OUTPUT
# ${PEAKDIR}/sexChr/MACS/All_peaks.narrowPeak.filt 
# ${PEAKDIR}/sexChr/MACS/All_peaks.narrowPeak.filt.saf
# ${PEAKDIR}/sexChr/EPIC2/All_peaks.broadPeak.filt 
# ${PEAKDIR}/sexChr/EPIC2/All_peaks.broadPeak.filt.saf
# ${PEAKCOUNTDIR}/sexChr/macs2.narrowPeak
# ${PEAKCOUNTDIR}/sexChr/epic2.broadPeak

cd ${ALIGNEDDIR}/sexChr/


PULLDOWN=$(ls *PC*.GRCh38.karyo.deduplicated.bam)
CONTROL=$(ls *IC*.GRCh38.karyo.deduplicated.bam)
mkdir -p ${PEAKDIR}/sexChr/MACS/
	
macs2 callpeak -t ${PULLDOWN} -c ${CONTROL} -f BAM -g hs -n All -q 0.01 --outdir ${PEAKDIR}/sexChr/MACS/


mkdir -p ${PEAKDIR}/sexChr/EPIC2/

epic2 -t ${PULLDOWN} -c ${CONTROL} --genome hg38 -fdr 0.01 -o ${PEAKDIR}/sexChr/EPIC2/All.txt

## exclude psuedo autosomal regions
bedtools intersect -a ${PEAKDIR}/sexChr/MACS/All_peaks.narrowPeak -b ${PAR} -v > ${PEAKDIR}/sexChr/MACS/All_peaks.narrowPeak.filt
bedtools intersect -a ${PEAKDIR}/sexChr/EPIC2/All.txt -b ${PAR} -v > ${PEAKDIR}/sexChr/EPIC2/All_peaks.broadPeak.filt

rm ${PEAKDIR}/sexChr/MACS/All_peaks.narrowPeak
rm ${PEAKDIR}/sexChr/MACS/All_peaks.xls
rm ${PEAKDIR}/sexChr/EPIC2/All.txt

## count reads in peaks
## reformat peaks for featureCounts input
awk 'BEGIN { OFS="\t";  print "GeneID\tChr\tStart\tEnd\tStrand"} {print $4,$1,$2,$3,"."}' ${PEAKDIR}/sexChr/MACS/All_peaks.narrowPeak.filt > ${PEAKDIR}/sexChr/MACS/All_peaks.narrowPeak.filt.saf

awk 'BEGIN { OFS="\t";  print "GeneID\tChr\tStart\tEnd\tStrand" } {print "EPIC2Peak_"NR,$1,$2,$3,"."}' ${PEAKDIR}/sexChr/EPIC2/All_peaks.broadPeak.filt > ${PEAKDIR}/sexChr/EPIC2/All_peaks.broadPeak.filt.saf

mkdir -p ${PEAKCOUNTDIR}/sexChr
featureCounts -F SAF -T 8 -p -B -a ${PEAKDIR}/sexChr/MACS/All_peaks.narrowPeak.filt.saf -o ${PEAKCOUNTDIR}/sexChr/macs2.narrowPeak ${PULLDOWN} ${CONTROL}

featureCounts -F SAF -T 8 -p -B -a ${PEAKDIR}/sexChr/EPIC2/All_peaks.broadPeak.filt.saf -o ${PEAKCOUNTDIR}/sexChr/epic2.broadPeak ${PULLDOWN} ${CONTROL}



