## Performs peak calling for sex chr using shifted tagAlign files from all samples

## EXECUTION
# sh ./ATACSeq/preprocessing/sexChrPeaks.sh
# where 
# <sampleName> is the file prefix
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, PEAKDIR

## REQUIRES the following software
# bedtools, macs2, htseq-count

## INPUT
# shifted tag align file
#

## OUTPUT
# ${PEAKDIR}/MACS/ShiftedTagAlign/*.broadPeak (and other macs output)
# ${PEAKDIR}/MACS/BAMPE/*.broadPeak
# ${PEAKDIR}/MACS/ShiftedTagAlign/${sampleName}.broadPeak.filt
# ${PEAKDIR}/MACS/BAMPE/${sampleName}.broadPeak.filt

mkdir -p ${PEAKDIR}/MACS/ShiftedTagAlign/sexChr


cd ${ALIGNEDDIR}

for chr in chrX chrY
do
    files=($(ls sexChr/*${chr}.tn5.tagAlign.gz))

    macs2 callpeak -t ${files[@]} --outdir ${PEAKDIR}/MACS/ShiftedTagAlign/sexChr -n ${chr} -f BED -g 2.9e9 -q 1e- --keep-dup all --shift 100 --extsize 200 --nomodel --broad
done

## X chr extract peaks nearest XIST and FIRRE

## Y chr exclude peaks overlapping psuedoautosomal regions
	bedtools intersect -v -a ${PEAKDIR}/MACS/ShiftedTagAlign/sexChr/${chr}_peaks.broadPeak -b ${PAR} \
	  | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
	  | grep -P 'chr[\dXY]+[ \t]' > ${PEAKDIR}/MACS/ShiftedTagAlign/sexChr/${chr}.broadPeak.filt

rm ${PEAKDIR}/MACS/ShiftedTagAlign/sexChr/${chr}_peaks.broadPeak


mkdir -p ${PEAKDIR}/MACS/ShiftedTagAlign/sexChr/GFF
cd ${PEAKDIR}/MACS/ShiftedTagAlign/sexChr

## create gff files
awk -v OFS='\t' '{print $1,"MACS2","BroadPeak",$2,$3,$5,$6,".","ID="$4";"}' ${PEAKDIR}/MACS/ShiftedTagAlign/sexChr/${chr}.broadPeak.filt > GFF/${chr}.gff

## count reads in peaks
htseq-count -f bam -t BroadPeak -i ID --stranded=no ${BAMFILES[@]} GFF/${chr}.gff > ${PEAKCOUNTS}/${chrX}

  