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
# shifted tag align files
#

## OUTPUT
# ${PEAKDIR}/MACS/ShiftedTagAlign/sexChr/GFF/chr*.gff
# ${PEAKDIR}/MACS/ShiftedTagAlign/sexChr/chr*.broadPeak.filt


mkdir -p ${PEAKDIR}/MACS/ShiftedTagAlign/sexChr
mkdir -p ${PEAKCOUNTS}/MACS/ShiftedTagAlign/sexChr/

cd ${PEAKDIR}/MACS/ShiftedTagAlign/sexChr

for chr in chrX chrY
do
	echo "Processing" ${chr}

    files=($(ls ${ALIGNEDDIR}/sexChr/*${chr}.tn5.tagAlign.gz))
	## only save really significant peaks
    macs2 callpeak -t ${files[@]} -n ${chr} -f BED -g 2.9e9 -q 1e-4 --keep-dup all --shift 100 --extsize 200 --nomodel --broad --broad-cutoff 1e-4

	if [ ${chr} == chrX ]
	then 
		## X chr extract peaks nearest XIST and FIRRE
		bedtools closest -D a -id -io -a ${XCHRBED} -b chrX_peaks.broadPeak > chrX.broadPeak.filt		
	else
		## Y chr exclude peaks overlapping psuedoautosomal regions
		bedtools intersect -v -a chrY_peaks.broadPeak -b ${PAR} \
		  | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
		 | grep -P 'chr[\dXY]+[ \t]' > chrY.broadPeak.filt
	fi
	
    rm *_peaks.*
	
	## count reads in peaks
	touch ${PEAKCOUNTS}/MACS/ShiftedTagAlign/sexChr/${chr}.peakcounts.txt
	for each in ${files[@]}
	do
		bedtools intersect -C -filenames -a ${chr}.broadPeak.filt -b ${each} | awk -v var=${each} '{ print $0, var}' >> ${PEAKCOUNTS}/MACS/ShiftedTagAlign/sexChr/${chr}.peakcounts.txt
    done
done

  