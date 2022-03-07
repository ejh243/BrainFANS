## peak calling with epic2 designed to call broader peaks
## run at relaxed q-value threshold to enable classification of null regions in the genome


type_input=$1


type_control=${type_input/Input/Control}
type=$(basename $type_input)
type=${type%.txt}
type=${type##PeakCallingInputFiles}
	
input_samples=$(cat ${type_input})
control_samples=$(cat ${type_control})

cd ${ALIGNEDDIR}
epic2 -t ${input_samples[@]} -c ${control_samples[@]} --genome hg38 -fdr 0.5 -o ${PEAKDIR}/EPIC2/${type}.txt

## filter out blacklist peaks

bedtools intersect -a ${PEAKDIR}/EPIC2/${type}.txt -b ${REFDIR}/ENCODE/Blacklist/hg38-blacklist.v2.bed -v > ${PEAKDIR}/EPIC2/${type}_peaks.broadPeak.filt