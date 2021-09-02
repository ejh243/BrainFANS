## peak calling with epic2 designed to call broader peaks
## run at relaxed q-value threshold to enable classification of null regions in the genome

cd ${ALIGNEDDIR}

mkdir -p ${PEAKDIR}/EPIC2

sampleSets=($(ls ../sampleLists/PeakCallingInputFiles*))

for type_input in ${sampleSets[@]}
do
	type_control=${type_input/Input/Control}
	type=$(basename $type_input)
	type=${type%.txt}
	type=${type##PeakCallingInputFiles}
	
	input_samples=$(cat ${type_input})
	control_samples=$(cat ${type_control})
	
	epic2 -t ${input_samples[@]} -c ${control_samples[@]} --genome hg38 -fdr 0.5 -o ${PEAKDIR}/EPIC2/${type}.txt

	
done

## run peaks across all smaples

input_samples=($(ls *PC*_L00.bml.GRCh38.karyo.deduplicated.bam))
control_samples=($(ls *IC*_L00.bml.GRCh38.karyo.deduplicated.bam))

epic2 -t ${input_samples[@]} -c ${control_samples[@]} --genome hg38 -fdr 0.5 -o ${PEAKDIR}/EPIC2/ALL.txt