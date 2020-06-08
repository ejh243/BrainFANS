## Written by Eilis


BAMFILES=($(ls ${ALIGNEDDIR}/*_depDup_q30.bam))
mkdir -p ${PEAKDIR}

echo "Number of bam files found for peak calling:" "${#BAMFILES[@]}"

for f in ${BAMFILES[@]};
do
  fileName=$(basename ${f})
  sampleName=${fileName/_*}
  macs2 callpeak -t ${f} --outdir ${PEAKDIR} -n ${sampleName} -g 2.9e9  -B --nomodel --shift 0
done

