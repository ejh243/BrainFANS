


sampleName=$1
sampleName=${sampleName%.selfSM}
vcfid=$2
vcfid="${vcfid%"${vcfid##*[![:space:]]}"}" 

cd ${ALIGNEDDIR}/

echo "processing" ${sampleName} "with" ${vcfid}

${VERIFYBAMID} --vcf ${GENODIR}/hg38/allchr_filt_rsq_maf.vcf.gz --bam baseRecalibrate/${sampleName}_baserecal.bam --out genotypeConcordance/${sampleName} --ignoreRG --smID ${vcfid} --best

