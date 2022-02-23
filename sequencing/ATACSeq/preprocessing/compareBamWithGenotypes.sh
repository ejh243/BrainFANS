# prepare bam file for comparision with verifyBamID

## EXECUTION
# sh ./ATACSeq/preprocessing/formatBamForVerifyBamID.sh <bam file>
# script needs to be executed from <git repo>/sequencing/

## REQUIRES the following variables in config file
# ALIGNEDDIR, GENOMEFASTA, KGREF

## REQUIRES the following software
# gatk, samtools, picard

## INPUT
# sorted bam file

## OUTPUT
#

## base recalibrate bam files for use in verifyBamID
cd ${ALIGNEDDIR}/

mkdir -p baseRecalibrate

sampleName=$1
vcfid=$2

if [[ "${vcfid}" != "#N/A" ]]
then 
    echo "processing" ${sampleName} "with" ${vcfid}
    bamfile=${sampleName}_sorted.bam

    re='^[0-9]{5}'
    if [[ $sampleName =~ $re ]] ;
    then
       projectID=${sampleName%%_*}
    else
       projectID="NA"
    fi

    ## mark duplicates only
    java -jar $EBROOTPICARD/picard.jar  MarkDuplicates INPUT=${bamfile} OUTPUT=baseRecalibrate/${sampleName}_dedup.bam METRICS_FILE=baseRecalibrate/${sampleName}_metrics.txt     


    ## add read group
    java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
           I=baseRecalibrate/${sampleName}_dedup.bam \
           O=baseRecalibrate/${sampleName}_dedup_rglabelled.bam \
           RGID=${projectID} \
           RGLB=lib1 \
           RGPL=ILLUMINA \
           RGPU=unit1 \
           RGSM=${sampleName}

    rm baseRecalibrate/${sampleName}_dedup.bam

    ##index
    samtools index baseRecalibrate/${sampleName}_dedup_rglabelled.bam

    # recalibrate bases in bam files
    gatk BaseRecalibrator \
        -R ${GENOMEFASTA} \
        -I baseRecalibrate/${sampleName}_dedup_rglabelled.bam \
        --known-sites ${KGREF}/1000G_omni2.5.hg38.vcf.gz \
        -O baseRecalibrate/${sampleName}_recal_data.table

    gatk ApplyBQSR \
       -R ${GENOMEFASTA} \
       -I baseRecalibrate/${sampleName}_dedup_rglabelled.bam \
       --bqsr-recal-file baseRecalibrate/${sampleName}_recal_data.table \
       -O baseRecalibrate/${sampleName}_baserecal.bam
       
    rm baseRecalibrate/${sampleName}_dedup_rglabelled.bam*
	rm baseRecalibrate/${sampleName}_recal_data.table

    mkdir -p genotypeConcordance

    ${VERIFYBAMID} --vcf ${GENODIR}/ImputationOutput/All/verifyBamID.vcf.gz --bam baseRecalibrate/${sampleName}_baserecal.bam --out genotypeConcordance/${sampleName} --verbose --ignoreRG --smID ${vcfid}

else

    echo "No genotype data available"

fi
 