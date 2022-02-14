#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=logFiles/ATAC/%u/sexCheck-%A_%a.o
#SBATCH --error=logFiles/ATAC/%u/sexCheck-%A_%a.e
#SBATCH --job-name=sexCheck-%A_%a.e

## call peaks for sex chromosomes

## read counts in sex chromosomes



## merge chr X variants

## create sample map to merge into single dataset
cd ${ALIGNEDDIR}


awk '{print $1,"SNPs/" $1 "_chrX.gvcf"}' ${METADIR}/matchedVCFIDs.txt > SNPs/cohort.sample_map

mkdir tmp 

gatk GenomicsDBImport \
       --genomicsdb-workspace-path SNPs/gatkDB \
       --batch-size 50 \
       -L chrX \
       --sample-name-map SNPs/cohort.sample_map \
       --tmp-dir tmp \
       --reader-threads 5

gatk GenotypeGVCFs \
    -R ${GENOMEFASTA} \
    -V gendb://SNPs/gatkDB \
    -O SNPs/mergedSamples.chrX.vcf	   

$PLINK/plink --vcf SNPs/mergedSamples.chrX.vcf --split-x b37 --make-bed --out SNPs/atacChrX --double-id
$PLINK/plink --bfile SNPs/atacChrX --check-sex --maf 0.05 --out atacChrX

