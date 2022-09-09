## filter vcf files to high quality snps only for use in verifyBamID

module load SAMtools

cd $1
for i in {1..22}
do
   ## identify variants with rsq > 0.8 to keep
   gzip -cd  chr$i.info.gz | awk '{if($7 > 0.8) print $1,$7}' > chr$i.keep
   vcftools --gzvcf chr$i.dose.vcf.gz --chr $i --remove-indels --snps chr$i.keep --recode --recode-INFO-all --out chr$i_SNPs

done

vcf-concat chr*.vcf | bgzip -c > verifyBamID.vcf.gz

rm chr$i_SNPs*.vcf*
rm chr$i.keep

## get list of sample IDs
vcf-query -l chr$i.dose.vcf.gz > vcfSampleIDs.txt

