## collate output of imputation server, convert vcf files to plink and filter variants on INFO score
## format files for use with Michegan Imputation Server

## EXECUTION
# sh SNPArray/preprocessing/7_liftoverhg38.sh <imputation output directory> <CHR>
# where 
# <imputation output directory> is the folder where the output vcf and dose files from imputation are located, with one file per chr
# <chr> is a numeric chromosome
# script needs to be executed from <git repo>/array/
# 

## REQUIRES the following variables in config file
# REFGENOME, UCSCUTILS

## REQUIRES the following software
# picard, BCFtools

## INPUT
#  vcf and dose files from imputation

## OUTPUT
# vcf files with variants in hg38 positions


#provide directory of output imputation files on command line
HG19DIR=$1
HG38DIR=${HG19DIR}/hg38

## take chromosome from command line
i=$2

mkdir -p ${HG38DIR}


cd ${HG38DIR}

# loop through all chr files with lift over 
echo Chr${i} started on:
date -u
 

echo unzipping chr${i}.dose.vcf.gz...
gunzip -c  ${HG19DIR}/chr${i}.dose.vcf.gz > chr${i}.dose.vcf

echo adding in chr notation...
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' chr${i}.dose.vcf > chr${i}.dose_withChr.vcf
rm chr${i}.dose.vcf # remove original unzipped file to save disk space

echo lifting over chr${i} from hg19 to hg38...

java -Xmx230G -jar $EBROOTPICARD/picard.jar LiftoverVcf \
  I=chr${i}.dose_withChr.vcf \
  O=chr${i}.dose.hg38.vcf \
  CHAIN=${UCSCUTILS}/hg19ToHg18.over.chain.gz \
  REJECT=chr${i}.rejected_variants.vcf \
  R=$REFGENOME \
  MAX_RECORDS_IN_RAM=100000 \
  WMC=TRUE



#check how many lines in rejected variants and output file
outVCF=$(grep -v "^#" chr${i}.dose.hg38.vcf | wc -l)
rejects=$(grep -v "^#" chr${i}.rejected_variants.vcf | wc -l)
inVCF=$(grep -v "^#" chr${i}.dose_withChr.vcf | wc -l)

echo chr${i} original imputation file contained ${inVCF} variants. 
echo lifted over VCF file contains ${outVCF} variants, with ${rejects} variants unmapped

rm chr${i}.dose_withChr.vcf # remove lift over input file to save disk space


# convert liftover output to compressed format and index
bcftools view -Oz -o chr${i}.dose.hg38.compressed.vcf.gz chr${i}.dose.hg38.vcf
bcftools index chr${i}.dose.hg38.compressed.vcf.gz

rm chr${i}.dose.hg38.vcf # remove uncompressed version to save disk space


#remove sites that have been assigned to a different chr
bcftools view chr${i}.dose.hg38.compressed.vcf.gz -t chr${i} -Oz -o chr${i}.dose.hg38.vcf.gz
correctChr=$(zgrep -v "^#" chr${i}.dose.hg38.vcf.gz | wc -l)
echo ${correctChr} variants were lifted over to chr${i} 

#create another file containing the sites that have been removed for future reference.
bcftools view chr${i}.dose.hg38.compressed.vcf.gz -t ^chr${i} -Oz -o chr${i}.dose_hg38_otherChrs.vcf.gz
otherChr=$(zgrep -v "^#" chr${i}.dose_hg38_otherChrs.vcf.gz | wc -l)
echo ${otherChr} variants were lifted over to a different chr

rm chr${i}.dose.hg38.compressed.vcf.gz # remove compressed total version to save disk space

echo Liftover of chr${i} complete! Job ended on:
date -u

