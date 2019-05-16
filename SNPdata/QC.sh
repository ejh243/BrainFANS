## Written by Eilis & Gemma
## This script takes genotype data outputted by genome studio and performs data QC
## including: removal of duplicate samples, update sample ids, perform sex check
## It cleans up intermediate files as it goes

####### 

## NOTE: Do not store confidential information in this file use the config file

######


cd ${DATADIR}/SNPdata/

## some samples were run twice (two different brain regions) identify these:
$KINGPATH/king -b SCZ2_gr38_binary.bed --duplicate

## calc missingness rate to identify those to exclude
cut -f 1,2 king.con > duplicateSamples.tmp
cut -f 3,4 king.con >> duplicateSamples.tmp
sort duplicateSamples.tmp | uniq > duplicateSamples.txt
rm duplicateSamples.tmp

${PLINK}/plink --bfile SCZ2_gr38_binary --missing --out duplicateSamples

## use python script to identify duplicated with greatest missingness
python ../scripts/SNPdata/ExcludeDuplicates.py king.con duplicateSamples.imiss dupsToExclude.txt

## remove duplicates
${PLINK}/plink --bfile SCZ2_gr38_binary --remove dupsToExclude.txt --make-bed --out SCZ2_gr38_update_1

## update sample ids
${PLINK}/plink --bfile SCZ2_gr38_update_1 --update-ids ID_sentrix_ID.txt --make-bed --out SCZ2_gr38_update_2

## update sex in fam file
${PLINK}/plink --bfile SCZ2_gr38_update_2 --update-sex MRC2_UpdateSex.txt --make-bed --out SCZ2_gr38_update_3 
 

## perform sex check on samples with enough data
## there are a couple of mismatches, and 1 which does not have a sex predicted. 
${PLINK}/plink --bfile SCZ2_gr38_update_3 --mind 0.02 --check-sex --out SexCheck

## exclude sample which does not have a sex predicted
## For time being sex mismatches have been kept in
awk '{if ($4 == 0) print $1,$2 }' SexCheck.sexcheck > NoSex.txt
${PLINK}/plink --bfile SCZ2_gr38_update_3 --remove NoSex.txt --make-bed --out SCZ2_gr38_update_4

## check for runs of homozygosity
awk '{if ($1 >= 1 && $1 <= 22) print $2}' SCZ2_gr38_update_3.bim > autosomalVariants.txt
${PLINK}/plink --bfile SCZ2_gr38_update_3 --extract autosomalVariants.txt --maf 0.01 --hwe 0.00001 --mind 0.02 --geno 0.05 --indep-pairwise 5000 1000 0.2 --out ld.auto
${PLINK}/plink --bfile SCZ2_gr38_update_3 --extract lad.auto.prune.in --het --out roh
## exclude anyone with |Fhet| > 0.2
awk '{if ($6 > 0.2 || $6 < -0.2) print $1,$2}' roh.het > excessHet.txt
${PLINK}/plink --bfile SCZ2_gr38_update_4 --remove excessHet.txt --make-bed --out SCZ2_gr38_update_5

## filter sample and variant missingness, HWE, rare variants
${PLINK}/plink --bfile SCZ2_gr38_update_5 --maf 0.001 --hwe 0.00001 --mind 0.02 --geno 0.05 --make-bed --out SCZ2_QCd

## write list of samples that passed QC for CNV calling
cut -f 1 SCZ2_QCd.fam > CNV/Samples.txt

## clean up intermediate files
rm SCZ2_gr38_update_*.* 



