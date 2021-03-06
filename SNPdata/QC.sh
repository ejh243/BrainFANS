## Written by Eilis & Gemma
## This script takes genotype data outputted by genome studio and performs data QC
## including: removal of duplicate samples, update sample ids, perform sex check
## It cleans up intermediate files as it goes

####### 

## NOTE: Do not store confidential information in this file use the config file

######


cd ${DATADIR}/SNPdata/Merged

## some samples were run twice (two different brain regions) identify these:
$KINGPATH/king -b SCZ_gr38_binary.bed --duplicate

## calc missingness rate to identify those to exclude
cut -f 1,2 king.con > duplicateSamples.tmp
cut -f 3,4 king.con >> duplicateSamples.tmp
sort duplicateSamples.tmp | uniq > duplicateSamples.txt
rm duplicateSamples.tmp

${PLINK}/plink --bfile SCZ_gr38_binary --missing --out duplicateSamples

## use python script to identify duplicated with greatest missingness
python ../../scripts/SNPdata/ExcludeDuplicates.py king.con duplicateSamples.imiss dupsToExclude.txt

## remove duplicates
${PLINK}/plink --bfile SCZ_gr38_binary --remove dupsToExclude.txt --make-bed --out SCZ2_gr38_update_1

## update sample ids
${PLINK}/plink --bfile SCZ2_gr38_update_1 --update-ids UpdateIDs.txt --make-bed --out SCZ2_gr38_update_2

## remove variants at the same position (i.e. triallelic)
awk '{if ($1 != 0) print $1":"$4}' SCZ2_gr38_update_2.bim > pos.tmp
sort pos.tmp | uniq -d > dupLocs.txt
#awk --delimiter=":" '{print $1,$2}' dupLocs.txt
awk -F ":" '{print $1,$2-1,$2,"set1", "set2"}' dupLocs.txt > positionsExclude.txt

${PLINK}/plink --bfile SCZ2_gr38_update_2 --exclude range positionsExclude.txt --make-bed --out SCZ2_gr38_update_3;
rm pos.tmp
rm dupLocs.txt

## update sex in fam file
${PLINK}/plink --bfile SCZ2_gr38_update_3 --update-sex UpdateSex.txt --make-bed --out SCZ2_gr38_update_4 

## perform sex check on samples with enough data
## there are a couple of mismatches, and 1 which does not have a sex predicted. 
${PLINK}/plink --bfile SCZ2_gr38_update_4 --mind 0.02 --check-sex --out SexCheck

## exclude sample which does not have a sex predicted
## exclude mismatched samples
## retain samples with missing sex info
awk '{if ($4 == 0) print $1,$2 }' SexCheck.sexcheck > sexErrors.txt
awk '{if ($4 != $3 && $3 != 0) print $1,$2 }' SexCheck.sexcheck >> sexErrors.txt
${PLINK}/plink --bfile SCZ2_gr38_update_4 --remove sexErrors.txt --make-bed --out SCZ2_gr38_update_5

## check for runs of homozygosity
awk '{if ($1 >= 1 && $1 <= 22) print $2}' SCZ2_gr38_update_5.bim > autosomalVariants.txt
${PLINK}/plink --bfile SCZ2_gr38_update_5 --extract autosomalVariants.txt --maf 0.01 --hwe 0.00001 --mind 0.02 --geno 0.05 --indep-pairwise 5000 1000 0.2 --out ld.auto
${PLINK}/plink --bfile SCZ2_gr38_update_5 --extract ld.auto.prune.in --het --out roh
## exclude anyone with |Fhet| > 0.2
awk '{if ($6 > 0.2 || $6 < -0.2) print $1,$2}' roh.het > excessHet.txt
${PLINK}/plink --bfile SCZ2_gr38_update_5 --remove excessHet.txt --make-bed --out SCZ2_gr38_update_6
rm autosomalVariants.txt


## filter sample and variant missingness, HWE, rare variants and exclude variants with no position
awk '{if ($1 == 0) print $2}' SCZ2_gr38_update_6.bim > noLocPos.tmp
${PLINK}/plink --bfile SCZ2_gr38_update_6 --exclude noLocPos.tmp --maf 0.001 --hwe 0.00001 --mind 0.02 --geno 0.05 --make-bed --out SCZ2_QCd


## write list of samples that passed QC for CNV calling

cut -f 1,2 --delimiter=" " SCZ_gr38_binary.fam > ${DATADIR}/SNPdata/CNV/ID_Map.txt
cut -f 2 --delimiter=" " SCZ2_QCd.fam > ${DATADIR}/SNPdata/CNV/Samples.txt

## clean up intermediate files but keep log files
rm SCZ2_gr38_update_*.b*
rm SCZ2_gr38_update_*.fam
 

## calc PCS within sample only
# LD prune
${PLINK}/plink --bfile SCZ2_QCd --indep 50 5 1.5 --out SCZ2_QCd.ld
${PLINK}/plink --bfile SCZ2_QCd --extract SCZ2_QCd.ld.prune.in --make-bed --out SCZ2_QCd.ld.prune

${GCTA}/gcta64 --bfile SCZ2_QCd.ld.prune --make-grm-bin --autosome --out SCZ2_QCd_GCTA
${GCTA}/gcta64 --grm SCZ2_QCd_GCTA --pca --out SCZ2_QCd.pca

rm SCZ2_QCd.ld.prune*

## extract SNP probes for comparision with DNAm data
${PLINK}/plink --bfile SCZ2_QCd --extract ../../References/EPICArray/RSprobes.txt --recodeA --out SCZ_59DNAmSNPs

