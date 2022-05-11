## collate output of imputation server, convert vcf files to plink and filter variants on INFO score
## format files for use with Michegan Imputation Server

## EXECUTION
# sh SNPArray/preprocessing/6_combineImputationOutput.sh <imputation output directory>
# where 
# <imputation output directory> is the folder where the output vcf and dose files from imputation are located, with one file per chr
# script needs to be executed from <git repo>/array/
# 

## REQUIRES the following variables in config file
# PROCESSDIR, IMPUTEDIR, FILEPREFIX

## REQUIRES the following software
# plink, vcftools

## INPUT
#  vcf and dose files from imputation

## OUTPUT
# plink files 


cd $1
for i in {1..22}
do
   ## identify variants with rsq > 0.3 to keep
   gzip -cd  chr$i.info.gz | awk '{if($7 > 0.3) print $1,$7}' > chr$i.keep
   vcftools --gzvcf chr$i.dose.vcf.gz --plink --chr $i --out chr$i
   
   ## convert to bim,bed,fam filter variants to MAF > 0.01
   ${PLINK}/plink --ped chr$i.ped --map chr$i.map --extract chr$i.keep --maf 0.01 --make-bed --out chr${i}_rsq0.3 
   rm chr$i.ped
   rm chr$i.map

done

# Merge them into one dataset

for i in {2..22}
do 
    echo "chr${i}_rsq0.3"
done > mergefile.txt

${PLINK}/plink --bfile chr1_rsq0.3 --merge-list mergefile.txt --make-bed --out ${FILEPREFIX}_rsq0.3

rm chr*_rsq0.3.bed
rm chr*_rsq0.3.bim
rm chr*_rsq0.3.fam

# Combine info files into a single file

head -n1 chr1.keep > AllVariantsrsq0.3.info

for i in {1..22}
do
    awk ' NR>1 {print $0}' < chr$i.keep |cat >> AllVariantsrsq0.3.info
done

rm chr*.keep


## update sex in fam file run hwe filter
awk '{split($1,a,"_"); print $1,$2,a[1]"_"a[2],a[3]}' ${FILEPREFIX}_rsq0.3.fam > UpdateIDs.txt
${PLINK}/plink --bfile ${FILEPREFIX}_rsq0.3 --update-ids UpdateIDs.txt --make-bed --out ${FILEPREFIX}_rsq0.3_QCd 

rm ${FILEPREFIX}_rsq0.3.b*
rm ${FILEPREFIX}_rsq0.3.fam


## recalculate PCs
# LD prune
${PLINK}/plink --bfile ${FILEPREFIX}_rsq0.3_QCd --maf 0.05 --indep 50 5 1.5 --out ${FILEPREFIX}_rsq0.3_QCd.ld
${PLINK}/plink --bfile ${FILEPREFIX}_rsq0.3_QCd --extract ${FILEPREFIX}_rsq0.3_QCd.ld.prune.in --make-bed --out ${FILEPREFIX}_rsq0.3_QCd.ld.prune


# use GCTA to calc PCs
${GCTA}/gcta64 --bfile ${FILEPREFIX}_rsq0.3_QCd.ld.prune --make-grm-bin --autosome --out ${FILEPREFIX}_rsq0.3_QCd
${GCTA}/gcta64 --grm ${FILEPREFIX}_rsq0.3_QCd --pca --out ${FILEPREFIX}_rsq0.3_QCd.pca

rm ${FILEPREFIX}_rsq0.3_QCd.ld.prune.b*
rm ${FILEPREFIX}_rsq0.3_QCd.ld.prune.fam


## plot PCs to identify outliers
Rscript ${SCRIPTDIR}/SNPArray/utilitys/plotPCs.r ${FILEPREFIX}_rsq0.3_QCd.pca.eigenvec 3
