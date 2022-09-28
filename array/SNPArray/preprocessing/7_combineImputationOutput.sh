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

   # filter on rsq > 0.3 and MAF > 0.01(https://www.biostars.org/p/205856/)
	bcftools view -i 'R2>.3 & MAF>0.01' -Ov chr${i}.dose.hg38.vcf.gz > chr${i}_filt.vcf
	vcftools --vcf chr${i}_filt.vcf --plink --chr chr$i --out chr$i 
	   
   ## convert to bim,bed,fam 
   ${PLINK}/plink --ped chr$i.ped --map chr$i.map --make-bed --out chr${i}_rsq0.3 
   rm chr$i.ped
   rm chr$i.map

done

# Merge vcf files into one

vcf-concat chr*_filt.vcf | gzip -c > allchr_filt_rsq_maf.vcf.gz
rm chr*_filt.vcf

# Merge them into one plink dataset

for i in {2..22}
do 
    echo "chr${i}_rsq0.3"
done > mergefile.txt


${PLINK}/plink --bfile chr1_rsq0.3 --merge-list mergefile.txt --make-bed --out ${FILEPREFIX}_rsq0.3

rm chr*_rsq0.3.bed
rm chr*_rsq0.3.bim
rm chr*_rsq0.3.fam


awk '{split($1,a,"_"); print $1,$2,a[1]"_"a[2],a[3]}' ${FILEPREFIX}_rsq0.3.fam > UpdateIDs.txt
${PLINK}/plink --bfile ${FILEPREFIX}_rsq0.3 --update-ids UpdateIDs.txt --make-bed --out ${FILEPREFIX}_rsq0.3_QCd 

rm ${FILEPREFIX}_rsq0.3.b*
rm ${FILEPREFIX}_rsq0.3.fam

## update snp names

awk '{print $2,$1":"$4}' ${FILEPREFIX}_rsq0.3_QCd.bim > updateVariantIDs.txt
${PLINK}/plink --bfile ${FILEPREFIX}_rsq0.3_QCd --update-name updateVariantIDs.txt --out ${FILEPREFIX}_rsq0.3_QCd_tmp --make-just-bim

mv ${FILEPREFIX}_rsq0.3_QCd_tmp.bim ${FILEPREFIX}_rsq0.3_QCd.bim 


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

## calculate MAF
${PLINK}/plink --bfile ${FILEPREFIX}_rsq0.3_QCd --freq --out ${FILEPREFIX}_rsq0.3_QCd_maf
