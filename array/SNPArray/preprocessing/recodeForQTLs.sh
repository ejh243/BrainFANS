## separate plink files by chr, filter snps and reformat to 0,1,2 for matrixQTL

## EXECUTION
# sh SNPArray/preprocessing/recodeForQTLs.sh <input directory> <output directory>
# where 
# <imputation output directory> is the folder combined imputation files are
# script needs to be executed from <git repo>/array/
# 

## REQUIRES the following variables in config file
# PLINK, FILEPREFIX

## REQUIRES the following software
# plink

## INPUT
#  plink binary files from imputation

## OUTPUT
# plink files 

cd $1
OUTDIR=$2

mkdir -p ${OUTDIR}/

for chr in {1..22}
do
	${PLINK}/plink --bfile ${FILEPREFIX}_rsq0.3_QCd --chr ${chr} --maf 0.05 --hwe 0.00001 --mind 0.02 --geno 0.05  --recode A-transpose --out ${OUTDIR}/${FILEPREFIX}_rsq0.3_${chr} 
	${PLINK}/plink --bfile ${FILEPREFIX}_rsq0.3_QCd --chr ${chr} --maf 0.05 --hwe 0.00001 --mind 0.02 --geno 0.05  --make-just-bim --out ${OUTDIR}/${FILEPREFIX}_rsq0.3_${chr}
done

