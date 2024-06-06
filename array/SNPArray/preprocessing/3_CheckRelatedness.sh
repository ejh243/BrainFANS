## because our sample is mixed ethnicity this looks funky return to after checking ethnicity.
## 

## EXECUTION
# sh 3_CheckRelatedness.sh
# where 

## REQUIRES the following variables in config file
# PROCESSDIR, FILEPREFIX, SCRIPTDIR

## REQUIRES the following software
# king, plink, R 

## INPUT
# ${FILEPREFIX}_QCd # binary plink files following prelim QC

## OUTPUT
# QCoutput/${FILEPREFIX}_${2}_QCd_king
# QCoutput/${FILEPREFIX}_${2}_QCd_ibd

cd ${PROCESSDIR}

${PLINK}/plink --bfile ${FILEPREFIX}_QCd --keep $1 --make-bed --out QCoutput/${FILEPREFIX}_${2}_QCd

## check for relatedness with other samples with KING
$KINGPATH/king -b QCoutput/${FILEPREFIX}_${2}_QCd.bed --kinship --prefix QCoutput/${FILEPREFIX}_${2}_QCd_king

## check for relatedness with other samples with plink
${PLINK}/plink --bfile QCoutput/${FILEPREFIX}_${2}_QCd --genome --min 0.1 --out QCoutput/${FILEPREFIX}_${2}_QCd_ibd

Rscript ${SCRIPTDIR}/utilitys/plotKinshipCoeff.r QCoutput/${FILEPREFIX}_${2}_QCd_king.kin0

rm QCoutput/${FILEPREFIX}_${2}_QCd.*


