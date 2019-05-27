## because our sample is mixed ethnicity this looks funky return to after checking ethnicity.

cd ${DATADIR}/SNPdata/

## check for relatedness with other samples first with KING
$KINGPATH/king -b SCZ2_QCd.bed --kinship

## check for relatedness with other samples second with plink
${PLINK}/plink --bfile SCZ2_QCd --genome --min 0.1 --out ibd
