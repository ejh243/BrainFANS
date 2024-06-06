## This script takes genotype data outputted by genome studio and performs data QC
## including: removal of duplicate samples, update sample ids, perform sex check

## EXECUTION
# sh 1_QC.sh 

## REQUIRES the following variables to be loaded config file
# RAWDATADIR, FILEPREFIX, METADIR, CNVDIR

## REQUIRES the following software
# king, plink, python, gcta

## INPUT
# ${METADIR}/UpdateIDs.txt # txt file mapping sentrix ids to sample ids
# ${METADIR}/UpdateSex.txt # txt file mapping with sexes for sex check

## OUTPUT
# QCoutput/ # summary of samples failung qc steps and related metrics
# ${FILEPREFIX}_QCd plink files # QC'd SNP data
# GCTA/ # gcta grm and genotype pcas 



cd ${PROCESSDIR}

mkdir -p QCoutput

## some samples were run twice (two different brain regions) identify these:
$KINGPATH/king -b ${RAWDATADIR}/${FILEPREFIX}.bed --duplicate --prefix QCoutput/king


cut -f 1,2 QCoutput/king.con > QCoutput/duplicateSamples.tmp
cut -f 3,4 QCoutput/king.con >> QCoutput/duplicateSamples.tmp
sort QCoutput/duplicateSamples.tmp | uniq > QCoutput/duplicateSamples.txt
rm QCoutput/duplicateSamples.tmp

##if any duplicate samples calc missingness rate to identify those to exclude
if [ -s QCoutput/duplicateSamples.txt ]
then 
    ${PLINK}/plink --bfile ${RAWDATADIR}/${FILEPREFIX} --missing --out duplicateSamples

	## use python script to identify duplicated with greatest missingness
	python ${SCRIPTDIR}/utilitys/ExcludeDuplicates.py QCoutput/king.con duplicateSamples.imiss QCoutput/dupsToExclude.txt

	## remove duplicates
	${PLINK}/plink --bfile ${RAWDATADIR}/${FILEPREFIX} --remove QCoutput/dupsToExclude.txt --make-bed --out ${FILEPREFIX}_update_1
else
    cp ${RAWDATADIR}/${FILEPREFIX}.bed ${FILEPREFIX}_update_1.bed
    cp ${RAWDATADIR}/${FILEPREFIX}.bim ${FILEPREFIX}_update_1.bim
    cp ${RAWDATADIR}/${FILEPREFIX}.fam ${FILEPREFIX}_update_1.fam
fi	

	
## update sample ids
${PLINK}/plink --bfile ${FILEPREFIX}_update_1 --update-ids ${METADIR}/UpdateIDs.txt --make-bed --out ${FILEPREFIX}_update_2

## remove variants at the same position (i.e. triallelic)
awk '{if ($1 != 0) print $1":"$4}' ${FILEPREFIX}_update_2.bim > pos.tmp
sort pos.tmp | uniq -d > dupLocs.txt
#awk --delimiter=":" '{print $1,$2}' dupLocs.txt
awk -F ":" '{print $1,$2-1,$2,"set1", "set2"}' dupLocs.txt > positionsExclude.txt

${PLINK}/plink --bfile ${FILEPREFIX}_update_2 --exclude range positionsExclude.txt --make-bed --out ${FILEPREFIX}_update_3;
rm pos.tmp
rm dupLocs.txt

## update sex in fam file
${PLINK}/plink --bfile ${FILEPREFIX}_update_3 --update-sex ${METADIR}/UpdateSex.txt --make-bed --out ${FILEPREFIX}_update_4 

## perform sex check on samples with enough data
${PLINK}/plink --bfile ${FILEPREFIX}_update_4 --mind 0.02 --check-sex --out QCoutput/SexCheck

## exclude samples which do not have sex predicted
## exclude mismatched samples
## retain samples with missing sex info
awk '{if ($4 == 0) print $1,$2 }' QCoutput/SexCheck.sexcheck > QCoutput/sexErrors.txt
awk '{if ($4 != $3 && $3 != 0) print $1,$2 }' QCoutput/SexCheck.sexcheck >> QCoutput/sexErrors.txt
${PLINK}/plink --bfile ${FILEPREFIX}_update_4 --remove QCoutput/sexErrors.txt --make-bed --out ${FILEPREFIX}_update_5

## check for runs of homozygosity
awk '{if ($1 >= 1 && $1 <= 22) print $2}' ${FILEPREFIX}_update_5.bim > autosomalVariants.txt
${PLINK}/plink --bfile ${FILEPREFIX}_update_5 --extract autosomalVariants.txt --maf 0.01 --hwe 0.00001 --mind 0.02 --geno 0.05 --indep-pairwise 5000 1000 0.2 --out ld.auto
${PLINK}/plink --bfile ${FILEPREFIX}_update_5 --extract ld.auto.prune.in --het --out QCoutput/roh
## exclude anyone with |Fhet| > 0.2
awk '{if ($6 > 0.2 || $6 < -0.2) print $1,$2}' QCoutput/roh.het > QCoutput/excessHet.txt
${PLINK}/plink --bfile ${FILEPREFIX}_update_5 --remove QCoutput/excessHet.txt --make-bed --out ${FILEPREFIX}_update_6
rm autosomalVariants.txt


## filter sample and variant missingness, HWE, rare variants and exclude variants with no position
awk '{if ($1 == 0) print $2}' ${FILEPREFIX}_update_6.bim > noLocPos.tmp
${PLINK}/plink --bfile ${FILEPREFIX}_update_6 --exclude noLocPos.tmp --maf 0.001 --hwe 0.00001 --mind 0.02 --geno 0.05 --make-bed --out ${FILEPREFIX}_QCd


## write list of samples that passed QC for CNV calling
cut -f 1,2 --delimiter=" " ${FILEPREFIX}.fam > ${CNVDIR}/ID_Map.txt
cut -f 2 --delimiter=" " ${FILEPREFIX}_QCd.fam > ${CNVDIR}/Samples.txt

## clean up intermediate files but keep log files
rm ${FILEPREFIX}_update_*.b*
rm ${FILEPREFIX}_update_*.fam
 

## calc PCS within sample only
# LD prune
${PLINK}/plink --bfile ${FILEPREFIX}_QCd --indep 50 5 1.5 --out ${FILEPREFIX}_QCd.ld
${PLINK}/plink --bfile ${FILEPREFIX}_QCd --extract ${FILEPREFIX}_QCd.ld.prune.in --make-bed --out ${FILEPREFIX}_QCd.ld.prune

mkdir -p GCTA

${GCTA}/gcta64 --bfile ${FILEPREFIX}_QCd.ld.prune --make-grm-bin --autosome --out GCTA/${FILEPREFIX}_QCd_GCTA
${GCTA}/gcta64 --grm GCTA/${FILEPREFIX}_QCd_GCTA --pca --out GCTA/${FILEPREFIX}_QCd.pca

rm ${FILEPREFIX}_QCd.ld.prune*

## extract SNP probes for comparison with DNAm data
${PLINK}/plink --bfile ${FILEPREFIX}_QCd --extract ${EPICREF}/RSprobes.txt --recodeA --out ${FILEPREFIX}_59DNAmSNPs

