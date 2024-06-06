## format files for use with Michegan Imputation Server

## EXECUTION
# sh 4_formatForimputation.sh <population> <SNP ref file>
# where 
# <population > is 3 letter code for super population state ALL for no subsetting by population
# <SNP ref file> is an input file of 


## REQUIRES the following variables in config file
# PROCESSDIR, IMPUTEDIR, FILEPREFIX

## REQUIRES the following software
# plink, perl,

## INPUT
#  # binary plink files following prelim QC

## OUTPUT
# vcf files split by chr for upload to michegan imputation server

population=$1
refFile=$2

cd ${IMPUTEDIR}/

mkdir -p ImputationInput

cd ImputationInput

mkdir -p ${population}



## use tool to check data prior to upload https://www.well.ox.ac.uk/~wrayner/tools/
## for All use 1000G
cd ${population}/

## subset samples
if [ $population != "ALL" ]
then
    ${PLINK}/plink --bfile ${PROCESSDIR}/${FILEPREFIX}_QCd --keep ${PROCESSDIR}/${population}Samples.txt --maf 0.05 --out ${FILEPREFIX}_QCd_${population} --make-bed
else
	cp ${PROCESSDIR}/${FILEPREFIX}_QCd.bim ${FILEPREFIX}_QCd_${population}.bim
	cp ${PROCESSDIR}/${FILEPREFIX}_QCd.bed ${FILEPREFIX}_QCd_${population}.bed
	cp ${PROCESSDIR}/${FILEPREFIX}_QCd.fam ${FILEPREFIX}_QCd_${population}.fam
fi

## liftover to hg19 for imputation
${PLINK}/plink --bfile ${FILEPREFIX}_QCd_${population} --update-map ${GSAREF}liftoverhg19.txt 3 --make-bed --out ${FILEPREFIX}_QCd_hg19

## for HRC check tool need freq file
${PLINK}/plink --bfile ${FILEPREFIX}_QCd_hg19 --freq --out ${FILEPREFIX}_QCd_hg19_freq

if [[ $(basename ${refFile}) == HRC* ]] ;
then
perl ${KINGPATH}/HRC-1000G-check-bim.pl -b ${FILEPREFIX}_QCd_hg19.bim -f ${FILEPREFIX}_QCd_hg19_freq.frq -r ${refFile} -g --hrc
else 
perl ${KINGPATH}/HRC-1000G-check-bim.pl -b ${FILEPREFIX}_QCd_hg19.bim -f ${FILEPREFIX}_QCd_hg19_freq.frq -r ${refFile} -g --1000g
fi


sed -i 's=plink=${PLINK}/plink=g' Run-plink.sh
sh Run-plink.sh
#for file in *.vcf; do awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' ${file} > with_chr_${file}; vcf-sort with_chr_${file} | bgzip -c > ${file}.gz;done
for file in *.vcf; do vcf-sort ${file} | bgzip -c > ${file}.gz;done
rm *.vcf
rm ${FILEPREFIX}_QCd*.*[^gz]

