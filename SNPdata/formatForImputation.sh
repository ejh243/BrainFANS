## format files for use with Michegan Imputation Server
cd ${DATADIR}/SNPdata/Merged

mkdir -p ImputationInput

cd ImputationInput

mkdir -p All
mkdir -p EUR


## use tool to check data prior to upload https://www.well.ox.ac.uk/~wrayner/tools/
## for all use 1000G
cd All/

cp ../../SCZ2_QCd.bim .
cp ../../SCZ2_QCd.bed .
cp ../../SCZ2_QCd.fam .

## liftover to hg19 for imputation

${PLINK}/plink --bfile SCZ2_QCd --update-map ${REF}GSAArray/liftoverhg19.txt 3 --make-bed --out SCZ2_QCd_hg19

## for HRC check tool need freq file
${PLINK}/plink --bfile SCZ2_QCd_hg19 --freq --out SCZ2_QCd_hg19_freq
perl ${KINGPATH}/HRC-1000G-check-bim.pl -b SCZ2_QCd_hg19.bim -f SCZ2_QCd_hg19_freq.frq -r ${KGG}/1000GP_Phase3_combined.legend -g
sed -i 's=plink=${PLINK}/plink=g' Run-plink.sh
sh Run-plink.sh
#for file in *.vcf; do awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' ${file} > with_chr_${file}; vcf-sort with_chr_${file} | bgzip -c > ${file}.gz;done
for file in *.vcf; do vcf-sort ${file} | bgzip -c > ${file}.gz;done
rm *.vcf
rm SCZ2_QCd*.*[^gz]



## for EUR use HRC
cd ../EUR
${PLINK}/plink --bfile ../../SCZ2_QCd --keep ../../EURSamples.txt --maf 0.05 --out SCZ2_QCd_EUR --make-bed
${PLINK}/plink --bfile SCZ2_QCd_EUR --update-map ${REF}GSAArray/liftoverhg19.txt 3 --make-bed --out SCZ2_QCd_EUR_hg19
${PLINK}/plink --bfile SCZ2_QCd_EUR_hg19 --freq --out SCZ2_QCd_hg19_EUR_freq

perl ${KINGPATH}/HRC-1000G-check-bim.pl -b SCZ2_QCd_EUR_hg19.bim -f SCZ2_QCd_hg19_EUR_freq.frq -r ${KGG}/../HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

sed -i 's=plink=${PLINK}/plink=g' Run-plink.sh
sh Run-plink.sh

#for file in *.vcf; do awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' ${file} > with_chr_${file}; vcf-sort with_chr_${file} | bgzip -c > ${file}.gz;done
for file in *.vcf; do vcf-sort ${file} | bgzip -c > ${file}.gz;done

rm *.vcf
rm SCZ2_QCd*.*[^gz]
