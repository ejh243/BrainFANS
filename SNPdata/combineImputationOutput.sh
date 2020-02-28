## collate output of imputation server, convert vcf files to plink and filter variants on INFO score
cd All
for i in {7..22}
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

${PLINK}/plink --bfile chr${i}_rsq0.3 --merge-list mergefile.txt --make-bed --out BrainFANS_1000GALL_rsq0.3

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
awk '{split($1,a,"_"); print $1,$2,a[1]"_"a[2],a[3]}' BrainFANS_1000GALL_rsq0.3.fam > UpdateIDs.txt
${PLINK}/plink --bfile BrainFANS_1000GALL_rsq0.3 --update-ids UpdateIDs.txt --make-bed --out BrainFANS_1000GALL_rsq0.3_1 
${PLINK}/plink --bfile BrainFANS_1000GALL_rsq0.3_1 --update-sex ../../UpdateSex.txt --hwe 0.00001 --make-bed --out BrainFANS_1000GALL_rsq0.3_QCd 

rm BrainFANS_1000GALL_rsq0.3.b*
rm BrainFANS_1000GALL_rsq0.3.fam
rm BrainFANS_1000GALL_rsq0.3_1.b*
rm BrainFANS_1000GALL_rsq0.3_1.fam

## recalculate PCs
# LD prune
${PLINK}/plink --bfile BrainFANS_1000GALL_rsq0.3_QCd --maf 0.05 --indep 50 5 1.5 --out BrainFANS_1000GALL_rsq0.3_QCd.ld
${PLINK}/plink --bfile BrainFANS_1000GALL_rsq0.3_QCd --extract BrainFANS_1000GALL_rsq0.3_QCd.ld.prune.in --make-bed --out BrainFANS_1000GALL_rsq0.3_QCd.ld.prune

# use GCTA to calc PCs
${GCTA}/gcta64 --bfile BrainFANS_1000GALL_rsq0.3_QCd.ld.prune --make-grm-bin --autosome --out BrainFANS_1000GALL_rsq0.3_QCd
${GCTA}/gcta64 --grm BrainFANS_1000GALL_rsq0.3_QCd --pca --out BrainFANS_1000GALL_rsq0.3_QCd.pca

rm BrainFANS_1000GALL_rsq0.3_QCd.ld.prune.b*
rm BrainFANS_1000GALL_rsq0.3_QCd.ld.prune.fam


## plot PCs to identify outliers
Rscript ../../../../scripts/SNPdata/plotPCs.r BrainFANS_1000GALL_rsq0.3_QCd.pca.eigenvec 3

cd ../EUR
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

${PLINK}/plink --bfile chr${i}_rsq0.3 --merge-list mergefile.txt --make-bed --out BrainFANS_HRC1.1EUR_rsq0.3

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

##update IDs
awk '{split($1,a,"_"); print $1,$2,a[1]"_"a[2],a[3]}' BrainFANS_HRC1.1EUR_rsq0.3.fam > UpdateIDs.txt
${PLINK}/plink --bfile BrainFANS_HRC1.1EUR_rsq0.3 --update-ids UpdateIDs.txt --make-bed --out BrainFANS_HRC1.1EUR_rsq0.3_1

## update sex in fam file run hwe filter
${PLINK}/plink --bfile BrainFANS_HRC1.1EUR_rsq0.3_1 --update-sex ../../UpdateSex.txt --hwe 0.00001 --make-bed --out BrainFANS_HRC1.1EUR_rsq0.3_QCd 

rm BrainFANS_HRC1.1EUR_rsq0.3.b*
rm BrainFANS_HRC1.1EUR_rsq0.3.fam

## recalculate PCs
# LD prune
${PLINK}/plink --bfile BrainFANS_HRC1.1EUR_rsq0.3_QCd --indep 50 5 1.5 --out BrainFANS_HRC1.1EUR_rsq0.3_QCd.ld
${PLINK}/plink --bfile BrainFANS_HRC1.1EUR_rsq0.3_QCd --extract BrainFANS_HRC1.1EUR_rsq0.3_QCd.ld.prune.in --make-bed --out BrainFANS_HRC1.1EUR_rsq0.3_QCd.ld.prune


# use GCTA to calc PCs
${GCTA}/gcta64 --bfile BrainFANS_HRC1.1EUR_rsq0.3_QCd.ld.prune --make-grm-bin --autosome --out BrainFANS_HRC1.1EUR_rsq0.3_QCd
${GCTA}/gcta64 --grm BrainFANS_HRC1.1EUR_rsq0.3_QCd --pca --out BrainFANS_HRC1.1EUR_rsq0.3_QCd.pca

rm BrainFANS_HRC1.1EUR_rsq0.3_QCd.ld.prune.b*
rm BrainFANS_HRC1.1EUR_rsq0.3_QCd.ld.prune.fam

## plot PCs to identify outliers
Rscript ../../../../scripts/SNPdata/plotPCs.r BrainFANS_HRC1.1EUR_rsq0.3_QCd.pca.eigenvec 3

## pull out outliers on PC1 
grep -w "Outlier PC1"  OutliersFromPC_3SDfromMean.txt | cut -f 2-3 --delim=" " > filterSamples.txt
cp BrainFANS_HRC1.1EUR_rsq0.3_QCd.bed BrainFANS_HRC1.1EUR_rsq0.3_QCd_tmp.bed
cp BrainFANS_HRC1.1EUR_rsq0.3_QCd.bim BrainFANS_HRC1.1EUR_rsq0.3_QCd_tmp.bim
cp BrainFANS_HRC1.1EUR_rsq0.3_QCd.fam BrainFANS_HRC1.1EUR_rsq0.3_QCd_tmp.fam

${PLINK}/plink --bfile BrainFANS_HRC1.1EUR_rsq0.3_QCd --update-sex ../../UpdateSex.txt --maf 0.01 --hwe 0.00001 --make-bed --out BrainFANS_HRC1.1EUR_rsq0.3_QCd_tmp

