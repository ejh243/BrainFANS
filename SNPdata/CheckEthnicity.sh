

## check sample ethnicity by comparing to 1000G
cd ${DATADIR}/SNPdata/

# change variant ids to chr:bp
awk '{if ($1 != 0) print $2,"chr"$1":"$4}' SCZ2_QCd.bim > updateTo1KGFormat.txt
${PLINK}/plink --bfile SCZ2_QCd --update-name updateTo1KGFormat.txt --make-bed --out SCZ2_QCd_1kgIDs

# first merge with 1000 genomes and filter variants to those in common
# need to test initially in case of error with triallelic variants
${PLINK}/plink --bfile SCZ2_QCd_1kgIDs --bmerge ${KGG}/1000G_gr38.bed ${KGG}/1000G_gr38.bim ${KGG}/1000G_gr38.fam --maf 0.1 --geno 0.1 --make-bed --out mergedw1000G_test

## issue with variants at same position but different alleles - exclude these
${PLINK}/plink --bfile SCZ2_QCd_1kgIDs --exclude mergedw1000G_test-merge.missnp --make-bed --out SCZ2_gr38_update_3_1kgIDs_forMerge

${PLINK}/plink --bfile SCZ2_gr38_update_3_1kgIDs_forMerge --bmerge ${KGG}/1000G_gr38.bed ${KGG}/1000G_gr38.bim ${KGG}/1000G_gr38.fam --maf 0.2 --geno 0.05 --make-bed --out mergedw1000G
# LD prune
${PLINK}/plink --bfile mergedw1000G --indep 50 5 1.5 --out mergedw1000G.ld
${PLINK}/plink --bfile mergedw1000G --extract mergedw1000G.ld.prune.in --make-bed --out mergedw1000G.ld.prune

rm SCZ2_gr38_update_3_1kgIDs_forMerge*
rm mergedw1000G_test*

# use GCTA to calc PCs
${GCTA}/gcta64 --bfile mergedw1000G.ld.prune --make-grm-bin --autosome --out mergedw1000G
${GCTA}/gcta64 --grm mergedw1000G --pca --out mergedw1000G.pca

