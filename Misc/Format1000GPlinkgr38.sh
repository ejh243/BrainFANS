## convert 1000 genomes data from vcf files to plink files 
## adapted tutorial: https://www.biostars.org/p/335605/ but using hg38 files

cd ${WD}

# Convert the 1000 Genomes files to BCF
## NB uses reference genome downloaded for alignment

for chr in {1..22}; do
   
	## need to add chr to chromosome notation
	gunzip -f ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz
	
	awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' ALL.chr${chr}_GRCh38.genotypes.20170504.vcf > ALL.chr${chr}_GRCh38.genotypes.20170504_with_chr.vcf
	
	bgzip -c ALL.chr${chr}_GRCh38.genotypes.20170504_with_chr.vcf > ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz
	tabix -p vcf ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz
	rm ALL.chr${chr}_GRCh38.genotypes.20170504_with_chr.vcf

	# Ensure that multi-allelic calls are split and that indels are left-aligned compared to reference genome (1st pipe)
	# Sets the ID field to a unique value: CHROM:POS:REF:ALT (2nd pipe)
	# Removes duplicates (3rd pipe) 
	# -I +'%CHROM:%POS:%REF:%ALT' means that unset IDs will be set to CHROM:POS:REF:ALT 
	# -x ID -I +'%CHROM:%POS:%REF:%ALT' first erases the current ID and then sets it to CHROM:POS:REF:ALT
    bcftools norm -m-any --check-ref w -f ${REFGENOME}genome.fa \
    ALL.chr"${chr}"_GRCh38.genotypes.20170504.vcf.gz | \

    bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' |

    bcftools norm -Ob --rm-dup both \
    > ALL.chr"${chr}"_GRCh38.genotypes.20170504.genotypes.bcf ;

    bcftools index ALL.chr"${chr}"_GRCh38.genotypes.20170504.genotypes.bcf ;
done

# Convert the BCF files to PLINK format
for chr in {1..22}; do
    ${PLINK}/plink --bcf ALL.chr"${chr}"_GRCh38.genotypes.20170504.genotypes.bcf  --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b38 no-fail --make-bed --out ALL.chr"${chr}"_GRCh38.genotypes.20170504.genotypes ;
done

# Exclude variants not on the coding strand
# NB - This step is only for microarray studies where the probes may only target one strand or the other (sense or non-sense)

# Prune variants from each chromosome
# --maf 0.10, only retain SNPs with MAF greater than 10%
# --indep [window size] [step size/variant count)] [Variance inflation factor (VIF) threshold]
# e.g. indep 50 5 1.5, Generates a list of markers in approx. linkage equilibrium - takes 50 SNPs at a time and then shifts by 5 for the window. VIF (1/(1-r^2)) is the cut-off for linkage disequilibrium

mkdir Pruned ;

for chr in {1..22}; do
    ${PLINK}/plink --noweb --bfile ALL.chr"${chr}"_GRCh38.genotypes.20170504.genotypes \
    --maf 0.10 --indep 50 5 1.5 \
    --out Pruned/ALL.chr"${chr}"_GRCh38.genotypes.20170504.genotypes ;

    ${PLINK}/plink --noweb --bfile ALL.chr"${chr}"_GRCh38.genotypes.20170504.genotypes \
    --extract Pruned/ALL.chr"${chr}"_GRCh38.genotypes.20170504.genotypes.prune.in --make-bed \
    --out Pruned/ALL.chr"${chr}"_GRCh38.genotypes.20170504.genotypes ;
done

#Get a list of all PLINK files
find . -name "*.bim" | grep -e "Pruned" > ForMerge.list ;

sed -i 's/.bim//g' ForMerge.list ;

# Merge all projects into a single PLINK file
${PLINK}/plink --merge-list ForMerge.list --out Merge ;