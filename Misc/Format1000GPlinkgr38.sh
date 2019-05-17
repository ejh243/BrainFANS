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
	rm ALL.chr${chr}_GRCh38.genotypes.20170504.vcf

	# Ensure that multi-allelic calls are split and that indels are left-aligned compared to reference genome (1st pipe)
	# Sets the ID field to a unique value: CHROM:POS:REF:ALT (2nd pipe)
	# Removes duplicates (3rd pipe) 
	# -I +'%CHROM:%POS:%REF:%ALT' means that unset IDs will be set to CHROM:POS:REF:ALT 
	# -x ID -I +'%CHROM:%POS:%REF:%ALT' first erases the current ID and then sets it to CHROM:POS
    bcftools norm -m-any --check-ref w -f ${REFGENOME}genome.fa \
    ALL.chr"${chr}"_GRCh38.genotypes.20170504.vcf.gz | \

    bcftools annotate -x ID -I +'%CHROM:%POS' |

    bcftools norm -Ob --rm-dup both \
    > ALL.chr"${chr}"_GRCh38.genotypes.20170504.genotypes.bcf ;
	
    bcftools index ALL.chr"${chr}"_GRCh38.genotypes.20170504.genotypes.bcf ;
	
	bcftools index ALL.chr"${chr}"_GRCh38.biallelic.bcf
	
done

# Convert the BCF files to PLINK format
for chr in {1..22}; do
    ${PLINK}/plink --bcf ALL.chr"${chr}"_GRCh38.genotypes.20170504.genotypes.bcf  --biallelic-only --vcf-idspace-to _ --const-fid --allow-extra-chr 0 --split-x b38 no-fail --make-bed --out ALL.chr"${chr}"_GRCh38.20170504.biallelic ;
done


#Get a list of all PLINK files
find . -name "ALL.chr*.bim" > ForMerge.list ;

sed -i 's/.bim//g' ForMerge.list ;

# Merge all projects into a single PLINK file
${PLINK}/plink --merge-list ForMerge.list --out 1000G_gr38;

## remove variants that are at the same position (i.e. triallelic) 
cut -f 2 1000G_gr38_2.bim | uniq -d > dupVariants.txt
${PLINK}/plink --bfile 1000G_gr38 --exclude dupVariants.txt --make-bed --out 1000G_gr38_biallelic;








