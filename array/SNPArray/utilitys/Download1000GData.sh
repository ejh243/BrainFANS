## download 1000 genomes data 
## adapted tutorial: https://www.biostars.org/p/335605/ but using hg38 files

cd ${WD}

# Download the files as VCF.gz (and tab-indices)
prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr" ;

suffix="_GRCh38.genotypes.20170504.vcf.gz" ;

for chr in {1..22}; do
    wget "${prefix}""${chr}""${suffix}" "${prefix}""${chr}""${suffix}".tbi ;
done

# Download 1000 Genomes PED file
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped ;