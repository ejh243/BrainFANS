all_pkgs <- c(
    "pander", "kableExtra", "gplots", "diptest", "corrplot", "pheatmap",
    "mixtools", "RColorBrewer", "rmarkdown", "e1071", "stringdist",
    "data.table", "matrixStats", "quadprog", "reshape2", "wateRmelon",
    "bigmelon", "genefilter", "minfi", "CETYGO", "cdegUtilities",
    "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylation450kmanifest",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
    "IlluminaHumanMethylationEPICmanifest",
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
    "IlluminaHumanMethylationEPICv2manifest"
)

'%ni%' <- Negate('%in%')
if (all(all_pkgs %in% rownames(installed.packages()))) {
    print("All packages successfully installed")
} else {
    absent <- all_pkgs[all_pkgs %ni% rownames(installed.packages())]
    stop("Failed installation of ", length(absent), " packages: ", absent)
}
