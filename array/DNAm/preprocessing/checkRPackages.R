required_packages <- c(
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

"%ni%" <- Negate("%in%")
if (!all(required_packages %in% rownames(installed.packages()))) {
    absent_packages <-
        required_packages[required_packages %ni% rownames(installed.packages())]
    number_of_missing_packages <- length(absent_packages)
    absent_packages <- paste(absent_packates, sep = ", ")
    stop(
        "Failed installation of ",
        number_of_missing_packages,
        " packages: ",
        absent_packages
    )
}
