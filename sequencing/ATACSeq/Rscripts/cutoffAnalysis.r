## ================================================================================================##
##                                 PIPPA STEP 9.4: Cutoff analysis                                 ##
## ================================================================================================##
## EXECUTION: Rscript ./Rscripts/cutoffAnalysis.r <project_config> <cf>                            ||
## - execute from scripts directory                                                                ||
##                                                                                                 ||
## DESCRIPTION: This scripts reads peaks called at relaxed threshold and performs segmented        ||
## analysis as an alternative to cutoff analysis to decide a better threshold for peak calling.    ||
##                                                                                                 ||
## INPUTS:                                                                                         ||
## - <project_config> : path to config.r file.                                                     ||
## - <cf> : cell type to perform segmented analysis on                                             ||
##                                                                                                 ||
## OUTPUTS:                                                                                        ||
## -  plot with segmented analysis results: "cutoffAnalysis_<cf>.pdf"                              ||
##                                                                                                 ||
## REQUIRES:                                                                                       ||
## - R version > 4.3, R libraries: rtracklayer, GenomicRanges, segmented,tidyr, ggplot2            ||
## - A peak file generated in STEP 9.3: "<cf>_<macs_thres_group>.sorted.narrowPeak.filt"           ||
## - Parameters in config.r file: col_vals, peaksAnalysis, macs_thres_group, peaksCutoff           ||
## ================================================================================================##

## ==========##
##   SET-UP  ##
## ==========##

suppressWarnings(suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(rtracklayer)
  library(segmented)
  library(tidyr)
  library(dplyr)
}))


args <- commandArgs(trailingOnly=TRUE)
configFile<-args[1]
cf<-args[2]
source(configFile)

plot_colour <- col_vals_plots[cf]

## ===============##
##    FUNCTIONS   ##
## ===============##

format_q_label <- function(q) {
  
  exponent <- round(log10(q))
  
  mantissa <- q / (10^exponent)
  
  nearest <- c(1, 5)[
    which.min(abs(c(1, 5) - mantissa))
  ]
  
  paste0(nearest, "e", exponent)
}

## ===============##
##   GATHER DATA  ##
## ===============##

## get peaks called with q5e-2, only blacklisted regions peaks filtered
peakFile <- paste0(peaksAnalysis,"/", cf,"/", cf,"_",macs_thres_group,".narrowPeak.filt")
peaks.gr<-makeGRangesFromDataFrame(read.table(file.path(peakFile)), keep.extra.columns=TRUE, seqnames.field="V1",
                                   start.field="V2",
                                   end.field="V3",
                                   ignore.strand=TRUE)

peaks.df <- as.data.frame(peaks.gr)
colnames(peaks.df) <- c("chr", "start", "end", "width", "strand", "peakID", "score", "strand2", "signal", "pvalue", "qvalue", "score2")

peaks.df$pvals <- 10^(-peaks.df$pvalue)
peaks.df$qvals<- 10^(-peaks.df$qvalue)

## Set of thresholds to be considered. From 0.05 to 1e-100
exponents <- -seq(2, 100, by = 2)

thresholds <- 10^exponents
order_levels.og <- as.character(thresholds)

## Number of peaks vs q-values
peaks.df$qvals <- 10^(-peaks.df$qvalue) 
counts_list <- sapply(
  thresholds,
  function(t) {
    sum(peaks.df$qvals < t, na.rm = TRUE)
  }
)
df <- data.frame(
    threshold = thresholds,
    logQ = -log10(thresholds),
    count = counts_list
)

## ======================##
##   SEGMENTED ANALYSIS  ##
## ======================##

#segmented analysis: fits two linear regimes and estimates a breakpoint.
model <- lm(count ~ logQ, data = df)
seg_model <- segmented(
    model,
    seg.Z = ~logQ
  )
  
seg_breaks <- abs(seg_model$psi[, "Est."])
ci_segmented <- confint(seg_model)
ci_df <- as.data.frame(ci_segmented)
ci_lower <- ci_df[, grep("low",colnames(ci_df),ignore.case = TRUE)]
ci_upper <- ci_df[, grep( "up",colnames(ci_df),ignore.case = TRUE)]

seg_labels <- format_q_label(10^(-seg_breaks))
ci_lower_labels <- format_q_label(10^(-ci_lower))
ci_upper_labels <- format_q_label(10^(-ci_upper))

segmented_breaks = seg_breaks
segmented_CI = ci_segmented


summary_df <- data.frame(Breakpoint_log10Q = seg_breaks, Breakpoint_Qvalue = signif(10^(-seg_breaks),3),
      Breakpoint_Label = seg_breaks, CI95_Lower_log10Q = ci_lower, CI95_Upper_log10Q = ci_upper,
      CI95_Lower_Qvalue = signif(10^(-ci_upper),3),
      CI95_Upper_Qvalue = signif(10^(-ci_lower),3),
      CI95_Lower_Label = ci_upper_labels,CI95_Upper_Label = ci_lower_labels)

ymax <- max(df$count)
  
seg_y <- seq(ymax * 0.95,ymax * 0.80,length.out = length(seg_breaks))
ci_low_y <- seq(ymax * 0.75,ymax * 0.65,length.out = length(ci_lower))
ci_up_y <- seq(ymax * 0.60,ymax * 0.50,length.out = length(ci_upper))

## ===============================##
##   GENERATE PLOT WITH ANALYSIS  ##
## ===============================##

pdf(file.path(paste0(peaksCutoff, "/cutoffAnalysis_", cf,".pdf")), width = 10, height = 5)

p <- ggplot(df,aes(x = logQ, y = count)) +
  geom_point(size = 1.3,color = plot_colour) +
  geom_line(size = 1.2, color = plot_colour) +
  geom_vline(xintercept = seg_breaks,color = "black",linewidth = 1) +
  annotate("text", x = seg_breaks + 0.15, y = ymax * 1.02,label = seg_labels, color = "black",vjust = 0,hjust = 0,size = 4.5,fontface = "bold") +
  geom_vline(xintercept = ci_lower,color = "darkblue",linewidth = 1,linetype = "dashed",alpha = 0.7) +
  annotate("text",x = ci_lower - 1,y = ci_low_y,label = ci_lower_labels,color = "black",hjust = 1,vjust = -0.3,size = 4) +
  geom_vline(xintercept = ci_upper, color = "darkblue", linetype = "dashed", linewidth = 1,alpha = 0.7) +
   annotate( "text", x = ci_upper + 1, y = ci_up_y,label = ci_upper_labels,color = "black",hjust = 0, vjust = -0.3,size = 4) +
  scale_y_continuous(labels = scales::label_number( suffix = "K",scale = 1e-3))+
  scale_x_continuous(name = expression(paste("Q value threshold (1e"^{-n}, ")")),breaks = seq(2, 100, by = 2),labels = seq(2, 100, by = 2)) +
  labs(y = "Number of peaks", title = paste0("Cutoff analysis for ",cf, " peaks")) + 
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,),
  plot.title = element_text(size = 12,face = "bold",hjust = 0.5))

print(p)
dev.off()