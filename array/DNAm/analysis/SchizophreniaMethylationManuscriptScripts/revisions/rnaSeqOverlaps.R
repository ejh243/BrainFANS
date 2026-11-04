##---------------------------------------------------------------------#
##
## Title: Overlap DMPs with RNAseq data from Ruzicka et al. 
##
##
## Input files: 
##
## Output files: 
##---------------------------------------------------------------------#



#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#



#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(ggplot2)



dataDir <- "Schizophrenia"
resFile <- file.path(dataDir, "4_analysis/EWAS/allCellType_LM_FDR.rdata")

outDir <- file.path(dataDir, "4_analysis/RNAseq/")
s4table <- file.path(outDir, "science.adg5136_data_s4.xlsx")

cellType <- "NeuN+"

setwd(outDir)

pThresh <- 9e-08
#pThresh <- 1e-04


#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

# load results files - annotated with gene info and snps/xhyb probes filtered out
load(resFile)
res <- res[[cellType]]
res <- res[res$chrm != "*",] # remove 30 probes which have '*' chr and start values of 0/-1

# subset to discovery DMPs
sig <- res[res$NullModel_SCZ_P < pThresh,]
sig <- sig[order(sig$NullModel_SCZ_P),]
sig$Probe_ID <- row.names(sig)


# get gene names
genes <- sig$GeneNames %>%
  str_split(";") %>%     
  unlist() %>%          
  str_trim() %>%        
  unique()              

genes <- genes[genes != "." & genes != ""]
# 15 for 9e-08

# remove these
genes <- genes[!genes %in% c("CTA-256D12.11", "AL163953.2", "AL163953.3")]

n_tests <- length(genes)  # number of genes tested

#----------------------------------------------------------------------#
# Get overlaps
#----------------------------------------------------------------------#

sheets <- excel_sheets(s4table)

results <- map_dfr(sheets, function(s) {
  read_excel(s4table, sheet = s) %>%
    dplyr::select(gene, Meta_logFC, Meta_P.Value, Meta_adj.P.Val) %>%
    mutate(gene = strsplit(gene, ";")) %>%
    unnest(gene) %>%
    mutate(gene = trimws(gene)) %>%
    dplyr::filter(gene %in% genes) %>%
    distinct() %>%
    mutate(cell_type = s, .before = 1)
})




sig_long <- sig %>%
  dplyr::select(Probe_ID, GeneNames, NullModel_SCZ_coeff, NullModel_SCZ_P) %>%
  mutate(gene = strsplit(GeneNames, ";")) %>%
  unnest(gene) %>%
  mutate(gene = trimws(gene)) %>%
  dplyr::filter(gene != ".", gene != "") %>%
  dplyr::select(-GeneNames) %>%
  distinct()

results_annot <- results %>%
  left_join(sig_long, by = "gene", relationship = "many-to-many")


write.csv(results_annot, paste0("rnaSeq", pThresh, "DMPs.csv"))


#----------------------------------------------------------------------#
# Tile Plot for genome-wide threhsold (coloured by logFC)
#----------------------------------------------------------------------#


cat_levels   <- c("Both negative", "Opposite directions")
group_levels <- c("Neurons", "Glia")


# function to format P value for plot
fmt_p <- function(p) {
  sapply(p, function(x) {
    s <- signif(x, 3)
    if (round(s, 3) == s) {
      formatC(s, format = "f", digits = 3, drop0trailing = TRUE)   
    } else {
      formatC(s, format = "e", digits = 2)                        
    }
  })
}


# create plot dataframe 
plot_df <- results %>%
  dplyr::filter(gene %in% genes) %>%
  group_by(cell_type, gene) %>%
  summarise(logFC = Meta_logFC[which.min(Meta_P.Value)],
            p_val = min(Meta_P.Value),
            .groups = "drop") %>%
  mutate(
    fill_val = logFC,
    sig       = p_val < 0.05 / n_tests,
    label    = if_else(sig, fmt_p(p_val), ""),
    group = case_when(
      str_detect(cell_type, regex("^\\s*ex[-_. ]", ignore_case = TRUE)) ~ "Neurons",
      str_detect(cell_type, regex("^\\s*in[-_. ]", ignore_case = TRUE)) ~ "Neurons",
      TRUE                                                              ~ "Glia"
    ),
    group     = factor(group, levels = group_levels),
    cell_type = factor(cell_type, levels = sheets),
    text_col  = if_else(abs(fill_val) > 0.6 * max(abs(fill_val)), "white", "black"),
    fontface  = if_else(sig, "bold", "plain")
  )

lim <- max(abs(plot_df$fill_val))





if(pThresh == 9e-08){
  
# Tile Plot

p <- ggplot(plot_df, aes(x = cell_type, y = gene, fill = fill_val)) +
  geom_tile(aes(alpha = sig), colour = "grey85", linewidth = 0.4) +
  geom_tile(data = dplyr::filter(plot_df, sig),
            fill = NA, colour = "black", linewidth = 0.8) +
  geom_text(aes(label = label, colour = text_col, fontface = fontface), size = 2) +
  scale_colour_identity() +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 1), guide = "none") +
  facet_grid(. ~ group, scales = "free_x", space = "free_x", drop = TRUE) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-lim, lim),
    name = "Meta logFC",
    guide = guide_colourbar(barheight = 10)
  ) +
  scale_x_discrete(expand = c(0, 0), drop = TRUE) +
  scale_y_discrete(limits = rev(genes)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank(),
        panel.spacing.x = unit(8, "pt"),
        strip.text = element_text(face = "bold", size = 14,
                                  margin = margin(b = 6)),
        axis.title = element_blank())




pdf(paste0("rnaSeq", pThresh, "DMPtilePlotColouredbylogFC.pdf"), height = 11, width = 15)
print(p)
dev.off()


} else {

# pThresh <- 1e-04


keep <- plot_df %>%
  group_by(gene) %>%
  summarise(n_sig = sum(sig), .groups = "drop") %>%
  dplyr::filter(n_sig > 0) %>%
  arrange(desc(n_sig))

keep <- keep %>% arrange(match(gene, genes))

plot_sub <- plot_df %>% dplyr::filter(gene %in% keep$gene)


p <- ggplot(plot_sub, aes(x = cell_type, y = gene, fill = fill_val)) +
  geom_tile(aes(alpha = sig), colour = "grey85", linewidth = 0.4) +
  geom_tile(data = dplyr::filter(plot_sub, sig),
            fill = NA, colour = "black", linewidth = 0.8) +
  geom_text(aes(label = label, colour = text_col, fontface = fontface), size = 2) +
  scale_colour_identity() +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 1), guide = "none") +
  facet_grid(. ~ group, scales = "free_x", space = "free_x", drop = TRUE) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-lim, lim),
    name = "Meta logFC",
    guide = guide_colourbar(barheight = 10)
  ) +
  scale_x_discrete(expand = c(0, 0), drop = TRUE) +
  scale_y_discrete(limits = rev(keep$gene)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank(),
        panel.spacing.x = unit(8, "pt"),
        strip.text = element_text(face = "bold", size = 14,
                                  margin = margin(b = 6)),
        axis.title = element_blank())



pdf(paste0("rnaSeq", pThresh, "DMPtilePlotColouredbylogFC.pdf"), height = 11, width = 15)
print(p)
dev.off()


}

#----------------------------------------------------------------------#
# Plot - coloured by DMPS / diff Expr genes are in same dir or not
#----------------------------------------------------------------------#

# RUN FOR BOTH pThresh <- 9e-08 and  pThresh <- 1e-04 

dmp_dir <- sig_long %>%
  group_by(gene) %>%
  summarise(dmp_dir = case_when(
    all(NullModel_SCZ_coeff < 0) ~ "Hypomethylated DMP",
    all(NullModel_SCZ_coeff > 0) ~ "Hypermethylated DMP",
    TRUE                         ~ "Mixed DMPs"
  ), .groups = "drop")

dir_levels <- c("Hypomethylated DMP", "Hypermethylated DMP", "Mixed DMPs")
cat_levels <- c("Same direction", "Opposite direction", "Mixed DMPs")
group_levels <- c("Neurons", "Glia")

plot_df <- results %>%
  dplyr::filter(gene %in% genes) %>%
  group_by(cell_type, gene) %>%
  summarise(logFC = Meta_logFC[which.min(Meta_P.Value)],
            p_val = min(Meta_P.Value),
            .groups = "drop") %>%
  left_join(dmp_dir, by = "gene") %>%          # this creates the dmp_dir column
  mutate(
    fill_val = logFC,
    sig      = p_val < 0.05 / n_tests,
    label    = if_else(sig, fmt_p(p_val), ""),
    dmp_dir  = factor(dmp_dir, levels = dir_levels),
    category = case_when(
      dmp_dir == "Mixed DMPs"                         ~ "Mixed DMPs",
      dmp_dir == "Hypomethylated DMP"  & logFC < 0    ~ "Same direction",
      dmp_dir == "Hypermethylated DMP" & logFC > 0    ~ "Same direction",
      TRUE                                            ~ "Opposite direction"
    ),
    category  = factor(category, levels = cat_levels),
    group     = if_else(str_detect(cell_type, regex("^\\s*(ex|in)[-_. ]", ignore_case = TRUE)),
                        "Neurons", "Glia"),
    group     = factor(group, levels = group_levels),
    cell_type = factor(cell_type, levels = sheets),
    text_col  = if_else(abs(fill_val) > 0.6 * max(abs(fill_val)), "white", "black")
  )




if(pThresh == 9e-08){
  

keep <- plot_df %>%
  group_by(gene) %>%
  summarise(n_sig = sum(sig), min_p = min(p_val), .groups = "drop") %>%
  dplyr::filter(n_sig > 0) %>%
  arrange(desc(n_sig), min_p)

plot_sub <- plot_df %>%
  dplyr::filter(gene %in% keep$gene) %>%
  mutate(gene = factor(gene, levels = rev(genes))) 


p_dir <- ggplot(plot_sub, aes(x = cell_type, y = gene, fill = category)) +
  geom_tile(aes(alpha = sig), colour = "grey85", linewidth = 0.4) +
  geom_tile(data = dplyr::filter(plot_sub, sig),
            fill = NA, colour = "black", linewidth = 0.8) +
  geom_text(aes(label = label), size = 2.5, fontface = "bold") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.5), guide = "none") +
  facet_grid(. ~ group, scales = "free_x", space = "free_x", drop = TRUE) +
  
  scale_fill_manual(
    values = c("Same direction"     = "#9ECAE1",   # pale blue
               "Opposite direction" = "#F4A6A6",   # pale red
               "Mixed DMPs"         = "grey80"),
    breaks = cat_levels, limits = cat_levels, name = NULL
  ) +
  guides(fill = guide_legend(override.aes = list(colour = "grey50"))) +
  scale_x_discrete(expand = c(0, 0), drop = TRUE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank(),
        panel.spacing = unit(8, "pt"),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_blank())


pdf(paste0("rnaSeq", pThresh, "DMPtilePlotColouredBySameOrOpposite.pdf"), height = 11, width = 15)
print(p_dir)
dev.off()



} else {
  

plot_sub <- plot_df %>%
  dplyr::filter(gene %in% keep$gene) %>%
  mutate(gene = factor(gene, levels = rev(keep$gene)))


p_dir <- ggplot(plot_sub, aes(x = cell_type, y = gene, fill = category)) +
  geom_tile(aes(alpha = sig), colour = "grey85", linewidth = 0.4) +
  geom_tile(data = dplyr::filter(plot_sub, sig),
            fill = NA, colour = "black", linewidth = 0.8) +
  geom_text(aes(label = label), size = 2.5, fontface = "bold") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.5), guide = "none") +
  facet_grid(dmp_dir ~ group, scales = "free", space = "free", drop = TRUE) +
  scale_fill_manual(
    values = c("Same direction"     = "#9ECAE1",   # pale blue
               "Opposite direction" = "#F4A6A6",   # pale red
               "Mixed DMPs"         = "grey80"),
    breaks = cat_levels, limits = cat_levels, name = NULL
  ) +
  guides(fill = guide_legend(override.aes = list(colour = "grey50"))) +
  scale_x_discrete(expand = c(0, 0), drop = TRUE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        panel.grid = element_blank(),
        panel.spacing = unit(8, "pt"),
        strip.text = element_text(face = "bold", size = 12),
        axis.title = element_blank())


pdf(paste0("rnaSeq", pThresh, "DMPtilePlotColouredBySameOrOpposite.pdf"), height = 11, width = 15)
print(p_dir)
dev.off()

}
