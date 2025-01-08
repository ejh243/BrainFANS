##---------------------------------------------------------------------#
##
## Title: Test age acceleration against disease status
##
## Purpose of script: compare performance of epigenetic clocks by cell type
## and against disease status
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# DEFINE FUNCTION FOR CALCULATING PREDICTED AGE
#----------------------------------------------------------------------#
# taken from Horvath code base
anti.trafo= function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }

# calculate RMSE
rmse <- function(var1, var2){
    round(sqrt(mean(var1-var2)^2),2)
}
rmseByGroup <- function(var1, var2, group){
    groupLabels <- unique(group)
    res <- rep(NA, length(groupLabels))
    names(res) <- groupLabels
    for(each in groupLabels){
        index<-which(group == each)
        res[each]<-round(sqrt(mean(var1[index]-var2[index])^2),2)
    }
    return(res)
}
corByGroup <- function(var1, var2, group){
    groupLabels <- unique(group)
    res <- rep(NA, length(groupLabels))
    names(res) <- groupLabels
    for(each in groupLabels){
        index<-which(group == each)
        res[each]<-cor(var1[index],var2[index])
    }
    return(res)
}

## test for cell type specific age cor
runAgeCTCorTests <- function(data, dependent_var, fixed_effects, random_effects_full, random_effects_null) {
  # Ensure dependent_var is treated as a symbol
  dependent_var <- rlang::ensym(dependent_var)
  
  # Create full and null model formulas
  full_model_formula <- as.formula(
    paste(rlang::as_label(dependent_var), "~", fixed_effects, "+", random_effects_full)
  )
  null_model_formula <- as.formula(
    paste(rlang::as_label(dependent_var), "~", fixed_effects, "+", random_effects_null)
  )
  
  # Fit the full and null models
  full_model <- lmer(full_model_formula, REML = FALSE, data = data)
  null_model <- lmer(null_model_formula, REML = FALSE, data = data)
  
  return(c(fixef(full_model)[1] ,
    as.numeric(t(coef(summary(full_model))[
    c("Age", "Cell.typeNeuN+", "Cell.typeSox10+", "Age:Cell.typeNeuN+", "Age:Cell.typeSox10+"),
    c("Estimate", "Pr(>|t|)")])),
    anova(full_model, null_model)["full_model","Pr(>Chisq)"]))
}

plotScatterAgainstAge <- function(data, value_column, group_column = "Cell.type"){
    # Ensure the value_column argument is evaluated correctly
    value_column <- rlang::ensym(value_column)
    group_column <- rlang::ensym(group_column)

    ggplot(data, aes(x = Age, y = !!value_column, colour = !!group_column))+
        geom_point() +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")
}

plotViolinWithMeans <- function(data, value_column, group_column = "Phenotype", 
    facet_column = "Cell.type") {
  # Ensure the value_column argument is evaluated correctly
  value_column <- rlang::ensym(value_column)
  group_column <- rlang::ensym(group_column)
  facet_column <- rlang::ensym(facet_column)

  # Calculate mean for each group
  mean_data <- data %>%
    group_by(!!facet_column, !!group_column) %>%
    summarize(mean_value = mean(!!value_column, na.rm = TRUE), .groups = "drop")

  # Create the violin plot
  ggplot(data, aes(x = !!group_column, y = !!value_column, fill = !!group_column)) +
    geom_violin(scale = "width", position = position_dodge(width = 0.9)) +
    facet_wrap(as.formula(paste("~", rlang::as_label(facet_column)))) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    geom_point(data = mean_data, aes(y = mean_value, group = !!group_column),
               color = "black", shape = 18, size = 2,
               position = position_dodge(width = 0.9))
}

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(ggplot2)
library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggpubr)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
refDir <- args[2]

normData<-file.path(dataDir, "/3_normalised/normalised.rdata")
resPath<-file.path(dataDir, "/4_analysis/EWAS/")

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)

## remove total samples and cell types with less than 20 samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type != "Total"),]
nSample<-table(QCmetrics$Cell.type)
QCmetrics<-QCmetrics[QCmetrics$Cell.type %in% names(nSample[which(nSample > 19)]),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

# filter to those with Age
QCmetrics<-QCmetrics[!is.na(QCmetrics$Age),]

celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]
cellTypes<-unique(QCmetrics$Cell.type)
celltypeNormbeta<-na.omit(celltypeNormbeta)

## read in the CorticalClockCoefficients from Shireby et al
CC_coef<-read.csv(file.path(refPath, "CortexClock/CorticalClockCoefficients.csv"), stringsAsFactors = FALSE)
## calc age with cortical clock & age acceleration
CCDNAmAge<-	data.frame("Basename" = QCmetrics$Basename,
    "CCDNAmAge" = anti.trafo(as.numeric(CC_coef[1,2] + 
    t(celltypeNormbeta[CC_coef[-1,1],])  %*% CC_coef[-1,2])),
    "Age" = QCmetrics$Age,
    "Cell.type" = QCmetrics$Cell.type, 
    "Phenotype" = QCmetrics$Phenotype, 
    "Sex" = QCmetrics$Sex, 
    "Individual.ID" = QCmetrics$Indidivual.ID, 
    "Tissue.Centre" = QCmetrics$Tissue.Centre)

CCDNAmAge<-na.omit(CCDNAmAge)

CCDNAmAge$AAbrain <- residuals(lm(CCDNAmAge~Age, data=CCDNAmAge))
CCDNAmAge$AAbrainByCT <- residuals(lm(CCDNAmAge~Age + Cell.type, data=CCDNAmAge))

## calculate Horvath clock
coefHorvath<- read.csv(file.path(refPath, "HorvathClock/HorvathClockCoefficients.csv"), stringsAsFactors = FALSE)
coefHorvath <- rbind(coefHorvath[1,], 
    coefHorvath[coefHorvath$CpGmarker %in% rownames(celltypeNormbeta),])
CCDNAmAge$DNAmAge<-anti.trafo(coefHorvath$CoefficientTraining[1] + 
    t(celltypeNormbeta[coefHorvath$CpGmarker[-1],]) %*% coefHorvath$CoefficientTraining[-1])
CCDNAmAge$AAhorvath <- residuals(lm(DNAmAge ~ Age, data=CCDNAmAge))
CCDNAmAge$AAhorvathByCT <- residuals(lm(DNAmAge ~ Age + Cell.type, data=CCDNAmAge))


#----------------------------------------------------------------------#
# CREATE PLOTS & STATISTICAL TESTS
#----------------------------------------------------------------------#

p1<-plotScatterAgainstAge(CCDNAmAge, CCDNAmAge) + 
    labs(colour = "Cell Type")
p2<-plotScatterAgainstAge(CCDNAmAge, DNAmAge) + 
    labs(colour = "Cell Type")
combinedPlot <- ggarrange(p1,p2, 
                           ncol = 2, nrow = 1, 
                           common.legend = TRUE, legend = "bottom")

pdf(file.path(resPath, "Plots", paste0("ScatterplotEpigeneticAgeVsChronologicalAge.pdf")), 
    width = 10, height = 5)
combinedPlot
dev.off()

rmsePerCT<-rbind("CCDNAmAge" = c(rmse(CCDNAmAge$Age, CCDNAmAge$CCDNAmAge),rmseByGroup(CCDNAmAge$Age, CCDNAmAge$CCDNAmAge, CCDNAmAge$Cell.type)),
"DNAmAge" = c(rmse(CCDNAmAge$Age, CCDNAmAge$DNAmAge),rmseByGroup(CCDNAmAge$Age, CCDNAmAge$DNAmAge, CCDNAmAge$Cell.type)))
corPerCT<-rbind("CCDNAmAge" = c(cor(CCDNAmAge$Age, CCDNAmAge$CCDNAmAge),corByGroup(CCDNAmAge$Age, CCDNAmAge$CCDNAmAge, CCDNAmAge$Cell.type)),
"DNAmAge" = c(cor(CCDNAmAge$Age, CCDNAmAge$DNAmAge),corByGroup(CCDNAmAge$Age, CCDNAmAge$DNAmAge, CCDNAmAge$Cell.type)))

write.csv(rmsePerCT,file.path(resPath, "Tables", "RMSEEpigeneticClockAgeByCelltypes.csv"))
write.csv(corPerCT, file.path(resPath, "Tables", "CorEpigeneticClockAgeByCelltypes.csv"))


p1 <- plotViolinWithMeans(CCDNAmAge, value_column = "AAbrain") + xlab("") + ylab("Age Acceleration (Cortical)")
p2 <- plotViolinWithMeans(CCDNAmAge, value_column = "AAbrainByCT") + xlab("") + ylab("Age Acceleration (Cortical)")
p3 <- plotViolinWithMeans(CCDNAmAge, value_column = "AAhorvath") + xlab("") + ylab("Age Acceleration (Horvath)")
p4 <- plotViolinWithMeans(CCDNAmAge, value_column = "AAhorvathByCT") + xlab("") + ylab("Age Acceleration (Horvath)")

combinedPlot <- ggarrange(p1,p3, 
                           ncol = 1, nrow = 2, 
                           common.legend = TRUE, legend = "bottom")

pdf(file.path(resPath, "Plots", paste0("ViolinPlotAgeAccelerationAltogetherByCaseControl.pdf")), 
    width = 10, height = 10)
combinedPlot
dev.off()

combinedPlot <- ggarrange(p2,p4, 
                           ncol = 1, nrow = 2, 
                           common.legend = TRUE, legend = "bottom")

pdf(file.path(resPath, "Plots", paste0("ViolinPlotAgeAccelerationByCTByCaseControl.pdf")), 
    width = 10, height = 10)
combinedPlot
dev.off()

ageCorP<-matrix(data = NA, nrow = 2, ncol = 12)
rownames(ageCorP)<-c("CCDNAmAge", "DNAmAge")
colnames(ageCorP)<-c("Intercept", apply(expand.grid(
    c("Coef", "P"), 
    c("Age", "Cell.typeNeuN+","Cell.typeSox10+","Age_Cell.typeNeuN+" ,"Age_Cell.typeSox10+")), 1, paste, collapse = "_"),
    "P_Individual")
for(each in c("CCDNAmAge", "DNAmAge")){
    ageCorP[each,]<-runAgeCTCorTests(
    data = CCDNAmAge,
    dependent_var = !!sym(each),
    fixed_effects = "Age * Cell.type + Sex",
    random_effects_full = "(1 | Tissue.Centre) + (1 | Individual.ID)",
    random_effects_null = "(1 | Tissue.Centre)"
    )
}

write.csv(ageCorP, file.path(resPath, "Tables", "MLMEpigeneticClockAgeByCelltypes.csv"))

sczCorP<-matrix(data = NA, nrow = 2, ncol = 12)
rownames(ageCorP)<-c("CCDNAmAge", "DNAmAge")
colnames(ageCorP)<-c("Intercept", apply(expand.grid(
    c("Coef", "P"), 
    c("Age", "Cell.typeNeuN+","Cell.typeSox10+","Age_Cell.typeNeuN+" ,"Age_Cell.typeSox10+")), 1, paste, collapse = "_"),
    "P_Individual")
for(each in c("CCDNAmAge", "DNAmAge")){
    ageCorP[each,]<-runAgeCTCorTests(
    data = CCDNAmAge,
    dependent_var = !!sym(each),
    fixed_effects = "Age * Cell.type + Sex",
    random_effects_full = "(1 | Tissue.Centre) + (1 | Individual.ID)",
    random_effects_null = "(1 | Tissue.Centre)"
    )
}

AABySCZ<-NULL
for(each in c("CCDNAmAge", "DNAmAge")){
    for(CT in cellTypes){
        dependent_var <- rlang::ensym(each)
        # Define the full model formula dynamically
        full_model_formula <- as.formula(paste(rlang::as_label(dependent_var), "~ Age*Phenotype + Sex"))
        full_model <- lm(full_model_formula, data = CCDNAmAge, subset = which(Cell.type == CT))
        residuals_reduced <- residuals(lm(CCDNAmAge ~ Age + Sex, data = CCDNAmAge, subset = which(Cell.type == CT)))
        ttest_result <- t.test( residuals_reduced ~ CCDNAmAge$Phenotype[which(CCDNAmAge$Cell.type == CT)])
        AABySCZ<-rbind(AABySCZ, (c(each, CT, 
        as.numeric(t(summary(full_model)$coefficients[c("Age", "PhenotypeSchizophrenia", "Age:PhenotypeSchizophrenia"), c("Estimate", "Pr(>|t|)")])),
        diff(ttest_result$estimate), ttest_result$p.value)))
    }
}
colnames(AABySCZ)<-c("Clock", "CT", apply(expand.grid(
    c("Coef", "P"), 
    c("Age", "SchizophreniaStatus","AgeBySchizophrenia")), 1, paste, collapse = "_"), "Ttest_MeanDiff", "Ttest_P")

write.csv(AABySCZ, file.path(resPath, "Tables", "EpigeneticClockAgeBySchizophrenia.csv"))