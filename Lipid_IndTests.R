
### In this script, I will test for effects of diapause and cold acclimation on 
# on the relative levels of individual phospholipid classes and fatty acyl methyl ester species using a planned set of five
# pairwise comparisons. For metabolites with normally-distributed residuals,
# I will do this using an lmer followed by pairwise t-tests using emmeans.
# For species without normally-distributed residuals, I will run non-parametric
# tests (Kruskal test followed by pairwise wilcox).
# I will use the results of these tests to generate annotation columns for a 
# heatmap that also shows z-scored values across each of the six total 
# treatment groups.

# load packages
library(dplyr)
library(reshape2)
library(stringr)
library(tibble)
library(viridis)
library(paletteer)
library(tidyr)
library(stats)
library(pheatmap)
library(car)
library(lme4)
library(ggplot2)
library("emmeans")

######## read in both original data-sets
### read in lipid head group data 
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Data Files for Submission")
lipdata <- read.csv(file="PhospholipidClass_RawData.csv")

## rename factors 
lipdata$diapause_treat[lipdata$diapause_treat == 'SD'] <- 'Diapause'
lipdata$diapause_treat[lipdata$diapause_treat == 'LD'] <- 'Non-diapause'
lipdata$acc_treat[lipdata$acc_treat == 'PREACC'] <- 'Before Acclimation'
lipdata$acc_treat[lipdata$acc_treat == 'WARM'] <- 'Warm Acclimated (15d)'
lipdata$acc_treat[lipdata$acc_treat == 'COLD'] <- 'Cold Acclimated (15d)'
lipdata$acc_treat <- factor(lipdata$acc_treat, levels=c('Before Acclimation','Warm Acclimated (15d)','Cold Acclimated (15d)'))
lipdata$diapause_treat <- factor(lipdata$diapause_treat, levels=c('Non-diapause','Diapause'))

lipdata$totaltreatment <- paste(lipdata$diapause_treat, lipdata$acc_treat)
lipdata$totaltreatment_short <- as.factor(paste(substr(lipdata$totaltreatment,1,1),
                                                word(lipdata$totaltreatment,2)))

# add parental cage info
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Input Files")
lipdata_parentalcageinfo <- read.csv("SampleInfo_Fall2024.csv",header = T)
lipdata_parentalcageinfo <- lipdata_parentalcageinfo[c(1:36),c(1,2)]
lipdata <- merge(lipdata,lipdata_parentalcageinfo,by="Label")
lipdata$parental_cage_num <- as.character(lipdata$parental_cage_num)

################### now read in fatty acyl data and compute summary variables
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Data Files for Submission")
fa_data <- read.csv(file="FattyAcyls_molarprops_RawData.csv")
head(fa_data)

## rename factors 
fa_data$diapause_treat[fa_data$diapause_treat == 'SD'] <- 'Diapause'
fa_data$diapause_treat[fa_data$diapause_treat == 'LD'] <- 'Non-diapause'
fa_data$acc_treat[fa_data$acc_treat == 'PREACC'] <- 'Before Acclimation'
fa_data$acc_treat[fa_data$acc_treat == 'WARM'] <- 'Warm Acclimated (15d)'
fa_data$acc_treat[fa_data$acc_treat == 'COLD'] <- 'Cold Acclimated (15d)'
fa_data$acc_treat <- factor(fa_data$acc_treat, levels=c('Before Acclimation','Warm Acclimated (15d)','Cold Acclimated (15d)'))
fa_data$diapause_treat <- factor(fa_data$diapause_treat, levels=c('Non-diapause','Diapause'))

## add parental cage info
fa_data <- merge(fa_data, lipdata_parentalcageinfo, by = "Label")
fa_data$parental_cage_num <- as.character(fa_data$parental_cage_num)

fa_data$totaltreatment <- paste(fa_data$diapause_treat, fa_data$acc_treat)
fa_data$totaltreatment_short <- as.factor(paste(substr(fa_data$totaltreatment,1,1),
                                                word(fa_data$totaltreatment,2)))

### run loop to choose best transformation for each PL class
outcome_vars <- colnames(lipdata[,c(6:15)])
# Initialize an empty data frame to store results
results_assumptions_lip <- data.frame(Outcome = character(),
                                  raw_shapiro_p = numeric(),
                                  log_shapiro_p = numeric(),
                                  sqrt_shapiro_p = numeric(),
                                  stringsAsFactors = FALSE)

for (outcome in outcome_vars) {
  # Fit the linear model
  outcome_log <- log(lipdata[outcome])
  outcome_sqrt <- sqrt(lipdata[outcome])
  
  formula_raw <- as.formula(paste(outcome, "~ totaltreatment + (1|parental_cage_num)"))
  model_raw <- lmer(formula_raw, data = lipdata)
  shapiro_p_raw <- shapiro.test(resid(model_raw))$p.value
  
  formula_log <- as.formula(paste(outcome_log, "~ totaltreatment + (1|parental_cage_num)"))
  model_log <- lmer(formula_log, data = lipdata)
  shapiro_p_log <- shapiro.test(resid(model_log))$p.value
  
  formula_sqrt <- as.formula(paste(outcome_sqrt, "~ totaltreatment + (1|parental_cage_num)"))
  model_sqrt <- lmer(formula_sqrt, data = lipdata)
  shapiro_p_sqrt <- shapiro.test(resid(model_sqrt))$p.value
  
  results_assumptions_lip <- rbind(results_assumptions_lip, data.frame(Outcome = outcome,
                                                               raw_shapiro_p = shapiro_p_raw,
                                                               log_shapiro_p = shapiro_p_log,
                                                               sqrt_shapiro_p = shapiro_p_sqrt))
  
}

results_assumptions_lip_melt <- reshape2::melt(results_assumptions_lip,
                                           id.vars = "Outcome", 
                                           variable.name = "Transformation", 
                                           value.name = "p_value")
View(results_assumptions_lip_melt)

## choose the transformation for each PL class that gives the largest p-value from shapiro test
# (i.e., the transformation that leads to the most normal residuals)
final_transformations_lip <- results_assumptions_lip_melt %>%
  group_by(Outcome) %>%
  top_n(1, p_value) %>%
  ungroup()
View(final_transformations_lip)


#####################################################################################

## make averages table for each PL class to include as supplementary table
lip_averages <- lipdata %>% 
  group_by(totaltreatment) %>% 
  summarise(across(PG:LPC, mean)) %>%
  as.data.frame()
lip_averages
View(lipdata)

# Define the standard error function
se <- function(x, na.rm = FALSE) {
  if (na.rm) {
    x <- na.omit(x)
  }
  2*(sd(x) / sqrt(length(x)))
}

lip_averages <- lipdata %>% 
  group_by(totaltreatment_short) %>%
  summarise(across(c(PG:LPC),
                   list(mean = mean, two_se = se))) %>%
  mutate(across(c("PG_mean":"LPC_two_se"), round, 2))

lip_averages <- lip_averages %>%
  unite("PG", PG_mean, PG_two_se, sep = "±") %>%
  unite("PI", PI_mean, PI_two_se, sep = "±") %>%
  unite("CL", CL_mean, CL_two_se, sep = "±") %>%
  unite("PE", PE_mean, PE_two_se, sep = "±") %>%
  unite("PS", PS_mean, PS_two_se, sep = "±") %>%
  unite("SM_PE", SM_PE_mean, SM_PE_two_se, sep = "±") %>%  
  unite("LPE", LPE_mean, LPE_two_se, sep = "±") %>%
  unite("PC", PC_mean, PC_two_se, sep = "±") %>%
  unite("SM", SM_mean, SM_two_se, sep = "±") %>%
  unite("LPC", LPC_mean, LPC_two_se, sep = "±") %>%
  t() %>%
  data.frame()
  
colnames(lip_averages) <- lip_averages[1,]
lip_averages <- lip_averages[c(2:11),]

View(lip_averages)

lip_order <- lipdata %>%
  summarise(across(c(PG:LPC), mean)) %>%
  t() %>%
  as.data.frame() %>%
  arrange(-V1)

lip_averages_ordered <- lip_averages[rownames(lip_order), ]
lip_averages_ordered <- lip_averages_ordered[, c(4,6,5,1,3,2)] 
lip_averages_ordered <- lip_averages_ordered[c(1,2,3,4,5,7,8,9,10,11,6),] 

## save to csv
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Output Files")
write.csv(file = "lipid_headgroup_averages.csv", lip_averages_ordered , row.names = T, col.names = T)



### now calculate log2FC of each PL class for each of the five pairwise contrasts
outcome_vars <- colnames(lipdata[,c(6:15)])
outcome_vars
# Initialize an empty data frame to store results
results_log2fc_lip <- data.frame(headgroup = as.character(),
                             Ncold_Nwarm_log2fc = numeric(), 
                             Dcold_Dwarm_log2fc = numeric(),
                             Dbefore_Nbefore_log2fc = numeric(),
                             Dwarm_Nwarm_log2fc = numeric(),
                             Dcold_Ncold_log2fc = numeric(),
                             stringsAsFactors = FALSE)

for (lip in outcome_vars) {
  mean_Ncold <- lipdata %>%
    filter(totaltreatment == "Non-diapause Cold Acclimated (15d)") %>%
    dplyr::select(lip) %>%
    colMeans()
  mean_Nwarm <- lipdata %>%
    filter(totaltreatment == "Non-diapause Warm Acclimated (15d)") %>%
    dplyr::select(lip) %>%
    colMeans()
  mean_Dcold <- lipdata %>%
    filter(totaltreatment == "Diapause Cold Acclimated (15d)") %>%
    dplyr::select(lip) %>%
    colMeans()
  mean_Dwarm <- lipdata %>%
    filter(totaltreatment == "Diapause Warm Acclimated (15d)") %>%
    dplyr::select(lip) %>%
    colMeans()
  mean_Dbefore <- lipdata %>%
    filter(totaltreatment == "Diapause Before Acclimation") %>%
    dplyr::select(lip) %>%
    colMeans()
  mean_Nbefore <- lipdata %>%
    filter(totaltreatment == "Non-diapause Before Acclimation") %>%
    dplyr::select(lip) %>%
    colMeans()
  results_log2fc_lip <- rbind(results_log2fc_lip, data.frame(
    headgroup = lip,
    Ncold_Nwarm_log2fc = log2(mean_Ncold/mean_Nwarm), 
    Dcold_Dwarm_log2fc = log2(mean_Dcold/mean_Dwarm),
    Dbefore_Nbefore_log2fc = log2(mean_Dbefore/mean_Nbefore),
    Dwarm_Nwarm_log2fc = log2(mean_Dwarm/mean_Nwarm),
    Dcold_Ncold_log2fc = log2(mean_Dcold/mean_Ncold)
  ))
}

## run dummy model to get the coefficient numbers for running contrasts of interest
model <- lmer(PG ~ totaltreatment + (1|parental_cage_num) , data = lipdata)
emm_model <- emmeans(model, "totaltreatment")
coef(pairs(emm_model))

## these are contrasts
con <- list(
  Ncold_Nwarm = c(0,0,0,0,1,-1), 
  Dcold_Dwarm = c(0,1,-1,0,0,0),
  
  Dbefore_Nbefore = c(1,0,0,-1,0,0),
  Dwarm_Nwarm = c(0,0,1,0,0,-1),
  Dcold_Ncold = c(0,1,0,0,-1,0)
)
emmeans(model, list(~ totaltreatment), contr = con, adjust = "BH")


#### model and extract p-values for pairwise comparisons of interest 
## using transformations that worked best to model full dataset
outcome_vars <- colnames(lipdata[,c(6:15)])
outcome_vars
# Initialize an empty data frame to store results
results_padj_lip <- data.frame(headgroup = character(),
                           overall_chisquared = numeric(),
                           overall_df = numeric(),
                           overall_pvalue = numeric(),
                           Ncold_Nwarm_padj = numeric(), 
                           Dcold_Dwarm_padj = numeric(),
                           Dbefore_Nbefore_padj = numeric(),
                           Dwarm_Nwarm_padj = numeric(),
                           Dcold_Ncold_padj = numeric(),
                           stringsAsFactors = FALSE)

# Loop over each outcome variable and fit the linear model
for (outcome in outcome_vars) {
  transformation <- subset(final_transformations_lip,final_transformations_lip$Outcome == outcome)$Transformation
  
  if (transformation == "raw_shapiro_p") {
    data_tr <- lipdata[outcome]
  }
  
  if (transformation == "log_shapiro_p") {
    data_tr <- log(lipdata[outcome])
  }
  
  if (transformation == "sqrt_shapiro_p") {
    data_tr <- sqrt(lipdata[outcome])
  }
  
  # Fit the linear model
  formula <- as.formula(paste(data_tr, "~ totaltreatment + (1|parental_cage_num)"))
  model <- lmer(formula, data = lipdata)
  
  # Perform planned contrasts to get the significance of each contrast
  contrast_list <- as.data.frame(emmeans(model, list(~ totaltreatment), contr = con, adjust = "BH")$` of totaltreatment`)
  results_padj_lip <- rbind(results_padj_lip, data.frame(
    headgroup = outcome,
    overall_chisquared = Anova(model)[1,1],
    overall_df = Anova(model)[1,2],
    overall_pvalue = Anova(model)[1,3],
    Ncold_Nwarm_padj = contrast_list[1,"p.value"], 
    Dcold_Dwarm_padj = contrast_list[2,"p.value"],
    Dbefore_Nbefore_padj = contrast_list[3,"p.value"],
    Dwarm_Nwarm_padj = contrast_list[4,"p.value"],
    Dcold_Ncold_padj = contrast_list[5,"p.value"]))
}

## get rid of metabolites that violated normality assumption
nonnormal_hg <- final_transformations_lip %>%
  filter(p_value < 0.05) %>%
  dplyr::select("Outcome") 
results_padj_lip <- results_padj_lip %>%
  filter(!headgroup %in% nonnormal_hg$Outcome)

### adjust the p-values for the many metabolites tested 
results_padj_lip_adj = results_padj_lip %>% 
  mutate(across(overall_pvalue:Dcold_Ncold_padj, ~ p.adjust(.x,"BH")))

########################## same thing but with non-parametric test for metabolites that violated
## normality assumption

## write a loop to do this 
outcome_vars <- nonnormal_hg$Outcome
# Initialize an empty data frame to store results
results_padj_wilcox_lip <- data.frame(headgroup = character(),
                                  overall_chisquared = numeric(),
                                  overall_df = numeric(),
                                  overall_pvalue = numeric(),
                                  overall_chisquared = numeric(),
                                  Ncold_Nwarm_padj = numeric(), 
                                  Dcold_Dwarm_padj = numeric(),
                                  Dbefore_Nbefore_padj = numeric(),
                                  Dwarm_Nwarm_padj = numeric(),
                                  Dcold_Ncold_padj = numeric(),
                                  stringsAsFactors = FALSE)

# Loop over each outcome variable and fit the linear model
for (outcome in outcome_vars) {
  formula <- as.formula(paste(outcome, "~totaltreatment"))
  
  kruskal_testresults <- kruskal.test(formula, data = lipdata)
  
  Ncold_Nwarm_padj <- wilcox.test(formula, 
                                  data = lipdata[lipdata$totaltreatment=="Non-diapause Cold Acclimated (15d)"|lipdata$totaltreatment=="Non-diapause Warm Acclimated (15d)",])$p.value
  Dcold_Dwarm_padj <- wilcox.test(formula, 
                                  data = lipdata[lipdata$totaltreatment=="Diapause Cold Acclimated (15d)"|lipdata$totaltreatment=="Diapause Warm Acclimated (15d)",])$p.value
  Dbefore_Nbefore_padj <- wilcox.test(formula, 
                                      data = lipdata[lipdata$totaltreatment=="Diapause Before Acclimation"|lipdata$totaltreatment=="Non-diapause Before Acclimation",])$p.value
  
  Dwarm_Nwarm_padj <- wilcox.test(formula, 
                                  data = lipdata[lipdata$totaltreatment=="Diapause Warm Acclimated (15d)"|lipdata$totaltreatment=="Non-diapause Warm Acclimated (15d)",])$p.value
  
  Dcold_Ncold_padj <- wilcox.test(formula, 
                                  data = lipdata[lipdata$totaltreatment=="Diapause Cold Acclimated (15d)"|lipdata$totaltreatment=="Non-diapause Cold Acclimated (15d)",])$p.value
  
  results_padj_wilcox_lip <- rbind(results_padj_wilcox_lip, data.frame(
    headgroup = outcome,
    overall_chisquared = kruskal_testresults$statistic,
    overall_df = kruskal_testresults$parameter,
    overall_pvalue = kruskal_testresults$p.value,
    Ncold_Nwarm_padj = Ncold_Nwarm_padj, 
    Dcold_Dwarm_padj = Dcold_Dwarm_padj,
    Dbefore_Nbefore_padj = Dbefore_Nbefore_padj,
    Dwarm_Nwarm_padj = Dwarm_Nwarm_padj,
    Dcold_Ncold_padj = Dcold_Ncold_padj))
}


### adjust the p-values for multiple metabolites tested
results_padj_wilcox_lip_adj = results_padj_wilcox_lip %>% 
  mutate(across(Ncold_Nwarm_padj:Dcold_Ncold_padj, ~ p.adjust(.x,"BH")))

## attach results for parametric and non-parametric tests
results_padj_adj_final_lip <- rbind(results_padj_lip_adj,results_padj_wilcox_lip_adj)

## add log2FC for each pairwise contrast for each metabolite
final_lip_data <- merge(results_log2fc_lip,results_padj_adj_final_lip,by="headgroup")

####. save for use in met log2FC plots script
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Output Files")
write.csv(file = "lip_headgroup_pairwise_contrasts_stats.csv", final_lip_data, row.names = F)

## now replace values with either 0 (for more than 0.05) or 1 (for less than 0.05)
## for use in annotation column in heatmap
final_lip_data_binarypvalues <- final_lip_data %>%
  mutate(Ncold_Nwarm_padj = case_when(
    Ncold_Nwarm_padj < 0.05 ~ "Yes",
    Ncold_Nwarm_padj > 0.05 ~ "No")) %>%
  mutate(Dcold_Dwarm_padj = case_when(
    Dcold_Dwarm_padj < 0.05 ~ "Yes",
    Dcold_Dwarm_padj > 0.05 ~ "No")) %>%
  mutate(Dbefore_Nbefore_padj = case_when(
    Dbefore_Nbefore_padj < 0.05 ~ "Yes",
    Dbefore_Nbefore_padj > 0.05 ~ "No")) %>%
  mutate(Dwarm_Nwarm_padj = case_when(
    Dwarm_Nwarm_padj < 0.05 ~ "Yes",
    Dwarm_Nwarm_padj > 0.05 ~ "No")) %>%
  mutate(Dcold_Ncold_padj = case_when(
    Dcold_Ncold_padj < 0.05 ~ "Yes",
    Dcold_Ncold_padj > 0.05 ~ "No"))


## change lipids to cleaned up version
final_lip_data_binarypvalues$headgroup <- final_lip_data_binarypvalues$headgroup %>%
  gsub("_","",.)

## only keep metabolite column and p-adjusted columns
final_lip_data_padj <- final_lip_data_binarypvalues[,c(1,10:14)]

colnames(final_lip_data_padj) <- c("Outcome","ND Cold - ND Warm","D Cold - D Warm",
                                   "D Before - ND Before", "D Warm - ND Warm",
                                   "D Cold - ND Cold")
## keep only p-value columns
final_lip_data_padj <- final_lip_data_padj[,c(2:6)]

# add row names that are the headgroup names
row.names(final_lip_data_padj) <- final_lip_data_binarypvalues$headgroup

## change order of columns 
final_lip_data_padj <- final_lip_data_padj[,c(5,4,3,2,1)]


## now z-score all data so that it is scaled to the mean for each PL then the mean is centered to 0
lipdata_zscore <- lipdata
lipdata_zscore[,c(6:15)] <- scale(lipdata_zscore[,c(6:15)], scale = TRUE, center = TRUE)
## make a version of the z-scored data that is group averaged for each of the six treatments for each PL
lipdata_zscore_grouped <- lipdata_zscore %>%
  group_by(totaltreatment) %>%
  summarise(across(!c(no.,Label, diapause_treat,acc_treat, totaltreatment_short, Mass, parental_cage_num), mean))

## change orientation from PL columns and treatment groups rows to
## PL rows and treatment group columns
lipdata_zscore_grouped_long <- lipdata_zscore_grouped %>%
  dplyr::select(2:11) %>%
  t()

## add shortened versions of the treatment groups as column names 
colnames(lipdata_zscore_grouped_long) <- paste(substr(lipdata_zscore_grouped$totaltreatment,1,1),
                                               word(lipdata_zscore_grouped$totaltreatment,2))
colnames(lipdata_zscore_grouped_long)

## change order of treatment group columns
lipdata_zscore_grouped_long <- lipdata_zscore_grouped_long[, c(4,6,5,1,3,2)] 

## convert from data_frame to matrix 
lip_matrix_grouped <- lipdata_zscore_grouped_long %>% 
  as.matrix()

## plot initial heatmap using matrix of zscores averaged across treatment groups
pheatmap(lip_matrix_grouped,cluster_cols = F,
         cutree_rows = 5)

## add abbreviated PL names to the PL matrix generated above
row.names(lip_matrix_grouped) <- row.names(lip_matrix_grouped) %>%
  gsub("_","",.) 

lip_matrix_grouped
## provide info on color pattern for annotation columns
my_colour_pairwise <- list(
  'ND Cold - ND Warm' = c(Yes = "black", No = "white"),
  'D Cold - D Warm' = c(Yes = "black", No = "white"),
  'D Before - ND Before' = c(Yes = "darkgrey", No = "white"),
  'D Warm - ND Warm' = c(Yes = "darkgrey", No = "white"),
  'D Cold - ND Cold' = c(Yes = "darkgrey", No = "white")
)

## plot heat map
heatmap_lip_pairwisepadj <- pheatmap(lip_matrix_grouped,
                                     cluster_cols = F,
                                     #color = myColor,
                                     cutree_rows = 4,
                                     annotation_row = final_lip_data_padj,
                                     annotation_colors = my_colour_pairwise,
                                     angle_col = 45,
                                     annotation_legend = F,
                                     fontsize = 10,
                                     fontsize_row = 12,
                                     fontsize_col = 12,
                                     cellwidth = 25,
                                     cellheight = 12,
                                     gaps_col = 3,
                                     border_color = "black")

setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Figures/Lipidomics/")

ggsave("heatmap_lip_headgroup_emmeans_pairwisepadj_nosummaryv.png",plot=heatmap_lip_pairwisepadj,dpi=600,units='in',
       width=9,height=4)

#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
## now same for all that but with FA dataset 
### run loop to choose best transformation for each metabolite
outcome_vars <- colnames(fa_data[,c(6:23)])
outcome_vars
# Initialize an empty data frame to store results
results_assumptions_fa <- data.frame(Outcome = character(),
                                      raw_shapiro_p = numeric(),
                                      log_shapiro_p = numeric(),
                                      sqrt_shapiro_p = numeric(),
                                      stringsAsFactors = FALSE)

for (outcome in outcome_vars) {
  # Fit the linear model
  outcome_log <- log(fa_data[outcome])
  outcome_sqrt <- sqrt(fa_data[outcome])
  
  formula_raw <- as.formula(paste(outcome, "~ totaltreatment + (1|parental_cage_num)"))
  model_raw <- lmer(formula_raw, data = fa_data)
  shapiro_p_raw <- shapiro.test(resid(model_raw))$p.value
  
  formula_log <- as.formula(paste(outcome_log, "~ totaltreatment + (1|parental_cage_num)"))
  model_log <- lmer(formula_log, data = fa_data)
  shapiro_p_log <- shapiro.test(resid(model_log))$p.value
  
  formula_sqrt <- as.formula(paste(outcome_sqrt, "~ totaltreatment + (1|parental_cage_num)"))
  model_sqrt <- lmer(formula_sqrt, data = fa_data)
  shapiro_p_sqrt <- shapiro.test(resid(model_sqrt))$p.value
  
  results_assumptions_fa <- rbind(results_assumptions_fa, data.frame(Outcome = outcome,
                                                                       raw_shapiro_p = shapiro_p_raw,
                                                                       log_shapiro_p = shapiro_p_log,
                                                                       sqrt_shapiro_p = shapiro_p_sqrt))
  
}

results_assumptions_fa_melt <- reshape2::melt(results_assumptions_fa,
                                               id.vars = "Outcome", 
                                               variable.name = "Transformation", 
                                               value.name = "p_value")

## choose the transformation for each metabolite that gives the largest p-value from shapiro test
# (i.e., the transformation that leads to the most normal residuals)
final_transformations_fa <- results_assumptions_fa_melt %>%
  group_by(Outcome) %>%
  top_n(1, p_value) %>%
  ungroup()


#####################################################################################
View(fa_data[,c(6:23)])
## calculate average +/- 2SE for each FAME in each of the six total treatment group combinations for supplementary table
fa_averages <- fa_data[,c(25,6:23)] %>% 
  group_by(totaltreatment) %>% 
  summarise(across(C12.0:C20.5n3, mean)) %>%
  as.data.frame()
fa_averages

fa_averages <- fa_data %>% 
  group_by(totaltreatment_short) %>%
  summarise(across(c(C12.0:C20.5n3),
                   list(mean = mean, two_se = se))) %>%
  mutate(across(c("C12.0_mean":"C20.5n3_two_se"), round, 2))

colnames(fa_averages)

fa_averages <- fa_averages %>%
  unite("C12.0", C12.0_mean, C12.0_two_se, sep = "±") %>%
  unite("C13.0", C13.0_mean, C13.0_two_se, sep = "±") %>%
  unite("C14.0", C14.0_mean, C14.0_two_se, sep = "±") %>%
  unite("C14.1n5", C14.1n5_mean, C14.1n5_two_se, sep = "±") %>%
  unite("C15.0", C15.0_mean, C15.0_two_se, sep = "±") %>%
  unite("C16.0", C16.0_mean, C16.0_two_se, sep = "±") %>%
  unite("C16.1n7", C16.1n7_mean, C16.1n7_two_se, sep = "±") %>%
  unite("C17.0", C17.0_mean, C17.0_two_se, sep = "±") %>%
  unite("C18.0", C18.0_mean, C18.0_two_se, sep = "±") %>%
  unite("C18.1n9c", C18.1n9c_mean, C18.1n9c_two_se, sep = "±") %>%
  unite("C18.2n6c", C18.2n6c_mean, C18.2n6c_two_se, sep = "±") %>%
  unite("C20.0", C20.0_mean, C20.0_two_se, sep = "±") %>%  
  unite("C18.3n6", C18.3n6_mean, C18.3n6_two_se, sep = "±") %>%  
  unite("C18.3n3", C18.3n3_mean, C18.3n3_two_se, sep = "±") %>%  
  unite("C22.0", C22.0_mean, C22.0_two_se, sep = "±") %>%  
  unite("C20.3n6", C20.3n6_mean, C20.3n6_two_se, sep = "±") %>%  
  unite("C20.4n6", C20.4n6_mean, C20.4n6_two_se, sep = "±") %>%  
  unite("C20.5n3", C20.5n3_mean, C20.5n3_two_se, sep = "±") %>% 
  t() %>%
  data.frame()

colnames(fa_averages) <- fa_averages[1,]
fa_averages <- fa_averages[c(2:23),]

fa_order <- fa_data %>%
  summarise(across(c(C12.0:C20.5n3), mean)) %>%
  t() %>%
  as.data.frame() %>%
  arrange(-V1)

fa_averages_ordered <- fa_averages[rownames(fa_order), ]
View(fa_averages_ordered)
fa_averages_ordered <- fa_averages_ordered[, c(4,6,5,1,3,2)] 

## save to csv
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Output Files")
write.csv(file = "lipid_fa_averages.csv", fa_averages_ordered , row.names = T, col.names = T)


############## now generate log2FC for each of the five planned pairwise contrasts
outcome_vars <- colnames(fa_data[,c(6:23)])
outcome_vars
# Initialize an empty data frame to store results
results_log2fc_fa <- data.frame(Outcome = as.character(),
                                 Ncold_Nwarm_log2fc = numeric(), 
                                 Dcold_Dwarm_log2fc = numeric(),
                                 Dbefore_Nbefore_log2fc = numeric(),
                                 Dwarm_Nwarm_log2fc = numeric(),
                                 Dcold_Ncold_log2fc = numeric(),
                                 stringsAsFactors = FALSE)

for (fa in outcome_vars) {
  mean_Ncold <- fa_data %>%
    filter(totaltreatment == "Non-diapause Cold Acclimated (15d)") %>%
    dplyr::select(fa) %>%
    colMeans()
  mean_Nwarm <- fa_data %>%
    filter(totaltreatment == "Non-diapause Warm Acclimated (15d)") %>%
    dplyr::select(fa) %>%
    colMeans()
  mean_Dcold <- fa_data %>%
    filter(totaltreatment == "Diapause Cold Acclimated (15d)") %>%
    dplyr::select(fa) %>%
    colMeans()
  mean_Dwarm <- fa_data %>%
    filter(totaltreatment == "Diapause Warm Acclimated (15d)") %>%
    dplyr::select(fa) %>%
    colMeans()
  mean_Dbefore <- fa_data %>%
    filter(totaltreatment == "Diapause Before Acclimation") %>%
    dplyr::select(fa) %>%
    colMeans()
  mean_Nbefore <- fa_data %>%
    filter(totaltreatment == "Non-diapause Before Acclimation") %>%
    dplyr::select(fa) %>%
    colMeans()
  results_log2fc_fa <- rbind(results_log2fc_fa, data.frame(
    Outcome = fa,
    Ncold_Nwarm_log2fc = log2(mean_Ncold/mean_Nwarm), 
    Dcold_Dwarm_log2fc = log2(mean_Dcold/mean_Dwarm),
    Dbefore_Nbefore_log2fc = log2(mean_Dbefore/mean_Nbefore),
    Dwarm_Nwarm_log2fc = log2(mean_Dwarm/mean_Nwarm),
    Dcold_Ncold_log2fc = log2(mean_Dcold/mean_Ncold)
  ))
}

## these are contrasts for the totaltreatment version of the model
con <- list(
  Ncold_Nwarm = c(0,0,0,0,1,-1), 
  Dcold_Dwarm = c(0,1,-1,0,0,0),
  
  Dbefore_Nbefore = c(1,0,0,-1,0,0),
  Dwarm_Nwarm = c(0,0,1,0,0,-1),
  Dcold_Ncold = c(0,1,0,0,-1,0)
)

#### model and extract p-values for pairwise comparisons of interest 
## using transformations that worked best to model full dataset
outcome_vars <- colnames(fa_data[,c(6:23)])
outcome_vars
# Initialize an empty data frame to store results
results_padj_fa <- data.frame(Outcome = character(),
                               overall_chisquared = numeric(),
                               overall_df = numeric(),
                               overall_pvalue = numeric(),
                               Ncold_Nwarm_padj = numeric(), 
                               Dcold_Dwarm_padj = numeric(),
                               Dbefore_Nbefore_padj = numeric(),
                               Dwarm_Nwarm_padj = numeric(),
                               Dcold_Ncold_padj = numeric(),
                               stringsAsFactors = FALSE)

# Loop over each outcome variable and fit the linear model
for (outcome in outcome_vars) {
  transformation <- subset(final_transformations_fa,final_transformations_fa$Outcome == outcome)$Transformation
  
  if (transformation == "raw_shapiro_p") {
    data_tr <- fa_data[outcome]
  }
  
  if (transformation == "log_shapiro_p") {
    data_tr <- log(fa_data[outcome])
  }
  
  if (transformation == "sqrt_shapiro_p") {
    data_tr <- sqrt(fa_data[outcome])
  }
  
  # Fit the linear model
  formula <- as.formula(paste(data_tr, "~ totaltreatment + (1|parental_cage_num)"))
  model <- lmer(formula, data = fa_data)
  
  # Perform planned contrasts to get the significance of each contrast
  contrast_list <- as.data.frame(emmeans(model, list(~ totaltreatment), contr = con, adjust = "BH")$` of totaltreatment`)
  results_padj_fa <- rbind(results_padj_fa, data.frame(
    Outcome = outcome,
    overall_chisquared = Anova(model)[1,1],
    overall_df = Anova(model)[1,2],
    overall_pvalue = Anova(model)[1,3],
    Ncold_Nwarm_padj = contrast_list[1,"p.value"], 
    Dcold_Dwarm_padj = contrast_list[2,"p.value"],
    Dbefore_Nbefore_padj = contrast_list[3,"p.value"],
    Dwarm_Nwarm_padj = contrast_list[4,"p.value"],
    Dcold_Ncold_padj = contrast_list[5,"p.value"]))
}

## get rid of species that violated normality assumption
nonnormal_fa <- final_transformations_fa %>%
  filter(p_value < 0.05) %>%
  dplyr::select("Outcome") 

results_padj_fa <- results_padj_fa %>%
  filter(!Outcome %in% nonnormal_fa$Outcome)

### adjust the p-values for the many species tested 
results_padj_fa_adj = results_padj_fa %>% 
  mutate(across(overall_pvalue:Dcold_Ncold_padj, ~ p.adjust(.x,"BH")))

########################## same thing but with non-parametric test for species that violated
## normality assumption

## write a loop to do this 
outcome_vars <- nonnormal_fa$Outcome
# Initialize an empty data frame to store results
results_padj_wilcox_fa <- data.frame(Outcome = character(),
                                      overall_chisquared = numeric(),
                                      overall_df = numeric(),
                                      overall_pvalue = numeric(),
                                      overall_chisquared = numeric(),
                                      Ncold_Nwarm_padj = numeric(), 
                                      Dcold_Dwarm_padj = numeric(),
                                      Dbefore_Nbefore_padj = numeric(),
                                      Dwarm_Nwarm_padj = numeric(),
                                      Dcold_Ncold_padj = numeric(),
                                      stringsAsFactors = FALSE)

# Loop over each outcome variable and fit the linear model
for (outcome in outcome_vars) {
  formula <- as.formula(paste(outcome, "~totaltreatment"))
  
  kruskal_testresults <- kruskal.test(formula, data = fa_data)
  
  Ncold_Nwarm_padj <- wilcox.test(formula, 
                                  data = fa_data[fa_data$totaltreatment=="Non-diapause Cold Acclimated (15d)"|fa_data$totaltreatment=="Non-diapause Warm Acclimated (15d)",])$p.value
  Dcold_Dwarm_padj <- wilcox.test(formula, 
                                  data = fa_data[fa_data$totaltreatment=="Diapause Cold Acclimated (15d)"|fa_data$totaltreatment=="Diapause Warm Acclimated (15d)",])$p.value
  Dbefore_Nbefore_padj <- wilcox.test(formula, 
                                      data = fa_data[fa_data$totaltreatment=="Diapause Before Acclimation"|fa_data$totaltreatment=="Non-diapause Before Acclimation",])$p.value
  
  Dwarm_Nwarm_padj <- wilcox.test(formula, 
                                  data = fa_data[fa_data$totaltreatment=="Diapause Warm Acclimated (15d)"|fa_data$totaltreatment=="Non-diapause Warm Acclimated (15d)",])$p.value
  
  Dcold_Ncold_padj <- wilcox.test(formula, 
                                  data = fa_data[fa_data$totaltreatment=="Diapause Cold Acclimated (15d)"|fa_data$totaltreatment=="Non-diapause Cold Acclimated (15d)",])$p.value
  
  results_padj_wilcox_fa <- rbind(results_padj_wilcox_fa, data.frame(
    Outcome = outcome,
    overall_chisquared = kruskal_testresults$statistic,
    overall_df = kruskal_testresults$parameter,
    overall_pvalue = kruskal_testresults$p.value,
    Ncold_Nwarm_padj = Ncold_Nwarm_padj, 
    Dcold_Dwarm_padj = Dcold_Dwarm_padj,
    Dbefore_Nbefore_padj = Dbefore_Nbefore_padj,
    Dwarm_Nwarm_padj = Dwarm_Nwarm_padj,
    Dcold_Ncold_padj = Dcold_Ncold_padj))
}

View(results_padj_wilcox_fa)

### adjust the p-values for multiple metabolites tested
results_padj_wilcox_fa_adj = results_padj_wilcox_fa %>% 
  mutate(across(Ncold_Nwarm_padj:Dcold_Ncold_padj, ~ p.adjust(.x,"BH")))

## attach results for parametric and non-parametric tests
results_padj_adj_final_fa <- rbind(results_padj_fa_adj,results_padj_wilcox_fa_adj)

## add log2FC for each pairwise contrast for each metabolite
final_fa_data <- merge(results_log2fc_fa,results_padj_adj_final_fa,by="Outcome")
View(final_fa_data)
####. save for use in met log2FC plots script
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Output Files")
write.csv(file = "lip_fa_pairwise_contrasts_stats.csv", final_fa_data, row.names = F)


## now replace values with either 0 (for more than 0.05) or 1 (for less than 0.05)
## for use in annotation column in heatmap
final_fa_data_binarypvalues <- final_fa_data %>%
  mutate(Ncold_Nwarm_padj = case_when(
    Ncold_Nwarm_padj < 0.05 ~ "Yes",
    Ncold_Nwarm_padj > 0.05 ~ "No")) %>%
  mutate(Dcold_Dwarm_padj = case_when(
    Dcold_Dwarm_padj < 0.05 ~ "Yes",
    Dcold_Dwarm_padj > 0.05 ~ "No")) %>%
  mutate(Dbefore_Nbefore_padj = case_when(
    Dbefore_Nbefore_padj < 0.05 ~ "Yes",
    Dbefore_Nbefore_padj > 0.05 ~ "No")) %>%
  mutate(Dwarm_Nwarm_padj = case_when(
    Dwarm_Nwarm_padj < 0.05 ~ "Yes",
    Dwarm_Nwarm_padj > 0.05 ~ "No")) %>%
  mutate(Dcold_Ncold_padj = case_when(
    Dcold_Ncold_padj < 0.05 ~ "Yes",
    Dcold_Ncold_padj > 0.05 ~ "No"))


## only keep outcome column and p-adjusted columns
View(final_fa_data_binarypvalues)
final_fa_data_padj <- final_fa_data_binarypvalues[,c(1,10:14)]

colnames(final_fa_data_padj) <- c("Outcome","ND Cold - ND Warm","D Cold - D Warm",
                                   "D Before - ND Before", "D Warm - ND Warm",
                                   "D Cold - ND Cold")


## keep only p-value columns
final_fa_data_padj <- final_fa_data_padj[,c(2:6)]

# add row names that are the headgroup names
row.names(final_fa_data_padj) <- final_fa_data_binarypvalues$Outcome

## change order of columns 
final_fa_data_padj <- final_fa_data_padj[,c(5,4,3,2,1)]
colnames(final_fa_data_padj)

## provide info on color pattern for annotation columns
my_colour_pairwise <- list(
  'ND Cold - ND Warm' = c(Yes = "black", No = "white"),
  'D Cold - D Warm' = c(Yes = "black", No = "white"),
  'D Before - ND Before' = c(Yes = "darkgrey", No = "white"),
  'D Warm - ND Warm' = c(Yes = "darkgrey", No = "white"),
  'D Cold - ND Cold' = c(Yes = "darkgrey", No = "white")
)


## now z-score all data so that it is scaled to the mean for each metabolite then the mean is centered to 0
fa_data_zscore <- fa_data
fa_data_zscore[,c(6:23)] <- scale(fa_data_zscore[,c(6:23)], scale = TRUE, center = TRUE)
View(fa_data_zscore)
## make a version of the z-scored data that is group averaged for each of the six treatments for each metabolite
fa_data_zscore_grouped <- fa_data_zscore %>%
  group_by(totaltreatment) %>%
  summarise(across(!c(no.,Label, diapause_treat,acc_treat, totaltreatment_short, Mass, parental_cage_num), mean))
View(fa_data_zscore_grouped)

## change orientation from metabolite columns and treatment groups rows to
## metabolite rows and treatment group columns
fa_data_zscore_grouped_long <- fa_data_zscore_grouped %>%
  dplyr::select(2:19) %>%
  t()

## add shortened versions of the treatment groups as column names 
colnames(fa_data_zscore_grouped_long) <- paste(substr(fa_data_zscore_grouped$totaltreatment,1,1),
                                               word(fa_data_zscore_grouped$totaltreatment,2))
colnames(fa_data_zscore_grouped_long)

## change order of treatment group columns
fa_data_zscore_grouped_long <- fa_data_zscore_grouped_long[, c(4,6,5,1,3,2)] 

## convert from data_frame to matrix 
fa_matrix_grouped <- fa_data_zscore_grouped_long %>% 
  as.matrix()

## plot heat map
heatmap_fa_pairwisepadj <- pheatmap(fa_matrix_grouped,
                                     cluster_cols = F,
                                     #color = myColor,
                                     cutree_rows = 7,
                                     annotation_row = final_fa_data_padj,
                                     annotation_colors = my_colour_pairwise,
                                     angle_col = 45,
                                     annotation_legend = F,
                                     fontsize = 10,
                                     fontsize_row = 12,
                                     fontsize_col = 12,
                                     cellwidth = 25,
                                     cellheight = 12,
                                     gaps_col = 3,
                                     border_color = "black")

setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Figures/Lipidomics/")
ggsave("heatmap_lip_fa_emmeans_pairwisepadj_nosummaryv.png",plot=heatmap_fa_pairwisepadj,dpi=600,units='in',
       width=9,height=9)


