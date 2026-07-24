
### In this script, I will test for effects of diapause and cold acclimation on 
# on the relative levels of individual metabolites using a planned set of five
# pairwise comparisons. For metabolites with normally-distributed residuals,
# I will do this using an lmer followed by pairwise t-tests using emmeans.
# For metabolites without normally-distributed residuals, I will run non-parametric
# tests (Kruskal test followed by pairwise wilcox).
# I will use the results of these tests to generate annotation columns for a 
# heatmap that also shows z-scored metabolite values across each of the six total 
# treatment groups.


## load all packages
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


## read in metabolite data and sample info (sample info has data on parental cage ID per sample)
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Data Files for Submission")
metdata <- read.csv(file="Metabolite_RawData.csv")
metdata_parentalcageinfo <- read.csv(file="SampleInfo_Fall2024.csv")


## remove columns where any of the samples have an NF
metdata <- metdata[ , colSums(is.na(metdata)) == 0]

## first need to scale each metabolite value to the starting mass of each sample
metdata <- metdata%>%mutate_at(vars(X2.Oxoglutaric_acid:Valine),list(MS=~./Mass..mg.))
#metdata <- metdata%>%mutate_at(vars(X2_and_3.phosphoglyceric_acid:Valine),list(MS=~./Mass)) ## when including NA value metabolites

## only keep mass-scaled data columns
metdata_scaled <- metdata[,c(1,2,3,4,51:95)]

metdata_scaled$totaltreatment <- paste(metdata_scaled$diapause_treat, metdata_scaled$acc_treat)

## add info about parental cage number to metdata_scaled dataframe
metdata_scaled <- merge(metdata_scaled,metdata_parentalcageinfo,by="Label")
table_parentalcage <- metdata_scaled %>% count(diapause_treat,parental_cage_num)
metdata_scaled$parental_cage_num <- as.character(metdata_scaled$parental_cage_num)


##############################################################################
##############################################################################
##############################################################################

### Below I used emmeans to run t-tests on the estimated marginal means
## but in cases where metabolites do not have normal residuals, I run non-parametric wilcox tests instead

### run loop to choose best transformation for each metabolite
outcome_vars <- colnames(metdata_scaled[,c(5:49)])
outcome_vars
# Initialize an empty data frame to store results
results_assumptions <- data.frame(Outcome = character(),
                      raw_shapiro_p = numeric(),
                      log_shapiro_p = numeric(),
                      sqrt_shapiro_p = numeric(),
                      stringsAsFactors = FALSE)

for (outcome in outcome_vars) {
  # Fit the linear model
  outcome_log <- log(metdata_scaled[outcome])
  outcome_sqrt <- sqrt(metdata_scaled[outcome])
  
  formula_raw <- as.formula(paste(outcome, "~ totaltreatment + (1|parental_cage_num)"))
  model_raw <- lmer(formula_raw, data = metdata_scaled)
  shapiro_p_raw <- shapiro.test(resid(model_raw))$p.value
  
  formula_log <- as.formula(paste(outcome_log, "~ totaltreatment + (1|parental_cage_num)"))
  model_log <- lmer(formula_log, data = metdata_scaled)
  shapiro_p_log <- shapiro.test(resid(model_log))$p.value
  
  formula_sqrt <- as.formula(paste(outcome_sqrt, "~ totaltreatment + (1|parental_cage_num)"))
  model_sqrt <- lmer(formula_sqrt, data = metdata_scaled)
  shapiro_p_sqrt <- shapiro.test(resid(model_sqrt))$p.value
  
  results_assumptions <- rbind(results_assumptions, data.frame(Outcome = outcome,
                                       raw_shapiro_p = shapiro_p_raw,
                                       log_shapiro_p = shapiro_p_log,
                                       sqrt_shapiro_p = shapiro_p_sqrt))
  
}

results_assumptions_melt <- reshape2::melt(results_assumptions,
                                 id.vars = "Outcome", 
                                 variable.name = "Transformation", 
                                 value.name = "p_value")
View(results_assumptions_melt)

## choose the transformation for each metabolite that gives the largest p-value from shapiro test
# (i.e., the transformation that leads to the most normal residuals)
final_transformations <- results_assumptions_melt %>%
  group_by(Outcome) %>%
  top_n(1, p_value) %>%
  ungroup()
View(final_transformations)


#####################################################################################
#####################################################################################
#### calculate log2FC values in case I want to put those in the annotation column

## make averages table for each metabolite across each treatment group to compute log2FC for contrasts 
met_averages <- metdata_scaled %>% 
  group_by(totaltreatment) %>% 
  summarise(across(X2.Oxoglutaric_acid_MS:Valine_MS, mean)) %>%
  as.data.frame()

outcome_vars <- colnames(metdata_scaled[,c(5:49)])
outcome_vars
# Initialize an empty data frame to store results
results_log2fc <- data.frame(metabolite = as.character(),
                             Ncold_Nwarm_log2fc = numeric(), 
                             Dcold_Dwarm_log2fc = numeric(),
                             Dbefore_Nbefore_log2fc = numeric(),
                             Dwarm_Nwarm_log2fc = numeric(),
                             Dcold_Ncold_log2fc = numeric(),
                             stringsAsFactors = FALSE)

for (met in outcome_vars) {
  mean_Ncold <- metdata_scaled %>%
    filter(totaltreatment == "Non-diapause Cold Acclimated (15d)") %>%
    dplyr::select(met) %>%
    colMeans()
  mean_Nwarm <- metdata_scaled %>%
    filter(totaltreatment == "Non-diapause Warm Acclimated (15d)") %>%
    dplyr::select(met) %>%
    colMeans()
  mean_Dcold <- metdata_scaled %>%
    filter(totaltreatment == "Diapause Cold Acclimated (15d)") %>%
    dplyr::select(met) %>%
    colMeans()
  mean_Dwarm <- metdata_scaled %>%
    filter(totaltreatment == "Diapause Warm Acclimated (15d)") %>%
    dplyr::select(met) %>%
    colMeans()
  mean_Dbefore <- metdata_scaled %>%
    filter(totaltreatment == "Diapause Before Acclimation") %>%
    dplyr::select(met) %>%
    colMeans()
  mean_Nbefore <- metdata_scaled %>%
    filter(totaltreatment == "Non-diapause Before Acclimation") %>%
    dplyr::select(met) %>%
    colMeans()
  results_log2fc <- rbind(results_log2fc, data.frame(
    metabolite = met,
    Ncold_Nwarm_log2fc = log2(mean_Ncold/mean_Nwarm), 
    Dcold_Dwarm_log2fc = log2(mean_Dcold/mean_Dwarm),
    Dbefore_Nbefore_log2fc = log2(mean_Dbefore/mean_Nbefore),
    Dwarm_Nwarm_log2fc = log2(mean_Dwarm/mean_Nwarm),
    Dcold_Ncold_log2fc = log2(mean_Dcold/mean_Ncold)
  ))
}
View(results_log2fc)

## set up dummy model to get the contrast coefficients
model <- lmer(Alanine_MS ~ totaltreatment + (1|parental_cage_num) , data = metdata_scaled)
Anova(model)
emm_model <- emmeans(model, "totaltreatment")
coef(pairs(emm_model))

## these are contrasts for the totaltreatment version of the model
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
outcome_vars <- colnames(metdata_scaled[,c(5:49)])
outcome_vars
# Initialize an empty data frame to store results
results_padj <- data.frame(metabolite = character(),
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
  transformation <- subset(final_transformations,final_transformations$Outcome == outcome)$Transformation
  
  if (transformation == "raw_shapiro_p") {
    data_tr <- metdata_scaled[outcome]
  }
  
  if (transformation == "log_shapiro_p") {
    data_tr <- log(metdata_scaled[outcome])
  }
  
  if (transformation == "sqrt_shapiro_p") {
    data_tr <- sqrt(metdata_scaled[outcome])
  }
  
  # Fit the linear model
  formula <- as.formula(paste(data_tr, "~ totaltreatment + (1|parental_cage_num)"))
  model <- lmer(formula, data = metdata_scaled)
  
  # Perform planned contrasts to get the significance of each contrast
  contrast_list <- as.data.frame(emmeans(model, list(~ totaltreatment), contr = con, adjust = "BH")$` of totaltreatment`)
  results_padj <- rbind(results_padj, data.frame(
    metabolite = outcome,
    overall_chisquared = Anova(model)[1,1],
    overall_df = Anova(model)[1,2],
    overall_pvalue = Anova(model)[1,3],
    Ncold_Nwarm_padj = contrast_list[1,"p.value"], 
    Dcold_Dwarm_padj = contrast_list[2,"p.value"],
    Dbefore_Nbefore_padj = contrast_list[3,"p.value"],
    Dwarm_Nwarm_padj = contrast_list[4,"p.value"],
    Dcold_Ncold_padj = contrast_list[5,"p.value"]))
}

## double check that the results match individual model being run 
model_sep <- lmer(log(Valine_MS) ~ totaltreatment + (1|parental_cage_num), data = metdata_scaled)
Anova(model_sep)

## get rid of metabolites that violated normality assumption
nonnormal_mets <- final_transformations %>%
  filter(p_value < 0.05) %>%
  dplyr::select("Outcome") 

results_padj <- results_padj %>%
  filter(!metabolite %in% nonnormal_mets$Outcome)

### adjust the p-values for the many metabolites tested 
results_padj_adj = results_padj %>% 
  mutate(across(overall_pvalue:Dcold_Ncold_padj, ~ p.adjust(.x,"BH")))

########################## same thing but with non-parametric test for metabolites that violated
## normality assumption

## write a loop to do this 
outcome_vars <- nonnormal_mets$Outcome
outcome_vars
# Initialize an empty data frame to store results
results_padj_wilcox <- data.frame(metabolite = character(),
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
  
  kruskal_testresults <- kruskal.test(formula, data = metdata_scaled)
  
  Ncold_Nwarm_padj <- wilcox.test(formula, 
                                  data = metdata_scaled[metdata_scaled$totaltreatment=="Non-diapause Cold Acclimated (15d)"|metdata_scaled$totaltreatment=="Non-diapause Warm Acclimated (15d)",])$p.value
  Dcold_Dwarm_padj <- wilcox.test(formula, 
                                  data = metdata_scaled[metdata_scaled$totaltreatment=="Diapause Cold Acclimated (15d)"|metdata_scaled$totaltreatment=="Diapause Warm Acclimated (15d)",])$p.value
  Dbefore_Nbefore_padj <- wilcox.test(formula, 
                                      data = metdata_scaled[metdata_scaled$totaltreatment=="Diapause Before Acclimation"|metdata_scaled$totaltreatment=="Non-diapause Before Acclimation",])$p.value
  
  Dwarm_Nwarm_padj <- wilcox.test(formula, 
                                  data = metdata_scaled[metdata_scaled$totaltreatment=="Diapause Warm Acclimated (15d)"|metdata_scaled$totaltreatment=="Non-diapause Warm Acclimated (15d)",])$p.value
  
  Dcold_Ncold_padj <- wilcox.test(formula, 
                                  data = metdata_scaled[metdata_scaled$totaltreatment=="Diapause Cold Acclimated (15d)"|metdata_scaled$totaltreatment=="Non-diapause Cold Acclimated (15d)",])$p.value
  
  results_padj_wilcox <- rbind(results_padj_wilcox, data.frame(
    metabolite = outcome,
    overall_chisquared = kruskal_testresults$statistic,
    overall_df = kruskal_testresults$parameter,
    overall_pvalue = kruskal_testresults$p.value,
    Ncold_Nwarm_padj = Ncold_Nwarm_padj, 
    Dcold_Dwarm_padj = Dcold_Dwarm_padj,
    Dbefore_Nbefore_padj = Dbefore_Nbefore_padj,
    Dwarm_Nwarm_padj = Dwarm_Nwarm_padj,
    Dcold_Ncold_padj = Dcold_Ncold_padj))
}

View(results_padj_wilcox)

### adjust the p-values for multiple metabolites tested
results_padj_wilcox_adj <- results_padj_wilcox %>% 
  mutate(across(Ncold_Nwarm_padj:Dcold_Ncold_padj, ~ p.adjust(.x,"BH")))

## attach results for parametric and non-parametric tests
results_padj_adj_final <- rbind(results_padj_adj,results_padj_wilcox_adj)
View(results_padj_adj_final)

## add log2FC for each pairwise contrast for each metabolite
final_met_data <- merge(results_log2fc,results_padj_adj_final,by="metabolite")
View(final_met_data)

####. save for use in met log2FC plots script
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Output Files")
write.csv(file = "met_pairwise_contrasts_stats.csv", final_met_data, row.names = F)

## now replace values with either 0 (for more than 0.05) or 1 (for less than 0.05)
## for use in annotation column in heatmap
final_met_data_binarypvalues <- final_met_data %>%
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

## only keep metabolite column and p-adjusted columns
final_met_data_padj <- final_met_data_binarypvalues[,c(1,10:14)]


colnames(final_met_data_padj) <- c("Outcome","ND Cold - ND Warm","D Cold - D Warm",
                              "D Before - ND Before", "D Warm - ND Warm",
                              "D Cold - ND Cold")

## change metabolites to cleaned up version (no _MS)
final_met_data_padj$Outcome <- final_met_data_padj$Outcome %>%
  gsub("_MS","",.)

### read in abbreviations for met variables 
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Input Files")
met_abbr <- read.csv(file="met_abbreviations.csv")
length(met_abbr$fullname)
met_abbr <- met_abbr %>%
  filter(met_abbr$fullname %in% final_met_data_padj$Outcome)
met_abbr$Outcome <- met_abbr$fullname
final_met_data_padj <- merge(final_met_data_padj, met_abbr, by = "Outcome")
row.names(final_met_data_padj) <- final_met_data_padj$abbr

## keep only p-value columns
final_met_data_padj <- final_met_data_padj[,c(2:6)]

## change order of columns 
final_met_data_padj <- final_met_data_padj[,c(5,4,3,2,1)]
colnames(final_met_data_padj)



## Prepare mass-scaled metbaolite data to make a heatmap of relative levels of all metabolites by z-scoring
## z-score all data so that it is scaled to the mean for each metabolite then the mean is centered to 0
metdata_zscore <- metdata_scaled[,c(1:49)]
metdata_zscore[,c(5:49)] <- scale(metdata_zscore[,c(5:49)], scale = TRUE, center = TRUE)

## make a version of the z-scored data that is group averaged for each of the six treatments for each metabolite
metdata_zscore$totaltreatment <- paste(metdata_zscore$diapause_treat, metdata_zscore$acc_treat)
metdata_zscore_grouped <- metdata_zscore %>%
  group_by(totaltreatment) %>%
  summarise(across(!c(no.,Label, diapause_treat,acc_treat), mean))
View(metdata_zscore_grouped)

## change orientation from metabolite columns and treatment groups rows to
## metabolite rows and treatment group columns
metdata_zscore_grouped_long <- metdata_zscore_grouped %>%
  dplyr::select(2:46) %>%
  t()

## add shortened versions of the treatment groups as column names 
colnames(metdata_zscore_grouped_long) <- paste(substr(metdata_zscore_grouped$totaltreatment,1,1),
                                               word(metdata_zscore_grouped$totaltreatment,2))
colnames(metdata_zscore_grouped_long)

## change order of treatment group columns
metdata_zscore_grouped_long <- metdata_zscore_grouped_long[, c(4,6,5,1,3,2)] 

## convert from data_frame to matrix 
met_matrix_grouped <- metdata_zscore_grouped_long %>% 
  as.matrix()

## add abbreviated metabolite names to the metabolite matrix generated above
row.names(met_matrix_grouped) <- row.names(met_matrix_grouped) %>%
  gsub("_MS","",.) 

### read in abbreviations for met variables 
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Input Files")
met_abbr <- read.csv(file="met_abbreviations.csv")
length(met_abbr$fullname)
met_abbr <- met_abbr %>%
  filter(met_abbr$fullname %in% row.names(met_matrix_grouped))
View(met_abbr)
row.names(met_matrix_grouped) <- met_abbr$abbr
View(met_matrix_grouped)

## provide info on color pattern for annotation columns
my_colour_pairwise <- list(
  'ND Cold - ND Warm' = c(Yes = "black", No = "white"),
  'D Cold - D Warm' = c(Yes = "black", No = "white"),
  'D Before - ND Before' = c(Yes = "darkgrey", No = "white"),
  'D Warm - ND Warm' = c(Yes = "darkgrey", No = "white"),
  'D Cold - ND Cold' = c(Yes = "darkgrey", No = "white")
)

## plot heat map
heatmap_met_pairwisepadj <- pheatmap::pheatmap(met_matrix_grouped,
                                 cluster_cols = F,
                                 #color = myColor,
                                 cutree_rows = 8,
                                 annotation_row = final_met_data_padj,
                                 annotation_colors = my_colour_pairwise,
                                 angle_col = 90,
                                 annotation_legend = F,
                                 fontsize = 10,
                                 fontsize_row = 12,
                                 fontsize_col = 12,
                                 cellwidth = 25,
                                 cellheight = 12,
                                 gaps_col = 3,
                                 border_color = "black")


setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Figures/Metabolomics")

ggsave("heatmap_met_emmeans_pairwisepadj_straightlabels.png",plot=heatmap_met_pairwisepadj,dpi=600,units='in',
       width=9,height=10)

