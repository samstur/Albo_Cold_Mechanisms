

## make summary plots and run stats for summary variables related to lipid stuff:
# PE:PC, average degree of unsaturation, proportion polyunsatuated, average length of fatty acyls
# for variables with normal distributions of residuals, run lmer with ANOVA on totaltreatment then all pairwise comparisons
# for variables with non-normal residuals, run Kruskal-wallis with totaltreatment, then do all pairwise comparisons

library(ggplot2)
library(patchwork)
library(stringr)
library(rcompanion)
library(multcompView)
library(broom)

## read in both original datasets
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

## add ratio of PE to PC
lipdata$PEPCratio <- lipdata$PE/lipdata$PC


################### now read in fatty acyl data and compute summary variables
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Data Files for Submission")
fa_data <- read.csv(file="FattyAcyls_molarprops_RawData.csv")

## rename factors 
fa_data$diapause_treat[fa_data$diapause_treat == 'SD'] <- 'Diapause'
fa_data$diapause_treat[fa_data$diapause_treat == 'LD'] <- 'Non-diapause'
fa_data$acc_treat[fa_data$acc_treat == 'PREACC'] <- 'Before Acclimation'
fa_data$acc_treat[fa_data$acc_treat == 'WARM'] <- 'Warm Acclimated (15d)'
fa_data$acc_treat[fa_data$acc_treat == 'COLD'] <- 'Cold Acclimated (15d)'
fa_data$acc_treat <- factor(fa_data$acc_treat, levels=c('Before Acclimation','Warm Acclimated (15d)','Cold Acclimated (15d)'))
fa_data$diapause_treat <- factor(fa_data$diapause_treat, levels=c('Non-diapause','Diapause'))

### percent unsaturated
fa_data$unsat <- fa_data$C14.1n5 + fa_data$C16.1n7 + fa_data$C18.1n9c +
  fa_data$C18.2n6c + fa_data$C18.3n6 + fa_data$C18.3n3 +
  fa_data$C20.3n6 + fa_data$C20.4n6 + fa_data$C20.5n3


### now look at the relative proportions across treatments of fatty acids with 
## no double bonds, 1 double bond, 3, 4 
fa_data$noDB <- fa_data$C12.0 + fa_data$C13.0 + fa_data$C14.0 + fa_data$C15.0 +
  fa_data$C16.0 + fa_data$C17.0 + fa_data$C18.0 + fa_data$C20.0 + fa_data$C22.0
fa_data$oneDB <- fa_data$C14.1n5 + fa_data$C16.1n7 + fa_data$C18.1n9c
fa_data$twoDB <-   fa_data$C18.2n6c
fa_data$threeDB <- fa_data$C18.3n6 + fa_data$C18.3n3 + fa_data$C20.3n6
fa_data$fourDB <- fa_data$C20.4n6
fa_data$fiveDB <- fa_data$C20.5n3

## check that adding up gets to 100 for each sample
fa_data$noDB+fa_data$oneDB+fa_data$twoDB+fa_data$threeDB+fa_data$fourDB+fa_data$fiveDB

### now look at changes in the average length of the fatty acids across samples
fa_data$totalC12 <- fa_data$C12.0
fa_data$totalC13 <- fa_data$C13.0
fa_data$totalC14 <- fa_data$C14.0 + fa_data$C14.1n5
fa_data$totalC15 <- fa_data$C15.0
fa_data$totalC16 <- fa_data$C16.0 + fa_data$C16.1n7
fa_data$totalC17 <- fa_data$C17.0
fa_data$totalC18 <- fa_data$C18.1n9c + fa_data$C18.2n6c + fa_data$C18.0 + fa_data$C18.3n6 + fa_data$C18.3n3
fa_data$totalC20 <- fa_data$C20.0 + fa_data$C20.3n6 + fa_data$C20.4n6 + fa_data$C20.5n3
fa_data$totalC22 <- fa_data$C22.0
## check that adding up gets to 100 for each sample
fa_data$totalC12 + fa_data$totalC13 + fa_data$totalC14 + fa_data$totalC15 + fa_data$totalC16 +
  fa_data$totalC17 + fa_data$totalC18 + fa_data$totalC20 + fa_data$totalC22

## calculate the average degree of unsaturation (average number of double bonds) => this is slightly
## different than the proportion of unsaturated FA because degree of unsaturation gives extra weight to fatty acyls that have multiple 
## double bonds while the prop of unsaturated FA treats a MUFA and PUFA the same 
fa_data$av_unsat <- ((fa_data$noDB*0) + (fa_data$oneDB*1) + (fa_data$twoDB*2) + (fa_data$threeDB*3) +
                       (fa_data$fourDB*4) + (fa_data$fiveDB*5))/100 ## divide by 100 because the DB columns are in percentages, not proportions (25 = 25% = 0.25 in prop)

## now calculate average length 
fa_data$av_length <- ((fa_data$totalC12*12) + (fa_data$totalC13*13) + (fa_data$totalC14*14) +
                        (fa_data$totalC15*15) + (fa_data$totalC16*16) + (fa_data$totalC17*17) + (fa_data$totalC18*18) +
                        (fa_data$totalC20*20) + (fa_data$totalC22*22))/100

## now calculate prop PUFA (poly-unsaturated)
fa_data$prop_PUFA <- fa_data$twoDB + fa_data$threeDB + fa_data$fourDB + fa_data$fiveDB

## ratio of unsaturated to saturated 
fa_data$ratio_UFA_SFA <- fa_data$unsat/fa_data$noDB



## now plot across all six groups 
PEPC_plot <- ggplot(lipdata,aes(x=acc_treat,y=PEPCratio,group=diapause_treat,shape=diapause_treat))+
  stat_summary(fun.data=mean_se,size=0.8,color="black",aes(colour=acc_treat, fill=acc_treat,shape=diapause_treat),
               position=position_dodge(width=0.4),fun.args = list(mult = 2))+
  xlab("")+
  scale_x_discrete(limits = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                   breaks = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                   labels = c("Before", "Warm", "Cold")) + 
  scale_fill_manual(limits = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                    breaks = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                    labels = c("Before", "Warm", "Cold"),
                    values=c("#747474", "#C00000","#4E95D9"),guide = 'none')+
  scale_colour_manual(limits = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                      breaks = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                      labels = c("Before", "Warm", "Cold"),
                      values=c("#747474", "#C00000","#4E95D9"),guide = 'none')+
  scale_shape_manual(values=c(22,21), name = "")+
  theme_bw() +
  labs(tag="A")+
  ylab("PE:PC Ratio")+
  theme(axis.title.x = element_text(size = 11,face="bold"),
        axis.text.x = element_text(size= 11),
        axis.title.y = element_text(size = 11,face="bold"),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size = 11),
        plot.tag = element_text(face = "bold")) +
  theme(legend.position = "top",
        plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
PEPC_plot

av_unsat_plot <- ggplot(fa_data,aes(x=acc_treat,y=av_unsat,group=diapause_treat,shape=diapause_treat))+
  stat_summary(fun.data=mean_se,size=0.8,color="black",aes(colour=acc_treat, fill=acc_treat,shape=diapause_treat),
               position=position_dodge(width=0.4),fun.args = list(mult = 2))+
  xlab("Acclimation Group")+
  scale_x_discrete(limits = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                   breaks = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                   labels = c("Before", "Warm", "Cold")) + 
  scale_fill_manual(limits = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                    breaks = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                    labels = c("Before", "Warm", "Cold"),
                    values=c("#747474", "#C00000","#4E95D9"),guide = 'none')+
  scale_colour_manual(limits = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                      breaks = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                      labels = c("Before", "Warm", "Cold"),
                      values=c("#747474", "#C00000","#4E95D9"),guide = 'none')+
  scale_shape_manual(values=c(22,21), name = "")+
  theme_bw() +
  labs(tag="C")+
  ylab("FA Unsaturation Index")+
  theme(axis.title.x = element_text(size = 11,face="bold"),
        axis.text.x = element_text(size= 11),
        axis.title.y = element_text(size = 11,face="bold"),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size = 11),
        plot.tag = element_text(face = "bold")) +
  theme(legend.position = "top",
        plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
av_unsat_plot 


prop_polyunsat_plot <- ggplot(fa_data,aes(x=acc_treat,y=prop_PUFA,group=diapause_treat,shape=diapause_treat))+
  stat_summary(fun.data=mean_se,size=0.8,color="black",aes(colour=acc_treat, fill=acc_treat,shape=diapause_treat),
               position=position_dodge(width=0.4),fun.args = list(mult = 2))+
  xlab("Acclimation Group")+
  scale_x_discrete(limits = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                   breaks = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                   labels = c("Before", "Warm", "Cold")) + 
  scale_fill_manual(limits = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                    breaks = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                    labels = c("Before", "Warm", "Cold"),
                    values=c("#747474", "#C00000","#4E95D9"),guide = 'none')+
  scale_colour_manual(limits = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                      breaks = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                      labels = c("Before", "Warm", "Cold"),
                      values=c("#747474", "#C00000","#4E95D9"),guide = 'none')+
  scale_shape_manual(values=c(22,21), name = "")+
  theme_bw() +
  labs(tag="D")+
  ylab("% Poly-unsaturated FA")+
  theme(axis.title.x = element_text(size = 11,face="bold"),
        axis.text.x = element_text(size= 11),
        axis.title.y = element_text(size = 11,face="bold"),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size = 11),
        plot.tag = element_text(face = "bold")) +
  theme(legend.position = "top",
        plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
prop_polyunsat_plot

av_length_plot <- ggplot(fa_data,aes(x=acc_treat,y=av_length,group=diapause_treat,shape=diapause_treat))+
  stat_summary(fun.data=mean_se,size=0.8,color="black",aes(colour=acc_treat, fill=acc_treat,shape=diapause_treat),
               position=position_dodge(width=0.4),fun.args = list(mult = 2))+
  xlab("")+
  scale_x_discrete(limits = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                   breaks = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                   labels = c("Before", "Warm", "Cold")) + 
  scale_fill_manual(limits = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                    breaks = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                    labels = c("Before", "Warm", "Cold"),
                    values=c("#747474", "#C00000","#4E95D9"),guide = 'none')+
  scale_colour_manual(limits = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                      breaks = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                      labels = c("Before", "Warm", "Cold"),
                      values=c("#747474", "#C00000","#4E95D9"),guide = 'none')+
  scale_shape_manual(values=c(22,21), name = "")+
  theme_bw() +
  labs(tag="B")+
  ylab("Average FA Length")+
  theme(axis.title.x = element_text(size = 11,face="bold"),
        axis.text.x = element_text(size= 11),
        axis.title.y = element_text(size = 11,face="bold"),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size = 11),
        plot.tag = element_text(face = "bold")) +
  theme(legend.position = "top",
        plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
av_length_plot


lipid_summary_plot <- guide_area() / (PEPC_plot + av_unsat_plot) / (prop_polyunsat_plot + av_length_plot) + plot_layout(heights = c(0.2,1,1)) + plot_layout(guides = "collect")
lipid_summary_plot

lipid_summary_plot_order2 <- guide_area() / (PEPC_plot + av_length_plot) / (av_unsat_plot + prop_polyunsat_plot) + plot_layout(heights = c(0.2,1,1)) + plot_layout(guides = "collect")
lipid_summary_plot_order2

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Figures/Lipidomics/")

ggsave("lipid_summary_plot_order2.png",plot=lipid_summary_plot_order2,
       dpi=600,units='in',width=7,height=7)



### now for each summary trait, run ANOVA across six treatment groups then posthoc contrasts

# add info on parental cage to both datasets 
## first add parental cage data 
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Input Files")
lipdata_parentalcageinfo <- read.csv("SampleInfo_Fall2024.csv",header = T)
lipdata_parentalcageinfo <- lipdata_parentalcageinfo[c(1:36),c(1,2)]
lipdata <- merge(lipdata,lipdata_parentalcageinfo,by="Label")
lipdata$parental_cage_num <- as.character(lipdata$parental_cage_num)
fa_data <- merge(fa_data, lipdata_parentalcageinfo, by = "Label")
fa_data$parental_cage_num <- as.character(fa_data$parental_cage_num)

# add total treatment to both datasets
lipdata$totaltreatment <- paste(lipdata$diapause_treat, lipdata$acc_treat)
fa_data$totaltreatment <- paste(fa_data$diapause_treat, fa_data$acc_treat)
lipdata$totaltreatment_short <- as.factor(paste(substr(lipdata$totaltreatment,1,1),
                                                    word(lipdata$totaltreatment,2)))
fa_data$totaltreatment_short <- as.factor(paste(substr(fa_data$totaltreatment,1,1),
                                                word(fa_data$totaltreatment,2)))


# first PE:PC ratio
PEPC_model <- lmer(sqrt(PEPCratio) ~ totaltreatment + (1|parental_cage_num), data = lipdata)
PEPC_model_anova <- Anova(PEPC_model, type = 2)  # Type II SS is most common
## sig effect of totaltreatment (p=0.033)
shapiro.test(resid(PEPC_model))$p.value ## p<0.05 and transformations don't help
hist(resid(PEPC_model))
qqnorm(resid(PEPC_model))
qqline(resid(PEPC_model))

# instead do non-parametric test
### so use non-parametric kruskal-wallis followed by pairwise wilcox test 
kruskal.test(PEPCratio ~ totaltreatment_short, data = lipdata) # p = 0.0119
wilcox_pepcratio <- pairwise.wilcox.test(lipdata$PEPCratio, lipdata$totaltreatment_short,
                                   p.adjust.method = "BH")


PT_PEPC = fullPTable(wilcox_pepcratio$p.value) ## converts lower triangle matrix to full square matrix needed for multcompLetters
exp_letters_PEPC <- multcompLetters(PT_PEPC)
exp_letters_PEPC
results_PEPC <- broom::tidy(wilcox_pepcratio)
results_PEPC <- as.data.frame(results_PEPC)
colnames(results_PEPC)[3] <- "p.value_PEPC"

# FA unsaturation index
unsat_model <- lmer(log(av_unsat) ~ totaltreatment + (1|parental_cage_num), data = fa_data)
unsat_model_anova <- Anova(unsat_model, type = 2)  # Type II SS is most common
## sig effect of totaltreatment
shapiro.test(resid(unsat_model))$p.value ## p<0.05 and transformations don't help
hist(resid(unsat_model))
qqnorm(resid(unsat_model))
qqline(resid(unsat_model))

# instead do non-parametric test
### so use non-parametric kruskal-wallis followed by pairwise wilcox test 
kruskal.test(av_unsat ~ totaltreatment_short, data = fa_data) # p = 0.001
wilcox_avunsat <- pairwise.wilcox.test(fa_data$av_unsat, fa_data$totaltreatment_short,
                                         p.adjust.method = "BH")

PT_avunsat = fullPTable(wilcox_avunsat$p.value) ## converts lower triangle matrix to full square matrix needed for multcompLetters
exp_letters_avunsat <- multcompLetters(PT_avunsat)
exp_letters_avunsat
results_avunsat <- broom::tidy(wilcox_avunsat)
results_avunsat <- as.data.frame(results_avunsat)
colnames(results_avunsat)[3] <- "p.value_avunsat"


# FA prop PUFA
pufa_model <- lmer(log(prop_PUFA) ~ totaltreatment + (1|parental_cage_num), data = fa_data)
pufa_model_anova <- Anova(pufa_model, type = 2)  # Type II SS is most common
## sig effect of totaltreatment
shapiro.test(resid(pufa_model))$p.value ## p<0.05 and transformations don't help
hist(resid(pufa_model))
qqnorm(resid(pufa_model))
qqline(resid(pufa_model))

# instead do non-parametric test
### so use non-parametric kruskal-wallis followed by pairwise wilcox test 
kruskal.test(prop_PUFA ~ totaltreatment_short, data = fa_data) # p = 0.02
wilcox_proppufa <- pairwise.wilcox.test(fa_data$prop_PUFA, fa_data$totaltreatment_short,
                                       p.adjust.method = "BH")
PT_proppufa = fullPTable(wilcox_proppufa$p.value) ## converts lower triangle matrix to full square matrix needed for multcompLetters
exp_letters_proppufa <- multcompLetters(PT_proppufa)
exp_letters_proppufa
results_proppufa <- broom::tidy(wilcox_proppufa)
results_proppufa <- as.data.frame(results_proppufa)
colnames(results_proppufa)[3] <- "p.value_proppufa"


# FA average length
avlength_model <- lmer(av_length ~ totaltreatment_short + (1|parental_cage_num), data = fa_data)
avlength_model_anova <- Anova(avlength_model, type = 2)  # Type II SS is most common
## sig effect of totaltreatment
shapiro.test(resid(avlength_model))$p.value ## p>0.05 
hist(resid(avlength_model))
qqnorm(resid(avlength_model))
qqline(resid(avlength_model))

kruskal.test(av_length ~ totaltreatment_short, data = fa_data) # p = 0.01
wilcox_avlength <- pairwise.wilcox.test(fa_data$av_length, fa_data$totaltreatment_short,
                                        p.adjust.method = "BH")

PT_avlength = fullPTable(wilcox_avlength$p.value) ## converts lower triangle matrix to full square matrix needed for multcompLetters
exp_letters_avlength <- multcompLetters(PT_avlength)
exp_letters_avlength
results_avlength <- broom::tidy(wilcox_avlength)
results_avlength <- as.data.frame(results_avlength)
colnames(results_avlength)[3] <- "p.value_avlength"

## now merge the posthoc test p-values to save to a file for supplementary table
finaltable <- merge(results_PEPC,results_avunsat, by = c("group1","group2"))
finaltable <- merge(finaltable,results_proppufa, by = c("group1","group2"))
finaltable <- merge(finaltable,results_avlength, by = c("group1","group2"))

finaltable$contrast <- paste(finaltable$group1, finaltable$group2, sep = " - ")
finaltable <- finaltable %>% dplyr::select(contrast, p.value_PEPC, p.value_avunsat, p.value_proppufa, p.value_avlength)

setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Output Files")
write.csv(finaltable, file = "Lipid_SummaryVariable_posthoctests.csv", row.names = F)
