

### install and load packages
install.packages("naniar")
install.packages("ggplot2")
install.packages("ggbeeswarm")
install.packages("reshape2")
install.packages("lme4")
install.packages("car")
install.packages("predictmeans")
install.packages("multcomp")
install.packages("ade4")
install.packages("factoextra")
install.packages("dplyr")
install.packages("marginaleffects")
install.packages("data.table")
install.packages("ade4")
install.packages("marginaleffects")
install.packages("data.table")
install.packages("factoextra")
devtools::install_github("thomasp85/patchwork")
install.packages("rcompanion")
install.packages("multcompView")

library("multcompView")
library(rcompanion)
library(dplyr)
library(naniar)
library(ggplot2)
library(ggbeeswarm)
library(reshape2)
library(lme4)
library(car)
library(predictmeans)
library(multcomp)
library(ade4)
library("factoextra")
library(gridExtra)
library(data.table)
library(marginaleffects)
library(tidyr)
library(tidyverse)
library("stringr")
library("devtools")
library("patchwork")

setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Data Files for Submission")
metdata <- read.csv(file="Metabolite_RawData.csv")
metdata_parentalcageinfo <- read.csv(file="SampleInfo_Fall2024.csv")
head(metdata)
View(metdata)

## remove columns where any of the samples have an NF
metdata <- metdata[ , colSums(is.na(metdata)) == 0]
head(metdata)

## first need to scale each metabolite value to the starting mass of each sample
metdata <- metdata%>%mutate_at(vars(X2.Oxoglutaric_acid:Valine),list(MS=~./Mass..mg.))
#metdata <- metdata%>%mutate_at(vars(X2_and_3.phosphoglyceric_acid:Valine),list(MS=~./Mass)) ## when including NA value metabolites
View(metdata)
metdata$totaltreatment <- paste(metdata$diapause_treat,metdata$acc_treat)


## subset data so it only contains the mass-scaled values
metdata_scaled <- metdata[,c(1,2,3,4,51:96)]

View(metdata_scaled)



##### visualize grouping with PCA
# first center the data
metdata_scaled_onlymet <- metdata_scaled[,c(5:49)] ## get rid of the ID and group columns
metdata_met_centered <- scale(metdata_scaled_onlymet, scale = FALSE, center = TRUE)
# second scale the data
metdata_met_centered_scaled <- scale(metdata_met_centered, center = FALSE, scale = TRUE)

colnames(metdata_met_centered_scaled) <- colnames(metdata_met_centered_scaled) %>%
  gsub("_MS","",.)

## read in abbreviations for met variables 
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Input Files")
met_abbr <- read.csv(file="met_abbreviations.csv")
length(met_abbr$fullname)
met_abbr <- met_abbr %>%
  filter(met_abbr$fullname %in% colnames(metdata_met_centered_scaled))
View(met_abbr)
colnames(metdata_met_centered_scaled) <- met_abbr$abbr

pca_data_met <- prcomp(metdata_met_centered_scaled)

groups <- as.factor(paste(substr(metdata$totaltreatment,1,1),
                word(metdata$totaltreatment,2)))

met_pca_acc <- fviz_pca_ind(pca_data_met,
                            col.ind = groups, # color by groups
                            addEllipses = TRUE, # Concentration ellipses
                            ellipse.type = "confidence",
                            ellipse.level=0.95,
                            legend.title = "Groups",
                            title = "",
                            repel = TRUE,
                            axes = c(1,2),
                            pointsize = 1.75,
                            alpha.ind = 0.65,
                            geom = "point") + 
#  labs(tag="A") +
  ggplot2::aes(linetype = groups)+
  theme_bw()+
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_text(size = 14,face="bold"),
        axis.text.y = element_text(size=13),
        strip.text = element_text(size = 13,face="bold"),
        strip.text.x = element_text(size = 13, margin = margin(t = 0.7, b = 0.7, unit = "pt"),hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = NA)) +
  theme(legend.position = c(0.75,0.15),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.box.background = element_rect(color="black", size=0.35))+
 # geom_point(size=0.5,alpha=0.5,aes(color=groups,fill=groups))+
  scale_shape_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     labels = c("ND Before", "D Before", "ND Warm", "D Warm", "ND Cold", "D Cold"),
                     values = c(15,16,15,16,15,16))+
  scale_fill_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     values = c("#747474","#747474","#C00000","#C00000","#4E95D9","#4E95D9"),
                    breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                    labels = c("ND Before", "D Before", "ND Warm", "D Warm", "ND Cold", "D Cold"))+
  scale_color_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     values = c("#747474","#747474","#C00000","#C00000","#4E95D9","#4E95D9"),
                     breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     labels = c("ND Before", "D Before", "ND Warm", "D Warm", "ND Cold", "D Cold"))+
  scale_linetype_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                        breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                          values = c("dashed","solid","dashed","solid","dashed","solid"),
                        labels = c("ND Before", "D Before", "ND Warm", "D Warm", "ND Cold", "D Cold"))+
  guides(color = guide_legend(nrow=3,ncol=2))+
  theme(plot.margin = unit(c(0,0,1,1), "cm"))+
  xlab("PC1 (35.9%)")+
  ylab("PC2 (18.2%)")
met_pca_acc  


################# run lmer on PC1 and PC2 to test for differences across treatment groups

pc_data <- as.data.frame(pca_data_met$x[,c(1:2)])
pc_data$diapause_treat <- as.factor(metdata_scaled$diapause_treat)
pc_data$acc_treat <- as.factor(metdata_scaled$acc_treat)
pc_data$Label <- metdata_scaled$Label
View(pc_data)

##### run separate anova's on PC1 and PC2
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Input Files")
metdata_parentalcageinfo <- read.csv("SampleInfo_Fall2024.csv",header = T)
metdata_parentalcageinfo <- metdata_parentalcageinfo[c(37:72),c(1,2)]
pc_data <- merge(pc_data,metdata_parentalcageinfo,by="Label")
pc_data$parental_cage_num <- as.character(pc_data$parental_cage_num)


## use non-parametric test because of non-normal residuals
### so use non-parametric kruskal-wallis followed by pairwise wilcox test 
pc_data$totaltreatment <- paste(pc_data$diapause_treat, pc_data$acc_treat)
pc_data$totaltreatment_short <- as.factor(paste(substr(pc_data$totaltreatment,1,1),
                                                         word(pc_data$totaltreatment,2)))
kruskal.test(PC2 ~ totaltreatment_short, data = pc_data)
wilcox_pc2 <- pairwise.wilcox.test(pc_data$PC2, pc_data$totaltreatment_short,
                     p.adjust.method = "BH")

PT2 = fullPTable(wilcox_pc2$p.value) ## converts lower triangle matrix to full square matrix needed for multcompLetters
exp_letters2 <- multcompLetters(PT2)
exp_letters2

### now non-parametric version for PC1
kruskal.test(PC1 ~ totaltreatment_short, data = pc_data)
wilcox_pc1 <- pairwise.wilcox.test(pc_data$PC1, pc_data$totaltreatment_short,
                                   p.adjust.method = "BH")
PT1 = fullPTable(wilcox_pc1$p.value) ## converts lower triangle matrix to full square matrix needed for multcompLetters
exp_letters1 <- multcompLetters(PT1)
exp_letters1
 

##### look for correlation between PCs and cold shock tolerance 

## read in CS data 
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Input Files")
CS_data <- read.csv(file = "coldshocksurvival_scaledaverages.csv")

CS_data$totaltreatment <- paste(CS_data$diapause_treat, CS_data$acc_treat)
CS_data <- CS_data %>%
  dplyr::select(totaltreatment, CS_survival_mean, CS_treat) %>%
  pivot_wider(
    names_from = CS_treat,
    values_from = CS_survival_mean) 

names(CS_data)[2] <- "CS_survival_minus9"
names(CS_data)[3] <- "CS_survival_minus11"

PC_coldshockdata <- merge(pc_data,CS_data, by = "totaltreatment")
PC_coldshockdata$totaltreatment_short <- as.factor(paste(substr(PC_coldshockdata$totaltreatment,1,1),
                          word(PC_coldshockdata$totaltreatment,2)))

PC2_minus9_CSsurvival_plot <-ggplot(PC_coldshockdata,aes(x=PC2,y=CS_survival_minus9))+
  geom_point(data=PC_coldshockdata,aes(x=PC2,y=CS_survival_minus9,
                                                            fill=totaltreatment_short,colour=totaltreatment_short,
                                       shape=totaltreatment_short), size =2)+
  theme_bw()+
  xlab("PC2")+
  ylab('-9°C Survival Rate')+
  coord_cartesian(ylim = c(0,1)) +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size=11),
        legend.text = element_text(size=9),
        legend.title= element_text(size=9))+
  scale_shape_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     labels = c("ND Before", "D Before", "ND Warm", "D Warm", "ND Cold", "D Cold"),
                     values = c(15,16,15,16,15,16))+
  scale_fill_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                    values = c("#747474","#747474","#C00000","#C00000","#4E95D9","#4E95D9"),
                    breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                    labels = c("ND Before", "D Before", "ND Warm", "D Warm", "ND Cold", "D Cold"))+
  scale_color_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     values = c("#747474","#747474","#C00000","#C00000","#4E95D9","#4E95D9"),
                     breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     labels = c("ND Before", "D Before", "ND Warm", "D Warm", "ND Cold", "D Cold"))+
  geom_line(stat="smooth",method = "lm",
            size = 1,
            linetype ="solid")+
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_text(size = 14,face="bold"),
        axis.text.y = element_text(size=13),
        strip.text = element_text(size = 13,face="bold"),
        strip.text.x = element_text(size = 13, margin = margin(t = 0.7, b = 0.7, unit = "pt"),hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = NA)) +
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
PC2_minus9_CSsurvival_plot

PC1_minus9_CSsurvival_plot <-ggplot(PC_coldshockdata,aes(x=PC1,y=CS_survival_minus9))+
  geom_point(data=PC_coldshockdata,aes(x=PC1,y=CS_survival_minus9,
                                       fill=totaltreatment_short,colour=totaltreatment_short,
                                       shape=totaltreatment_short), size=2)+
  theme_bw()+
  xlab("PC1")+
  ylab('-9°C Survival Rate')+
  coord_cartesian(ylim = c(0,1)) +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size=11),
        legend.text = element_text(size=9),
        legend.title= element_text(size=9))+
  scale_shape_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     labels = c("ND Before", "D Before", "ND Warm", "D Warm", "ND Cold", "D Cold"),
                     values = c(15,16,15,16,15,16))+
  scale_fill_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                    values = c("#747474","#747474","#C00000","#C00000","#4E95D9","#4E95D9"),
                    breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                    labels = c("ND Before", "D Before", "ND Warm", "D Warm", "ND Cold", "D Cold"))+
  scale_color_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     values = c("#747474","#747474","#C00000","#C00000","#4E95D9","#4E95D9"),
                     breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     labels = c("ND Before", "D Before", "ND Warm", "D Warm", "ND Cold", "D Cold"))+
  geom_line(stat="smooth",method = "lm",
            size = 1,
            linetype ="solid")+
  theme(axis.title.x = element_text(size = 14,face="bold"),
        axis.text.x = element_text(size = 13),
        axis.title.y = element_text(size = 14,face="bold"),
        axis.text.y = element_text(size=13),
        strip.text = element_text(size = 13,face="bold"),
        strip.text.x = element_text(size = 13, margin = margin(t = 0.7, b = 0.7, unit = "pt"),hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = NA)) +
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
PC1_minus9_CSsurvival_plot

cor.test(PC_coldshockdata$PC1, PC_coldshockdata$CS_survival_minus9, method = "pearson")
## p = 0.8519, corr = -0.032
cor.test(PC_coldshockdata$PC2, PC_coldshockdata$CS_survival_minus9, method = "kendall") ## Kendall is robust to outliers 
## p < 0.001, tau = 0.542


max(PC_coldshockdata$PC1)
PC1_plot <- ggplot(PC_coldshockdata,aes(x=acc_treat,y=PC1,group=diapause_treat,shape=diapause_treat))+
  stat_summary(fun.data=mean_se,size=0.7,color="black",aes(colour=acc_treat, fill=acc_treat,shape=diapause_treat),
               position=position_dodge(width=0.4),fun.args = list(mult = 2))+
  xlab("")+
  scale_x_discrete(limits = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                   breaks = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                   labels = c("Before Acc.", "Warm Acc.", "Cold Acc.")) + 
  scale_fill_manual(values=c("#747474","#4E95D9", "#C00000"),guide = 'none')+
  scale_colour_manual(values=c("#747474","#4E95D9", "#C00000"),guide = 'none')+
  scale_shape_manual(values=c(21,22))+
  coord_cartesian(ylim = c(-6, 6.5)) +
  theme_bw() +
  labs(tag = "A")+
  theme(axis.title.x = element_text(size = 11,face="bold"),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 11,face="bold"),
        axis.text.y = element_text(size=11),
        strip.text = element_text(size = 11,face="bold"),
        strip.text.x = element_text(size = 11, margin = margin(t = 0.7, b = 0.7, unit = "pt"),hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(face="bold")) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
PC1_plot

min(PC_coldshockdata$PC2)
max(PC_coldshockdata$PC2)
PC2_plot <- ggplot(PC_coldshockdata,aes(x=acc_treat,y=PC2,group=diapause_treat,shape=diapause_treat))+
  stat_summary(fun.data=mean_se,size=0.7,color="black",aes(colour=acc_treat, fill=acc_treat,shape=diapause_treat),
               position=position_dodge(width=0.4),fun.args = list(mult = 2))+
  xlab("Acclimation Group")+
  scale_x_discrete(limits = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                   breaks = c("Before Acclimation","Warm Acclimated (15d)", "Cold Acclimated (15d)"),
                   labels = c("Before Acc.", "Warm Acc.", "Cold Acc.")) + 
  scale_fill_manual(values=c("#747474","#4E95D9", "#C00000"),guide = 'none')+
  scale_colour_manual(values=c("#747474","#4E95D9", "#C00000"),guide = 'none')+
  scale_shape_manual(values=c(21,22))+
  theme_bw() +
  labs(tag = "B")+
  coord_cartesian(ylim = c(-5, 5)) +
  theme(axis.title.x = element_text(size = 11,face="bold"),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 11,face="bold"),
        axis.text.y = element_text(size=11),
        strip.text = element_text(size = 11,face="bold"),
        strip.text.x = element_text(size = 11, margin = margin(t = 0.7, b = 0.7, unit = "pt"),hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(face="bold")) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
PC2_plot

### Now add PCA loadings plots
met_loadings <- as.data.frame(pca_data_met$rotation)
met_loadings$met <- row.names(met_loadings)

met_loadings_PC1 <- met_loadings %>%
  arrange(desc(PC1))
# lock in factor level order
met_loadings_PC1$met <- factor(met_loadings_PC1$met, levels = met_loadings_PC1$met)


met_loadings_PC1_plot <- ggplot(data=met_loadings_PC1, aes(x = met, y = PC1)) +
  geom_col(aes(fill = PC1>0)) + 
  coord_flip() +
  scale_fill_manual(values = c("#90aba3","#c86d51"))+
  theme_bw() +
  scale_y_continuous(limits = c(-0.2,0.24),
                     breaks =c(-0.2,0,0.2)) +
  ylab("PC1 Loading")+
  theme(legend.position = "",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size=10,face = c("plain")),
        strip.text = element_text(size = 11,face="bold"),
        plot.title = element_text(size=11,face="bold"))+
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
met_loadings_PC1_plot

met_loadings_PC2 <- met_loadings %>%
  arrange(desc(PC2))
# lock in factor level order
met_loadings_PC2$met <- factor(met_loadings_PC2$met, levels = met_loadings_PC2$met)

met_loadings_PC2_plot <- ggplot(data=met_loadings_PC2, aes(x = met, y = PC2)) +
  geom_col(aes(fill = PC2>0)) + 
  coord_flip() +
  scale_fill_manual(values = c("#90aba3","#c86d51"))+
  theme_bw() +
  scale_y_continuous(limits = c(-0.33,0.25),
                     breaks =c(-0.25,0,0.25)) +
  ylab("PC2 Loading")+
  theme(legend.position = "",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size=10,face = c("plain")),
        strip.text = element_text(size = 11,face="bold"),
        plot.title = element_text(size=11,face="bold"))+
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
met_loadings_PC2_plot


layout_design <- "
BD
AC
"

met_cssurv_plot <-  (PC1_minus9_CSsurvival_plot / PC2_minus9_CSsurvival_plot) 


met_pca_cssurv_plot_loadings <- met_pca_acc  + (PC1_minus9_CSsurvival_plot / PC2_minus9_CSsurvival_plot) +  (met_loadings_PC1_plot + met_loadings_PC2_plot) + plot_layout(widths = c(2, 1, 1))


met_pca_plot_loadings <- met_pca_acc +  (met_loadings_PC1_plot + met_loadings_PC2_plot) + plot_layout(widths = c(1.5, 1))


met_pca_values <- guide_area() / PC1_plot / PC2_plot + plot_layout(guides = "collect") + plot_layout(heights = c(0.2,1,1)) 





setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Figures/Metabolomics/")

ggsave("met_pca_CSsurvival_plot.png",plot=met_pca_cssurv_plot,
       dpi=600,units='in',width=15,height=7)

ggsave("met_pca_loadings_plot.png",plot=met_pca_plot_loadings,
       dpi=600,units='in',width=13,height=7)

ggsave("met_pca_CSsurvival_loadings_plot.png",plot=met_pca_cssurv_plot_loadings,
       dpi=600,units='in',width=16,height=7)


ggsave("met_pca_wrappedloadings_plot.png",plot=wrapped_pca,
       dpi=600,units='in',width=10,height=9)
ggsave("met_pca_CSsurvival_plot.png",plot=met_cssurv_plot ,
       dpi=600,units='in',width=5,height=9)

ggsave("met_pca_pcvalues_plot.png",plot=met_pca_values,
       dpi=600,units='in',width=4,height=6)

