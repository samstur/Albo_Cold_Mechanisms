
###### In this script, I perform a separate PCA analysis on the phospholipid headgroup
## classes as well as the fatty acyl methyl ester species

# load packages
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
library("stringr")
library("tidyr")
library("dplyr")
library(rcompanion)


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

## need to only keep columns with lipid data in it
lipdata_onlylip<- lipdata[,c(6:15)] ## get rid of the ID and group columns

##### visualize grouping with PCA
# first center the data
lipdata_lip_centered <- scale(lipdata_onlylip, scale = FALSE, center = TRUE)
# second scale the data
lipdata_lip_centered_scaled <- scale(lipdata_lip_centered, center = FALSE, scale = TRUE)

pca_data_lip <- prcomp(lipdata_lip_centered_scaled)

lipdata$totaltreatment <- paste(lipdata$diapause_treat, lipdata$acc_treat)
groups <- as.factor(paste(substr(lipdata$totaltreatment,1,1),
                          word(lipdata$totaltreatment,2)))
lip_headgroup_pca <- fviz_pca_ind(pca_data_lip,
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
  ggplot2::aes(linetype = groups)+
  theme_bw()+
  theme(legend.position = c(0.25,0.89),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.box.background = element_rect(color="black", size=0.3))+
  # geom_point(size=0.5,alpha=0.5,aes(color=groups,fill=groups))+
  scale_shape_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     values = c(15,16,15,16,15,16))+
  scale_fill_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                    values = c("#747474","#747474","#C00000","#C00000","#4E95D9","#4E95D9"),
                    breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"))+
  scale_color_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     values = c("#747474","#747474","#C00000","#C00000","#4E95D9","#4E95D9"),
                     breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"))+
  scale_linetype_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                        breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                        values = c("dashed","solid","dashed","solid","dashed","solid"))+
  guides(color = guide_legend(nrow=3,ncol=2))+
  labs(tag = "A") +
  theme(axis.title.x = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size= 11),
        axis.title.y = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=11),
        strip.text = element_text(size = 11,face="bold"),
        plot.tag = element_text(size=14, face = "bold"))
lip_headgroup_pca


################## Now do stats and plots to compare headgroup PCA across the 6 treatment groups
##### run separate anova's on PC1 and PC2
lip_pc_data <- as.data.frame(pca_data_lip$x[,c(1:2)])
lip_pc_data$diapause_treat <- as.factor(lipdata$diapause_treat)
lip_pc_data$acc_treat <- as.factor(lipdata$acc_treat)
lip_pc_data$Label <- lipdata$Label

## first add parental cage data 
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Input Files")
lipdata_parentalcageinfo <- read.csv("SampleInfo_Fall2024.csv",header = T)

lipdata_parentalcageinfo <- lipdata_parentalcageinfo[c(1:36),c(1,2)]
lip_pc_data <- merge(lip_pc_data,lipdata_parentalcageinfo,by="Label")
lip_pc_data$parental_cage_num <- as.character(lip_pc_data$parental_cage_num)

lip_pc_data$totaltreatment <- paste(lip_pc_data$diapause_treat, lip_pc_data$acc_treat)
## instead use non-parametric test
### so use non-parametric kruskal-wallis followed by pairwise wilcox test 
lip_pc_data$totaltreatment_short <- as.factor(paste(substr(lip_pc_data$totaltreatment,1,1),
                                                word(lip_pc_data$totaltreatment,2)))
kruskal.test(PC1 ~ totaltreatment_short, data = lip_pc_data)
wilcox_pc1 <- pairwise.wilcox.test(lip_pc_data$PC1, lip_pc_data$totaltreatment_short,
                                   p.adjust.method = "BH")

PT1 = fullPTable(wilcox_pc1$p.value) ## converts lower triangle matrix to full square matrix needed for multcompLetters
exp_letters1 <- multcompLetters(PT1)
exp_letters1

### same but with PC2 
## instead use non-parametric test
### so use non-parametric kruskal-wallis followed by pairwise wilcox test
kruskal.test(PC2 ~ totaltreatment_short, data = lip_pc_data)
wilcox_pc2 <- pairwise.wilcox.test(lip_pc_data$PC2, lip_pc_data$totaltreatment_short,
                                   p.adjust.method = "BH")

PT2 = fullPTable(wilcox_pc2$p.value) ## converts lower triangle matrix to full square matrix needed for multcompLetters
exp_letters2 <- multcompLetters(PT2)
exp_letters2

## save lipid head group pc data to csv
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Output Files")
write.csv(file = "lipid_headgroup_PCA_values.csv", lip_pc_data, row.names = F)
 

## now plot the PC across all six groups 
lip_PC1_plot <- ggplot(lip_pc_data,aes(x=acc_treat,y=PC1,group=diapause_treat,shape=diapause_treat))+
  stat_summary(fun.data=mean_se,size=0.7,color="black",aes(colour=acc_treat, fill=acc_treat,shape=diapause_treat),
               position=position_dodge(width=0.3),fun.args = list(mult = 2))+
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
  scale_shape_manual(values=c(22,21))+
  labs(tag = "B") +
  theme_bw() +
  coord_cartesian(ylim = c(-4,4))+
  theme(axis.title.x = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size= 11),
        axis.title.y = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=11),
        strip.text = element_text(size = 11,face="bold"),
        strip.text.x = element_text(size = 11, margin = margin(t = 0.7, b = 0.7, unit = "pt"),hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(size = 14, face="bold")) +
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
lip_PC1_plot

lip_PC2_plot <- ggplot(lip_pc_data,aes(x=acc_treat,y=PC2,group=diapause_treat,shape=diapause_treat))+
  stat_summary(fun.data=mean_se,size=0.7,color="black",aes(colour=acc_treat, fill=acc_treat,shape=diapause_treat),
               position=position_dodge(width=0.3),fun.args = list(mult = 2))+
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
  scale_shape_manual(values=c(22,21))+
  theme_bw() +
  coord_cartesian(ylim = c(-3,4)) +
  labs(tag = "C") +
  theme(axis.title.x = element_text(size = 11,face="bold"),
        axis.text.x = element_text(size= 11),
        axis.title.y = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=11),
        strip.text = element_text(size = 11,face="bold"),
        strip.text.x = element_text(size = 11, margin = margin(t = 0.7, b = 0.7, unit = "pt"),hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(size = 14, face="bold")) +
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
lip_PC2_plot

########################## now look at correlations with CS survival for each PC 
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

lip_PC_coldshockdata <- merge(lip_pc_data,CS_data, by = "totaltreatment")
lip_PC_coldshockdata$totaltreatment_short <- as.factor(paste(substr(lip_PC_coldshockdata$totaltreatment,1,1),
                                                         word(lip_PC_coldshockdata$totaltreatment,2)))
View(lip_PC_coldshockdata)
lip_PC1_minus9_CSsurvival_plot <-ggplot(lip_PC_coldshockdata,aes(x=PC1,y=CS_survival_minus9))+
  geom_point(data=lip_PC_coldshockdata,aes(x=PC1,y=CS_survival_minus9,
                                       fill=totaltreatment_short,colour=totaltreatment_short,
                                       shape=totaltreatment_short), size =2)+
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
            linetype ="solid") +
  labs(tag = "D") +
  theme(axis.title.x = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=11),
        strip.text = element_text(size = 11,face="bold"),
        strip.text.x = element_text(size = 11, margin = margin(t = 0.7, b = 0.7, unit = "pt"),hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(size=14, face = "bold")) +
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
lip_PC1_minus9_CSsurvival_plot

lip_PC2_minus9_CSsurvival_plot <-ggplot(lip_PC_coldshockdata,aes(x=PC2,y=CS_survival_minus9))+
  geom_point(data=lip_PC_coldshockdata,aes(x=PC2,y=CS_survival_minus9,
                                           fill=totaltreatment_short,colour=totaltreatment_short,
                                           shape=totaltreatment_short), size =2)+
  theme_bw()+
  xlab("PC2")+
  ylab('-9°C Survival Rate')+
  coord_cartesian(ylim = c(0,1)) +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size=11))+
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
  labs(tag = "E") +
  theme(axis.title.x = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=11),
        strip.text = element_text(size = 11,face="bold"),
        strip.text.x = element_text(size = 11, margin = margin(t = 0.7, b = 0.7, unit = "pt"),hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(size=14, face = "bold")) +
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
lip_PC2_minus9_CSsurvival_plot

### calculate correlations

## both PC's are normally distributed
shapiro.test(lip_PC_coldshockdata$PC2)
shapiro.test(lip_PC_coldshockdata$PC1)

cor.test(lip_PC_coldshockdata$PC1, lip_PC_coldshockdata$CS_survival_minus9, method = "pearson")
## p = 0.4502, corr = -0.1299
cor.test(lip_PC_coldshockdata$PC2, lip_PC_coldshockdata$CS_survival_minus9, method = "pearson")
## p < 0.0001, corr = 0.6254


## add loadings plot to bottom of this 
## https://cran.r-project.org/web/packages/LearnPCA/vignettes/Vig_04_Scores_Loadings.pdf
## loadings = Since individual loadings are defined by a cosine function, they are limited to the range –1 (an angle of
 #rotation of 180◦) to +1 (an angle of rotation of 0◦).5 The sign of the loading indicates how a variable
 #contributes to the principal component. A positive loading indicates that a variable contributes to some
 #degree to the principal component, and a negative loading indicates that its absence contributes to some
 #degree to the principal component. The larger a loading’s relative magnitude, the more important is its
 #presence or absence to the principal component
lip_loadings <- as.data.frame(pca_data_lip$rotation)
lip_loadings$lip_hg <- row.names(lip_loadings)

lip_loadings_PC1 <- lip_loadings %>%
  arrange(desc(PC1))
# lock in factor level order
lip_loadings_PC1$lip_hg <- factor(lip_loadings_PC1$lip_hg, levels = lip_loadings_PC1$lip_hg)

lip_hg_loadings_PC1 <- ggplot(data=lip_loadings_PC1, aes(x = lip_hg, y = PC1)) +
  geom_col(aes(fill = PC1>0)) + 
  coord_flip() +
  scale_fill_manual(values = c("darkblue","darkred"))+
  theme_bw() +
  scale_y_continuous(breaks =c(-0.3,0,0.3)) +
  ylab("PC1 Loading")+
  labs(tag = "F")+
  theme(legend.position = "",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size=11,face = c("bold")),
        strip.text = element_text(size = 11,face="bold"),
        plot.title = element_text(size=11,face="bold"),
        plot.tag = element_text(size = 14, face = "bold"))+
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
lip_hg_loadings_PC1

lip_loadings_PC2 <- lip_loadings %>%
  arrange(desc(PC2))
# lock in factor level order
lip_loadings_PC2$lip_hg <- factor(lip_loadings_PC2$lip_hg, levels = lip_loadings_PC2$lip_hg)

lip_hg_loadings_PC2 <- ggplot(data=lip_loadings_PC2, aes(x = lip_hg, y = PC2)) +
  geom_col(aes(fill = PC2>0)) + 
  coord_flip() +
  scale_fill_manual(values = c("darkblue","darkred"))+
  theme_bw() +
  scale_y_continuous(breaks =c(-0.3,0,0.3,0.6)) +
  ylab("PC2 Loading")+
  labs(tag = "G")+
  theme(legend.position = "",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size=11,face = c("bold")),
        strip.text = element_text(size = 11,face="bold"),
        plot.title = element_text(size=11,face="bold"),
        plot.tag = element_text(size = 14, face = "bold"))+
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
lip_hg_loadings_PC2

lip_pca_cssurv_plot <- lip_headgroup_pca  + (lip_PC1_plot / lip_PC2_plot) + (lip_PC1_minus9_CSsurvival_plot / lip_PC2_minus9_CSsurvival_plot) + lip_hg_loadings_PC1 + lip_hg_loadings_PC2 + plot_layout(widths = c(2, 1, 1, 0.5,0.5))

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Figures/Lipidomics/")

ggsave("lip_hg_pca_CSsurvival_plot.png",plot=lip_pca_cssurv_plot,
       dpi=600,units='in',width=17.5,height=7)


################################################################################################
################################################################################################
# Now look at fatty acyl data

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

## need to only keep columns with fatty acyl data in it
fa_data_onlylip<- fa_data[,c(6:23)] ## get rid of the ID and group columns


##### visualize grouping with PCA
# first center the data
fa_data_onlylip_centered <- scale(fa_data_onlylip, scale = FALSE, center = TRUE)
# second scale the data
fa_data_onlylip_centered_scaled <- scale(fa_data_onlylip_centered, center = FALSE, scale = TRUE)
pca_data_fa <- prcomp(fa_data_onlylip_centered_scaled)

fa_data$totaltreatment <- paste(fa_data$diapause_treat, fa_data$acc_treat)
groups <- as.factor(paste(substr(fa_data$totaltreatment,1,1),
                          word(fa_data$totaltreatment,2)))
lip_fa_pca <- fviz_pca_ind(pca_data_fa,
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
  ggplot2::aes(linetype = groups)+
  theme_bw()+
  #theme(legend.position = c(0.76,0.89),
  theme(legend.position = c(0.72,0.89),
        legend.direction = "horizontal",
        legend.title=element_blank(),
        legend.box.background = element_rect(color="black", size=0.3))+
  coord_cartesian(xlim=c(-6,10))+
  labs(tag = "A") +
  # geom_point(size=0.5,alpha=0.5,aes(color=groups,fill=groups))+
  scale_shape_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     values = c(15,16,15,16,15,16))+
  scale_fill_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                    values = c("#747474","#747474","#C00000","#C00000","#4E95D9","#4E95D9"),
                    breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"))+
  scale_color_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                     values = c("#747474","#747474","#C00000","#C00000","#4E95D9","#4E95D9"),
                     breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"))+
  scale_linetype_manual(limits = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                        breaks = c("N Before","D Before","N Warm","D Warm","N Cold","D Cold"),
                        values = c("dashed","solid","dashed","solid","dashed","solid"))+
  guides(color = guide_legend(nrow=3,ncol=2))+
  theme(axis.title.x = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size= 11),
        axis.title.y = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=11),
        strip.text = element_text(size = 12,face="bold"),
        plot.tag = element_text(size = 14, face = "bold"))
lip_fa_pca


##### run separate anova's on PC1 and PC2
fa_pc_data <- as.data.frame(pca_data_fa$x[,c(1:2)])
fa_pc_data$diapause_treat <- as.factor(fa_data$diapause_treat)
fa_pc_data$acc_treat <- as.factor(fa_data$acc_treat)
fa_pc_data$Label <- fa_data$Label

## first add parental cage data 
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Input Files")
lipdata_parentalcageinfo <- read.csv("SampleInfo_Fall2024.csv",header = T)
View(lipdata_parentalcageinfo)
lipdata_parentalcageinfo <- lipdata_parentalcageinfo[c(1:36),c(1,2)]
fa_pc_data <- merge(fa_pc_data,lipdata_parentalcageinfo,by="Label")
fa_pc_data$parental_cage_num <- as.character(fa_pc_data$parental_cage_num)

fa_pc_data$totaltreatment <- paste(fa_pc_data$diapause_treat, fa_pc_data$acc_treat)

## save lipid head group pc data to csv
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Output Files")
write.csv(file = "lipid_fa_PCA_values.csv", fa_pc_data, row.names = F)


fa_pc_data$totaltreatment_short <- as.factor(paste(substr(fa_pc_data$totaltreatment,1,1),
                                                         word(fa_pc_data$totaltreatment,2)))

## instead use non-parametric versions for both PCs
### now non-parametric version for PC1
head(fa_pc_data)
kruskal.test(PC1 ~ totaltreatment_short, data = fa_pc_data)
wilcox_pc1 <- pairwise.wilcox.test(fa_pc_data$PC1, fa_pc_data$totaltreatment_short,
                                   p.adjust.method = "BH")
PT1 = fullPTable(wilcox_pc1$p.value) ## converts lower triangle matrix to full square matrix needed for multcompLetters
exp_letters1 <- multcompLetters(PT1)
exp_letters1

kruskal.test(PC2 ~ totaltreatment_short, data = fa_pc_data) ## Not significant 


## now plot the PC across all six groups 
fa_PC1_plot <- ggplot(fa_pc_data,aes(x=acc_treat,y=PC1,group=diapause_treat,shape=diapause_treat))+
  stat_summary(fun.data=mean_se,size=0.7,color="black",aes(colour=acc_treat, fill=acc_treat,shape=diapause_treat),
               position=position_dodge(width=0.3),fun.args = list(mult = 2))+
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
  scale_shape_manual(values=c(22,21))+
  theme_bw() +
  coord_cartesian(ylim = c(-5,6)) +
  labs(tag = "B") +
  theme(axis.title.x = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size= 11),
        axis.title.y = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=11),
        strip.text = element_text(size = 11,face="bold"),
        strip.text.x = element_text(size = 11, margin = margin(t = 0.7, b = 0.7, unit = "pt"),hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(size = 14, face = "bold")) +
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
fa_PC1_plot

fa_PC2_plot <- ggplot(fa_pc_data,aes(x=acc_treat,y=PC2,group=diapause_treat,shape=diapause_treat))+
  stat_summary(fun.data=mean_se,size=0.7,color="black",aes(colour=acc_treat, fill=acc_treat,shape=diapause_treat),
               position=position_dodge(width=0.3),fun.args = list(mult = 2))+
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
  scale_shape_manual(values=c(22,21))+
  theme_bw() +
  labs(tag = "C") +
  coord_cartesian(ylim = c(-2,2.5)) + 
  theme(axis.title.x = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size= 11),
        axis.title.y = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=11),
        strip.text = element_text(size = 11,face="bold"),
        strip.text.x = element_text(size = 11, margin = margin(t = 0.7, b = 0.7, unit = "pt"),hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(size = 14, face = "bold")) +
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
fa_PC2_plot

########################## now look at correlations with CS survival for each PC 
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

fa_PC_coldshockdata <- merge(fa_pc_data,CS_data, by = "totaltreatment")
fa_PC_coldshockdata$totaltreatment_short <- as.factor(paste(substr(fa_PC_coldshockdata$totaltreatment,1,1),
                                                             word(fa_PC_coldshockdata$totaltreatment,2)))

fa_PC1_minus9_CSsurvival_plot <-ggplot(fa_PC_coldshockdata,aes(x=PC1,y=CS_survival_minus9))+
  geom_point(data=fa_PC_coldshockdata,aes(x=PC1,y=CS_survival_minus9,
                                           fill=totaltreatment_short,colour=totaltreatment_short,
                                           shape=totaltreatment_short), size =2)+
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
            linetype ="solid") +
  labs(tag = "D")+
  theme(axis.title.x = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=11),
        strip.text = element_text(size = 11,face="bold"),
        strip.text.x = element_text(size = 11, margin = margin(t = 0.7, b = 0.7, unit = "pt"),hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(size = 14, face = "bold")) +
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
fa_PC1_minus9_CSsurvival_plot

fa_PC2_minus9_CSsurvival_plot <-ggplot(fa_PC_coldshockdata,aes(x=PC2,y=CS_survival_minus9))+
  geom_point(data=fa_PC_coldshockdata,aes(x=PC2,y=CS_survival_minus9,
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
  labs(tag = "E") +
  theme(axis.title.x = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size = 11),
        axis.title.y = element_text(size = 12,face="bold"),
        axis.text.y = element_text(size=11),
        strip.text = element_text(size = 11,face="bold"),
        strip.text.x = element_text(size = 11, margin = margin(t = 0.7, b = 0.7, unit = "pt"),hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = NA),
        plot.tag = element_text(size = 14, face = "bold")) +
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
fa_PC2_minus9_CSsurvival_plot

## both PC's are normally distributed
shapiro.test(fa_PC_coldshockdata$PC2)
shapiro.test(fa_PC_coldshockdata$PC1)

cor.test(fa_PC_coldshockdata$PC1, fa_PC_coldshockdata$CS_survival_minus9, method = "pearson")
## p < 0.001, corr = 0.5339
cor.test(fa_PC_coldshockdata$PC2, fa_PC_coldshockdata$CS_survival_minus9, method = "pearson")
## p = 0.6077, corr = 0.089

fa_loadings <- as.data.frame(pca_data_fa$rotation)
fa_loadings$fa_species <- row.names(fa_loadings)

fa_loadings_PC1 <- fa_loadings %>%
  arrange(desc(PC1))
# lock in factor level order
fa_loadings_PC1$fa_species <- factor(fa_loadings_PC1$fa_species, levels = fa_loadings_PC1$fa_species)

fa_loadings_PC1_plot <- ggplot(data=fa_loadings_PC1, aes(x = fa_species, y = PC1)) +
  geom_col(aes(fill = PC1>0)) + 
  coord_flip() +
  scale_fill_manual(values = c("darkblue","darkred"))+
  theme_bw() +
  scale_y_continuous(limits = c(-0.3,0.3),
                     breaks =c(-0.3,0,0.3)) +
  ylab("PC1 Loading") +
  labs(tag = "F") +
  theme(legend.position = "",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size=11,face = c("bold")),
        strip.text = element_text(size = 11,face="bold"),
        plot.title = element_text(size=11,face="bold"),
        plot.tag = element_text(size = 14, face = "bold"))+
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
fa_loadings_PC1_plot

fa_loadings_PC2 <- fa_loadings %>%
  arrange(desc(PC2))
# lock in factor level order
fa_loadings_PC2$fa_species <- factor(fa_loadings_PC2$fa_species, levels = fa_loadings_PC2$fa_species)

fa_loadings_PC2_plot <- ggplot(data=fa_loadings_PC2, aes(x = fa_species, y = PC2)) +
  geom_col(aes(fill = PC2>0)) + 
  coord_flip() +
  scale_fill_manual(values = c("darkblue","darkred"))+
  theme_bw() +
  scale_y_continuous(breaks =c(-0.6,-0.3,0,0.3,0.6)) +
  ylab("PC2 Loading")+
  labs(tag = "G") +
  theme(legend.position = "",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12,face="bold"),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size=11,face = c("bold")),
        strip.text = element_text(size = 11,face="bold"),
        plot.title = element_text(size=11,face="bold"),
        plot.tag = element_text(size = 14, face = "bold"))+
  theme(legend.position = "",
        plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
fa_loadings_PC2_plot

fa_pca_cssurv_plot <- lip_fa_pca  + (fa_PC1_plot / fa_PC2_plot) + (fa_PC1_minus9_CSsurvival_plot / fa_PC2_minus9_CSsurvival_plot) + fa_loadings_PC1_plot + fa_loadings_PC2_plot + plot_layout(widths = c(2, 1, 1, 0.5,0.5))

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Figures/Lipidomics/")

ggsave("fa_pca_CSsurvival_plot.png",plot=fa_pca_cssurv_plot,
       dpi=600,units='in',width=17.5,height=7)


