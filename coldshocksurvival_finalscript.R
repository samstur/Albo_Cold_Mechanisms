

### open data file containing cold shock survival results
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Data Files for Submission")
data <- read.csv(file="coldshocksurvival_finaldata.csv",header=TRUE) 

### load the necessary packages, installing first if necessary
library(ggplot2)
library("ggbeeswarm")
library("multcomp")
library(lme4)
library(stats)
library("plotrix")
library(dplyr)
library(PsychHelperFunctions)
library("glmmTMB")
library("car")
library("emmeans")

### rename cold shock treatment factors 
data$CS_treat[data$CS_treat == '-9'] <- '-9°C for 20hrs'
data$CS_treat[data$CS_treat == '-11'] <- '-11°C for 20hrs'

### change order of factor levels for visualization purposes
data$acc_treat <- factor(data$acc_treat, levels=c('Before Acclimation','Warm Acclimated (15d)','Cold Acclimated (15d)','Cold Acclimated (22d)','De-acclimated'))
data$diapause_treat <- factor(data$diapause_treat, levels=c('Non-diapause','Diapause'))
data$CS_treat<- factor(data$CS_treat, levels = c("Control", "-9°C for 20hrs","-11°C for 20hrs"))

#### get counts per group for sample size table
table <- data %>% count(diapause_treat,acc_treat, CS_treat)
colnames(table) <- c("Diapause Treatment", "Acclimation Treatment", "Cold Shock Treatment", "N")
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Output Files")
write.csv(table, file = "coldshock_samplesizes.csv", row.names=F)

### calculate average live hatch rates +/- 2 standard errors for each group (not yet scaled to control rates)
data_averages <- data %>% 
  group_by(diapause_treat,acc_treat, CS_treat) %>% 
  summarise(mean = mean(livehatchrate),
            std = 2*(std.error(livehatchrate)))
## compare control average hatch rates
data_averages_control <- subset(data_averages, data_averages$CS_treat=="Control")
mean(data_averages_control$mean) ## 0.947
min(data_averages_control$mean) ## 0.918
max(data_averages_control$mean) ## 0.970

## get table of reps per parental cage and split across treatments
table_parentalcages <- data %>% count(parental_cage_num,diapause_treat,acc_treat, CS_treat)
View(table_parentalcages)
colnames(table_parentalcages) <- c("Parental Cage ID","Diapause Treatment", "Acclimation Treatment", "Cold Shock Treatment", "N")
length(unique(table_parentalcages$`Parental Cage ID`)) #14


####### plot just control hatch rates for all treatment groups 
plot_CS_survival_control <- ggplot(data[data$CS_treat=="Control",],
                           aes(x=acc_treat,y=livehatchrate,color=acc_treat,group=acc_treat))+
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               fun.args = list(mult = 2), width = 0.15,size=0.6,colour="black")+
  stat_summary(fun.data=mean_se,size=0.6,color="black",aes(fill=acc_treat),shape=21)+
  geom_quasirandom(alpha=0.5,size=0.7)+
  scale_fill_manual(values=c("#747474","#C00000","#4E95D9","#0331E5","#A981E1"))+
  scale_colour_manual(values=c("#747474","#C00000","#4E95D9","#0331E5","#A981E1"))+
  xlab("Acclimation Treatment")+
  ylab("Hatch Rate")+
  theme_bw()+
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))+
  theme(axis.title.x = element_text(size = 11,face="bold"),
        axis.text.x = element_text(size = 10,angle=45,vjust=1,hjust=1),
        axis.title.y = element_text(size = 11,face="bold"),
        axis.text.y = element_text(size=10),
        strip.text = element_text(size = 11,face="bold"),
        strip.text.x = element_text(size = 11, margin = margin(t = 0.7, b = 0.7, unit = "pt"),hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = NA),
        legend.position = "")+
  scale_x_discrete(limits=c("Before Acclimation","Warm Acclimated (15d)","Cold Acclimated (15d)",
                            'Cold Acclimated (22d)','De-acclimated'),
                   labels=c("Before Acc.","Warm Acc. 15-d","Cold Acc. 15-d",
                            'Cold Acc. 22-d','De-acc.22-d'))+
  facet_grid(cols = vars(diapause_treat))
plot_CS_survival_control

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Figures/CS_survival/")
ggsave("CS_survival_control.png",plot=plot_CS_survival_control,dpi=600,units='in',
       width=5,height=6)

######## next need to scale all CS hatch rates to the relevant control group 

## compute averages for each diapause and acclimation combo
control_averages <- data %>%
  filter(CS_treat=="Control") %>%
  group_by(diapause_treat, acc_treat) %>%
  summarise(livehatchrate_mean=mean(livehatchrate))
View(control_averages)

## so now I have the control averages for each diapause_treat x acc_treat combination
## now I need to go through each row of the original dataframe, then make it pick the right average value and divide the live hatch rate by that
## then save that value to a new column
data$scaled_hatchrate <- NA
data$control_hatchrate <- NA
for (i in 1:nrow(data)) {
  diap_id <- data[i,"diapause_treat"]
  acc_id <- data[i,"acc_treat"]
  av <- subset(control_averages,control_averages$diapause_treat==diap_id & control_averages$acc_treat==acc_id)$livehatchrate_mean
  scaled <- data[i,"livehatchrate"] / av
  data$scaled_hatchrate[i] <- scaled
  data$control_hatchrate[i] <- av
}

View(data)

####### get counts and averages per group for scaled hatch rates 
data_averages_scaled <- data %>% 
  filter(CS_treat != "Control") %>%
  group_by(diapause_treat,acc_treat, CS_treat) %>% 
  summarise(CS_survival_mean = mean(scaled_hatchrate),
            CS_survival_std = 2*(std.error(scaled_hatchrate)))
View(data_averages_scaled)
mean(subset(data_averages_scaled, data_averages_scaled$CS_treat=="-9°C for 20hrs")$CS_survival_mean)

## compute averages for -9 and -11 exposures
data_minus9_minus11_averages_scaled <- data %>% 
  filter(CS_treat != "Control") %>%
  group_by(CS_treat) %>% 
  summarise(CS_survival_mean = mean(scaled_hatchrate),
            CS_survival_std = 2*(std.error(scaled_hatchrate)))
View(data_minus9_minus11_averages_scaled)

### save this for use in the metabolomics and lipidomics scripts 
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Input Files")
write.csv(data_averages_scaled, file = "coldshocksurvival_scaledaverages.csv", row.names=F)


### plot scaled averages
plot_CS_survival_scaled <- ggplot(data[data$CS_treat!="Control",],
                           aes(x=acc_treat,y=scaled_hatchrate,color=acc_treat,group=acc_treat))+
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               fun.args = list(mult = 2), width = 0.2,size=0.8,colour="black")+
  stat_summary(fun.data=mean_se,size=0.7,color="black",aes(fill=acc_treat),shape=21)+
  geom_quasirandom(alpha=0.5,size=0.8)+
  scale_fill_manual(values=c("#747474","#C00000","#4E95D9","#0331E5","#A981E1"))+
  scale_colour_manual(values=c("#747474","#C00000","#4E95D9","#0331E5","#A981E1"))+
  xlab("Acclimation Treatment")+
  ylab("Survival (scaled)")+
  theme_bw()+
  scale_y_continuous(limits=c(0,1.1),breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))+
  theme(axis.title.x = element_text(size = 11,face="bold"),
        axis.text.x = element_text(size = 10,angle=45,vjust=1,hjust=1),
        axis.title.y = element_text(size = 11,face="bold"),
        axis.text.y = element_text(size=10),
        strip.text = element_text(size = 11,face="bold"),
        strip.text.x = element_text(size = 11, margin = margin(t = 0.7, b = 0.7, unit = "pt"),hjust = 0.5),
        strip.background = element_rect(fill = "white", colour = NA),
        legend.position = "")+
  scale_x_discrete(limits=c("Before Acclimation","Warm Acclimated (15d)","Cold Acclimated (15d)",
                            'Cold Acclimated (22d)','De-acclimated'),
                   labels=c("Before Acc.","Warm Acc. 15-d","Cold Acc. 15-d",
                            'Cold Acc. 22-d','De-acc.22-d'))+
  facet_grid(cols = vars(diapause_treat),rows = vars(CS_treat))

plot_CS_survival_scaled

setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Figures/CS_survival/")
ggsave("CS_survival_scaled.png",plot=plot_CS_survival_scaled,dpi=600,units='in',
       width=5.5,height=7)



###### stat analysis
## first perform analysis of -9C dataset
data_minus9 <- subset(data,data$CS_treat == "-9°C for 20hrs")
data_minus9$scaled_hatchrate <- ifelse(data_minus9$scaled_hatchrate>1,1,data_minus9$scaled_hatchrate)
model.prop_minus9 <- glmmTMB(scaled_hatchrate ~ diapause_treat*acc_treat + days_postoviposition
                             + (1|parental_cage_num),
                             data = data_minus9,
                             family=betabinomial(link="logit"),
                             weights = egg_count)
summary(model.prop_minus9)
Anova(model.prop_minus9) ## sig effects of diapause, acclimation, and the interaction
df.residual(model.prop_minus9) ## 234 residual df
posthoc_minus9 <- emmeans(model.prop_minus9, list(pairwise ~ diapause_treat*acc_treat), adjust = "BH")
posthoc_minus9_df <- as.data.frame(posthoc_minus9$`pairwise differences of diapause_treat, acc_treat`)
View(posthoc_minus9_df)
colnames(posthoc_minus9_df)[1] <- "contrast"
posthoc_minus9_df <- posthoc_minus9_df %>%
  dplyr::select(contrast, estimate, SE, z.ratio, p.value)
## save stats for posthoc tests
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Output Files")
write.csv(posthoc_minus9_df, file = "CS_minus9_posthocstats.csv", row.names=F)


## then perform analysis of -11C dataset
### issue where ND Before Acc group is all 0s. This causes complete separation which
## messes with the model. So, I will first look for effects of diapause and cold acclimation
## and the interaction without the before acc groups (either ND or D Before Acc).
## Then I will create another model that is based on total treatment combination, and
## excludes just the ND Before Acc group. I will use this model for posthoc contrasts
data_minus11 <- subset(data,data$CS_treat == "-11°C for 20hrs")
data_minus11$scaled_hatchrate <- ifelse(data_minus11$scaled_hatchrate>1,1,data_minus11$scaled_hatchrate)
#model.prop_minus11 <- glmmTMB(scaled_hatchrate ~ diapause_treat*acc_treat + days_postoviposition
#                             + (1|parental_cage_num),
#                             data = data_minus11,
#                             family=betabinomial(link="logit"),
#                             weights = egg_count)
#summary(model.prop_minus11)
#Anova(model.prop_minus11) ## significant effect of diapause, acclimation, days_postoviposition, but not the interaction of diap:acc (p=0.13468)
#emmeans(model.prop_minus11, list(pairwise ~ diapause_treat*acc_treat), adjust = "BH")

model.prop_minus11 <- glmmTMB(scaled_hatchrate ~ diapause_treat*acc_treat + days_postoviposition
                              + (1|parental_cage_num),
                              data = data_minus11[data_minus11$acc_treat!="Before Acclimation",],
                              family=betabinomial(link="logit"),
                              weights = egg_count)
summary(model.prop_minus11)
Anova(model.prop_minus11) ## sig effect of diapause, acc, and dayspostov, but not diapause:acc (p = 0.076)
df.residual(model.prop_minus11) ## 173 residual df
View(data_minus11)

data_minus11$totaltreatment <- paste(data_minus11$diapause_treat, data_minus11$acc_treat)
model.prop_minus11_totaltreatment <- glmmTMB(scaled_hatchrate ~ totaltreatment + days_postoviposition
                              + (1|parental_cage_num),
                              data = data_minus11[data_minus11$totaltreatment!="Non-diapause Before Acclimation",],
                              family=betabinomial(link="logit"),
                              weights = egg_count)
summary(model.prop_minus11_totaltreatment)
Anova(model.prop_minus11_totaltreatment) ## sig effect of totaltreatment

posthoc_minus11 <- emmeans(model.prop_minus11_totaltreatment, list(pairwise ~ totaltreatment), adjust = "BH")
posthoc_minus11_df <- as.data.frame(posthoc_minus11$`pairwise differences of totaltreatment`)
colnames(posthoc_minus11_df)[1] <- "contrast"
posthoc_minus11_df <- posthoc_minus11_df %>%
  dplyr::select(contrast, estimate, SE, z.ratio, p.value)
## save stats for posthoc tests
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Output Files")
write.csv(posthoc_minus11_df, file = "CS_minus11_posthocstats.csv", row.names=F)




