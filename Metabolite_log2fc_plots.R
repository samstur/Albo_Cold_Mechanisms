
### In this script, I want to create lollipop plots that show the magnitude of
# the log2FC and the size of the p-value for each metabolite that was significantly
# different in each pairwise contrast.

# load packages
library(naniar)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
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

### read in pairwise stats from metabolite_clustering script
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Output Files")
met_contrasts <- read.csv(file = "met_pairwise_contrasts_stats.csv")

## remove _MS from metabolite names 
met_contrasts$metabolite <- met_contrasts$metabolite %>%
  gsub("_MS","",.)

## add in abbreviated met names 
setwd("/Users/samsturiale/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Input Files")
met_abbr <- read.csv(file="met_abbreviations.csv")
length(met_abbr$fullname)
met_abbr <- met_abbr %>%
  filter(met_abbr$fullname %in% met_contrasts$metabolite)
met_abbr <- rename(met_abbr, metabolite = fullname, met_abbr = abbr)
met_contrasts <- merge(met_contrasts,met_abbr, by = "metabolite")

## only keep metabolites significantly different between NCold and NWarm
sig_Ncold_Nwarm <- met_contrasts %>%
  filter(Ncold_Nwarm_padj < 0.05)

## only keep metabolites significantly different between DCold and DWarm
sig_Dcold_Dwarm <- met_contrasts %>%
  filter(Dcold_Dwarm_padj < 0.05)

## get the list of metabolites that were significantly affected by cold acclimation in both 
# diapause and non-diapause eggs. These will be bolded in the figure
sig_coldacc <- met_contrasts %>%
  filter(metabolite %in% intersect(sig_Ncold_Nwarm$metabolite,sig_Dcold_Dwarm$metabolite))

sig_Ncold_Nwarm$shared_acc <- ifelse(sig_Ncold_Nwarm$metabolite %in% sig_coldacc$metabolite, "yes", "no")
sig_Dcold_Dwarm$shared_acc <- ifelse(sig_Dcold_Dwarm$metabolite %in% sig_coldacc$metabolite, "yes", "no")

## arrange metabolites in descending order by Log2FC
sig_Ncold_Nwarm <- sig_Ncold_Nwarm %>%
  arrange(desc(Ncold_Nwarm_log2fc))

# lock in factor level order
sig_Ncold_Nwarm$met_abbr <- factor(sig_Ncold_Nwarm$met_abbr, levels = sig_Ncold_Nwarm$met_abbr)

### same but for other contrast with diapause eggs
sig_Dcold_Dwarm <- sig_Dcold_Dwarm %>%
  arrange(desc(Dcold_Dwarm_log2fc))

# lock in factor level order
sig_Dcold_Dwarm$met_abbr <- factor(sig_Dcold_Dwarm$met_abbr, levels = sig_Dcold_Dwarm$met_abbr)


########## now diapause-related contrasts 

sig_Dbefore_Nbefore <- met_contrasts %>%
  filter(Dbefore_Nbefore_padj < 0.05)
sig_Dwarm_Nwarm <- met_contrasts %>%
  filter(Dwarm_Nwarm_padj < 0.05)
sig_Dcold_Ncold <- met_contrasts %>%
  filter(Dcold_Ncold_padj < 0.05)

overlap_list <- intersect(sig_Dbefore_Nbefore$metabolite, sig_Dwarm_Nwarm$metabolite)
overlap_list <- intersect(overlap_list,sig_Dcold_Ncold$metabolite)

sig_diapause <- met_contrasts %>%
  filter(metabolite %in% overlap_list)

sig_Dbefore_Nbefore$shared_diapause <- ifelse(sig_Dbefore_Nbefore$metabolite %in% sig_diapause$metabolite, "yes", "no")
sig_Dwarm_Nwarm$shared_diapause <- ifelse(sig_Dwarm_Nwarm$metabolite %in% sig_diapause$metabolite, "yes", "no")
sig_Dcold_Ncold$shared_diapause <- ifelse(sig_Dcold_Ncold$metabolite %in% sig_diapause$metabolite, "yes", "no")

sig_Dbefore_Nbefore <- sig_Dbefore_Nbefore %>%
  arrange(desc(Dbefore_Nbefore_log2fc))

# lock in factor level order
sig_Dbefore_Nbefore$met_abbr <- factor(sig_Dbefore_Nbefore$met_abbr, levels = sig_Dbefore_Nbefore$met_abbr)


sig_Dwarm_Nwarm <- sig_Dwarm_Nwarm %>%
  arrange(desc(Dwarm_Nwarm_log2fc))
# lock in factor level order
sig_Dwarm_Nwarm$met_abbr <- factor(sig_Dwarm_Nwarm$met_abbr, levels = sig_Dwarm_Nwarm$met_abbr)

sig_Dcold_Ncold <- sig_Dcold_Ncold %>%
  arrange(desc(Dcold_Ncold_log2fc))
# lock in factor level order
sig_Dcold_Ncold$met_abbr <- factor(sig_Dcold_Ncold$met_abbr, levels = sig_Dcold_Ncold$met_abbr)


############## Now make lollipop plots

sig_Ncold_Nwarm$Ncold_Nwarm_padj_log10 <- -log10(sig_Ncold_Nwarm$Ncold_Nwarm_padj)
lp_Ncold_Nwarm <- ggplot(sig_Ncold_Nwarm, aes(x = met_abbr, y = Ncold_Nwarm_log2fc)) +
  geom_segment(aes(x = met_abbr, xend = met_abbr, y = 0, yend = Ncold_Nwarm_log2fc)) +
  geom_point(pch = 21, color = "black", aes(fill = Ncold_Nwarm_log2fc>0, size = Ncold_Nwarm_padj_log10)) +
 # geom_hline(yintercept = c(-0.58, 0.58), col = "darkgray", linetype = 'dashed') +
  coord_flip() +
  theme_bw() +
  ggtitle("(1) ND Cold - ND Warm")+
  scale_fill_manual(values = c("#7CBFAC","#7C7FBF"), guide = "none")+
  scale_size_continuous(breaks = c(1,5,10,15),
                        limits = c(1,15),
                        range = c(1,7))+
  ylab(bquote(Log[2]*FC))+
  labs(size = bquote(-Log[10]*(p[adj])))+
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 11,face="bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size=11,face = c("plain","plain","plain","plain","plain",
                                                    "bold","bold","plain","bold","plain",
                                                    "plain","plain","plain","plain","bold",
                                                    "bold","bold","plain","bold","plain",
                                                    "bold","plain","bold","plain","bold")),
        strip.text = element_text(size = 11,face="bold"),
        plot.title = element_text(size=11,face="bold"))
  
lp_Ncold_Nwarm

sig_Dcold_Dwarm$Dcold_Dwarm_padj_log10 <- -log10(sig_Dcold_Dwarm$Dcold_Dwarm_padj)
lp_Dcold_Dwarm <- ggplot(sig_Dcold_Dwarm, aes(x = met_abbr, y = Dcold_Dwarm_log2fc)) +
  geom_segment(aes(x = met_abbr, xend = met_abbr, y = 0, yend = Dcold_Dwarm_log2fc)) +
  geom_point(pch = 21, color = "black", aes(fill = Dcold_Dwarm_log2fc>0, size = Dcold_Dwarm_padj_log10)) +
  # geom_hline(yintercept = c(-0.58, 0.58), col = "darkgray", linetype = 'dashed') +
  coord_flip() +
  theme_bw() +
  ggtitle("(2) D Cold - D Warm")+
  scale_fill_manual(values = c("#7CBFAC","#7C7FBF"), guide = "none")+
  scale_size_continuous(breaks = c(1,5,10,15),
                        limits = c(1,15),
                        range = c(1,7))+
  ylab(bquote(Log[2]*FC))+
  labs(size = bquote(-Log[10]*(p[adj])))+
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 11,face="bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size=11,face = c("bold","bold","plain","plain","plain",
                                                    "plain","plain","bold","bold","plain",
                                                    "bold","bold","bold","bold","bold",
                                                    "bold")),
        strip.text = element_text(size = 11,face="bold"),
        plot.title = element_text(size=11,face="bold"))
lp_Dcold_Dwarm


sig_Dbefore_Nbefore$Dbefore_Nbefore_padj_log10 <- -log10(sig_Dbefore_Nbefore$Dbefore_Nbefore_padj)
lp_Dbefore_Nbefore <- ggplot(sig_Dbefore_Nbefore, aes(x = met_abbr, y = Dbefore_Nbefore_log2fc)) +
  geom_segment(aes(x = met_abbr, xend = met_abbr, y = 0, yend = Dbefore_Nbefore_log2fc)) +
  geom_point(pch = 21, color = "black", aes(fill = Dbefore_Nbefore_log2fc>0, size = Dbefore_Nbefore_padj_log10)) +
  # geom_hline(yintercept = c(-0.58, 0.58), col = "darkgray", linetype = 'dashed') +
  coord_flip() +
  theme_bw() +
  ggtitle("(3) D Before - ND Before")+
  scale_fill_manual(values = c("#7CBFAC","#7C7FBF"), guide = "none")+
  scale_size_continuous(breaks = c(1,5,10,15),
                        limits = c(1,15),
                        range = c(1,7))+
  ylab(bquote(Log[2]*FC))+
  labs(size = bquote(-Log[10]*(p[adj])))+
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 11,face="bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size=11,face = c("plain","plain","bold","plain","plain",
                                                    "plain","bold","plain","plain","plain",
                                                    "plain","plain","bold","plain","plain",
                                                    "plain","plain","plain","bold","bold",
                                                    "bold","plain","plain","plain","bold",
                                                    "bold","plain","bold","plain","plain")),
        strip.text = element_text(size = 11,face="bold"),
        plot.title = element_text(size=11,face="bold"))
lp_Dbefore_Nbefore

sig_Dwarm_Nwarm$Dwarm_Nwarm_padj_log10 <- -log10(sig_Dwarm_Nwarm$Dwarm_Nwarm_padj)
lp_Dwarm_Nwarm <- ggplot(sig_Dwarm_Nwarm, aes(x = met_abbr, y = Dwarm_Nwarm_log2fc)) +
  geom_segment(aes(x = met_abbr, xend = met_abbr, y = 0, yend = Dwarm_Nwarm_log2fc)) +
  geom_point(pch = 21, color = "black", aes(fill = Dwarm_Nwarm_log2fc>0, size = Dwarm_Nwarm_padj_log10)) +
  # geom_hline(yintercept = c(-0.58, 0.58), col = "darkgray", linetype = 'dashed') +
  coord_flip() +
  theme_bw() +
  ggtitle("(4) D Warm - ND Warm")+
  scale_fill_manual(values = c("#7CBFAC","#7C7FBF"), guide = "none")+
  scale_size_continuous(breaks = c(1,5,10,15),
                        limits = c(1,15),
                        range = c(1,7))+
  ylab(bquote(Log[2]*FC))+
  labs(size = bquote(-Log[10]*(p[adj])))+
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 11,face="bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size=11,face = c("bold","bold","plain","plain","bold",
                                                    "plain","plain","bold","plain","bold",
                                                    "bold","bold","bold","bold")),
        strip.text = element_text(size = 11,face="bold"),
        plot.title = element_text(size=11,face="bold"))
lp_Dwarm_Nwarm

sig_Dcold_Ncold$Dcold_Ncold_padj_log10 <- -log10(sig_Dcold_Ncold$Dcold_Ncold_padj)
lp_Dcold_Ncold <- ggplot(sig_Dcold_Ncold, aes(x = met_abbr, y = Dcold_Ncold_log2fc)) +
  geom_segment(aes(x = met_abbr, xend = met_abbr, y = 0, yend = Dcold_Ncold_log2fc)) +
  geom_point(pch = 21, color = "black", aes(fill = Dcold_Ncold_log2fc>0, size = Dcold_Ncold_padj_log10)) +
  # geom_hline(yintercept = c(-0.58, 0.58), col = "darkgray", linetype = 'dashed') +
  coord_flip() +
  theme_bw()+
  ggtitle("(5) D Cold - ND Cold")+
  scale_fill_manual(values = c("#7CBFAC","#7C7FBF"), guide = "none")+
  scale_size_continuous(breaks = c(1,5,10,15),
                        limits = c(1,15),
                        range = c(1,7))+
  ylab(bquote(Log[2]*FC))+
  labs(size = bquote(-Log[10]*(p[adj])))+
  theme(legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 11,face="bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size=11,face = c("plain","plain","plain","bold","bold",
                                                    "plain","plain","plain","plain","bold",
                                                    "plain","bold","bold","bold","plain",
                                                    "plain","bold","plain","plain","bold",
                                                    "bold")),
        strip.text = element_text(size = 11,face="bold"),
        plot.title = element_text(size=11,face="bold"))
lp_Dcold_Ncold

design <- "
  123
  456
"
lp_contrasts <- lp_Ncold_Nwarm + lp_Dcold_Dwarm + guide_area() + lp_Dbefore_Nbefore + lp_Dwarm_Nwarm + lp_Dcold_Ncold + plot_layout(design = design) + plot_layout(guides = "collect")
lp_contrasts

lp_contrasts_acc <- lp_Ncold_Nwarm + lp_Dcold_Dwarm + plot_layout(guides = "collect")
lp_contrasts_diap <- lp_Dbefore_Nbefore + lp_Dwarm_Nwarm + lp_Dcold_Ncold

setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Figures/Metabolomics")
ggsave("lp_contrasts.png",plot=lp_contrasts,
       dpi=600,units='in',width=11,height=11)

setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Figures/Metabolomics")
ggsave("lp_contrasts_acc.png",plot=lp_contrasts_acc,
       dpi=600,units='in',width=8,height=6)
ggsave("lp_contrasts_diap.png",plot=lp_contrasts_diap,
       dpi=600,units='in',width=11,height=6)

