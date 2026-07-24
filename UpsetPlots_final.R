

## script to make two upset plots that show the degree of overlap in both
## metabolites and lipid species across the five pairwise comparisons of interest
## -- one plot for upregulated things and one for downregulated things 
## tutorial = https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html

install.packages("UpSetR")
library(UpSetR)
BiocManager::install("ComplexHeatmap")
library("ComplexHeatmap")
install.packages("ComplexUpset")
library("ComplexUpset")

## first read in files of the stats / log2FC of the pairwise comparisons generated from
## the clustering scripts of metabolite and lipid stuff. 
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Output Files")
met_contrasts <- read.csv(file = "met_pairwise_contrasts_stats.csv")
head(met_contrasts)

lip_hg_contrasts <- read.csv(file = "lip_headgroup_pairwise_contrasts_stats.csv")
head(lip_hg_contrasts)

lip_fa_contrasts <- read.csv(file = "lip_fa_pairwise_contrasts_stats.csv")
head(lip_fa_contrasts)

## need to get the data organized with a column for the metabolites/lipids, then five binary (0,1)
## value columns -- one for each of the five pairwise comparisons. 
## but first need to sort by things higher or lower in response to each pairwise comparison 

## first rename the first column of each of these datasets so they can be rbind -ed together
colnames(met_contrasts)[1] <- "Outcome"
colnames(lip_hg_contrasts)[1] <- "Outcome"
colnames(lip_fa_contrasts)[1] <- "Outcome"

all_contrasts <- rbind(met_contrasts, lip_hg_contrasts)
all_contrasts <- rbind(all_contrasts, lip_fa_contrasts)

View(all_contrasts)

## need to add five new columns that say whether each outcome is upregulated significantly in response 
## to each pairwise comparisons (needs to have log2FC > 0 and p-value < 0.05 for that comparison)
## then repeat for downregulation for each comparison 

all_contrasts <- all_contrasts %>%
  mutate(Ncold_Nwarm_up_binary = case_when(
    Ncold_Nwarm_log2fc>0 & Ncold_Nwarm_padj<0.05 ~ 1,
    Ncold_Nwarm_log2fc<0 | Ncold_Nwarm_padj>0.05 ~ 0)) %>%
  mutate(Dcold_Dwarm_up_binary = case_when(
    Dcold_Dwarm_log2fc>0 & Dcold_Dwarm_padj<0.05 ~ 1,
    Dcold_Dwarm_log2fc<0 | Dcold_Dwarm_padj>0.05 ~ 0)) %>%
  mutate(Dbefore_Nbefore_up_binary = case_when(
    Dbefore_Nbefore_log2fc>0 & Dbefore_Nbefore_padj<0.05 ~ 1,
    Dbefore_Nbefore_log2fc<0 | Dbefore_Nbefore_padj>0.05 ~ 0)) %>%
  mutate(Dwarm_Nwarm_up_binary = case_when(
    Dwarm_Nwarm_log2fc>0 & Dwarm_Nwarm_padj<0.05 ~ 1,
    Dwarm_Nwarm_log2fc<0 | Dwarm_Nwarm_padj>0.05 ~ 0)) %>%
  mutate(Dcold_Ncold_up_binary = case_when(
    Dcold_Ncold_log2fc>0 & Dcold_Ncold_padj<0.05 ~ 1,
    Dcold_Ncold_log2fc<0 | Dcold_Ncold_padj>0.05 ~ 0)) %>%
  mutate(Ncold_Nwarm_down_binary = case_when(
    Ncold_Nwarm_log2fc<0 & Ncold_Nwarm_padj<0.05 ~ 1,
    Ncold_Nwarm_log2fc>0 | Ncold_Nwarm_padj>0.05 ~ 0)) %>%
  mutate(Dcold_Dwarm_down_binary = case_when(
    Dcold_Dwarm_log2fc<0 & Dcold_Dwarm_padj<0.05 ~ 1,
    Dcold_Dwarm_log2fc>0 | Dcold_Dwarm_padj>0.05 ~ 0)) %>%
  mutate(Dbefore_Nbefore_down_binary = case_when(
    Dbefore_Nbefore_log2fc<0 & Dbefore_Nbefore_padj<0.05 ~ 1,
    Dbefore_Nbefore_log2fc>0 | Dbefore_Nbefore_padj>0.05 ~ 0)) %>%
  mutate(Dwarm_Nwarm_down_binary = case_when(
    Dwarm_Nwarm_log2fc<0 & Dwarm_Nwarm_padj<0.05 ~ 1,
    Dwarm_Nwarm_log2fc>0 | Dwarm_Nwarm_padj>0.05 ~ 0)) %>%
  mutate(Dcold_Ncold_down_binary = case_when(
    Dcold_Ncold_log2fc<0 & Dcold_Ncold_padj<0.05 ~ 1,
    Dcold_Ncold_log2fc>0 | Dcold_Ncold_padj>0.05 ~ 0))

## now split by upregulation and downregulation
all_contrasts_up <- all_contrasts[,c(1,15:19)]
all_contrasts_down <- all_contrasts[,c(1,20:24)]

## rename columns to be better for figures
colnames(all_contrasts_up) <- c("Outcome", "ND Cold - ND Warm", "D Cold - D Warm",
                                "D Before - ND Before", "D Warm - ND Warm", "D Cold - ND Cold")
colnames(all_contrasts_down) <- c("Outcome", "ND Cold - ND Warm", "D Cold - D Warm",
                                "D Before - ND Before", "D Warm - ND Warm", "D Cold - ND Cold")

## add another column in each dataframe that says whether it's a metabolite or lipid

all_contrasts_up <- all_contrasts_up %>%
  mutate(responsetype = case_when(
    Outcome %in% met_contrasts$Outcome ~ "Metabolite",
    Outcome %in% lip_hg_contrasts$Outcome ~ "Lipid",
    Outcome %in% lip_fa_contrasts$Outcome ~ "Lipid"
  ))

all_contrasts_down <- all_contrasts_down %>%
  mutate(responsetype = case_when(
    Outcome %in% met_contrasts$Outcome ~ "Metabolite",
    Outcome %in% lip_hg_contrasts$Outcome ~ "Lipid",
    Outcome %in% lip_fa_contrasts$Outcome ~ "Lipid"
  ))

upset_plot_down <- UpSetR::upset(all_contrasts_down, 
              #order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE),
              order.by = "degree", decreasing = FALSE,
              sets = c("D Cold - ND Cold", "D Warm - ND Warm", "D Before - ND Before",
                       "D Cold - D Warm", "ND Cold - ND Warm"),
              queries = list(list(query=elements, params = list("responsetype","Lipid"), color="#79a48c" , active=T)),
              keep.order = TRUE,
              main.bar.color = "#a47991",
              sets.bar.color = "#a47991",
              mainbar.y.label = "Intersection Size",
              sets.x.label = "Set Size",
              point.size = 3,
              line.size = 1,
              mb.ratio = c(0.50,0.50),
              text.scale = c(1.3, 1, 1.3, 1, 1.3, 1),
              mainbar.y.max = 20,
              set_size.scale_max = 30)

upset_plot_down

upset_plot_up <- UpSetR::upset(all_contrasts_up, 
                                 #order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE),
                                 order.by = "degree", decreasing = FALSE,
                                 sets = c("D Cold - ND Cold", "D Warm - ND Warm", "D Before - ND Before",
                                          "D Cold - D Warm", "ND Cold - ND Warm"),
                                 queries = list(list(query=elements, params = list("responsetype","Lipid"), color="#79a48c" , active=T)),
                                 keep.order = TRUE,
                                 main.bar.color="#a47991",
                                 sets.bar.color = "#a47991",
                                 mainbar.y.label = "Intersection Size",
                                 sets.x.label = "Set Size",
                                 point.size = 3,
                                 line.size = 1,
                                 mb.ratio = c(0.50,0.50),
                                 text.scale = c(1.3, 1, 1.3, 1, 1.3, 1),
                                 mainbar.y.max = 20,
                                 set_size.scale_max = 30)

upset_plot_up$New_data
View(upset_plot_up)
## check plot splitting
View(subset(all_contrasts_up,all_contrasts_up$`ND Cold - ND Warm`>0 &
         all_contrasts_up$`D Cold - D Warm`==0 &
         all_contrasts_up$`D Before - ND Before`==0 &
         all_contrasts_up$`D Warm - ND Warm`==0 &
         all_contrasts_up$`D Cold - ND Cold`==0))

View(subset(all_contrasts_down,all_contrasts_down$`ND Cold - ND Warm`>0 &
              all_contrasts_down$`D Cold - D Warm`==0 &
              all_contrasts_down$`D Before - ND Before`==0 &
              all_contrasts_down$`D Warm - ND Warm`==0 &
              all_contrasts_down$`D Cold - ND Cold`==0))

## get split counts for each set size
View(subset(all_contrasts_up,all_contrasts_up$`ND Cold - ND Warm`>0)) ## 10 met, 18 lip
View(subset(all_contrasts_up,all_contrasts_up$`D Cold - D Warm`>0)) ## 6 met, 1 lip
View(subset(all_contrasts_up,all_contrasts_up$`D Before - ND Before`>0)) ## 11 met, 5 lip
View(subset(all_contrasts_up,all_contrasts_up$`D Warm - ND Warm`>0)) ## 2 met, 7 lip
View(subset(all_contrasts_up,all_contrasts_up$`D Cold - ND Cold`>0)) ## 8 met, 1 lip

View(subset(all_contrasts_down,all_contrasts_down$`ND Cold - ND Warm`>0)) ## 15 met, 4 lip
View(subset(all_contrasts_down,all_contrasts_down$`D Cold - D Warm`>0)) ## 11 met, 0 lip
View(subset(all_contrasts_down,all_contrasts_down$`D Before - ND Before`>0)) ## 19 met, 0 lip
View(subset(all_contrasts_down,all_contrasts_down$`D Warm - ND Warm`>0)) ## 12 met, 1 lip
View(subset(all_contrasts_down,all_contrasts_down$`D Cold - ND Cold`>0)) ## 13 met, 2 lip

## to calculate size of lipid part of bar: =(X*22.227)/30 where X is number of lipid responses

## save plots 
png("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Figures/upset_plot_up.png",res=600,width=6,height=4.5,units='in')
upset_plot_up
dev.off()

png("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Figures/upset_plot_down.png",res=600,width=6,height=4.5,units='in')
upset_plot_down
dev.off()

## note that the formatting of the final figure was modified manually using Inkscape
