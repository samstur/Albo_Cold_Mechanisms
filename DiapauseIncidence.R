#### Script to calculate average and standard error diapause incidence for 
## egg replicates produced by long-day (16hr:8hr light:dark) exposed parental 
## cages versus by short-day (8hr:16hr light:dark) exposed parental cages.

### open data file containing diapause incidence results
setwd("~/Documents/Georgetown/Armbruster Lab/Cold Tolerance/Data Files for Submission")
data <- read.csv(file="DiapauseIncidenceCheck_finaldata.csv",header=TRUE) 
head(data)


### calculate diapause incidence +/- 1 standard errors for each group 
data_averages <- data %>% 
  group_by(diapause_treat) %>% 
  summarise(mean = mean(DI),
            std = (std.error(DI)))
View(data_averages)
