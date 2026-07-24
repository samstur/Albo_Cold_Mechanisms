
## Cold shock survival
Aedes albopictus eggs were exposed to two different diapause treatments (diapause and nondiapause) and five different acclimation treatments (before acclimation, warm acclimation-15d, cold acclimation-15d, de-acclimation-22d, cold acclimation-22d) in a full factorial manner. Then they were exposed to either no cold shock (control), a cold shock of -9C for 20hours, or a cold shock of -11C for 20 hours. Cold shock survival was measured as the proportion of each replicate (i.e., group of ~35 eggs) that hatched successfully.

Diapause incidence was measured in a separate subset of eggs to ensure that diapause treatments worked as expected: [link](https://github.com/samstur/Albo_Cold_Mechanisms/blob/main/Raw_Data_Files/DiapauseIncidenceCheck_finaldata.csv)

The raw data file used for this analysis is here: [link](https://github.com/samstur/Albo_Cold_Mechanisms/blob/main/Raw_Data_Files/coldshocksurvival_finaldata.csv)

The script used for statistical analyses, data summary, and data visualization is here: [link](https://github.com/samstur/Albo_Cold_Mechanisms/blob/main/Raw_Data_Files/coldshocksurvival_finalscript.R)

This data file contains the mean survival after -9C and -11C exposure for each diapause and acclimation treatment group (after scaling for control hatch rates for each group) : [link](https://github.com/samstur/Albo_Cold_Mechanisms/blob/main/Raw_Data_Files/coldshocksurvival_scaledaverages.csv)

## Metabolite profile
Metabolite profiles (i.e., relative levels of 45 metabolites of a priori interest) were generated for six treatment groups: ND Before, ND Warm, ND Cold, D Before, D Warm, D Cold by Vladimír Košťál and Martin Moos.

The raw data file used for this analysis is here: [link](https://github.com/samstur/Albo_Cold_Mechanisms/blob/main/Raw_Data_Files/Metabolite_RawData.csv)

The script used for conducting a principal components analysis and follow-up statistical analyses and data visualization: [link](https://github.com/samstur/Albo_Cold_Mechanisms/blob/main/Metabolite_Analysis_PCA_finalscript.R)

The script used for testing whether individual metabolites are significantly affected by diapause and cold acclimation using a planned set of five pairwise contrasts. Generates a datasheet with p-values for each pairwise contrast for each metabolite as well as Log2FC values. Also generates the metabolite heatmap: [link](https://github.com/samstur/Albo_Cold_Mechanisms/blob/main/Metabolite_IndMetTests.R)

This is the list of metabolite abbreviations used in the script for data visualization: [link](https://github.com/samstur/Albo_Cold_Mechanisms/blob/main/Raw_Data_Files/met_abbreviations.csv)

This file contains information on the sample IDs, tube labels, diapause and acclimation treatments, and parental cage numbers for each metabolite and lipid sample: [link](https://github.com/samstur/Albo_Cold_Mechanisms/blob/main/Raw_Data_Files/SampleInfo_Fall2024.csv)

## Lipid profile
Phospholipid classes were quantified relatively and fatty acyl membrane ester species were quantified absolutely for six treatment groups: ND Before, ND Warm, ND Cold, D Before, D Warm, D Cold by Vladimír Košťál and Martin Moos.

The data file used for this analysis is here:

The script used for statistical analyses, data summary, and data visualization is here:
