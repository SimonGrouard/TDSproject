
# packages ----------------------------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(impute))

#loading the data--------------------------
telomere_length <- data.frame(readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Genomics_data_recoded.rds"))
exposures <- data.frame(readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final_onehot.rds"))

## Adding telomere length variable to dataframe
data <- exposures %>%
  mutate(AdjTSRatio = telomere_length$AdjTSRatio.0.0)
rownames(data) <- rownames(telomere_length)

## discard participants without telomere length
data <- data %>% filter(!is.na(AdjTSRatio))

# select numeric variable for imputation: NA <= 20% in each columns ----------------------------------------------

## NA columns
NA_num <- as.data.frame(colSums(is.na(data)))
colnames(NA_num) <- c('NAnumber_num')
NA_num$names <- rownames(NA_num)
saveRDS(NA_num, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/ImputeResPath/NA_num.rds")

NA_ratio <- as.data.frame(round(NA_num$NAnumber_num/dim(data)[1], digits = 2))
colnames(NA_ratio) <- c('NAratio_num')
rownames(NA_ratio) <- NA_num$names
NA_ratio$names <- NA_num$names
NA_ratio <- arrange(NA_ratio, desc(NA_ratio$NAratio_num)) # rank by NA%
saveRDS(NA_ratio, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/ImputeResPath/NA_ratio.rds")

SelVarName <- as.data.frame(subset(NA_ratio, NA_ratio$NAratio_num <= 0.2))
subdata <- data[,colnames(data)%in%SelVarName$names]
ExpUnimp <- na.omit(subdata)

saveRDS(ExpUnimp, "/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/final/Exposures_unimputed.rds")

