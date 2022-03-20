
# packages ----------------------------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(impute))

# data cleaning of categorical variable ------------------------------------------

## load the data needed
imputedNum <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/final/Exposures_numeric_imputed.rds")
telomere_length <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Genomics_data_recoded.rds")
exposures <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final_onehot.rds") 

## Adding telomere length variable to dataframe
data <- exposures %>%
  mutate(AdjTSRatio = telomere_length$AdjTSRatio.0.0)

## discard participants without telomere length
data <- data %>% filter(!is.na(AdjTSRatio))

cat <- data[,119:342]
num <- data[,-(119:342)]

# count NA numbers in categorical data
NA_cat <- as.data.frame(colSums(is.na(cat)))
colnames(NA_cat) <- c("NAnumber")

# how many variables will remain in our categorical dataset after excluding variables with > 20% NA
cat_remain <- as.data.frame(subset(NA_cat, NAnumber<= 0.2*(nrow(cat)))) # 119/224 variables remained in our categorical dataset (sample size =__)
print(paste0(nrow(cat_remain), "/244 categorical variables will be kept under the threshold of 20% NA."))

# categorical dataset with < 20% NA
cat_selected <- cat %>% select(colnames = rownames(cat_remain))
colnames(cat_selected) <- rownames(cat_remain) 

# combine imputed numeric data with selected categorical data -----------
data <- cbind(imputedNum, cat_selected)

# get the final dataset
exposure_imputed_full <- na.omit(data)
print(paste0(nrow(exposure_imputed_full), " participants will be kept in our final exposure dataset."))
saveRDS(exposure_imputed_full, "/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/final/Exposures_imputed_full.rds")

