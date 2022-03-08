# Imputation of exposure variables
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(mice))

## Parameters

args=commandArgs(trailingOnly=TRUE)
data_path=toString(args[1])
nchunks=as.numeric(args[2])
ichunk=as.numeric(args[3])


# Loading data after onehot encoding
data <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final_onehot.rds")
df <- data[, which(colMeans(!is.na(data)) < 0.8)] # drop variables with NA proportions > 80%

# imputation
imputedDf <- impute.knn(df, k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)

# Saving the imputed dataset
saveRDS(imputedDf, "/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Exposures_covariates_onehot_imputed.rds")