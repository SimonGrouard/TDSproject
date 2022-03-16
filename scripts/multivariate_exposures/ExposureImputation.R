

# packages ----------------------------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(bnstruct))


# parameters --------------------------------------------------------------
args=commandArgs(trailingOnly=TRUE)
data_path=toString(args[1])
nchunks=as.numeric(args[2])
ichunk=as.numeric(args[3])

#loading the data--------------------------
telomere_length <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Genomics_data_recoded.rds")
exposures <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final_onehot.rds") 

## Adding telomere length variable to dataframe
data <- exposures %>%
  mutate(AdjTSRatio = telomere_length$AdjTSRatio.0.0)

cat <- data[,119:342]
num <- data[,1:118] %>%
  mutate(AdjTSRatio = telomere_length$AdjTSRatio.0.0)


# select data with <= 20% NA ----------------------------------------------
NA_num <- as.data.frame(colSums(is.na(num)))
colnames(NA_num) <- c('NAnumber_num')

NA_ratio <- as.data.frame(round(NA_num/502461, digits = 2))
colnames(NA_ratio) <- c('NAratio_num')
NA_ratio <- arrange(NA_ratio, desc(NA_ratio$NAratio_num)) # rank by NA%

imputeVarName <- subset(NA_ratio, NA_ratio$NAratio_num <= 0.2)
numToimpute <- num[, colnames(num) %in% rownames(imputeVarName)]

## numToimpute <- numToimpute[sample(1:nrow(numToimpute), 1000), ]
numToimpute <- as.matrix(numToimpute)


# Imputation-------------------------
impData <- knn.impute(
  numToimpute,
  k = 10,
  cat.var = NULL,
  to.impute = 1:nrow(numToimpute),
  using = 1:nrow(numToimpute)
)
impData <- data.frame(impData)
colnames(impData) <- colnames(numToimpute)

saveRDS(impData,("/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/exposure_imputed.rds"))

