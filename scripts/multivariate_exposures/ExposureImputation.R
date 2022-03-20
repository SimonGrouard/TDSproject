
# packages ----------------------------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))
suppressMessages(library(impute))


# parameters --------------------------------------------------------------
args=commandArgs(trailingOnly=TRUE)
data_path=toString(args[1])
nchunks=as.numeric(args[2])
ichunk=as.numeric(args[3])

#loading the data--------------------------
telomere_length <- data.frame(readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Genomics_data_recoded.rds"))
exposures <- data.frame(readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final_onehot.rds"))

## Adding telomere length variable to dataframe
data <- exposures %>%
  mutate(AdjTSRatio = telomere_length$AdjTSRatio.0.0)
rownames(data) <- rownames(telomere_length)


## discard participants without telomere length
data <- data %>% filter(!is.na(AdjTSRatio))

cat <- data[,119:342]
num <- data[,-(119:342)]

# select numeric variable for imputation: NA <= 20% in each columns ----------------------------------------------

## NA columns
NA_num <- as.data.frame(colSums(is.na(num)))
colnames(NA_num) <- c('NAnumber_num')
NA_num$names <- rownames(NA_num)
saveRDS(NA_num, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/ImputeResPath/NA_num.rds")

NA_ratio <- as.data.frame(round(NA_num$NAnumber_num/dim(num)[1], digits = 2))
colnames(NA_ratio) <- c('NAratio_num')
rownames(NA_ratio) <- NA_num$names
NA_ratio$names <- NA_num$names
NA_ratio <- arrange(NA_ratio, desc(NA_ratio$NAratio_num)) # rank by NA%
saveRDS(NA_ratio, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/ImputeResPath/NA_ratio.rds")

imputeVarName <- as.data.frame(subset(NA_ratio, NA_ratio$NAratio_num <= 0.2))
rownames(imputeVarName) <- imputeVarName$names
saveRDS(imputeVarName, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/ImputeResPath/imputeVarName.rds")

numToimpute <- num[, colnames(num) %in% imputeVarName$names]
colnames(numToimpute) <- imputeVarName$names
saveRDS(numToimpute, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/ImputeResPath/numToimpute.rds")

print(paste0(dim(numToimpute)[2], " numeric exposure variables will be imputed."))

# numToimpute <- numToimpute[sample(1:nrow(numToimpute), 1000), ]
numToimpute <- as.matrix(numToimpute)

# Imputation-------------------------
## using impute package
## rowmax = 1: no limitation for NA% in each row
## maxp: The largest block of genes imputed using the knn algorithm inside impute.knn (default 1500); larger blocks are divided by two-means clustering (recursively) prior to imputation. If maxp=p, only knn imputation is done.
imputedNum <- impute.knn(numToimpute, k = 10, rowmax = 1, colmax = 0.21, maxp = dim(numToimpute)[1], rng.seed=362436069)
imputedNum <-imputedNum$data

# saveRDS(imputedNum, "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/ImputeResPath/Exposures_numeric_imputed.rds")
saveRDS(imputedNum, "/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/final/Exposures_numeric_imputed.rds")