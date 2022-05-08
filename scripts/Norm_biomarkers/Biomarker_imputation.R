#packages-------------------
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(bnstruct))


## Parameters----------------------------------
args=commandArgs(trailingOnly=TRUE)
data_path=toString(args[1])
nchunks=as.numeric(args[2])
ichunk=as.numeric(args[3])


#loading the data--------------------------
Genomics_data_recoded <- readRDS('../../extraction_and_recording/outputs/final/Genomics_data_final.rds')
gen <- Genomics_data_recoded

data <- readRDS('../../extraction_and_recording/outputs/final/biomarker_nNightingale_final.rds')
data <- data.frame(data)

data$tlen <- gen$AdjTSRatio.0.0


#data filtering based on a 20% NA threshold------------------------------

#filtering rows that have more than 20% missingness
#data_fil <- data %>% filter(rowSums(is.na(data))< 0.2*(ncol(data)))

#filtering out those with no telomere length
data_fil <- data %>% filter(!is.na(tlen))

##filtering out columns with more than 20% NA

df_fin <- data_fil %>% select(-contains("LPA"))
df_fin <- df_fin %>% select(-contains("E2"))
df_fin <- df_fin %>% select(-contains("RF"))
df_fin <- df_fin %>% select(-contains("MiAlbUr"))

#df_fin <- df_fin[sample(1:nrow(df_fin), 1000), ]
df_fin <- as.matrix(df_fin)
#Imputation-------------------------
tempData <- knn.impute(
  df_fin,
  k = 10,
  cat.var = NULL,
  to.impute = 1:nrow(df_fin),
  using = 1:nrow(df_fin)
)
#tempData <- mice(df_fin,m=5,maxit=50,meth='pmm',seed=500)

saveRDS(tempData,("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/final/bio_imputed_2.rds"))
library(dplyr)

biomarker_nNightingale_final$tlen <- telomere_length$AdjTSRatio.0.0
biomarker_nNightingale_final_2 <- biomarker_nNightingale_final[!is.na(biomarker_nNightingale_final$tlen),]

#bio_imputed <- data.frame(bio_imputed)

#exposures <- exposures %>% filter(!is.na(tlen))
row.names(bio) <- row.names(exposures_2)
saveRDS(biomarker_nNightingale_final_2, 'bio_not_impute.rds')
