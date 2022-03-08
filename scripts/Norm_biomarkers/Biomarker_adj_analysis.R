###SCRIPT TO RUN UNIVARIATE ANALYSIS ADJUSTED FOR AGE AND SEX


## Loading packages and data

suppressMessages(library(tidyverse))

data <- readRDS('/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/final/biomarker_nNightingale_final.rds')
data <- data.frame(data)

Genomics_data_recoded <- readRDS('/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/final/Genomics_data_final.rds')
gen <- Genomics_data_recoded

Exposures_covariates_recoded_combined_final <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final.rds")

data$tlen <- gen$AdjTSRatio.0.0
data$sex <- gen$GeneticSex.0.0
data$age <- Exposures_covariates_recoded_combined_final$AgeAssess

#running only on count data
data <- data %>% select(-contains("Count"))

#--------------------------------------------------------------
#univariate adjusted analysis

pval <- rep(0,ncol(data)) 

for (j in 1:ncol(data)){
  model <- lm(data$tlen~. + data$age + data$sex ,data=data[,c(j,ncol(data))], na.action = na.omit) ##for some reason the only way I could get it to 
  coefs <- summary(model)$coefficients
  pval[j] <- coefs[2,'Pr(>|t|)']
}

#saving results
final <- cbind(colnames(data), pval)
ifelse(dir.exists("Results_biomarker_univ_3"),"",dir.create("Results_biomarker_univ_3"))
saveRDS(final,("Results_biomarker_univ_3/univ_pval_adj_biomarker_perc.rds"))