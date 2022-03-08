##SCRIPT TO RUN UNIVARIATE ANALYSIS- HAVE RESULTS FROM THIS DON'T NEED TO RUN AGAIN

suppressMessages(library(tidyverse))

## Parameters

args=commandArgs(trailingOnly=TRUE)
data_path=toString(args[1])
nchunks=as.numeric(args[2])
ichunk=as.numeric(args[3])


## Loading packages and data
data <- readRDS('../extraction_and_recording/outputs/final/biomarker_nNightingale_final.rds')
data <- data.frame(data)

Genomics_data_recoded <- readRDS('../extraction_and_recording/outputs/final/Genomics_data_final.rds')
gen <- Genomics_data_recoded

data$tlen <- gen$AdjTSRatio.0.0

#--------------------------------------------------------------
#univariate analysis and plotting
pval <- rep(0,ncol(data)) 
for (j in 1:ncol(data)){
  model <- lm(tlen~.,data=data[,c(j,ncol(data))], na.action = na.omit)
  coefs <- summary(model)$coefficients
  pval[j] <- coefs[2,'Pr(>|t|)']
}

pvalues <- cbind(colnames(data), pval)
ifelse(dir.exists("Results_biomarker_univ_2"),"",dir.create("Results_biomarker_univ_2"))
saveRDS(pvalues, paste0("Results_biomarker_univ_2/univ_pval_biomarker_", ichunk, ".rds"))
