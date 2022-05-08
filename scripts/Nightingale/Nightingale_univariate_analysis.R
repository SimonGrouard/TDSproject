suppressMessages(library(tidyverse))

## Parameters

args=commandArgs(trailingOnly=TRUE)
data_path=toString(args[1])
nchunks=as.numeric(args[2])
ichunk=as.numeric(args[3])


## Loading packages and data

nightingale <- readRDS('../extraction_and_recording/outputs/final/Nightingale_metabolomics_final.rds')
Genomics_data_recoded <- readRDS('../extraction_and_recording/outputs/final/Genomics_data_final.rds')

genomics <- data.frame(Genomics_data_recoded)

nightingale$TSR <- genomics$AdjTSRatio.0.0


##filtering out all those who did not take part in the metabolomics study
df_fil <- nightingale %>%
  filter_at(vars(TotChol.0.0), all_vars(!is.na(.)))  

##filtering out any columns relating to technical covariates
df_new <- df_fil %>% select(-contains("QCF"))
df_new <- df_new %>% select(-contains("Spectr"))
df_new <- df_new %>% select(-contains("MQF"))
df_new <- df_new %>% select(-contains("Ship"))
df_new <- df_new %>% select(-contains("HLactate"))
df_new <- df_new %>% select(-contains("HPyruvate"))
df_new <- df_new %>% select(-contains("LGlucose"))
df_new <- df_new %>% select(-contains("LProtein"))
#removing further NA values
df_new <- na.omit(df_new)

##running univariate analysis 
pval <- rep(0,ncol(df_new)) 
for (j in 1:ncol(df_new)){
  model <- lm(TSR~.,data=df_new[,c(j,ncol(df_new))])
  coefs <- summary(model)$coefficients
  pval[j] <- coefs[2,'Pr(>|t|)']
}

pvalues <- cbind(colnames(df_new), pval)
ifelse(dir.exists("Results_nightingale_univ"),"",dir.create("Results_nightingale_univ"))
saveRDS(pvalues, paste0("Results_nightingale_univ/univ_pval_nightingale_", ichunk, ".rds"))
