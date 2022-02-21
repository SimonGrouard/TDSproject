args=commandArgs(trailingOnly=TRUE)
nchunks=as.numeric(args[1])
ichunk=as.numeric(args[2])

rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

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
df_new <- na.omit(df_new)

##running univariate analysis 
get_pvalues = function(X) {
  model0 <- lm(X ~ 1, data = df_new)
  model1 <- lm(X ~ exposure, data = df_new)
  pval <- anova(model0, model1)$`Pr(>F)`[2]
  return(pval)
}

ids=as.character(cut(1:ncol(df_new), breaks = nchunks, labels = 1:nchunks))

t0=Sys.time()
pvalues = apply(df_new[,ids==ichunk], 2, FUN = get_pvalues)
t1=Sys.time()
print(t1-t0)

ifelse(dir.exists("Results_nightingale_univ"),"",dir.create("Results_nightingale_univ"))
saveRDS(pvalues, paste0("Results_nightingale_univ/univ_pval_nightingale_", ichunk, ".rds"))
