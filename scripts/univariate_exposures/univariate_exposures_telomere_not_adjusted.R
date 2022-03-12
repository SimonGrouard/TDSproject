suppressMessages(library(tidyverse))

## Parameters

args=commandArgs(trailingOnly=TRUE)
data_path=toString(args[1])
nchunks=as.numeric(args[2])
ichunk=as.numeric(args[3])


## Loading packages and data

expo_path <- paste0(data_path,"Exposures_covariates_recoded_combined_final_onehot.rds")
output_path <- paste0(data_path,"Genomics_data_recoded.rds")

expo<-readRDS(expo_path)
tel_length<-readRDS(output_path)["AdjTSRatio.0.0"]

## Running univariate models
get_pvalues = function(X) {
  df <- data.frame(exposure = X, AdjTSRatio = tel_length$AdjTSRatio.0.0) %>% drop_na()
  
  ## To remove the outliers
  Q <- quantile(df$AdjTSRatio, probs=c(.25, .75), na.rm = T)
  iqr <- IQR(df$AdjTSRatio, na.rm = T)
  df <- df %>% filter(AdjTSRatio > (Q[1] - 1.5*iqr) & 
                      AdjTSRatio < (Q[2] + 1.5*iqr))
  
  if (dim(df)[1] < 1000 | length(table(df[['exposure']])) == 1){ 
    return(1) # too many NAs or only one value in the whole exposure column
  }
  
  model0 <- lm(AdjTSRatio ~ 1, data = df)
  model1 <- lm(AdjTSRatio ~ exposure, data = df)
  pval <- anova(model0, model1)$`Pr(>F)`[2]
  return(pval)
}

ids=as.character(cut(1:ncol(expo), breaks = nchunks, labels = 1:nchunks))
t0=Sys.time()
pvalues = apply(expo[,ids==ichunk], 2, FUN = get_pvalues)
t1=Sys.time()
print(t1-t0)

ifelse(dir.exists("../../Results_univariate_exposures/hotencoded_notadjusted"),"",dir.create("../../Results_univariate_exposures/hotencoded_notadjusted"))
saveRDS(pvalues, paste0("../../Results_univariate_exposures/hotencoded_notadjusted/univ_exposures_no_outliers_", ichunk, ".rds"))

