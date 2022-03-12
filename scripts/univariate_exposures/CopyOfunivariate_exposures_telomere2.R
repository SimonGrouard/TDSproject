suppressMessages(library(tidyverse))

## Loading packages and data
data_path <- '/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/'
expo_path <- paste0(data_path,"Exposures_covariates_recoded_combined_final_onehot.rds")
output_path <- paste0(data_path,"Genomics_data_recoded.rds")

expo<-readRDS(expo_path)
tel_length<-readRDS(output_path)["AdjTSRatio.0.0"]

####
get_factors <- function(name){
  if (is.factor(exposure[[name]])){
    return(name)
  }
}
exp_factors <- lapply(colnames(exposure), get_factors)
exp_factors <- exposure[,unlist(exp_factors)]
high_cat <- function(name){
  if (length(table(exp_factors[[name]]))>10){
    return(name)
  }
}
high_factors <- lapply(colnames(exp_factors), high_cat)
high_factors <- exp_factors[,unlist(high_factors)]
####

tel_length <- data.frame(age = expo$AgeAssess, sex = expo$Sex_Female, AdjTSRatio = tel_length$AdjTSRatio.0.0)
expo <- expo %>% select(-one_of('AgeAssess', 'Sex_Female')) # we remove AgeAssess and Sex_female so they are used as variate 
expo_trunc <- expo %>% select(one_of(colnames(expo)[1:15]))
colnames(expo_trunc)

## Running univariate models
get_pvalues = function(X) {
  df <- data.frame(exposure = X, tel_length) %>% drop_na()
  print(length(table(df[['exposure']])))
  ## To remove the outliers
  Q <- quantile(df$AdjTSRatio, probs=c(.25, .75), na.rm = T)
  iqr <- IQR(df$AdjTSRatio, na.rm = T)
  df <- df %>% filter(AdjTSRatio > (Q[1] - 1.5*iqr) & 
                        AdjTSRatio < (Q[2] + 1.5*iqr))
  
  if (dim(df)[1] < 1000 | length(table(df[['exposure']])) == 1){ 
    return(1) # too many NAs
  }

  model0 <- lm(AdjTSRatio ~ age + sex, data = df)
  model1 <- lm(AdjTSRatio ~ exposure + age + sex, data = df)
  pval <- anova(model0, model1)$`Pr(>F)`[2]
  #return(pval)
  return(list(pvalue = pval, beta_coef = model1$coefficients[[2]]))
}

#t0=Sys.time()
pvalues = apply(expo_trunc, 2, FUN = get_pvalues)
p_value <- NULL
for (i in 1:length(colnames(expo_trunc))){
  print(colnames(expo_trunc)[i])
  pval <- get_pvalues(expo_trunc[[colnames(expo_trunc)[i]]])
  p_value <- c(p_value, pval)
}
#t1=Sys.time()
#print(t1-t0)
names(p_value) <- colnames(expo_trunc)
p_value2 <- p_value[1:as.integer(length(p_value)/2)]
p_value3 <- p_value[(as.integer(length(p_value)/2)+1):length(p_value)]

ifelse(dir.exists("../Results_univariate_exposures"),"",dir.create("../Results_univariate_exposures"))
saveRDS(pvalues, paste0("/rds/general/project/hda_21-22/live/TDS/Group_6/Results_univariate_exposures/hotencoded_adjusted_volcano/univ_exposures_no_outliers_1.rds"))
#saveRDS(p_value2, paste0("../Results_univariate_exposures/univ_exposures_no_outliers_2.rds"))
#saveRDS(p_value3, paste0("../Results_univariate_exposures/univ_exposures_no_outliers_3.rds"))
