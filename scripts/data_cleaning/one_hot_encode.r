suppressMessages(library(tidyverse))
suppressMessages(library(mltools))
suppressMessages(library(data.table))

expo_path <- "/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/"
final_path <- paste(expo_path,"Exposures_covariates_recoded_combined_final.rds",sep="")
expo<-readRDS(final_path)
print(dim(expo))

# take the subset of columns which are factors
is_factor <- sapply(expo, is.factor)
expo_factors <- expo[,which(is_factor)]
print(dim(expo_factors))

# divide in two parts: equal to two levels, more than two levels
multilevel_factors <- sapply(expo_factors, function(x) {length(levels(x))>2})
expo_factors_superior2 <- expo_factors[,which(multilevel_factors)] # more than two levels
expo_factors_2 <- expo_factors[,which(!multilevel_factors)] # equal to two levels
print(dim(expo_factors_superior2))
print(dim(expo_factors_2))

# deal with factors which have two levels
change_colname <- function(x){
  if (levels(expo_factors_2[[x]])[1] != '0'){
    x <- paste0(x,'_',levels(expo_factors_2[[x]])[1])
  }
  return(x)
}
colnames(expo_factors_2) <- sapply(colnames(expo_factors_2), change_colname)

change_column <- function(x){
  x <- as.factor(x)
  if (levels(x)[1] != '0'){
    levels(x) <- c(1,0)
  }
  return(as.integer(levels(x))[x])
}
expo_factors_2 <- data.frame(apply(expo_factors_2, 2, change_column))

# deal with factors which have more than two levels
expo_factors_superior2 <- one_hot(as.data.table(expo_factors_superior2))

# put them back in expo
expo <- data.frame(expo, expo_factors_2, expo_factors_superior2)
expo <- expo %>% select(!colnames(expo_factors))
print(dim(expo))

saveRDS(expo, paste0(expo_path, "Exposures_covariates_recoded_combined_final_onehot.rds"))
