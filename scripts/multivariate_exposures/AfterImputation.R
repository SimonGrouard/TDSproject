# get the results from imputation -------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(tidyverse))

results_path <- "/rds/general/project/hda_21-22/live/TDS/Group_6/Results/Results_multivariate_exposures/ImputeResPath/"

## results from the analysis
get_results <- function(path, results){
  results <- c(results, readRDS(paste0(results_path,path)))
  return(results)
}

results <- readRDS(paste0(results_path, list.files(results_path)[1]))
for (i in 2:length(list.files(results_path))){ 
  results <- get_results(list.files(results_path)[i], results)
}

length(results)

# data cleaning after imputation ------------------------------------------
imputedNum <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/.rds")

data <- cbind(imputedNum, cat)

NA_Exposure <- as.data.frame(colSums(is.na(data)))
NA_Exposure <- as.data.frame(t(NA_Exposure))
rownames(NA_Exposure) <- c("NAnumber")

tNA_Exposure <- as.data.frame(t(NA_Exposure))
df_lessNA_var <- as.data.frame(subset(tNA_Exposure, NAnumber<= 0.2*(nrow(data)))) # only 171 of the 344 variables remained in our dataset (sample size = 150649)
df_lessNA <- data %>% select(colnames = rownames(df_lessNA_var))
colnames(df_lessNA) <- rownames(df_lessNA_var)
lessNA_Exposure <- t(as.data.frame(colSums(is.na(df_lessNA))))
rownames(lessNA_Exposure) <- c("NA_counts")
df <- na.omit(df_lessNA)
