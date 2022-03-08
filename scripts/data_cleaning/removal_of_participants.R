participants_drop_out <- read.csv(here::here("w19266_20220222.csv"))
exposures_covariates <- readRDS(here::here("/extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined.rds"))
genomics_covariates <- readRDS(here::here("extraction_and_recording/outputs/recoded/Genomics_data_recoded.rds"))
non_nightingale <- readRDS(here::here("extraction_and_recording/outputs/recoded/biomarker_nNightingale_recoded.rds"))
nightingale <- readRDS(here::here("extraction_and_recording/outputs/recoded/Nightingale_metabolomics_recoded.rds"))
technical <- readRDS(here::here("extraction_and_recording/outputs/recoded/Technical_covariates_recoded.rds"))

removing_participants <- function(x) {
  subset(x, !rownames(x) %in% participants_drop_out$X1031783)
}

# list <- c(1000016, 1000021)
# list <- data.frame(list)
# list <- as.integer(list$list)
# list <- as.integer(participants_drop_out$X1031783)
table(rownames(technical) %in% participants_drop_out$X1031783)
  # To check if there are really no common values
technical_removed <- removing_participants(technical)
nightingale_removed <- removing_participants(nightingale)
