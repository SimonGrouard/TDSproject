
# Libraries ---------------------------------------------------------------



# Data --------------------------------------------------------------------

telomere_length <- readRDS(here::here("extraction_and_recording/outputs/recoded/Genomics_data_recoded.rds"))
exposures <- readRDS(("../../extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final.rds"))
bio <- readRDS(here::here("extraction_and_recording/outputs/final/bio_imputed.rds"))
