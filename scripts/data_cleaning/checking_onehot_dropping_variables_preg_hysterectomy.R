exposures_onehot <- readRDS(here::here("extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final_onehot.rds"))
exposures_banded <- readRDS(here::here("extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_banded.rds"))

sym_diff <- function(a,b) setdiff(union(a,b), intersect(a,b))

col_names_onehot <- colnames(exposures_onehot)
col_names_banded <- colnames(exposures_banded)

sym_diff(col_names_onehot, col_names_banded)
setdiff(names(exposures_onehot), names(exposures_banded))

exposures_onehot$Preg_No[exposures_onehot$Preg_Unsure == 1] <- NA
exposures_onehot$Hysterectomy_No[exposures_onehot$Hysterectomy_Not.sure == 1] <- NA

drops <- c("Hysterectomy_Yes", "Hysterectomy_Not.sure", "Preg_Yes", "Preg_Unsure")
exposures_onehot_dropped <- exposures_onehot[ , !(names(exposures_onehot) %in% drops)]
saveRDS(exposures_onehot_dropped, here::here("extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final_onehot.rds"))
