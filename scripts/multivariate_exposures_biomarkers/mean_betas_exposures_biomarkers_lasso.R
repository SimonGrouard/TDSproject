library(dplyr)

mean_betas <- readRDS(here::here("Results/exposures_biomarker_lasso/mean_betas.rds")) # Your mean betas
exposures_bio_selprop <- readRDS(here::here("Results/exposures_biomarker_lasso/exposures_bio_stability_selection_proportions.rds")) # Your selection proportions
lasso_betas <- read.csv(here::here("Results/exposures_biomarker_lasso/variables_selected_with_betas_direction_lambdamin.csv"))

# To extract the betas for those selected >= 90% of the time
sel_prop_90 <- exposures_bio_selprop[exposures_bio_selprop >= 0.9]
mean_betas_sel_prop_90 <- mean_betas[names(mean_betas) %in% names(sel_prop_90)]
mean_betas_sel_prop_90 <- as.data.frame(mean_betas_sel_prop_90)

# To compare lasso betas and selection proportion betas
rownames(lasso_betas) <- lasso_betas[,1]
lasso_betas[,1] <- NULL
lasso_selection_combined <- merge(mean_betas_sel_prop_90, lasso_betas, by = 0)
