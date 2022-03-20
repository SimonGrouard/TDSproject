
# Libraries ---------------------------------------------------------------
library(focus)
library(pheatmap)
library(glmnet)
library(igraph)


# Function ----------------------------------------------------------------
check_positive_negative <- function(x){
  if (x < 0){
    x = "Negative"
  } else if (x > 0){
    x = "Positive"
  }
}


# Data --------------------------------------------------------------------

telomere_length <- readRDS(here::here("extraction_and_recording/outputs/recoded/Genomics_data_recoded.rds"))
# exposures_imputed <- readRDS(here::here("extraction_and_recording/outputs/final/Exposures_imputed_full.rds"))
exposures_unimputed <- readRDS(here::here("extraction_and_recording/outputs/final/Exposures_unimputed.rds"))
bio_imputed <- readRDS(here::here("extraction_and_recording/outputs/final/bio_imputed_2.rds"))
bio_imputed <- data.frame(bio_imputed)

# Data Cleaning -----------------------------------------------------------
# Identifying rownames that are common
rown <- rownames(bio_imputed)
rowb <- rownames(exposures_unimputed)
exposures <- exposures_unimputed[rownames(exposures_unimputed) %in% rown, ] 
telomere_length <- telomere_length[rownames(telomere_length) %in% rowb, ]
bio_imputed_reduced <- bio_imputed[rownames(bio_imputed) %in% rowb, ]

# Combining dataframes
bio_exposures <- cbind(exposures, bio_imputed_reduced)

# To remove telomere length variable from bio_exposures
bio_exposures_wo_TS <- subset(bio_exposures, select = -c(AdjTSRatio, tlen))

# Analysis ----------------------------------------------------------------
bio_exposures_test <- head(bio_exposures, 1000)

# Setting penalty factor for Age and Sex, where if they are 0 they will always be selected
penalty_factor <-
  names(bio_exposures_wo_TS) == "AgeAssess" | names(bio_exposures_wo_TS) == "Sex_Female"

penalty_factor <-
  as.integer(!penalty_factor)

x = as.matrix(bio_exposures_wo_TS)
y = as.matrix(bio_exposures$AdjTSRatio)
set.seed(1)
lasso_bio_exposures <- cv.glmnet(x = x, y = y, penalty.factor = penalty_factor)

betas <- coef(lasso_bio_exposures, s='lambda.min')[-1]
names(betas) <- rownames(coef(lasso_bio_exposures))[-1]
variables_selected <- betas[betas != 0] # 115 variables selected
names_variables_selected <- names(variables_selected)

# Checking beta direction of variables selected
variables_selected_beta_direction <- 
  t(data.frame(lapply(variables_selected, check_positive_negative)))


# Including exact beta values together with the beta direction
variables_selected_direction_value <-
  cbind(variables_selected, variables_selected_beta_direction)

write.csv(variables_selected_direction_value, here::here("Results/exposures_biomarker_lasso/variables_selected_with_betas_direction_lambdamin.csv"))

# Extracting the names of variables selected
biomarkers_selected <- names_variables_selected[which(names(variables_selected) %in% names(bio_imputed_reduced))]
exposures_selected <- names_variables_selected[which(names(variables_selected) %in% names(exposures_unimputed))]


lasso_bio_exposures_plot_MSE <- plot(lasso_bio_exposures)
