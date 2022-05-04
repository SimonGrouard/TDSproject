
# Libraries ---------------------------------------------------------------
library(focus)
library(pheatmap)
library(glmnet)
library(igraph)

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
# bio_exposures_test <- head(bio_exposures_wo_TS, 1000)
# bio_exposures_test_TS <- head(bio_exposures, 1000)

# Setting penalty factor for Age and Sex, where if they are 0 they will always be selected
penalty_factor <-
  names(bio_exposures_wo_TS) == "AgeAssess" | names(bio_exposures_wo_TS) == "Sex_Female"

penalty_factor <-
  as.integer(!penalty_factor)

x = as.matrix(bio_exposures_wo_TS)
y = as.matrix(bio_exposures$AdjTSRatio)

# Stability selection -------------------------

stability_selection_lasso <- 
  VariableSelection(xdata = x, ydata = y, verbose=FALSE, 
                    penalty.factor= penalty_factor, 
                    family="gaussian")

pdf(here::here('Results/exposures_biomarker_lasso/exposures_bio_stability_selection_calibration.pdf'))
CalibrationPlot(stability_selection_lasso)
dev.off()

# Calibrated selection proportions 
selprop <- SelectionProportions(stability_selection_lasso)

saveRDS(selprop, here::here('Results/exposures_biomarker_lasso/exposures_bio_stability_selection_proportions.rds'))

# Calibrated parameters
hat_params <- Argmax(stability_selection_lasso)
print(hat_params)
saveRDS(hat_params, here::here('Results/exposures_biomarker_lasso/exposures_bio_stability_hat_params.rds'))

# Visualisation of selection proportions

pdf(here::here('Results/exposures_biomarker_lasso/exposures_bio_stability_selection_proportions.pdf'))
par(mar = c(10, 5, 5, 1))

plot(selprop, type = "h", lwd = 3, las = 1, 
     xlab = "", ylab = "Selection Proportions", xaxt = "n", 
     col = ifelse(selprop >= hat_params[2], yes = "blue", no = "grey"),
     cex.lab = 1.5)

abline(h = hat_params[2], lty = 2, col = "darkred")

for (i in 1:length(selprop)){
  axis(side = 1, at = i, labels = names(selprop)[i], las = 2, 
       cex.axis = 0.25,
       col = ifelse(selprop >= hat_params[2], yes = "blue", no = "grey"), 
       col.axis = ifelse(selprop[i] >=hat_params[2], yes = "blue", no = "grey"))
}
dev.off()




