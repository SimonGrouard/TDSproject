# Libraries ---------------------------------------------------------------
library(dplyr)
library(glmnet)
library(igraph)
library(devtools)
library(focus)
library(pheatmap)

# Data --------------------------------------------------------------------
#rm(list=ls())
telomere_length <- readRDS(("../../extraction_and_recording/outputs/recoded/Genomics_data_recoded.rds"))
exposures <- readRDS(("../../extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final.rds"))
bio <- readRDS(("../../extraction_and_recording/outputs/final/bio_imputed.rds"))
rown <- rownames(bio)
exposures <- exposures[rownames(exposures) %in% rown, ] 
telomere_length <- telomere_length[rownames(telomere_length) %in% rown, ] 
#test_telomere <- head(telomere_length, 100)
#test_biomarker <- head(bio, 100)
exposures$Sex <- ifelse(exposures$Sex == 'Female', 0, 1)


# Data Cleaning -----------------------------------------------------------
## Adding telomere length variable to biomarker dataframe
bio_telomere <- bio %>%
  mutate(Age = exposures$AgeAssess,
         Sex = exposures$Sex,
         AdjTSRatio = telomere_length$AdjTSRatio.0.0
         )


bio_telomere <- bio_telomere %>% select(-contains("tlen"))

#test_bio_telomere <- head(bio_telomere, 1000)
## Removal of NAs

#filtering rows that have more than 20% missingness
#bio_fil <- test_bio_telomere %>% filter(rowSums(is.na(test_bio_telomere))<= 0.2*(ncol(test_bio_telomere)))
#bio_fil <- bio_telomere %>% filter(rowSums(is.na(bio_telomere))<= 0.2*(ncol(bio_telomere)))


#filtering out those with no telomere length
#bio_fil <- bio_fil %>% filter(!is.na(AdjTSRatio))

##filtering out columns with more than 20% NA

#bio_fin <- bio_fil %>% select(-contains("LPA"))
#bio_fin <- bio_fin %>% select(-contains("E2"))
#bio_fin <- bio_fin %>% select(-contains("RF"))
#bio_fin <- bio_fin %>% select(-contains("MiAlbUr"))
#bio_fin = na.omit(bio_fin)

# Analysis ----------------------------------------------------------------
## LASSO Models
x = as.matrix(bio_telomere[,1:48])
y= bio_telomere$AdjTSRatio
set.seed(1)
model.lasso <- cv.glmnet(x = x, y = y, alpha = 1)

betas=coef(model.lasso, s='lambda.1se')[-1]
names(betas)=rownames(coef(model.lasso))[-1]

#pdf('lasso_plot_imp_2.pdf')
#plot(model.lasso)
#dev.off()

best_lam <- model.lasso$lambda.1se
best_model <- glmnet(x = x, 
                     y = y, lambda = best_lam)

table(coef(model.lasso, s='lambda.1se')[-1]!=0)
betas_2=coef(best_model, s='lambda.1se')[-1]
names(betas_2)=rownames(coef(best_model, s='lambda.1se'))[-1]

#pdf('lasso_plot_best_imp_2.pdf')
#plot(best_model)
#dev.off

#saveRDS(betas, '/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/lasso_betas_imp_2.rds')
#saveRDS(betas_2, '/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/lasso_best_betas_imp_2.rds')


#stability selection 2-------------------------

out=VariableSelection(xdata=x, ydata=y, verbose=FALSE, 
                      penalty.factor=c(rep(1,ncol(x))), 
                      family="gaussian")

pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/stab_selection_bio_2.pdf')
CalibrationPlot(out)
dev.off()

# Calibrated selection proportions 
selprop=SelectionProportions(out)
selprop <- data.frame(selprop)

saveRDS(selprop, '/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/selection_proportions_bio_2.rds')

# Calibrated parameters
hat_params=Argmax(out)
print(hat_params)
saveRDS(hat_params, '/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/hat_params_bio_2.rds')

# Visualisation of selection proportions
pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/vis_sel_prop_2.pdf')
par(mar=c(10,5,1,1))
plot(selprop, type="h", lwd=3, las=1, xlab="", ylab="Selection Proportion", xaxt="n",
     col=ifelse(selprop>=hat_params[2], yes="red", no="grey"), cex.lab=1.5)
abline(h=hat_params[2], lty=2, col="darkred")
for (i in 1:length(selprop)){
  axis(side=1, at=i, labels=names(selprop)[i], las=2, 
       col=ifelse(selprop[i]>=hat_params[2], yes="red", no="grey"),
       col.axis=ifelse(selprop[i]>=hat_params[2], yes="red", no="grey"))
}
dev.off()

