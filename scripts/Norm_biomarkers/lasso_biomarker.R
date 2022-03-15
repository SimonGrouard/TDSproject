# Libraries ---------------------------------------------------------------
library(dplyr)
library(glmnet)
library(igraph)
library(focus)
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
         AdjTSRatio = telomere_length$AdjTSRatio.0.0,
         Sex = exposures$Sex)
#test_bio_telomere <- head(bio_telomere, 1000)

bio_telomere <- bio_telomere %>% select(-contains("tlen"))


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
x = as.matrix(bio_telomere[,1:47])
y= bio_telomere$AdjTSRatio
set.seed(1)
model.lasso <- cv.glmnet(x = x, y = y, alpha = 1)

betas=coef(model.lasso, s='lambda.1se')[-1]
names(betas)=rownames(coef(model.lasso))[-1]

pdf('lasso_plot_imp_2.pdf')
plot(model.lasso)
dev.off()

best_lam <- model.lasso$lambda.1se
best_model <- glmnet(x = x, 
                     y = y, lambda = best_lam)

table(coef(model.lasso, s='lambda.1se')[-1]!=0)
betas_2=coef(best_model, s='lambda.1se')[-1]
names(betas_2)=rownames(coef(best_model, s='lambda.1se'))[-1]

pdf('lasso_plot_best_imp_2.pdf')
plot(best_model)
dev.off

saveRDS(betas, '/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/lasso_betas_imp_2.rds')
saveRDS(betas_2, '/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/lasso_best_betas_imp_2.rds')

#Stability analysis-----------------------
#LassoSub = function(k=1, Xdata, Ydata, family="gaussian", penalty.factor=NULL) {
 # if (is.null(penalty.factor)){
#    penalty.factor=rep(1,ncol(Xdata))
 # }
  #set.seed(k)
#  s = sample(nrow(Xdata), size = 0.8 * nrow(Xdata))
 # Xsub = Xdata[s, ]
  #Ysub = Ydata[s]
  #model.sub = cv.glmnet(x = Xsub, y = Ysub, alpha = 1, family=family, penalty.factor=penalty.factor)
  #coef.sub = coef(model.sub, s = "lambda.1se")[-1]
  #return(coef.sub)
#}

#niter=100
#lasso.stab=sapply(1:niter, FUN=LassoSub, Xdata=x, Ydata=y)

#apply(lasso.stab, 2, FUN=function(x){sum(x!=0)})

#lasso.prop=apply(lasso.stab, 1, FUN=function(x){sum(x!=0)/length(x)})
#names(lasso.prop)=colnames(x)

#lasso.prop=sort(lasso.prop, decreasing = TRUE)

#pdf('/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/plot_prop_imp.pdf')
#plot(lasso.prop[lasso.prop>0.2], type = 'h', col='navy', lwd=3, xaxt='n', 
 #    xlab='', ylab=expression(beta), ylim=c(0,1.2), las=1)
#text(lasso.prop[lasso.prop>0.2]+0.07, labels = names(lasso.prop[lasso.prop>0.2]), 
#     pos=3, srt=90, cex=0.7)
#dev.off()

#saveRDS(lasso.prop, '/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/Norm_biomarkers/lasso_prop_imp.rds')

#stability selection 2-------------------------







