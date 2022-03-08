Genomics_data_recoded <- readRDS('../extraction_and_recording/outputs/final/Genomics_data_final.rds')
gen <- Genomics_data_recoded


Q <- quantile(gen$AdjTSRatio.0.0, probs=c(.25, .75))

iqr <- IQR(gen$AdjTSRatio.0.0)

gen2 <- gen %>% filter(AdjTSRatio.0.0 > (Q[1] - 1.5*iqr) & 
                         AdjTSRatio.0.0 < (Q[2] + 1.5*iqr))  
hist(gen2$AdjTSRatio.0.0)