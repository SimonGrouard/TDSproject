Genomics_data_recoded <- readRDS('../extraction_and_recording/outputs/final/Genomics_data_final.rds')
gen <- Genomics_data_recoded

  #distribution of Telomere length variable
  summary(df_new$TSR)
  ##appears there are some high outliers as third quartile is 0.9077 but max is 5.36
  (table(df_new$TSR>1))
  #1 person with TSR bigger than 5, 1 person between 4-5, 3 people between 3-4, 7 people between 2-3, over 1000 between 1-2

  ggplot(df_new, aes(x = TSR )) +
    geom_density()

  fil <- df_new %>% filter(TSR >0.9077)
  plot(density(fil$TSR))

  fil2 <- df_new %>% filter(TSR > 2)
  plot(density(fil2$TSR))


library(tidyverse)

#filtering NAs
gen <- gen %>% filter(!is.na(gen$AdjTSRatio.0.0))

#filtering TLs more than 1.5*IQR

Q <- quantile(gen$AdjTSRatio.0.0, probs=c(.25, .75))

iqr <- IQR(gen$AdjTSRatio.0.0)

gen2 <- gen %>% filter(AdjTSRatio.0.0 > (Q[1] - 1.5*iqr) & 
                       AdjTSRatio.0.0 < (Q[2] + 1.5*iqr))  
hist(gen2$AdjTSRatio.0.0)

#filtering those with genetic sex as NA
gen2 <- gen2 %>% filter(!is.na(gen2$GeneticSex.0.0))

saveRDS(gen2, file = "../extraction_and_recording/outputs/final/Genomics_data_filtered.rds")

