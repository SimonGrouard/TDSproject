# Data --------------------------------------------------------------------
telomere_length <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Genomics_data_recoded.rds")
exposures <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final_onehot.rds") 
## Adding telomere length variable to exposure dataframe
data <- exposures %>%
  mutate(AdjTSRatio = telomere_length$AdjTSRatio.0.0,
         ZAdjTSRatio = telomere_length$ZAdjTSRatio.0.0)
test_data <- head(data, 1000)
readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/dataForimpute.rds")

## handle NAs
NA_Exposure <- as.data.frame(colSums(is.na(data)))
NA_Exposure <- as.data.frame(t(NA_Exposure))
rownames(NA_Exposure) <- c("NAnumber")

### variables with n(NA) < 20% 
tNA_Exposure <- as.data.frame(t(NA_Exposure))
df_lessNA_var <- as.data.frame(subset(tNA_Exposure, NAnumber<= 0.2*(nrow(data)))) # only 171 of the 344 variables remained in our dataset (sample size = 150649)
df_lessNA <- data %>% select(colnames = rownames(df_lessNA_var))
colnames(df_lessNA) <- rownames(df_lessNA_var)
lessNA_Exposure <- t(as.data.frame(colSums(is.na(df_lessNA))))
rownames(lessNA_Exposure) <- c("NA_counts")
df <- na.omit(df_lessNA)
saveRDS(df, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/dataForlasso.rds')

### table of n(variables) kept and acceptance of ratio(NA), along with sample size.
tab <- matrix(c(291,271,183,180,174,173,171,153,
                0,0,0,15556,74365,123369,150649,213459), ncol=8, byrow=TRUE)
colnames(tab) <- c('NA < 80%','NA < 70%','NA < 60%','NA < 50%',
                   'NA < 40%','NA < 30%','NA < 20%','NA < 10%')
rownames(tab) <- c('n(variables)','n(participants)')
#### n(variables)
var <- t(data.frame(tab[1,]))
rownames(var) <- c("n(variables)")
#### n(participants)
participants <- t(data.frame(tab[2,]))
rownames(participants) <- c("n(participants)")
#### visualisation
barplot(as.matrix(var), names.arg = names(var), las = 2,
        cex.axis = 1, cex.names = 0.5,
        width = 1, horiz = F, 
        xlab = "proportion of NA", ylab = "n(variables)",
        main = "n(variables) kept after excluding exposures with high NA%")
barplot(as.matrix(participants), names.arg = names(participants), las = 2,
        cex.axis = 0.8, cex.names = 0.5,
        width = 1, horiz = F, 
        xlab = "NA numbers", ylab = "n(participants)",
        main = "n(participants) kept after excluding all rows with NA")

## NA pattern visualisation: failed (Error: cannot allocate vector of size 940.5 Gb)
suppressPackageStartupMessages(library(naniar))
png("NA_Exposure_pattern.png", 30, 21)
vis_miss(data.frame(df_lessNA), cluster= T, sort_miss = T, warn_large_data = F)
dev.off()

saveRDS(df, '/rds/general/project/hda_21-22/live/TDS/Group_6/Results_multivariate_exposures/dataForlasso.rds')

## test with first 1000 rows
test_df <- as.matrix(head(df, 1000))
X = test_df[,1:169]
Y = test_df[,170]