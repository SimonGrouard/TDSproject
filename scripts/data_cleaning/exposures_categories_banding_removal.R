
# Data --------------------------------------------------------------------
exposures <- readRDS(here::here("extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final.rds"))
exposures_checkpoint <- readRDS(here::here("extraction_and_recording/outputs/recoded/checkpoint_rebanding.rds"))

# Functions ---------------------------------------------------------------
## Looking for the strings and replacing them with NAs
strings_to_na <- c("None of the above", "Do not know", "Prefer not to answer")
removing_na <- function(x){
  replace(x, grepl(paste(strings_to_na, collapse = "|"), x), NA)
}

## Looking for strings and replacing them with specified value
replacement_fun <- function(x, y, z){
  # x = List to replace
  # y = Vector of strings to replace
  # z = String to replace with
  x = as.character(x)
  x = replace(x, grepl(paste(y, collapse = "|"), x, ignore.case = TRUE), z)
  x = as.factor(x)
}

# Function to check if the NAs correspond to original
code_check <- function(x, y){
  ifelse(table(is.na(x)) == table(is.na(y)), print("OKAY"), print("NOT OKAY"))
}
#removing variables------------------------------------------

exposures_rm <- exposures[ , ! names(exposures) %in% c("MOB", "BirthWKnow",'LTFUReason', 
                                                       
                                                       'FreqWlkPle4Wk', 
                                                       
                                                       'FreqStrenSport4WK', 
                                                       
                                                       'Nap', 
                                                       
                                                       'Insomnia', 
                                                       
                                                       'DaySleep', 
                                                       
                                                       'BirthCountryUK', 
                                                       
                                                       'SkinCol', 
                                                       
                                                       'Adopted', 
                                                       
                                                       'FreqDepress2Wk', 
                                                       
                                                       'FreqDisint2Wk', 
                                                       
                                                       'FreqTense2Wk', 
                                                       
                                                       'FreqTired2Wk', 
                                                       
                                                       'FreqOthExc4Wk', 
                                                       
                                                       'BirthCountryNonUK', 
                                                       
                                                       'JobWalkStand', 
                                                       
                                                       'JobHvyPhys', 
                                                       
                                                       'JobShift', 
                                                       
                                                       'TobacTypePrevSmok', 
                                                       
                                                       'EverSmokStop6MoPlus', 
                                                       
                                                       'JobNight', 
                                                       
                                                       'GestDiab', 
                                                       
                                                       'PrivateHc', 
                                                       
                                                       'PrevSmokAllDays', 
                                                       
                                                       'GPPhysRestrcHeart', 
                                                       
                                                       'FamilyIllness', 
                                                       
                                                       
                                                       'Bipolar', 
                                                       
                                                       'MajDepressEp', 
                                                       
                                                       'MajDepressModRecr', 
                                                       
                                                       'MajDepressSevRecr', 
                                                       
                                                       'EverAlcDepend', 
                                                       
                                                       'EverAlcAddict', 
                                                       
                                                       'EverAlcConcern', 
                                                       
                                                       'AlcInjure', 
                                                       
                                                       'MajOp', 
                                                       
                                                       'FreqAlc', 
                                                       
                                                       'CurrAlcAddict', 
                                                       
                                                       'Freq6UnitAlc', 
                                                       
                                                       'EverAddictIllicitDrug', 
                                                       
                                                       'EverSuicideAttempt', 
                                                       
                                                       'Portion_size',
                                                       
                                                       'Water2', 
                                                       
                                                       'Coffee', 
                                                       
                                                       'AlcConsume', 
                                                       
                                                       'VitSuppl', 
                                                       
                                                       'VigPhysTime', 
                                                       
                                                       'ModPhysTime', 
                                                       
                                                       'LightPhysTime')]


#converting 'none of the above', 'do not know' and 'prefer not to answer' to NA-----------------------

# exposures_rm[exposures_rm == "None of the above"] <- NA
# exposures_rm[exposures_rm == "Do not know"] <- NA
# exposures_rm[exposures_rm == "Prefer not to answer"] <- NA

# Using apply with removing_na function to replace strings with NAs as well as remove them as factors
exposures_rm_incl_na <- lapply(exposures_rm, removing_na)
exposures_rm_incl_na <- data.frame(exposures_rm_incl_na)
exposures_rm_incl_na <- lapply(exposures_rm_incl_na_test, function(x) {if(is.factor(x)) factor(x) else x})
exposures_rm_incl_na <- data.frame(exposures_rm_incl_na)

# exposures_rm_incl_na_test <- exposures_rm_incl_na
# # attr(exposures_rm_incl_na_test$AccomType, "ATT") <- NULL
# # 
# # 
# # exposures_rm_incl_na_test <- apply(exposures_rm_incl_na_test, 2, function(x) { attr(x, "ATT") <- NULL })
# # data.frame(exposures_rm_incl_na_test)
# # attr(exposures_rm_incl_na, "names") <- NULL
# exposures_rm_incl_na_test <- lapply(exposures_rm_incl_na_test, function(x) {if(is.factor(x)) factor(x) else x})
# exposures_rm_incl_na_test <- data.frame(exposures_rm_incl_na_test)

saveRDS(exposures_rm_incl_na, here::here("extraction_and_recording/outputs/recoded/checkpoint_rebanding.rds"))


# Special Cases ------------------------------------------
## Special case: LightSmoke convert NAs to Yes
# exposures_rm_incl_na_short <- head(exposures_checkpoint, 5000)
exposures_checkpoint$LightSmoke <-
  replace(exposures_checkpoint$LightSmoke, grepl("Smoked on most or all days", exposures_checkpoint$SmokePast), "Yes")

## EverCannabis banding for Yes
cannabis_strings <- c("Yes, 1-2 times", "Yes, 3-10 times")
cannabis_strings_11more <- c("Yes, 11-100 times", "Yes, more than 100 times")
exposures_checkpoint_cannabis <- exposures_checkpoint

### Changing to character first for replacement as can't simply replace factors
exposures_checkpoint_cannabis$EverCannabis <- as.character(exposures_checkpoint_cannabis$EverCannabis)

### Banding Yes, 1-2 times and Yes, 3-10 times as Yes, 1-10 times
exposures_checkpoint_cannabis$EverCannabis <-
  replace(exposures_checkpoint_cannabis$EverCannabis,  grepl(paste(cannabis_strings, collapse = "|"), 
                                                             exposures_checkpoint_cannabis$EverCannabis), "Yes, 1-10 times")

### Banding for Yes, 11 or more times
exposures_checkpoint_cannabis$EverCannabis <-
  replace(exposures_checkpoint_cannabis$EverCannabis,  grepl(paste(cannabis_strings_11more, collapse = "|"), 
                                                             exposures_checkpoint_cannabis$EverCannabis), "Yes, 11 or more times")

### Converting back to factor
exposures_checkpoint_cannabis$EverCannabis <- as.factor(exposures_checkpoint_cannabis$EverCannabis)

### Checkpoint 2
exposures_checkpoint2 <- exposures_checkpoint_cannabis

# Checking to make sure that NAs are the same as original (In case of coding mistake)
table(is.na(exposures_checkpoint$EverCannabis)) == table(is.na(exposures_checkpoint2$EverCannabis))

#banding variables-------------------------------------
## Banding Nutrition
### Banding once a week
strings_once_week <- c("Less than once a week", "Once a week")
columns_apply <- c("OilyFish1", "NonOilyFish", "ProcMeat", "Poultry1", "Beef1", "Lamb1", "Pork1", "Cheese")
exposures_checkpoint2_nutrition <- exposures_checkpoint2
exposures_checkpoint2_nutrition[columns_apply] <- apply(exposures_checkpoint2_nutrition[columns_apply], 2,
                                                        FUN = replacement_fun, y = strings_once_week, z = "Up to once a week")

### Banding 2-6 times a week
strings_6_week <- c("2-4 times a week", "5-6 times a week")
exposures_checkpoint2_nutrition[columns_apply] <- apply(exposures_checkpoint2_nutrition[columns_apply], 2,
                                                        FUN = replacement_fun, y = strings_6_week, z = "Between 2-6 times a week")

exposures_checkpoint3 <- exposures_checkpoint2_nutrition

## SmokePast
smokepast_strings <- c("I have never smoked", "Just tried once or twice")
exposures_checkpoint3$SmokePast <- replacement_fun(exposures_checkpoint3$SmokePast, smokepast_strings, "Never smoked")

## HouseSmokers
housesmokers_strings <- c("Yes, more than one household member smokes", "Yes, one household member smokes")
exposures_checkpoint3$HouseSmokers <- replacement_fun(exposures_checkpoint3$HouseSmokers, housesmokers_strings, 
                                                      "Yes, one or more household member smokes")

## AlcFreq
### Banding daily and 3 or four times a week
alcfreq_strings_daily <- c("Daily or almost daily", "Three or four times a week")
exposures_checkpoint3$AlcFreq <- replacement_fun(exposures_checkpoint3$AlcFreq, alcfreq_strings_daily,
                                                 "Often")
alcfreq_strings_rarely <- c("Special occasions only", "One to three times a month")
exposures_checkpoint3$AlcFreq <- replacement_fun(exposures_checkpoint3$AlcFreq, alcfreq_strings_rarely,
                                                 "Rarely")

## HomeUrbanRural
### Postcode not linkable to NAs
exposures_checkpoint3$HomeUrbanRural <- replacement_fun(exposures_checkpoint3$HomeUrbanRural, "Postcode not linkable", NA)
### Banding categories to sparse and less sparse
#### Sparse
homeurbanrural_sparse <- c("- sparse$", "Accessible Rural$", "Accessible Small Town$", "Remote Small Town$", "Remote Rural$")
exposures_checkpoint3$HomeUrbanRural <- replacement_fun(exposures_checkpoint3$HomeUrbanRural, homeurbanrural_sparse,
                                                 "Sparse")

#### Less sparse
homeurbanrural_less_sparse <- c("less sparse$", "Urban Area$")
exposures_checkpoint3$HomeUrbanRural <- replacement_fun(exposures_checkpoint3$HomeUrbanRural, homeurbanrural_less_sparse,
                                                        "Non-Sparse")

## BipolarDepressStat
exposures_checkpoint4 <- exposures_checkpoint3
### Banding Bipolar disorders together
bipolardepressstat_bipolar <- c("^Bipolar")
exposures_checkpoint4$BipolarDepressStat <- replacement_fun(exposures_checkpoint4$BipolarDepressStat, bipolardepressstat_bipolar,
                                                            "Have Bipolar")
### Banding depression together
bipolardepressstat_depress <- c("major depression")
exposures_checkpoint4$BipolarDepressStat <- replacement_fun(exposures_checkpoint4$BipolarDepressStat, bipolardepressstat_depress,
                                                            "Have Depression")

## AvgAlcDay
### 1 to 4
one_four <- c("1 or 2", "3 or 4")
exposures_checkpoint4$AvgAlcDay <- replacement_fun(exposures_checkpoint4$AvgAlcDay, one_four, "1 to 4")

### 5 to 9
five_nine <- c("5 or 6", "8 or 9")
exposures_checkpoint4$AvgAlcDay <- replacement_fun(exposures_checkpoint4$AvgAlcDay, five_nine, "5 to 9")

## Happy variables
### Very happy
very_happy <- c("Extremely happy", "Very Happy")
columns_happy <- c("GenHappy", "GenHealthHappy")
exposures_checkpoint4[columns_happy] <- apply(exposures_checkpoint4[columns_happy], 2, replacement_fun,
                                             y = very_happy, z = "Very Happy")

### Very Unhappy
very_unhappy <- c("Extremely unhappy", "Very unhappy")
exposures_checkpoint4[columns_happy] <- apply(exposures_checkpoint4[columns_happy], 2, replacement_fun,
                                              y = very_unhappy, z = "Very Unhappy")

## Ethnicity
### Banding white
white <- c("British", "Irish", "Any other white background")
exposures_checkpoint4$Ethnicity <- replacement_fun(exposures_checkpoint4$Ethnicity, white, "White")

### Banding Mixed
mixed <- c("^White and", "Any other mixed background")
exposures_checkpoint4$Ethnicity <- replacement_fun(exposures_checkpoint4$Ethnicity, mixed, "Mixed")

### Banding Asian or Asian British
asian <- c("Indian", "Pakistani", "Bangladeshi", "Any other Asian background")
exposures_checkpoint4$Ethnicity <- replacement_fun(exposures_checkpoint4$Ethnicity, asian, "Asian")

### Banding Black or Black British
black <- c("Caribbean", "African", "Any other Black background")
exposures_checkpoint4$Ethnicity <- replacement_fun(exposures_checkpoint4$Ethnicity, black, "Black")

# Output data
exposures_final <- exposures_checkpoint4
saveRDS(exposures_final, here::here("extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_banded.rds"))
