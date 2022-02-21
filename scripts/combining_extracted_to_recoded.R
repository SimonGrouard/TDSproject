
exposures_recoded <- readRDS(here::here('extraction_and_recording/outputs/recoded/Exposures_covariates_recoded.rds'))
exposures_extracted <- readRDS(here::here('extraction_and_recording/outputs/Exposures_covariates_extracted.rds'))
test_recoded <- exposures_recoded
test_extracted <- exposures_extracted
variables <- c('CancerFirst', 'NonCancerFirst', 'OpFirst', 'FTEduAge', 'NDay10minWlk', 'NDayModPhys10min',
               'NDayVigPhys10min', 'SleepDur', 'HomeSmokeExp', 'OutHomeSmokeExp',
               'MumAge', 'NChildFath', 'NLiveBirth', 'FirstChldBW', 'FirstLiveBirthAge', 'LastLiveBirthAge',
               'FormSmokAgeStart', 'NCigPrevSmokeDay', 'AgeSmokEnd', 'NUnsuccSmokStop', 'FatherAge', 'AgeStartSmokCurr',
               'NCigCurrSmok', 'YrImmigratUK', 'StillBirths', 'SponMiscarr', 'PregTermination', 'DepressLongest', 'DepressEpi',
               'CurrentSmokeAmount', 'CookedVeg', 'RawVeg', 'FreshFruit', 'DriedFruit1', 'Tea1', 'Water1')
all_variable_names <- names(test_recoded[, grepl(paste(variables, collapse='|'), # Extracting all variable names that start with the names seen in "variables"
                                                 colnames(test_recoded))])

variables_to_combine_recoded <- test_recoded[, grepl(paste(variables, collapse='|'), # Extracting columns containing variables in "all_variable_names"
                                                     colnames(test_recoded))]
variables_to_combine_extracted <- test_extracted[, grepl(paste(variables, collapse='|'),
                                                         colnames(test_extracted))]

variables_to_combine_recoded[is.na(variables_to_combine_recoded)] <- variables_to_combine_extracted[is.na(variables_to_combine_recoded)] 
  # Replacing only NA values in recoded with the values in extracted

test_recoded[, all_variable_names] <- variables_to_combine_recoded[, all_variable_names]

saveRDS(test_recoded, "../outputs/recoded/Exposures_covariates_recoded_combined.rds")


# combined_recoded <- test_recoded[, !grepl(paste(variables, collapse='|'),
#                                          colnames(test_recoded))]


# combined_recoded <- cbind(combined_recoded, variables_to_combine_recoded)



# merge(variables_to_combine_extracted, variables_to_combine_recoded, by = all_variable_names)

# ifelse(is.na(variables_to_combine_recoded), variables_to_combine_extracted, variables_to_combine_recoded)

