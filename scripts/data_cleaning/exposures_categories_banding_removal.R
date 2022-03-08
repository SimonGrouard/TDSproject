exposures <- readRDS(here::here("extraction_and_recording/outputs/recoded/Exposures_covariates_recoded_combined_final.rds"))

apply(exposures, 2, table)

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

exposures_rm[exposures_rm == "None of the above"] <- NA
exposures_rm[exposures_rm == "Do not know"] <- NA
exposures_rm[exposures_rm == "Prefer not to answer"] <- NA


#banding variables-------------------------------------



