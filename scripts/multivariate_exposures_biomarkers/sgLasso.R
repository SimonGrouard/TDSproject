# 1- Packages and imports ------------------------------------------------------
#suppressPackageStartupMessages(library(sgPLS))
#suppressPackageStartupMessages(library(sglOptim))
suppressPackageStartupMessages(library(gglasso))
suppressPackageStartupMessages(library(tidyverse))

#args=commandArgs(trailingOnly=TRUE)
#data_path=toString(args[1])

#biomarkers <- readRDS(paste0(data_path, "final/bio_imputed_2.rds"))
#exposures <- readRDS(paste0(data_path,"final/Exposures_unimputed.rds"))
biomarkers <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/final/bio_imputed_2.rds")
exposures <- readRDS("/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/final/Exposures_unimputed.rds")

# 2- Groups of variables  ------------------------------------------------------

### exposure groups

#Cigarette impacts 
cigarette<-str_split("SmokeCurr SmokePast HouseSmokers LightSmoke TobacTypePrevSmok EverSmokStop6MoPlus SmokStatus PrevSmokAllDays EverSmok TobSmok PackYrsSmok AvgNPackYrLife HomeSmokeExp OutHomeSmokeExp FormSmokAgeStart NCigPrevSmokeDay AgeSmokEnd NUnsuccSmokStop AgeStopSmok CurrentSmokeAmount", " ")[[1]]

#Alcohol impacts 
alcohol<-str_split("AlcFreq FormAlcDrink AlcStatus AvgAlcDay EverAlcDepend EverAlcConcern EverAlcAddict AlcInjure FreqAlc CurrAlcAddict Freq6UnitAlc AlcConsume ALC", " ")[[1]]

#Drug impacts (could be combined with alcohol) 

drug<-str_split("EverAddict EverCannabis EverAddictIllictDrug", " ")[[1]]

#Socioeconomic status 

socioeconomic<-str_split("AccomType JobShift AvgHHIBefTax OwnRentAccom JobNight PrivateHc PayRentMorgage DepIndexRecr IdxMultDeprvE IncomeScrE EmployScrE EduScrE HousScrE CrimScrE ServW HousScrW LivEnviroScrW ComSafeScrW IdxMultDeprvW IdxMultDeprvS IncomeScrS EduScrS HousScrS ServS CrimScrS FTEduAge NJH Qualifications FinancialDifficulties MotorVehicle PublicTransport", " ")[[1]]

#Violence based stress 

violence<-str_split("EverWar PhysAbusePartner SexAbusePartner PhysViolentCrime WitnessViolentDeath VictimSexAssult", " ")[[1]]

#Health stress experienced (for self or family) 

stress<-str_split("Nap Insomnia DaySleep MoodSwing Miserableness Irritibility Sensitivity FedUp Nervous Anxious Tense WorryAftEmbarass Nerves Loneliness Guilty RiskTaking FreqDepress2Wk FreqDisint2Wk FreqTense2Wk FreqTired2Wk SeeGPMH SeePyschMH OverallHealth LongIllDis MajOp DiabetesDiag CancerDiag BoneBreak5Yr SeriousCondDiag EverEndPregDuring EverOralContr Hysterectomy GestDiab DepressWW EUDWW EMH2Days EHIA2Days ShortBreathLG GPPhysRestrcHeart Bipolar MajDepressEp MajDepressModRecr MajDepressSevRecr BipolarDepressStat EverAnxiousMonth EverWorryMoreAvg EverProlongedNoInterest EverProlongedSad GenHappy GenHealthHappy EverSelfHarm EverSuicideAttempt EverMania EverXtremeIrritable LifeThreatenAccident CriticalIllness HealthScrE HealthScrS PregTermination DepressLongest DepressEpi NeuroticismScore FamilyIllness SelfSeriousIllness FamilySeriousIllness FamilyDeath MaritalSeparation", " ")[[1]]

#Level of physical activity 

physactivity<-str_split("JobHvyPhys JobWalkStand FreqWlkPle4Wk FreqStrenSport4WK FreqOthExc4Wk VigPhysTime ModPhysTime LightPhysTime NDay10minWlk NDayModPhys10min NDayVigPhys10min SleepDur Walk Cycle", " ")[[1]]

#Gender, Ethnicity and other demographic factors 

demographic<-str_split("Sex BirthWKnow LTFUReason BirthCountryUK SkinCol Adopted Preg BirthCountryNonUK Ethnicity WeightManual BMI Weight DeathAge YOB NChildFath NLiveBirth FirstChldBW FirstLiveBirthAge LastLiveBirthAge AgeAssess AgeRecr MOB", " ")[[1]]

#Environment factors (e.g., pollution) 

environmental<-str_split("NoisyWork LoudMusEF HomeUrbanRural CMR NO2 NO PM.0 PM2.5 PM2.5absorb PM2.5_.0 ID_NR ID_NMR SRLMR1.0 NDAP2005 NDAP2006 NDAP2007 PMAP2007 ADSL_NP AESL_NP ANSL_NL GP.0.0 GP3.0 NEP.0.0 NEP3.0 DCoast TI_NR TI_NMR TTL_MR GasSolidFuelCookHeat WPchemical WPsmok", " ")[[1]]

#Diet 

diet<-str_split("OilyFish1 NonOilyFish ProcMeat Poultry1 Beef1 Lamb1 Pork1 Cheese SaltAdd Portion_size Water2 Coffee VitSuppl FW Energy Protein Fat CHO SF PUF sugar EDF Iron VB6 VB12 VB9 VC K Mg VA1 E160A VD Starch Ca2 VE CookedVeg RawVeg FreshFruit DriedFruit1 Tea1 Water1 VitMinSupp DietarySupp", " ")[[1]]

#Neo-natal/early life influences 

natal<-str_split("Breastfed MatSmokeBirth HatedAsChild PhysAbuseAsChild LovedAChild MolestedAsChild BirthWeight MumAge FatherAge", " ")[[1]]

## organise by group membership in exposures

exposure_variables <- data.frame(name = colnames(exposures))
exposure_variables$name <- as.character(exposure_variables$name)

exposure_variables <- exposure_variables %>% 
  mutate(group = ifelse(name %in% cigarette, "cigarette", NA),
         group = ifelse(name %in% alcohol, "alcohol", group),
         group = ifelse(name %in% drug, "drug", group),
         group = ifelse(name %in% socioeconomic, "socioeconomic", group),
         group = ifelse(name %in% violence, "violence", group),
         group = ifelse(name %in% stress, "stress", group),
         group = ifelse(name %in% physactivity, "physactivity", group),
         group = ifelse(name %in% demographic, "demographic", group),
         group = ifelse(name %in% environmental, "environmental", group),
         group = ifelse(name %in% diet, "diet", group),
         group = ifelse(name %in% natal, "natal", group))

for (i in 1:nrow(exposure_variables)){
  if (is.na(exposure_variables[i,]$group)){
    if (any(startsWith(exposure_variables[i,]$name, cigarette))){
      exposure_variables[i,]$group <- "cigarette"
    }else if (any(startsWith(exposure_variables[i,]$name, alcohol))){
      exposure_variables[i,]$group <- "alcohol"
    }else if (any(startsWith(exposure_variables[i,]$name, drug))){
      exposure_variables[i,]$group <- "drug"
    }else if (any(startsWith(exposure_variables[i,]$name, socioeconomic))){
      exposure_variables[i,]$group <- "socioeconomic"
    }else if (any(startsWith(exposure_variables[i,]$name, violence))){
      exposure_variables[i,]$group <- "violence"
    }else if (any(startsWith(exposure_variables[i,]$name, stress))){
      exposure_variables[i,]$group <- "stress"
    }else if (any(startsWith(exposure_variables[i,]$name, physactivity))){
      exposure_variables[i,]$group <- "physactivity"
    }else if (any(startsWith(exposure_variables[i,]$name, demographic))){
      exposure_variables[i,]$group <- "demographic"
    }else if (any(startsWith(exposure_variables[i,]$name, environmental))){
      exposure_variables[i,]$group <- "environmental"
    }else if (any(startsWith(exposure_variables[i,]$name, diet))){
      exposure_variables[i,]$group <- "diet"
    }else if (any(startsWith(exposure_variables[i,]$name, natal))){
      exposure_variables[i,]$group <- "natal"
    }else{
      exposure_variables[i,]$group <- "other"
    }
  }
}

exposure_variables <- arrange(exposure_variables, group)

exposures <- exposures[, exposure_variables$name] # the exposures are now regrouped

### biomarkers groups

# blood cells

blood_cells <- str_split("WBCCount PlateletCount LymphCount MonocyteCount NeutCount EosCount BasoCount NuclRBCCount LymphPerc MonoPerc NeuPerc EosPerc BasoPerc NucRBCPerc ReticPerc ReticCount", " ")[[1]]

biomarker_variables <- data.frame(name = colnames(biomarkers))
biomarker_variables$name <- as.character(biomarker_variables$name)

biomarker_variables <- biomarker_variables %>% mutate(group = NA)

for (i in 1:nrow(biomarker_variables)){
  if (is.na(biomarker_variables[i,]$group)){
    if (any(startsWith(biomarker_variables[i,]$name, blood_cells))){
      biomarker_variables[i,]$group <- "blood_cells"
    } else {
      biomarker_variables[i,]$group <- "non_blood_cells"
    }
  }
}

biomarker_variables <- arrange(biomarker_variables, group)

biomarkers <- biomarkers[, biomarker_variables$name] # the biomarkers are now regrouped


# 3- Combine exposures and biomarkers ------------------------------------------

# Identifying rownames that are common
rown <- rownames(biomarkers)
rowb <- rownames(exposures)
exposures <- exposures[rownames(exposures) %in% rown, ] 
biomarkers <- biomarkers[rownames(biomarkers) %in% rowb, ]

# Combining dataframes
bio_exposures <- cbind(exposures, biomarkers)

# To remove telomere length variable from bio_exposures
bio_exposures_wo_TS <- subset(bio_exposures, select = -c(AdjTSRatio, tlen))

# Create the group indexes 
exposure_biomarker_variables <- rbind(exposure_variables, biomarker_variables)
exposure_biomarker_variables <- exposure_biomarker_variables[!(exposure_biomarker_variables$name %in% c("AdjTSRatio", "tlen")),]

group_index <- which(exposure_biomarker_variables$group[-1] != exposure_biomarker_variables$group[-dim(exposure_biomarker_variables)[1]])
group_index <- c(1, group_index, dim(exposure_biomarker_variables)[1])
group_indexes <- NULL
for (i in 1:(length(group_index)-1)){
  group_indexes <- c(group_indexes, rep(i, group_index[i+1]-group_index[i]))
}
group_indexes <- c(1, group_indexes)

exposure_biomarker_variables <- cbind(exposure_biomarker_variables, group_indexes)

# 4- Analysis! Finally ---------------------------------------------------------

# if doesn't work: take subset of dataframe
# bio_exposures_subset <- bio_exposures[sample(nrow(bio_exposures), 215), ] #214 is okay, 215 takes a pretty long time, and more takes absurd amount of time...
# bio_exposures_wo_TS_subset <- subset(bio_exposures_subset, select = -c(AdjTSRatio, tlen))
bio_exposures_wo_TS <- subset(bio_exposures, select = -c(AdjTSRatio, tlen))

# actual analysis
x = as.matrix(bio_exposures_wo_TS)
y = bio_exposures$AdjTSRatio
set.seed(1)
glasso_bio_exposures <- cv.gglasso(x = x, y = y, group = group_indexes)

pre = coef(glasso_bio_exposures$gglasso.fit, s = glasso_bio_exposures$lambda.min) # all coef are equal to 0 is 215 rows or less...
write.csv(pre, "/rds/general/project/hda_21-22/live/TDS/Group_6/scripts/multivariate_exposures_biomarkers/sLasso_results.csv")





