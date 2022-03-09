## load packages
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggplot2))

## Parameters
args=commandArgs(trailingOnly=TRUE)
data_path=toString(args[1])
nchunks=as.numeric(args[2])
ichunk=as.numeric(args[3])

# get results of p_values -------------------------------------------------

results_path <- "/rds/general/project/hda_21-22/live/TDS/Group_6/Results_univariate_exposures/"

## results from the analysis
get_results <- function(path, results){
  results <- c(results, readRDS(paste0(results_path,path)))
  return(results)
}

results <- readRDS(paste0(results_path, list.files(results_path)[1])) # change 25 by 1 if with outliers
for (i in 2:24){ # change 26 by 2, and length by 24 if with outliers
  results <- get_results(list.files(results_path)[i], results)
}

length(results)

results <- data.frame(pval = results) %>% 
  mutate(bonf = ifelse(pval <= 0.05/dim(.)[1],1,0))

print(paste0(round(sum(results$bonf)*100/length(results$bonf)), "%")) # proportion significative
print(sum(results$bonf))


## all groups of exposures

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


## construct graph

results <- results %>% 
  mutate(log_pval = ifelse(-log(pval)<75, -log(pval),75),
         name = rownames(results),
         group = ifelse(name %in% cigarette, "cigarette", NA),
         group = ifelse(name %in% alcohol, "alcohol", group),
         group = ifelse(name %in% drug, "drug", group),
         group = ifelse(name %in% socioeconomic, "socioeconomic", group),
         group = ifelse(name %in% violence, "violence", group),
         group = ifelse(name %in% stress, "stress", group),
         group = ifelse(name %in% physactivity, "physactivity", group),
         group = ifelse(name %in% demographic, "demographic", group),
         group = ifelse(name %in% environmental, "environmental", group),
         group = ifelse(name %in% diet, "diet", group),
         group = ifelse(name %in% natal, "natal", group),
         group = ifelse(is.na(group), "other", group))

results <- results %>%
  arrange(group)

#
# get results of x (beta coefs) ------------------------------------------------




# Volcano plot ------------------------------------------------------------
## no one-hot encoding_unadjusted
## no one-hot encoding_withoutadjusted
## one-hot encoding_unadjusted
## one-hot encoding_withoutadjusted


## x: beta coefs, y: -log10(p_val)
results %>%
  ggplot(aes(x= beta, y=-log10(pvalue), col=group)) +
  geom_point() + 
  geom_vline(xintercept=c(-0.5, 0.5), col="red") +
  geom_hline(yintercept=-log10(0.05/dim(results)[1]), col="red") +
  ggtitle("Volcano plot of all p_values, coloured by group, with a Bonferroni correction: analysis without the outliers")


