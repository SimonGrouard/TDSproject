suppressPackageStartupMessages(library(sgPLS))

args=commandArgs(trailingOnly=TRUE)
data_path=toString(args[1])

biomarkers <- readRDS(paste0(data_path, "final/bio_imputed.rds"))
exposures <- readRDS()


## Functions


CalibratesgPLS=function(dataX, dataY, ncomp=1, Nrepeat=100, name="1", Xgroups, Alpha=seq(0.9, 0, by=-0.1)){
  TmpSummary=NULL
  SelectedGroups=NULL
  SelectedAlpha=NULL
  comp_error=NULL
  
  all_ids=1:length(dataY)
  for (comp in 1:ncomp){
    for (NGroups in 1:(length(Xgroups)+1)){
      print(NGroups)
      for (alpha in Alpha){
        error=NULL
        TmpGroups=c(SelectedGroups, NGroups)
        TmpAlpha=c(SelectedAlpha, alpha)
        TmpsgPLS <- sgPLS(dataX, dataY, ncomp = comp, ind.block.x = Xgroups, keepX = TmpGroups, alpha.x = TmpAlpha)
        pb=txtProgressBar(style=3)
        for (i in 1:Nrepeat){
          folds=MergeLists(lapply(split(1:length(dataY), f=dataY), FUN=function(x){split(x, rep(1:5, length.out=length(x)))}))
          for (l in 1:length(folds)){
            test_ids=folds[[l]]
            train_ids=all_ids[!all_ids%in%test_ids]
            dataX_train=dataX[train_ids,]
            dataY_train=dataY[train_ids]
            dataX_test=dataX[test_ids,]
            dataY_test=dataY[test_ids]
            
            sgPLS_train <- sgPLS(dataX_train, dataY_train, ncomp=comp, 
                                   ind.block.x=Xgroups, keepX=TmpGroups, alpha.x=TmpAlpha)
            predicted=predict(sgPLS_train, newdata=dataX_test, dist="max.dist")
            mytable=table(dataY_test, predicted$class$max.dist)
            misclassif_rate=1-sum(diag(mytable))/sum(mytable)
            error=c(error, misclassif_rate)
          }
          
          setTxtProgressBar(pb, i/Nrepeat)
        }
        print(alpha)
        TmpSummary=rbind(TmpSummary, c(comp, NGroups, alpha, mean(error)))
        TmpSummary=as.data.frame(TmpSummary, stringsAsFactors = FALSE)
        colnames(TmpSummary)=c("Ncomp", "NGroups", "alpha", "error")
      }
    }
    SelectedGroups=c(SelectedGroups, TmpSummary$NGroups[TmpSummary$Ncomp==comp][which.min(TmpSummary$error[TmpSummary$Ncomp==comp])])
    SelectedAlpha=c(SelectedAlpha, TmpSummary$alpha[TmpSummary$Ncomp==comp][which.min(TmpSummary$error[TmpSummary$Ncomp==comp])])
    comp_error=c(comp_error, min(TmpSummary$error[TmpSummary$Ncomp==comp]))
  }
  
  if (ncomp>1){
    names(comp_error)=paste0("Comp", seq(1,ncomp))
    res=list(Summary=TmpSummary, CompMinError=comp_error, 
             NComp=which.min(comp_error), 
             NGroup=SelectedGroups[1:which.min(comp_error)], alpha=SelectedAlpha[1:which.min(comp_error)])
  } else {
    res=list(Summary=TmpSummary, MinError=comp_error, 
             NGroup=TmpSummary$NGroups[which.min(TmpSummary$error)], alpha=TmpSummary$alpha[which.min(TmpSummary$error)])
  }
  return(res)
}

PlotCalib=function(res, ncomp_selected=1, main=NULL){ # only for sparse group PLS
    par(mar=c(5,7,3,1))
    plot(1:sum(res$Summary$Ncomp==ncomp_selected),
         res$Summary$error[res$Summary$Ncomp==ncomp_selected],
         type="l", xaxt="n", xlab="", ylab="Misclassification rate",
         main="sgPLS calibration", cex.lab=1.5)
    axis(side=1, at=1:sum(res$Summary$Ncomp==ncomp_selected),
         labels=res$Summary$alpha[res$Summary$Ncomp==ncomp_selected])
    tmp=c(which(!duplicated(res$Summary$NGroups[res$Summary$Ncomp==ncomp_selected]))-0.5, sum(res$Summary$Ncomp==ncomp_selected)+0.5)
    abline(v=tmp, lty=2, col="grey")
    axis(side=1, at=tmp, labels=NA, line=2.5)
    axis(side=1, at=apply(rbind(tmp[-1], tmp[-length(tmp)]),2,mean),
         labels=unique(res$Summary$NGroups), line=2.5, tick=FALSE)
    points(which.min(res$Summary$error[res$Summary$Ncomp==ncomp_selected]),
           res$Summary$error[res$Summary$Ncomp==ncomp_selected][which.min(res$Summary$error[res$Summary$Ncomp==ncomp_selected])],
           pch=19, col="red")
    mtext("Penalty", side=1, line=1, at=0, adj=1)
    mtext("Number of groups", side=1, line=3.5, at=0, adj=1)
}


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


## organise by group membership in exposures


exposures_variables <- data.frame(name = colnames(exposures))

exposures_variables <- exposure_variables %>% 
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

exposure_variables <- exposure_variables %>% group_by(group)

exposures <- exposures[, exposure_variables$name] # the exposures are now regrouped

# we compare the shifted vector to itself to see which index are associated with which groups:
group_index <- which(exposure_variables$group[-1] != exposure_variables[-length(exposure_variables)])


## calibrate sgPLS

# We keep only one component in the model
# Indeed, we will not reuse the data in each component for further analyses 
# (like we would do in dimensionality reduction), we just use PLS to infer which variables
# and which groups are important, and this can be seen directly in the first component


set.seed(1)
res_sgpls = CalibratesgPLS(dataX = exposures, dataY = biomarkers, Xgroups = group_index,
                               ncomp = 1, Nrepeat = 5)

ifelse(dir.exists("../../Results/Results_pls_analysis"),"",dir.create("../../Results/Results_pls_analysis"))
saveRDS(res_sgpls, "../../Results/Results_pls_analysis/results_calibration_sgpls.rds")

png("../../Results/Results_pls_analysis/calibration_sgpls.png")
PlotCalib(res = res_sgpls)
dev.off()


## run sgPLS on calibrated model


MysgPLS <- sgPLS(exposures, biomarkers, ncomp = 1, ind.block.x = group_index, 
                 keepX = res_sgpls$NGroups, alpha.x = res_sgpls$alpha)

saveRDS(MysgPLS, "../../Results/Results_pls_analysis/results_sgPLS.rds")







