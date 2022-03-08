suppressMessages(library(tidyverse))

expo_path <- "/rds/general/project/hda_21-22/live/TDS/Group_6/extraction_and_recording/outputs/recoded/"
final_path <- paste(expo_path,"Exposures_covariates_recoded_combined.rds",sep="")
expo<-readRDS(final_path)

print(dim(expo))

all_arrays <- function(aim_expo, all_expo){
  # inputs: all_expo = all the colnames of expo / aim_expo = the CodingName of the variable
  # output: all the arrays related to the aim_expo exposure
  return(all_expo[startsWith(all_expo, aim_expo)])
}

### Highest level

highest_level <- function(expo_truncated){
  level <- NA
  for (i in 1:length(expo_truncated)){
    # could do it quicker with max(level, expo_truncated[[i]]) but would have to recode low, intermediate, high is numbers (and NAs would fuck it up)
    if(!is.na(expo_truncated[[i]])){
      if (is.na(level)){ # then the ith value is better than level
        level <- expo_truncated[[i]]
      }else if (level == "Low" && !is.na(expo_truncated[[i]])){ # then the ith value is better or equal to Low
        level <- expo_truncated[[i]]
      }else if (level == "Intermediate" && expo_truncated[[i]] == "High"){ # then the ith value is better than Intermediate
        level <- expo_truncated[[i]]
      }
    }
  }  
  return(level)
}

qualif_arrays <- all_arrays("Qualifications", colnames(expo)) # 6 variables
all_highest <- apply(expo[,qualif_arrays], 1, highest_level) 
expo <- data.frame(expo, Qualifications = all_highest) # add new column in expo dataframe
expo <- expo[,!(names(expo) %in% qualif_arrays)] # remove previous qualification columns

### Pure binary 

pure_binary <- function(expo_truncated){
  level <- NA
  for (i in 1:length(expo_truncated)){
    if (!is.na(expo_truncated[[i]])){
      if (expo_truncated[[i]] == "None of the above"){
        level <- 0
      } else if(expo_truncated[[i]] %in% c("Prefer not to answer", "Do not know")){
        level <- NA
      } else if(!is.na(expo_truncated[[i]])){
        level <- 1
      }
    }
  }
  return(level)
}

to_binary <- c("GasSolidFuelCookHeat", "VitMinSupp", "DietarySupp", "FatherIllness", "MotherIllness", "SibIllness", "WPchemical", "WPsmok")
to_binary <- sapply(to_binary, all_arrays, all_expo = colnames(expo)) # 129 variables
for (i in 1:length(to_binary)){
  name <- names(to_binary)[i]
  all_names <- to_binary[[name]]
  assign(name, apply(expo[,all_names], 1, pure_binary)) # assign will put create the variable contained as a string in name
}
FamilyIllness <- ifelse(SibIllness == 1 | MotherIllness == 1 | FatherIllness == 1, 1, 0) # create the FamilyIllness variable
expo <- data.frame(expo, GasSolidFuelCookHeat = GasSolidFuelCookHeat, VitMinSupp = VitMinSupp, DietarySupp = DietarySupp, WPchemical = WPchemical, WPsmok = WPsmok, FamilyIllness = FamilyIllness)
expo <- expo[,!(names(expo) %in% unlist(to_binary))]

### One hot encoding, warning: terrible non optimised code incoming

to_hotencode <- c("IllMournStress2Yrs", "TransportType")
to_hotencode <- sapply(to_hotencode, all_arrays, all_expo = colnames(expo)) # 10 variables

# for IllMournStress2Yrs
SelfSeriousIllness <- rep(0, dim(expo)[1]) #initialise with all values being 0
FamilySeriousIllness <- rep(0, dim(expo)[1])
FamilyDeath <- rep(0, dim(expo)[1])
MaritalSeparation <- rep(0, dim(expo)[1])
FinancialDifficulties <- rep(0, dim(expo)[1])
for (i in 1:dim(expo[,to_hotencode[[1]]])[1]){
  elem <- expo[i,to_hotencode[[1]]]
  for (value in elem){
    if (!is.na(value)){
      if (value == "Death of a close relative" | value == "Death of a spouse or partner"){
        FamilyDeath[i] <- 1
      }else if (value == "Marital separation/divorce"){
        MaritalSeparation[i] <- 1
      }else if (value == "Financial difficulties"){
        FinancialDifficulties[i] <- 1
      }else if (value == "Serious illness, injury or assault of a close relative"){
        FamilySeriousIllness[i] <- 1
      }else if (value == "Serious illness, injury or assault to yourself"){
        SelfSeriousIllness[i] <- 1
      }
    }
  }
}
SelfSeriousIllness <- as.factor(SelfSeriousIllness)
FamilySeriousIllness <- as.factor(FamilySeriousIllness)
FamilyDeath <- as.factor(FamilyDeath)
MaritalSeparation <- as.factor(MaritalSeparation)
FinancialDifficulties <- as.factor(FinancialDifficulties)

# for TransportType
MotorVehicle <- rep(0, dim(expo)[1]) #initialise with all values being 0
Walk <- rep(0, dim(expo)[1])
PublicTransport <- rep(0, dim(expo)[1])
Cycle <- rep(0, dim(expo)[1])
for (i in 1:dim(expo[,to_hotencode[[1]]])[1]){
  elem <- expo[i,to_hotencode[[2]]]
  for (value in elem){
    if (!is.na(value)){
      if (value == "Car/motor vehicle"){
        MotorVehicle[i] <- 1
      }else if (value == "Walk"){
        Walk[i] <- 1
      }else if (value == "Public transport"){
        PublicTransport[i] <- 1
      }else if (value == "Cycle"){
        Cycle[i] <- 1
      }
    }
  }
}
MotorVehicle <- as.factor(MotorVehicle)
Walk <- as.factor(Walk)
PublicTransport <- as.factor(PublicTransport)
Cycle <- as.factor(Cycle)

expo <- data.frame(expo, SelfSeriousIllness = SelfSeriousIllness, FamilySeriousIllness = FamilySeriousIllness, FamilyDeath = FamilyDeath, MaritalSeparation = MaritalSeparation, FinancialDifficulties = FinancialDifficulties, MotorVehicle = MotorVehicle, Walk = Walk, PublicTransport = PublicTransport, Cycle = Cycle)
expo <- expo[,!(names(expo) %in% unlist(to_hotencode))]

### Remove the columns we didn't highlight

remove_col <- c("HomeHeatType", "CurrEmployment", "JobCommuteType", "VascHeartProbDiag", "ManicSymp", "PainLastMonth", "LeisureActivity", "PhysType4Wk", "WH_LC", "WPnoisy", "WPcold", "WPhot", "WPdusty", "ShiftWork", "DaySW", "WH_exact", "CancerFirst", "NonCancerFirst", "OpFirst")
remove_col <- unlist(sapply(remove_col, all_arrays, all_expo = colnames(expo))) # 419 variables!
expo <- expo[,!(names(expo) %in% remove_col)]
colnames(expo) <- sub(".0.0", "", colnames(expo)) # remove all ".0.0" in the column names of expo

# save the result
print(dim(expo))
saveRDS(expo, paste0(expo_path, "Exposures_covariates_recoded_combined_final.rds"))
