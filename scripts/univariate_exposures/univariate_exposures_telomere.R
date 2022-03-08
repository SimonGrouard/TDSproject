suppressMessages(library(tidyverse))

## Parameters

args=commandArgs(trailingOnly=TRUE)
data_path=toString(args[1])
nchunks=as.numeric(args[2])
ichunk=as.numeric(args[3])


## Loading packages and data

expo_path <- paste(data_path,"Exposures_covariates_recoded_combined_final.rds",sep="")
output_path <- paste(data_path,"Genomics_data_recoded.rds",sep="")

expo<-readRDS(expo_path)
tel_length<-readRDS(output_path)["AdjTSRatio.0.0"]


## Running univariate models

contains_letters <- function(x) {
  new <- sub("\\.","",x) # replace the first occurence of "." in a string by ""
  # so all doubles will be like integers
  return(grepl("\\D", new)) # true if contains any character that is not a number
                            # including "." so that is why we removed them 
                            # (otherwise doubles would have been TRUE)
}

convert_type <- function(df, colname){ 
  # convert the colname variable to the right class and type, takes a dataframe as input
  # it filters rows containing letters if the class is supposed to be integer or numeric
  # it converts to factors if it is not integer or numeric (may be redundant)
  
  col <- names(table(df[[colname]])) # all the different categories of the variable X
  cat_letters <- col[sapply(col, FUN = contains_letters)] # all the categories containing letters
  
  if (sum(length(cat_letters)) > 0){ # ie some categories contain letters
    
    if (sum(length(cat_letters)) < length(col)){ # ie it is not supposed to be a categorical variable
      
      if (is.factor(df[[colname]])){ # ie the current class if factor
        
        without_letters <- df[!(df[[colname]] %in% cat_letters),][[colname]][1] # the first element of the column (after removing all rows containing letters)
        if (!grepl("\\.", without_letters)) { # if the true class is integer
            # convert categories with letters to NA and the rest to integers
            df[[colname]] <- suppressWarnings(as.integer(levels(df[[colname]]))[df[[colname]]]) 
            df <- drop_na(df)
          }
        else if(grepl("\\.", without_letters)){ # if the true class is numeric
            # then convert categories with letters to NA and the rest to numeric
            df[[colname]] <- suppressWarnings(as.numeric(levels(df[[colname]]))[df[[colname]]])
            df <- drop_na(df)
          }
        }
    }
    else{ # then it is rightfully a categorical variable
      df[[colname]] <- as.factor(df[[colname]]) # may be redundant
    }
  }
  return(df)
}


tel_length <- data.frame(age = expo$AgeAssess, sex = expo$Sex, AdjTSRatio = tel_length$AdjTSRatio.0.0)
expo <- expo %>% select(-one_of('AgeAssess', 'Sex')) # we remove AgeAssess and Sex_female so they are used as variate
get_pvalues = function(X) {
  df <- data.frame(exposure = X, tel_length) %>% drop_na()
  # data.frame converts columns which contain characters to factors
  # this allows us to then convert them to numeric or integers if needed, in the following convert_type function
  
  ## To remove the outliers
  Q <- quantile(df$AdjTSRatio, probs=c(.25, .75), na.rm = T)
  iqr <- IQR(df$AdjTSRatio, na.rm = T)
  df <- df %>% filter(AdjTSRatio > (Q[1] - 1.5*iqr) & 
                      AdjTSRatio < (Q[2] + 1.5*iqr))  
  
  if (dim(df)[1] < 100 | (is.factor(df$exposure) && length(levels(droplevels(df$exposure))) < 2) | length(table(df[['exposure']])) == 1){ 
    # && doesn't read the second part of the statement is the first one is not fulfilled, contrary to &
    # we need that because droplevels outputs an error if not factor
    return(1) # too many NAs or not enough categories, so we just return a non significant pval
  }
  df <- convert_type(df, "exposure")

  if (length(table(droplevels(df)$sex)) < 2){ # if, after converting and dropping na, only one gender is represented
    return(1)
  }
  
  model0 <- lm(AdjTSRatio ~ age + sex, data = df)
  model1 <- lm(AdjTSRatio ~ exposure + age + sex, data = df)
  pval <- anova(model0, model1)$`Pr(>F)`[2]
  return(pval)
}

ids=as.character(cut(1:ncol(expo), breaks = nchunks, labels = 1:nchunks))
t0=Sys.time()
pvalues = apply(expo[,ids==ichunk], 2, FUN = get_pvalues)
t1=Sys.time()
print(t1-t0)

ifelse(dir.exists("../../Results_univariate_exposures/not_onehotencode_adjusted"),"",dir.create("../../Results_univariate_exposures/not_onehotencode_adjusted"))
saveRDS(pvalues, paste0("../../Results_univariate_exposures/not_onehotencode_adjusted/univ_exposures_no_outliers_", ichunk, ".rds"))

