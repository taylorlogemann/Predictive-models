# This function will handle installing necessary packages
import_package = function(x){
  if(x %in% rownames(installed.packages()) == FALSE) {
    install.packages(x)
    require(x, character.only = TRUE)
  }
  else {
    require(x, character.only = TRUE)
  }
}

elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(strptime(end_date, format = "%Y-%m-%d"))
  sd <- as.POSIXlt(strptime(start_date, format = "%Y-%m-%d"))
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}

get_pred_count = function(pred_obj, threshold){
  count = 0
  for (i in pred_obj){
    if(i <= threshold){ #.32 on rfsrc 200 trees -> 432
      count = count + 1
    }
  }
  return(count)
}

########## import packages ##########
import_package("missRanger")
import_package("survival")
import_package("randomForestSRC")
import_package("groupdata2")
import_package("ggRandomForests")
import_package("ggplot2")
import_package("pec")
import_package("ROSE")
import_package("ranger")
import_package("caret")
import_package("DMwR")
import_package("pROC")
import_package("Hmisc")
import_package("dplyr")
import_package("tidyr")
import_package("pROC")
import_package("ROCR")
import_package("survAUC")
import_package("boot")
import_package("parallel")
import_package("missForest")
import_package("resample")
import_package("cvTools")

########## special lines for use of survcomp package ##########
source("https://bioconductor.org/biocLite.R")
biocLite("survcomp")
# run the above lines, hit "n", "enter", then continue

library(survcomp)

########### set working directory ##########
setwd("C:/Users/Taylor.Logemann/Desktop/ANDA")

########### import data ##########
data = read.csv("most_recent.csv")
#write.csv(data,"less_recent.csv")
data$X = NULL
data$Close_Date.5 = NULL

# convert dates in dataset to date format for later subsetting
for(col in names(data)){
  if(grepl("_Date", col, fixed = TRUE)){
    data[,paste0(col)] = as.Date((data[,paste0(col)]), format = "%m/%d/%Y")
  }
}

# for applications which have not come back, set sponsor time as
# today's date - last action date (Close_Date.x)
for(x in 1:4){
  for(name in names(data)){
    for(row in 1:nrow(data)){
      if(grepl(paste0("Init_Date.", x), name, fixed = TRUE)){
        value = data[row,paste0("Init_Date.", x)]
        if(is.na(value)){
          start = data[row,paste0("Close_Date.", x-1)]
          end = Sys.Date()
          value = elapsed_months(end, start)
          data[row, paste0("Sponsor_Time.", x)] = value
        }
      }
    }
  }
}

survival_data.3 = data[data$Close_Event_Cont.3 == "Complete Response",]

# more IMS data
survival_data.3= survival_data.3[,c("Application_Number", "Censor_WasSubmitted", "Sponsor_Time.4", "Sponsor_Time.3", "Sponsor_Time.2", "Dosage_Form_Numeric", "DEA_Schedule_Numeric", "Route_Numeric", "Therapeutic_Class_Code_Numeric", "hasGuidance", "ComplexDF", "MR_oral", "Product_Type_Numeric", "Pharm_Class_Numeric", "Form_TLC2_Numeric", "Corp_Numeric", "Prod_Form_Numeric", "Bioequivalence_Numeric", "Biopharmaceutics_Numeric", "Clinical_Numeric", "Drug.Product_Numeric", "Drug.Substance_Numeric", "Facility_Numeric", "Labeling_Numeric", "Manufacturing.Process.and.Controls_Numeric", "Micro_Numeric", "Overall.Manufacturing.Inspection_Numeric", "Statistical_Numeric", "ATC3_Numeric", "ATC_class_Numeric","Sales_Avg", "Eaches", "Competitors", "Patent_Numeric", "Priority_Numeric", "Device", "Complex_Active")]

# more IMS data
survival_formula3 <- formula(paste('Surv(', 'Sponsor_Time.4', ',', 'Censor_WasSubmitted', ') ~ ',
                                   '+Dosage_Form_Numeric+DEA_Schedule_Numeric+Route_Numeric+Therapeutic_Class_Code_Numeric+Route_Numeric+hasGuidance',
                                   '+ComplexDF+MR_oral+Device+Product_Type_Numeric+Pharm_Class_Numeric+Drug.Substance_Numeric+Facility_Numeric+Labeling_Numeric',
                                   '+Manufacturing.Process.and.Controls_Numeric+Micro_Numeric+Overall.Manufacturing.Inspection_Numeric+Statistical_Numeric',
                                   '+ATC3_Numeric+ATC_class_Numeric+Form_TLC2_Numeric+Corp_Numeric+Prod_Form_Numeric+Bioequivalence_Numeric+Biopharmaceutics_Numeric',
                                   '+Clinical_Numeric+Drug.Product_Numeric+Sponsor_Time.2+Sponsor_Time.3+Sales_Avg+Eaches+Competitors+Patent_Numeric+Priority_Numeric+Complex_Active'))

set.seed(124)
splits <- partition(survival_data.3, p = 0.3, id_col = "Application_Number")
validate <- splits[[1]]
train <- splits[[2]]
validate$Application_Number = NULL
train$Application_Number = NULL

  # working on true data #


k = 3
folds = cvFolds(nrow(survival_data.3), K=k)
survival_data.3$holdoutpred = rep(0, nrow(survival_data.3))

for(i in 1:k){
  train = survival_data.3[folds$subsets[folds$which != i],]
  validate = survival_data.3[folds$subsets[folds$which == i],]
  
set.seed(124)
rsf = rfsrc(survival_formula3, data = train, ntree = 200, mtry = 6)
rsf_pred <- predict(rsf, 
                    validate

)

c = get_pred_count(rsf_pred, .32)

survival_data.3[folds$subsets[folds$which == i],]$holdoutpred = c
}

set.seed(124)
rsf = rfsrc(survival_formula3, data = train, ntree = 200, mtry = 6)
rsf_pred <- predict(rsf, 
                    validate
                    
)

boot_func = function(data){
  rsf = rfsrc(survival_formula3, data = data, ntree = 200, mtry = 6)
  rsf_pred <- predict(rsf, 
                      validate
                      
  )
  return(rsf_pred$survival)
}

# bootstrap testing
booty = bootstrap(data = train,
                  boot_func(train),
                  R = 3)
booty.ci = CI.percentile(booty, probs = c(0.025, 0.975))

# THIS C-INDEX FUNCTION WORKS
surv.test = Surv(validate$Sponsor_Time.4, validate$Censor_WasSubmitted)
survConcordance(surv.test ~ rsf_pred$survival[,1], data = validate)
# C-index 49.3

# THIS C-INDEX FUNCTION WORKS TOO
p1 = predictSurvProb(rsf, validate, times = 10)
harrelC1 = rcorr.cens(p1, with(validate, Surv(Sponsor_Time.4, Censor_WasSubmitted)))
# C-index 49.4

# THIS FUNCTION WORKS SOMETIMES
library(survcomp)
cindex = concordance.index(rsf_pred$survival[,2], surv.time = validate$Sponsor_Time.4,
                           surv.event = validate$Censor_WasSubmitted, method = "noether")
# C-index 51.6

  # run data through missRanger and repeat #
    # ranger#  

# impute ranger data (can't handle missing values)
ranger_train = missRanger(train, pmm.k = 0, seed = 124, num.trees = 100)
ranger_test = missRanger(validate, pmm.k = 0, seed = 124, num.trees = 100)

splitrules = c("logrank", "C", "extratrees", "maxstat")
for (split in splitrules){
  rsf_ranger = ranger(survival_formula3, data = ranger_train,
                      num.trees = 100,
                      mtry = 6,
                      #case.weights = read_smote_weights$x,
                      importance = "none",
                      write.forest = TRUE)
  rsf_ranger_pred = predict(rsf_ranger, ranger_test,
                            outcome = "test")
  print(split)
  surv.test = Surv(ranger_test$Sponsor_Time.4, ranger_test$Censor_WasSubmitted)
  c = survConcordance(surv.test ~ rsf_ranger_pred$survival[,1], data = ranger_test)
  print(c)
}
set.seed(124)
rsf_ranger = ranger(survival_formula3, data = ranger_train,
                    num.trees = 100,
                    mtry = 6,
                    #case.weights = read_smote_weights$x,
                    importance = "none",
                    write.forest = TRUE)
rsf_ranger

rsf_ranger_pred = predict(rsf_ranger, ranger_test,
                          outcome = "test")

# THIS C-INDEX FUNCTION WORKS
surv.test = Surv(ranger_test$Sponsor_Time.4, ranger_test$Censor_WasSubmitted)
survConcordance(surv.test ~ rsf_ranger_pred$survival[,1], data = ranger_test)
# C-index 37.1

# THIS FUNCTION WORKS SOMETIMES
library(survcomp)
cindex = concordance.index(rsf_ranger_pred$survival[,2], surv.time = validate$Sponsor_Time.4,
                           surv.event = validate$Censor_WasSubmitted, method = "noether")
# C-index 31

      # rfsrc #

set.seed(124)
rsf = rfsrc(survival_formula3, data = ranger_train, ntree = 200, mtry = 6)
rsf_pred <- predict(rsf, 
                    ranger_test
)

# THIS C-INDEX FUNCTION WORKS
surv.test = Surv(ranger_test$Sponsor_Time.4, ranger_test$Censor_WasSubmitted)
survConcordance(surv.test ~ rsf_pred$survival[,1], data = ranger_test)
# C-index 38.7

# THIS C-INDEX FUNCTION WORKS TOO
p1 = predictSurvProb(rsf, ranger_test, times = 10)
harrelC1 = rcorr.cens(p1, with(ranger_test, Surv(Sponsor_Time.4, Censor_WasSubmitted)))
# C-index 70.3

# THIS FUNCTION *DOES* WORK, BUT NOT WITH SMOTE DATA, MISSING DATA, OR RANGER
pec::cindex(rsf, survival_formula3, ranger_test, cause = 1)
# C-index 65.1

# THIS FUNCTION WORKS SOMETIMES
library(survcomp)
cindex = concordance.index(rsf_pred$survival[,2], surv.time = validate$Sponsor_Time.4,
                           surv.event = validate$Censor_WasSubmitted, method = "noether")
# C-index 41.5

  # run data through SMOTE and repeat #

train$Censor_WasSubmitted = as.factor(train$Censor_WasSubmitted)
smote_train = SMOTE(Censor_WasSubmitted ~., train, perc.over = 500, perc.under = 100)
smote_train$Censor_WasSubmitted = as.numeric(smote_train$Censor_WasSubmitted)
train$Censor_WasSubmitted = as.numeric(train$Censor_WasSubmitted)
# smote ranger train data (no missing smote data)
smote_ranger_train = missRanger(smote_train, pmm.k = 0, seed = 124, num.trees = 100)

# smote ranger test data
validate$Censor_WasSubmitted = as.factor(validate$Censor_WasSubmitted)
smote_validate = SMOTE(Censor_WasSubmitted ~., validate, perc.over = 400, perc.under = 100)
smote_validate$Censor_WasSubmitted = as.numeric(smote_validate$Censor_WasSubmitted)
validate$Censor_WasSubmitted = as.numeric(validate$Censor_WasSubmitted)
# smote ranger test data (no missing smote data)
smote_ranger_validate = missRanger(smote_validate, pmm.k = 0, seed = 124, num.trees = 100)

    # ranger #

set.seed(124)
smote_rsf_ranger = ranger(survival_formula3, data = smote_ranger_train,
                    num.trees = 100,
                    mtry = 6,
                    #case.weights = read_smote_weights$x,
                    importance = "none",
                    write.forest = TRUE)

smote_rsf_ranger_pred = predict(smote_rsf_ranger, ranger_test,
                          outcome = "test")

# THIS C-INDEX FUNCTION WORKS
surv.test = Surv(ranger_test$Sponsor_Time.4, ranger_test$Censor_WasSubmitted)
survConcordance(surv.test ~ smote_rsf_ranger_pred$survival[,1], data = ranger_test)
# C-index 37.07

    # rfsrc #

set.seed(124)
rsf = rfsrc(survival_formula3, data = smote_ranger_train, ntree = 200, mtry = 6)
rsf_pred <- predict(rsf, 
                    ranger_test
                    
)

# THIS C-INDEX FUNCTION WORKS
surv.test = Surv(ranger_test$Sponsor_Time.4, ranger_test$Censor_WasSubmitted)
survConcordance(surv.test ~ rsf_pred$predicted[,1], data = ranger_test)
# C-index 43.63

# THIS FUNCTION WORKS SOMETIMES
#library(survcomp)
cindex = concordance.index(rsf_pred$predicted[,2], surv.time = ranger_test$Sponsor_Time.4,
                           surv.event = ranger_test$Censor_WasSubmitted, method = "noether")
# C-index 68.4

# THIS FUNCTION *DOES* WORK, BUT NOT WITH SMOTE DATA, MISSING DATA, OR RANGER
pec::cindex(rsf, survival_formula3, ranger_test, cause = 1)
# C-index 51.8

#########################################################################################
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#* CYCLE 1 #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#########################################################################################

# subset columns needed for cycle 1 survival model
survival_data.1 = data[,c("Application_Number", "Censor_WasSubmitted", "Sponsor_Time.2", "Dosage_Form_Numeric", "DEA_Schedule_Numeric", "Route_Numeric", "Therapeutic_Class_Code_Numeric", "hasGuidance", "ComplexDF", "MR_oral", "Device", "Product_Type_Numeric", "Pharm_Class_Numeric", "Form_TLC2_Numeric", "Corp_Numeric", "Prod_Form_Numeric", "Bioequivalence_Numeric", "Biopharmaceutics_Numeric", "Clinical_Numeric", "Drug.Product_Numeric", "Drug.Substance_Numeric", "Facility_Numeric", "Labeling_Numeric", "Manufacturing.Process.and.Controls_Numeric", "Micro_Numeric", "Overall.Manufacturing.Inspection_Numeric", "Statistical_Numeric", "ATC3_Numeric", "ATC_class_Numeric", "Eaches", "Competitors", "Patent_Numeric", "Priority_Numeric", "Complex_Active", "Sales_Avg")]

# create cycle 1 survival formula
survival_formula <- formula(paste('Surv(', 'Sponsor_Time.2', ',', 'Censor_WasSubmitted', ') ~ ',
                                  '+Dosage_Form_Numeric+DEA_Schedule_Numeric+Route_Numeric+Therapeutic_Class_Code_Numeric+Route_Numeric+hasGuidance',
                                  '+ComplexDF+MR_oral+Device+Product_Type_Numeric+Pharm_Class_Numeric+Drug.Substance_Numeric+Facility_Numeric+Labeling_Numeric',
                                  '+Manufacturing.Process.and.Controls_Numeric+Micro_Numeric+Overall.Manufacturing.Inspection_Numeric+Statistical_Numeric',
                                  '+ATC3_Numeric+ATC_class_Numeric+Form_TLC2_Numeric+Corp_Numeric+Prod_Form_Numeric+Bioequivalence_Numeric+Biopharmaceutics_Numeric',
                                  '+Clinical_Numeric+Drug.Product_Numeric+Sales_Avg+Eaches+Competitors+Patent_Numeric+Priority_Numeric+Complex_Active'))


# set seed and split data
set.seed(124)
splits <- partition(survival_data.1, p = 0.3, id_col = "Application_Number")
validate <- splits[[1]]
train <- splits[[2]]
validate$Application_Number = NULL
train$Application_Number = NULL

# impute ranger data (can't handle missing values)
ranger_train = missRanger(train, pmm.k = 0, seed = 124, num.trees = 100)
ranger_test = missRanger(validate, pmm.k = 0, seed = 124, num.trees = 100)

# load rfsrc object from HPC
load("data.RData")

rsf_pred <- predict(rsf, 
                    ranger_test
                    
)

# THIS C-INDEX FUNCTION WORKS TOO
p1 = predictSurvProb(rsf, ranger_test, times = 10)
harrelC1 = rcorr.cens(p1, with(ranger_test, Surv(Sponsor_Time.2, Censor_WasSubmitted)))
# C-index 80

# THIS FUNCTION *DOES* WORK, BUT NOT WITH SMOTE DATA, MISSING DATA, OR RANGER
pec::cindex(rsf, survival_formula, ranger_test, cause = 1)
# C-index 80

#########################################################################################
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#* CYCLE 2 #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#########################################################################################

survival_data.2 = data[data$Close_Event_Cont.2 == "Complete Response",]

# more IMS data
survival_data.2 = survival_data.2[,c("Application_Number", "Censor_WasSubmitted", "Sponsor_Time.3", "Sponsor_Time.2", "Dosage_Form_Numeric", "DEA_Schedule_Numeric", "Route_Numeric", "Therapeutic_Class_Code_Numeric", "hasGuidance", "ComplexDF", "MR_oral", "Device", "Product_Type_Numeric", "Pharm_Class_Numeric", "Form_TLC2_Numeric", "Corp_Numeric", "Prod_Form_Numeric", "Bioequivalence_Numeric", "Biopharmaceutics_Numeric", "Clinical_Numeric", "Drug.Product_Numeric", "Drug.Substance_Numeric", "Facility_Numeric", "Labeling_Numeric", "Manufacturing.Process.and.Controls_Numeric", "Micro_Numeric", "Overall.Manufacturing.Inspection_Numeric", "Statistical_Numeric", "ATC3_Numeric", "ATC_class_Numeric","Sales_Avg", "Eaches", "Competitors", "Patent_Numeric", "Priority_Numeric", "Complex_Active")]

# create cycle 2 survival formula
survival_formula2 <- formula(paste('Surv(', 'Sponsor_Time.3', ',', 'Censor_WasSubmitted', ') ~ ',
                                   '+Dosage_Form_Numeric+DEA_Schedule_Numeric+Route_Numeric+Therapeutic_Class_Code_Numeric+Route_Numeric+hasGuidance',
                                   '+ComplexDF+MR_oral+Device+Product_Type_Numeric+Pharm_Class_Numeric+Drug.Substance_Numeric+Facility_Numeric+Labeling_Numeric',
                                   '+Manufacturing.Process.and.Controls_Numeric+Micro_Numeric+Overall.Manufacturing.Inspection_Numeric+Statistical_Numeric',
                                   '+ATC3_Numeric+ATC_class_Numeric+Form_TLC2_Numeric+Corp_Numeric+Prod_Form_Numeric+Bioequivalence_Numeric+Biopharmaceutics_Numeric',
                                   '+Clinical_Numeric+Drug.Product_Numeric+Sponsor_Time.2+Sales_Avg+Eaches+Competitors+Patent_Numeric+Priority_Numeric+Complex_Active'))

set.seed(124)
splits <- partition(survival_data.2, p = 0.3, id_col = "Application_Number")
validate <- splits[[1]]
train <- splits[[2]]
validate$Application_Number = NULL
train$Application_Number = NULL

# impute ranger data (can't handle missing values)
ranger_train = missRanger(train, pmm.k = 0, seed = 124, num.trees = 100)
ranger_test = missRanger(validate, pmm.k = 0, seed = 124, num.trees = 100)

# load rfsrc object from HPC
load("data2.RData")

rsf_pred <- predict(rsf, 
                    ranger_test
                    
)

# THIS C-INDEX FUNCTION WORKS TOO
p1 = predictSurvProb(rsf, ranger_test, times = 10)
harrelC1 = rcorr.cens(p1, with(ranger_test, Surv(Sponsor_Time.3, Censor_WasSubmitted)))
# C-index 76

# THIS FUNCTION *DOES* WORK, BUT NOT WITH SMOTE DATA, MISSING DATA, OR RANGER
pec::cindex(rsf, survival_formula2, ranger_test, cause = 1)
# C-index 66.1