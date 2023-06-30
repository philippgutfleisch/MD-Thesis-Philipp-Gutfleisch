# Cox Regression


#load packages
library(tidyverse)
library(readxl)
library(survminer)
library(survival)
library(WriteXLS)

# import dataset
df <- HEV_Auswertung

# define variables
time <- df$DSS_5y_months
event <- df$DSS_5y_event

# univariate cox regression
Cox <- coxph(Surv(time, event) ~ HEV_groups, 
                 data=df)

summary(Cox) 


# multivariate Cox Regression (alle Variablen, die bei univariaten signifikant waren, werden in multivariaten unabhängig voneinander geprüft)
# achte auf gleichen klin. Endpunkt wie in univariater und auf signifikante Variablen
df_multivariat <- df%>%
        select(Patient_ID,
               OS_5y_event,OS_5y_months,
               Gender,
               #Alcohol,
               HPV,
               #ALI, 
               #PNI,
               Subsite_groups,
               T_Status_groups, 
               N_Status_groups,
               Resection,
               Radiation, 
               ICS_risk)
# für alle Variablen muss Wert gegeben sein (NAs rausschmeißen)
sum(is.na.data.frame(df_multivariat))
df_multivariat <- df_multivariat%>%
        drop_na()
sum(is.na.data.frame(df_multivariat))

#variables
time_multivariat <- df_multivariat$OS_5y_months
event_multivariat <- df_multivariat$OS_5y_event
Gender_multivariat <- ifelse(df_multivariat$Gender == "Male",1,0)
#Smoking_multivariat <- ifelse(df_multivariat$Smoking == "Yes",1,0)
#Alcohol_multivariat <- ifelse(df_multivariat$Alcohol == "Yes",1,0)
HPV_multivariat <- ifelse(df_multivariat$HPV== "Positive",1,0)
#ALI_multivariat <- ifelse (df_multivariat$ALI == "Yes",1,0)
#PNI_multivariat <- ifelse(df_multivariat$PNI == "Yes",1,0)
Subsite_groups_multivariat <- ifelse(df_multivariat$Subsite_groups == "OPSCC",1,0)
T_Status_multivariat <- ifelse(df_multivariat$T_Status_groups == "T3-4",1,0)
N_Status_multivariat <- ifelse(df_multivariat$N_Status_groups == "N+",1,0)
Resection_multivariat <- ifelse(df_multivariat$Resection == "R1",1,0)
Radiation_multivariat <- ifelse(df_multivariat$Radiation == "Yes",1,0)
ICS_risk_multivariat <- ifelse(df_multivariat$ICS_risk == "low",1,0)

#Regression

CoxCD20_multivariat <- coxph(Surv(time_multivariat, event_multivariat) ~
                                     Gender_multivariat +
                                     #Alcohol_multivariat +
                                     HPV_multivariat +
                                     #ALI_multivariat + 
                                     #PNI_multivariat +
                                     Subsite_groups_multivariat +
                                     T_Status_multivariat + 
                                     N_Status_multivariat + 
                                     Resection_multivariat +
                                     Radiation_multivariat + 
                                     ICS_risk_multivariat,  
                             data=df_multivariat)

summary(CoxCD20_multivariat) 

# save df you worked on
write.xlsx(df,"Cox_regression.xlsx")





