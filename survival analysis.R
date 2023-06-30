# Survival Analysis

# Cox Regression
#load packages
library(tidyverse)
library(readxl)
library(survminer)
library(survival)
library(WriteXLS)

# import dataset
df <- df

# define variables
time <- df$OS_5y_time_days
event <- df$OS_5y_event
group <- df$ICS_risk

# univariate cox regression
Cox <- coxph(Surv(time, event) ~ ICS_risk, 
             data=df)

summary(Cox) 


#KM Plot

plottitle <- "ICS riskmodel"

#KMSurvival
KM_Curve <- survfit(Surv(time,event) ~ group)

ggsurvplot(KM_Curve, data = df,
           title = plottitle, xlab = "Survival Time", ylab = "Survival Probability",
           legend = c("right"),
           #legend.title = c("ICS_risk"),
           pval = T, pval.method = T,
           #log.rank.weights = "sqrtN",
           risk.table  = T)



summary(KM_Curve)

survdiff(Surv(time,event)~group)


WriteXLS(ACC,"df.xlsx")
# Save Plot 



