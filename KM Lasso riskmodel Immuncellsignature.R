# Survival Analysis

#install packages (tidyverse, survival, maxstat, survminer, ggfortify)
library(tidyverse)
library(readxl)
library(maxstat)
library(survminer)
library(survival)


#load data

df <- TCGA_clinical_Data

# define variables
time <- df$DSS_5y_months
event <- df$DSS_5y_event
group <- df$HEV_groups

plottitle <- "DSS 5y"


#KMSurvival
KM_Curve <- survfit(Surv(time,event) ~ group)

ggsurvplot(KM_Curve, data = df,
           title = plottitle,
           #xlab = "Survival Time (Years)", ylab = "Survival Rate",
           legend = c("right"),
           #legend.title = c("risk"), legend.labs = c("high", "low"),
           pval = T, pval.method = T,
           risk.table = T)

summary(KM_Curve)

survdiff(Surv(time,event)~group)

# Save Plot 




