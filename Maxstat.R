#Maxstat

library(tidyverse)
library(readxl)
library(maxstat)
library(survminer)
library(survival)
library(xlsx)

df <- df

stat <- maxstat.test(Surv(df$`5y_OS_time`,df$`5y_OS_event`)~df$ICS_riskscore, data = df, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)


ICS_risk <- as.vector(ifelse(df$ICS_riskscore > cutoff,"high","low"))
df <- cbind(df,ICS_risk)
