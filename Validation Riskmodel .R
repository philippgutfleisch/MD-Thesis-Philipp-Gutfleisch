#ICS riskmodel validation script

#clear environment and reload df you want to work with

library(tidyverse)
library(WriteXLS)
library(readxl)

#run xCell to get enrichment scores###############################
#xCell Analysis
# load xCell package
BiocManager::install("devtools")
devtools::install_github('dviraran/xCell')
library(xCell)

df <- expression_primary

# Gene Symbols have to be rownames
df <- df%>%
        column_to_rownames("Gene_ID") # if duplicated and NAs open script "duplicates and NAs"

xCell_result <- xCellAnalysis(df,
                              rnaseq = F # RNAseq Data?
)

# transpose xCell results and extract 5 distinct celltypes
xCell_result <- as.data.frame(t(xCell_result))
xCell_result <- xCell_result%>%
        rownames_to_column("Patient_ID")%>%
        dplyr::select(Patient_ID,`Class-switched memory B-cells`, `naive B-cells`, `CD8+ Tem`, `CD4+ memory T-cells`, `CD4+ naive T-cells`)

# run linear model with the distinct coefficients to get a riskscore
# Coefficients of the riskmodel
memB <- -6.71274322668653 #class switched memory b cells
naiveB <- -4.23725423930418 # naive b cells
CD8 <- 5.39842069512591 # CD8+ Tem
memCD4 <- 0.165010754294171 # CD4+ memory cells
naiveCD4 <- -3.27467162725497 # CD4+ naive cells

# calculate ICS riskscore
#calculate celltype score (coefficient*celltype) for each patient
memB_ICS <- memB*xCell_result[,2]
naiveB_ICS <- naiveB*xCell_result[,3]
CD8_ICS <- CD8*xCell_result[,4]
memCD4_ICS <- memCD4*xCell_result[,5]
naiveCD4_ICS <- naiveCD4*xCell_result[,6]
# calculate sum of all celltype scores to get ICS_riskscore
by_celltype <- as.data.frame(cbind(memB_ICS,naiveB_ICS,CD8_ICS,memCD4_ICS,naiveCD4_ICS))
ICS_riskscore <- rowSums(by_celltype)
xCell_result <- cbind(xCell_result,ICS_riskscore)


WriteXLS(xCell_result, "GSE65858 xCell results.xlsx", row.names = T)


# combine survival data with ICS riskscore#######################
#select only patient ID and ICS riskcore
ICS_riskscore_df <- xCell_result%>%
        select(Patient_ID, ICS_riskscore)

#load survival data and extract HPV positive cases (page:clinical data)
clinical_HPVneg <- GSE65858%>%
        filter(HPV == "Negative")%>%
        select(Patient_ID, OS_5y_time, OS_5y_event)

# merge HPVneg survival data with ICS_riskscore 
df_all <- merge(ICS_riskscore_df, clinical_HPVneg)

# calculate best cutoff with Maxstat ####################################
#Maxstat
library(maxstat)
library(survminer)
library(survival)

stat <- maxstat.test(Surv(df_all$OS_5y_time,df_all$OS_5y_event)~df_all$ICS_riskscore, data = df_all, smethod="LogRank",pmethod="exactGauss", abseps=0.01)
cutoff<-stat$estimate
plot(stat)


ICS_risk <- as.vector(ifelse(df_all$ICS_riskscore > cutoff,"high","low"))
df_result_without_HPV <- cbind(df_all,ICS_risk)


# note cutoff and safe maxstat plot 


# take HPVpos cases in analysis and calculate ICS risk for all cases
ICS_risk <- as.vector(ifelse(ICS_riskscore_df$ICS_riskscore > cutoff,"high","low"))
df_result_with_HPV <- cbind(ICS_riskscore_df, ICS_risk)    
df_result_with_HPV <- merge(GSE65858, df_result_with_HPV)


# start with survival analysis scripts ###########################

#load packages
library(tidyverse)
library(readxl)
library(survminer)
library(survival)
library(WriteXLS)

df <- df_result_without_HPV

# define variables
time <- df$OS_5y_time
event <- df$OS_5y_event
group <- df$ICS_risk

# univariate cox regression
Cox <- coxph(Surv(time, event) ~ ICS_risk, 
             data=df)
summary(Cox) 


# KM Plots
KM_Curve <- survfit(Surv(time,event) ~ group)
plottitle <- "Riskmodel Immuncellsignature"
ggsurvplot(KM_Curve, data = df,
           title = plottitle,
           #xlab = "Survival Time (Years)", ylab = "Survival Rate",
           legend = c("right"),
           #legend.title = c("risk"), legend.labs = c("high", "low"),
           pval = T, pval.method = T,
           risk.table = T)

WriteXLS(df, "E-MTAB-8588.xlsx")

