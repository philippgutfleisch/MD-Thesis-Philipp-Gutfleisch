# Crosstab

df <- cluster


crosstab <- table(df$ICS_risk , df$`results[[3]][["consensusClass"]]`)
crosstab
chisq.test(crosstab)

fisher.test(crosstab) # if warning: Chi-squared approximation may be incorrect



#on whole df
CHIS <- lapply(df[,-2], function(x) chisq.test(data[,2], x)); CHIS


