# Forrest Plot


library(tidyverse)
library(openxlsx)
library(xlsx)

#create df

df_without_HPV16 <- data.frame(riskmodel_and_endpoint = c("!HPV16 & 5y OS & 5y OS","!HPV16 & 5y OS & 5y PFS", "!HPV16 & 5y OS & 5y DSS", "!HPV16 & 5y DSS & 5y OS", "!HPV16 & 5y DSS & 5y PFS", "!HPV16 & 5y DSS & 5y DSS"),
                 index = 7:12,
                 HR = c (0.5657, 0.6842, 0.5753, 0.5424, 0.5148, 0.4480 ),
                 lower = c(0.4237, 0.5019, 0.3953, 0.4008, 0.3736, 0.3052),
                 upper = c(0.7554, 0.9327, 0.8372, 0.7341, 0.7095, 0.6575))

df <- df_without_HPV

df <- rbind(df_without_HPV,df_without_HPV16)

# create plot
ggplot(data=df, aes(y=index, x=HR, xmin=lower, xmax=upper)) +
        geom_point() + 
        geom_errorbarh(height=.1) +
        scale_y_continuous(name = "", breaks=1:nrow(df), labels=df$riskmodel_and_endpoint)+
        labs(title='ICS riskmodels without HPV positive cases', x='HR', y = "riskmodel") +
        geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
        theme_minimal()

