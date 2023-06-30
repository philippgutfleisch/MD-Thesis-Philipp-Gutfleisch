library(openxlsx)       
library(tidyverse)
library(ggpubr)
library(forcats)
library(psych)
library(WriteXLS)
library(car)
library(lsr)

#load data

df_violin <- GSVA_HIPO_TLS

y <- df_violin$CD20_groups
x <- df_violin$GSVA_score


describeBy(x,y) # auf sd schauen und schauen ob in etwa gleich. 
                # Levene Test berechnen

#auf Normalverteilung prüfen, wenn eine Gruppe n<30
shapiro.test(x) # wenn <0,05 -> keine Normalverteilung
zx <- scale(x)
qqnorm(zx)
qqline(zx)
hist(zx) 
# -> wenn normalverteilt: weiter mit levene Test
# -> wenn nicht normalverteilt: weiter mit Wilcoxon Test (siehe weiter unten)

#------------------------
#Levene Test:
leveneTest(x,y) # wenn p<0.05 -> Varianz verschieden -> var.equal = F
                # wenn p>0.05 -> Varianz gleich -> var.equal = T

# t-test:
t.test(x ~ y, var.equal = T)

# Cohens D für Effektstärke:
cohensD(x~y)    # cohensD ab 0.2 -> kleiner Effekt
                # cohensD ab 0.5 -> mittlerer Effekt
                # cohensD ab 0.8 -> großer Effekt

#------------------------
# Wilcoxon Test
wilcox.test(x ~ y)


#violinplot
ggplot(data = subset(df_violin, !is.na(CD20_groups)), 
       aes(x = CD20_groups , fill = CD20_groups ,
           y = GSVA_score))+
        geom_violin()+
        geom_boxplot(width=0.3)+ 
        theme_classic()+
        ggtitle("HIPO_TLS_genesignature")+
        stat_compare_means(method = "t.test", 
                           method.args = list(var.equal = T))+
        ylim(-1,1)+
        scale_fill_manual(values = c("red","blue", "yellow"))




