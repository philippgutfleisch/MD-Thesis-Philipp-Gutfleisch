#Anova Test



#load data

df <- Kruskal_Wallis_Test_TCGA_Immunstatus_MS4A1_ESR1

x <- df$ESR1_related_riskscore
y <- df$Immunstatus

# one way Anova Test

anova.test <- aov(x ~ y)
summary(anova.test)

#Normalverteilung der Residuen überprüfen 

hist(rstandard(anova.test))
plot(anova.test,2)
# -> wenn nicht normalverteilt zu Punkt: Kruskal Wallis Test

# post-hoc test
pairwise.t.test(x, y, p.adjust.method = "bonferroni")


# Kruskal Wallis Test
kruskal.test(x ~ y)

#post-hoc test
pairwise.wilcox.test(x , y , paired = F, p.adjust.method = "bonferroni")


#violinplot
ggplot(data = df_violin, aes(x = Immunstatus ,fill = Immunstatus , y = ESR1_related_riskscore))+
        geom_violin()+
        geom_boxplot(width=0.3)+ 
        theme_classic()+
        ylim(0,2)+
        scale_fill_manual(values = c("red","blue", "yellow"))
