# forestplot for meta analysis


library(forestplot)
library(tidyverse)

# WICHTIG: Feature, Group, n, p müssen character vectors sein!!!!
# WICHTIG: mean, lower, upper müssen numeric sein
# in Excel so vorbereiten

fplot <- TCGA

tiff("forestplot.tiff", width = 2000, height = 1500, res = 200) #für Auflösung

fplot%>%
        forestplot(labeltext = c (feature, n, p_value),
                   is.summary = c(T,rep(F,35)), # macht erste Reihe dick
                   hrzl_lines = list ( "2" = gpar(lty = 1, columns = c(1:3))), # fügt vertikale Linie unter erste Reihe von col1-4
                   xlog = T, # logarithmiert x Achse
                   boxsize = 0.2, 
                   xlab = "HR [95% CI]", # beschriftet x Achse
                   xticks = c(0.1,1,10), # gibt Zahlenbereich der x Achse vor
                   graphwidth = unit(10,"cm"), # Breite des plots
                   col = fpColors(box = "black", # Farbe des plots
                                  line = "black"))

dev.copy(tiff, "forestplot.tiff")
dev.off()
dev.off()

