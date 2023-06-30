# create VennDiagram
library(VennDiagram)
grid.newpage()
venn_plot_up <- draw.triple.venn(area1 = 385, #deseq2
                              area2 = 359, #edgeR
                              area3 = 140, # limma
                              n12 = 265, 
                              n23 = 106,
                              n13 = 113,
                              n123 = 106,
                              fill = c("red", "blue", "green"),
                              category = c("DESeq2 down", "edgeR down", "limma down"),
                              cat.col = c("red", "blue", "green"),
                              cat.dist = c(0.09,0.09,0.09),
                              cat.pos = c(220,140,0) # in degrees!
                              
)

