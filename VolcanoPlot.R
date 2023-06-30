#Volcano plots codes :

library(limma)
library(edgeR)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)


# select and rename columns: "log2FoldChange" as "log2FC", "adjusted p Value" as "pvalue"
# for DESeq2
DEGs <- DEG_DESeq2%>%
        dplyr::select(2,6)%>%
        dplyr::rename(c(log2FC = log2FoldChange, pvalue = padj))
# for edgeR
DEGs <- DEG_edgeR%>%
        dplyr::select(1,5)%>%
        dplyr::rename(c(log2FC = logFC, pvalue = FDR))
# for limma
DEGs <- DEG_limma_voom%>%
        dplyr::select(1,5)%>%
        dplyr::rename(c(log2FC = logFC, pvalue = adj.P.Val))

#put your DEGs here, to color the significant genes different 
#1 was my cutoff for log2fc
volcanoplot <- DEGs %>% 
        mutate(class = case_when(log2FC >= 1 & pvalue < 0.05 ~ "Up",
                                 log2FC <= -1 & pvalue < 0.05 ~ "Down",
                                 (log2FC < 1 & log2FC > -1) | pvalue >= 0.05 ~ "Not Sig"))


DEG_Volcanoplot <- ggplot(data = volcanoplot, mapping = aes(x=log2FC, y=-log(pvalue), colour = class)) + # defining which variables on which axis
        geom_point() +
        ggtitle ("ICS risk: DEGs (limma voom)") +
        scale_colour_manual(values = c(Up = "red", Down = "blue", `Not Sig` = "#999999"), name = "Regulation") + #setting the color for the points, you can also just write "green" "red etc.
        geom_vline(xintercept = 1, color ="gray") + # generate a vertical line
        geom_vline(xintercept = -1, color ="gray") +
        geom_hline(yintercept = -log(0.05), color = "gray") + #if you want a horizontal line
        coord_cartesian(ylim = c(0,50),xlim = c(-5, 5)) + #to set the window size, limits of x and y axis to show If your Plot looks good you dont need it
        xlab("log2 FC") + #labelling the x axis
        ylab("-log10 [adj. pValue]") +  #labelling the y axis
        theme_classic() +
        theme(title = element_text(size = 18),
              axis.title.y = element_text(size = 18), # adjusting the text sizes...
              axis.title.x = element_text(size = 18),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14)) 

DEG_Volcanoplot
