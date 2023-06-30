# Pathway Activity GSVA Analysis

#download GeneSets from MSigDB & Run GSVA Analysis
#to get GSVA Score for each GeneSet on each Patient

#estimate GSVA scores
library(tidyverse)
library(GSVA)
library(GSEABase)
library(limma)
library(WriteXLS)

geneset <- getGmt("/Users/philippgutfleisch/Desktop/Data/MSigDB/h.all.v2023.1.Hs.symbols.gmt")
geneset <- as.matrix(geneset)

expression_data <- as.matrix(TCGA_HNSC_TPM_proteincoding_HPVneg)

gsva_scores <- gsva(expression_data, geneset, min.sz=10, max.sz=500, verbose=TRUE)

write.table(gsva_scores, "gsva_scores_pathway_activity.txt", quote = F, sep = "\t", row.names = T)

# add variable ICS_risk for T-Tests
#Patient ID in columns and Genesets as row names

temp <- TCGA_clinical_Data
        
temp <- as.data.frame(colnames(gsva_scores)) %>%
        left_join(.,temp, by = c("colnames(gsva_scores)" = "Patient_ID")) %>%
        na.omit(ICS_risk)

gsva_scores <- gsva_scores %>%
        as.data.frame(.) %>%
        dplyr::select(temp$`colnames(gsva_scores)`)

# run limma to highlight most different genesets
design <- model.matrix(~0+temp$ICS_risk)
colnames(design) <- gsub("temp$ICS_risk", "", colnames(design), fixed = T)
design

contr.matrix <- makeContrasts(
        UvsO = high-low, 
        levels = colnames(design))
contr.matrix

vfit <- lmFit(gsva_scores, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)

limma_results <- topTable(efit, n=Inf, genelist = rownames(efit), sort.by = "logFC")
limma_results <- limma_results%>%
        arrange(desc(logFC))

# create Volanoplot
limma_results <- limma_results %>% #put here your gsva scores, to color the significant genes different 
        mutate(class = case_when(logFC > 0 & adj.P.Val < 0.05 ~ "Up",
                                 logFC < 0 & adj.P.Val < 0.05 ~ "Down",
                                 adj.P.Val >= 0.05 ~ "Not Sig"))
WriteXLS(limma_results, "limma_gsva_pathway_activity_scores.xlsx",row.names=FALSE)

DEG_Volcanoplot <- ggplot(data = limma_results, mapping = aes(x=logFC, y=-log(adj.P.Val), colour = class)) + # defining which variables on which axis
        geom_point() +
        scale_colour_manual(values = c(Up = "red", Down = "blue", `Not Sig` = "#999999"), name = "Regulation") + #setting the color for the points, you can also just write "green" "red etc.
        geom_hline(yintercept = -log(0.05), color = "gray") + #if you want a horizontal line
        coord_cartesian(ylim = c(0,30),xlim = c(-1, 1)) + #to set the window size, limits of x and y axis to show If your Plot looks good you dont need it
        xlab("log2 Fold Change") + #labelling the x axis
        ylab("-log10 adj.P") +
        theme(legend.position = "none") + # I have deleted the legend, but you can change this
        theme_classic() +
        theme(axis.title.y = element_text(size = 18), # adjusting the text sizes...
              axis.title.x = element_text(size = 18),
              legend.text=element_text(size=16),
              legend.title=element_text(size=16),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14)) 

plot(DEG_Volcanoplot)
ggsave("MSigDB_C7_ImmunSig_Volcanoplot.pdf", plot = DEG_Volcanoplot, width = 14, height= 14, units = "cm")



