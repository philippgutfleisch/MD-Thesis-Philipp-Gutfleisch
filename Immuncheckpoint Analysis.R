# Immuncheckpoint Analysis

ICs <- as.vector(Liste_Immuncheckpoint_Proteine$GeneID)

counts <- round(TCGA_HNSC_counts_HPVneg_all_genes)
counts <- counts[which(row.names(counts) %in% ICs),]

colData <- TCGA_clinical_Data



# LIMMA VOOM ----------------------------------------------------------------------------------------
library(limma)
library(edgeR)

#create dge object
design <- model.matrix(~ 0 + factor(colData$ICS_risk))
colnames(design) <- levels(factor(colData$ICS_risk))
rownames(design) <- colnames(counts)
design

dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = T, prior.count = 3)

# fit model
v <- voom(dge, design, plot = T, normalize = "quantile")
fit <- lmFit(v,design)

contr <- makeContrasts(contrasts = c("high-low"), levels = design)
fit2 <- contrasts.fit(fit,contr)
fit2 <- eBayes(fit2)

DEG_limma_voom <- topTable(fit2, coef='high-low', n=Inf)
DEG_limma_voom <- na.omit(DEG_limma_voom)
head(DEG_limma_voom)

WriteXLS(DEG_limma_voom, "Immuncheckpoints_limma_voom.xlsx", row.names = T)

limma_sig <- DEG_limma_voom%>%        
        filter(logFC > 1 | logFC < (-1))%>%
        filter(adj.P.Val <= 0.05)%>%
        dplyr::select(logFC, adj.P.Val)%>%
        dplyr::rename(c(log2FC=logFC, pvalue=adj.P.Val))
WriteXLS(limma_sig,"Immuncheckpoints_limma_sign.xlsx", row.names = T)

# Volcanoplot
DEGs <- DEG_limma_voom%>%
        dplyr::select(1,5)%>%
        dplyr::rename(c(log2FC = logFC, pvalue = adj.P.Val))

volcanoplot <- DEGs %>% 
        mutate(class = case_when(log2FC >= 0.2 & pvalue < 0.05 ~ "Up",
                                 log2FC <= -0.2 & pvalue < 0.05 ~ "Down",
                                 (log2FC < 0.2 & log2FC > -0.2) | pvalue >= 0.05 ~ "Not Sig"))

DEG_Volcanoplot <- ggplot(data = volcanoplot, mapping = aes(x=log2FC, y=-log(pvalue), colour = class)) + # defining which variables on which axis
        geom_point() +
        ggtitle ("Immuncheckpoints limma voom") +
        scale_colour_manual(values = c(Up = "red", Down = "blue", `Not Sig` = "#999999"), name = "Regulation") + #setting the color for the points, you can also just write "green" "red etc.
        geom_vline(xintercept = 0.2, color ="gray") + # generate a vertical line
        geom_vline(xintercept = -0.2, color ="gray") +
        geom_hline(yintercept = -log(0.05), color = "gray") + #if you want a horizontal line
        coord_cartesian(ylim = c(0,10),xlim = c(-1, 1)) + #to set the window size, limits of x and y axis to show If your Plot looks good you dont need it
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
