# DEG Analysis
# we use 3 Packages for DEG Analysis: DESeq2, limma voom, edgeR and compare the results of
# the analysis with a venndiagram

library(tidyverse)
library(WriteXLS)


# DESEQ2 --------------------------------------------------------------------------------------------------------
library(DESeq2)

# STEP1: PREPARE COUNTS AND CLINICAL DATA
# read in counts data

counts <- round(TCGA_HNSC_counts_proteincoding_HPVneg) #round: converts double to integer (counts are integres but are importet as doubles...)

# read in clinical data & select the variable that defines the subgroups (ICS_risk)
colData <- TCGA_clinical_Data%>%
        column_to_rownames("Patient_ID")

# make sure, that there is clinical data and count data for each case
all(colnames(counts) %in% rownames(colData))
#if false: 
cases <- rownames(colData)
counts <- counts%>%
        dplyr::select(any_of(cases))
#recheck:
all(colnames(counts) %in% rownames(colData))

# get them in the same order
all(colnames(counts) == rownames(colData))

#STEP2: CONSTRUCT DESEQDATASET (dds) OBJECT
dds <- DESeqDataSetFromMatrix(countData = counts,  
                       colData = colData,
                       design = ~ ICS_risk ) # defines the subgroups
#pre-filtering: remove rows with low gene count
#keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds
# set factor level
dds$ICS_risk <- relevel(dds$ICS_risk, ref = "low")

#STEP3: RUN DESEQ
dds <- DESeq(dds)
res_deseq2 <- results(dds)
res_deseq2 <- res_deseq2[order(res_deseq2$padj),]
DEG_deseq2 <- as.data.frame(res_deseq2)
#DEG_deseq2 <- DEG_deseq2[,c(2,6)]
#colnames(DEG_deseq2) <- c("log2FoldChange", "pvalue")
#draw(counts, DEG_deseq2, `DEseq2`, ICS_risk,1)

write.table(DEG_deseq2, "DEG_DESeq2.txt", sep = "\t", quote = F)
#---------------------------------------------------------------------------------------------------------

#EDGER ------------------------------------------------------------------------------------------------------
library(edgeR)

# read in counts data

# read in clinical data & select the variable that defines the subgroups (ICS_risk)

#create DGE Object
d0 <- DGEList(counts = counts,
              group = factor(colData$ICS_risk))
keep <- rowSums(cpm(d0)>1) >=2
table(keep)
d0 <- d0[keep, keep.lib.sizes = F]
d0 <- calcNormFactors(d0)
d0$samples
dge <- d0
design <- model.matrix(~ 0 + factor(colData$ICS_risk))
rownames(design) <- colnames(dge)
colnames(design) <- levels(factor(colData$ICS_risk))
dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge,design)
dge <- estimateGLMTagwiseDisp(dge,design)

# fit dge object to model & calculate DEGs
fit <- glmFit(dge,design)
lrt <- glmLRT(fit, contrast = c(1,-1))
DEG_edgeR <- topTags(lrt, n = nrow(dge))
DEG_edgeR <- as.data.frame(DEG_edgeR)
head(DEG_edgeR)

write.table(DEG_edgeR, "DEG_edgeR.txt", sep = "\t", quote = F)
#-------------------------------------------------------------------------------------------------------

# LIMMA VOOM ----------------------------------------------------------------------------------------
library(limma)

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

write.table(DEG_limma_voom, "DEG_limma_voom.txt", sep = "\t", quote = F)

#-------------------------------------------------------------------------------------------------------
# select the common Genes between all three algorithms

# cutoff for logFC is -1,1
# cutoff for adj. pValue is 0.05

DESeq2_sig <- DEG_deseq2%>%
        filter(log2FoldChange > 1 | log2FoldChange < (-1))%>%
        filter(padj <= 0.05)%>%
        dplyr::select(log2FoldChange, padj)%>%
        dplyr::rename(c(pvalue = padj, log2FC = log2FoldChange))
WriteXLS(DESeq2_sig,"DEG_deseq2_sign.xlsx", row.names = T)

edgeR_sig <- DEG_edgeR%>%
        filter(logFC > 1 | logFC < (-1))%>%
        filter(FDR <= 0.05)%>%
        dplyr::select(logFC, FDR)%>%
        dplyr::rename(c(log2FC=logFC, pvalue=FDR))
WriteXLS(edgeR_sig,"DEG_edger_sign.xlsx", row.names = T)

limma_sig <- DEG_limma_voom%>%        
        filter(logFC > 1 | logFC < (-1))%>%
        filter(adj.P.Val <= 0.05)%>%
        dplyr::select(logFC, adj.P.Val)%>%
        dplyr::rename(c(log2FC=logFC, pvalue=adj.P.Val))
WriteXLS(limma_sig,"DEG_limma_sign.xlsx", row.names = T)

# create a signifikant gene list of each algorithm
list_DESeq <- as.data.frame(rownames(DESeq2_sig))
list_edgeR <- as.data.frame(rownames(edgeR_sig))
list_limma <- as.data.frame(rownames(limma_sig))

# merge all 3 gene lists to get common DEGs
DEG_risk <- merge(list_DESeq,list_edgeR,
                  by.x = "rownames(DESeq2_sig)", by.y = "rownames(edgeR_sig)")
DEG_risk <- merge(DEG_risk, list_limma,
                  by.x = "rownames(DESeq2_sig)", by.y = "rownames(limma_sig)")
DEG_risk <- DEG_risk%>%
        dplyr::rename("DEG_ICS_risk" = "rownames(DESeq2_sig)")

WriteXLS(DEG_risk,"ICS_risk_DEGs.xlsx")


# create VennDiagram
library(VennDiagram)
grid.newpage()
venn.plot <- draw.triple.venn(area1 = 364, #deseq2
                              area2 = 343, #edgeR
                              area3 = 138, # limma
                              n12 = 242, 
                              n23 = 98,
                              n13 = 109,
                              n123 = 98,
                              fill = c("red", "blue", "green"),
                              category = c("DESeq2", "edgeR", "limma"),
                              cat.col = c("red", "blue", "green")
)

