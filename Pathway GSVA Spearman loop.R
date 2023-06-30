# Pathway Activity GSVA Analysis

#download GeneSets from MSigDB & Run GSVA Analysis
#to get GSVA Score for each GeneSet on each Patient

#estimate GSVA scores
library(tidyverse)
library(GSVA)
library(GSEABase)
library(WriteXLS)

#load required dfs
geneset <- getGmt("/Users/philippgutfleisch/Desktop/Data/MSigDB/c6.all.v2023.1.Hs.symbols.gmt")
geneset <- as.matrix(geneset)
expression_data <- as.matrix(TCGA_HNSC_TPM_proteincoding_HPVneg)
clinical <- TCGA_clinical_Data%>% # clinical data mit . 
        dplyr::select(Patient_ID, ICS_riskscore)

#calculate gsva scores
gsva_scores <- gsva(expression_data, geneset, min.sz=10, max.sz=500, verbose=TRUE)
write.table(gsva_scores, "gsva_scores_pathway_activity.txt", quote = F, sep = "\t", row.names = T)

# save gsva_score_df for downstream analysis to save only the significant 
# pathways for heatmap

sig_pathways <- gsva_scores

# data tidying to create df on which spearman correlation can be performed
gsva_scores <- t(gsva_scores)
gsva_scores <- as.data.frame(gsva_scores)
gsva_scores <- gsva_scores%>%
        rownames_to_column("Patient_ID")
spearman_df <- merge(gsva_scores, clinical)
spearman_df<- spearman_df%>%
        column_to_rownames("Patient_ID")

#calculate spearman correlation on whole df

cor_df <- do.call(rbind,lapply(1:(ncol(spearman_df)-1),function(x) { #ncol-1:weil ICS_riskscore ist immer an letzter Stelle im df
        cor.result <- cor.test(spearman_df$ICS_riskscore, spearman_df[,x], method = "spearman")
        estimate <- cor.result$estimate
        pvalue <- cor.result$p.value
        return(data.frame(estimate = estimate,pvalue = pvalue))
}))
pathways <- colnames(spearman_df[1:(ncol(spearman_df)-1)])
spearman_results <- cbind(pathways, cor_df)
rownames(spearman_results) <- NULL
spearman_results <- spearman_results%>%
        column_to_rownames("pathways")

#sort by pValue
spearman_results <- spearman_results%>%
        arrange(pvalue)

#adjust pValue
adj.p.Values <- p.adjust(spearman_results$pvalue, method = "bonferroni", n = length(spearman_results$pvalue))
spearman_results <- cbind(spearman_results, adj.p.Values)
spearman_results <- spearman_results%>%
        arrange(adj.p.Values)

WriteXLS(spearman_results, "GSVA_Pathway_Activity_results.xlsx", row.names = T)

#write txt file with only significant pathways for heatmap
sig_pathway_names <- spearman_results%>%
        filter(adj.p.Values < 0.05)
sig_pathway_names <- rownames(sig_pathway_names)

sig_pathways <- sig_pathways[which(row.names(sig_pathways) %in% sig_pathway_names),]
sig_pathways <- as.data.frame(sig_pathways)

write.table(sig_pathways, "significant_gsva_scores_pathways.txt", quote = F,
            sep = "\t", row.names = T)

