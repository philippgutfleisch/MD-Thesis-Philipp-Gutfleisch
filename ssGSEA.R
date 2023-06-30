# ssGSEA

BiocManager::install("GSVA", force = T)

library(GSVA)

patients <- TCGA_HNSC_TPM_proteincoding

gene_set <- HEV_ssGSEA_Analysis

patients <- as.matrix(patients)
gene_set <- as.matrix(gene_set)


enrichment_scores <- gsva(patients, list(gene_set), method="ssgsea", mx.diff=1, parallel.sz =1)

enrichment_scores <- as.data.frame(t(enrichment_scores))
enrichment_scores <- enrichment_scores%>%
        dplyr::rename("ssGSEA_scores" = "V1")
enrichment_scores <- enrichment_scores%>%
        rownames_to_column("Patient_ID")



write.xlsx(enrichment_scores,"ssGSEA_results.xlsx")
