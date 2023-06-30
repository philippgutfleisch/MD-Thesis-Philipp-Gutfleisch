if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
BiocManager::install(version = "3.17")

library(TCGAbiolinks)
library(tidyverse)

##### TPM Data  #########################################################################
query <- GDCquery(project = "TCGA-HNSC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "STAR - Counts",
)

GDCdownload(query = query,
            method = "api")


STAR <- GDCprepare(query)

rn <- c(STAR@rowRanges@elementMetadata@listData[["gene_name"]])
cn <- c(STAR@colData@rownames)

tpm <- as.data.frame(STAR@assays@data@listData[["unstranded"]], row.names = rn, col.names = cn)
colnames(tpm) <- cn


#create Variable "Patient_ID" with only first 12 spaces of the TCGA ID Code
tumor_index <- which(substr(colnames(tpm),14,15) == '01')
tumor_tpm <- tpm[,tumor_index]
Patient_ID <- substr(colnames(tumor_tpm),1,12)
tumor_tpm <- data.frame(t(tumor_tpm))
TPM <- cbind(Patient_ID, tumor_tpm)

# create df in form: Genes as rownames, Patient_ID as colnames
TPM <- TPM%>%
        rownames_to_column("TCGA_ID")%>%
        dplyr::select(!TCGA_ID)%>%
        column_to_rownames("Patient_ID")
TPM <- data.frame(t(TPM))

write.table(TPM,"TCGA_HNSC_counts_all_cases_all_genes.txt", sep="\t", quote = F)




# filter only Proteincoding genes

library(biomaRt)
ensembl=useMart(biomart='ensembl', dataset='hsapiens_gene_ensembl')
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
bbh = getBM(attributes = c('ensembl_gene_id',
                           'hgnc_symbol', 'gene_biotype'), mart = ensembl)
head(bbh)


library(dplyr)
PRC<- bbh[,c("ensembl_gene_id","hgnc_symbol","gene_biotype")]
Protein_coding_Genes<- PRC[bbh$gene_biotype == 'protein_coding',]
Protein_coding_Genes<- bbh %>%
        dplyr::select('ensembl_gene_id', 'hgnc_symbol', 'gene_biotype') %>%
        filter(gene_biotype == "protein_coding")


genes <- Protein_coding_Genes%>%
        dplyr::select(hgnc_symbol)
genes <- genes%>%
        dplyr::rename("Gene_ID" = "hgnc_symbol")%>%
        drop_na()%>%
        distinct() # there little over 3000 Gene_IDs duplicated, so we drop the duplicates 


HNSC <- TPM%>%
        rownames_to_column("Gene_ID")
HNSC_proteincoding <- merge(HNSC,genes)   

HNSC_proteincoding <- HNSC_proteincoding%>%
        column_to_rownames("Gene_ID")

write.table(TPM_HPVneg,"TCGA_HNSC_all_counts_HPVneg.txt", sep = "\t", quote = F)

TPM_HPVneg <- TPM%>%
        dplyr::select(any_of(HPVneg))














