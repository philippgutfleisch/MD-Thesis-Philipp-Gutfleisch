library(org.Hs.eg.db)

library(clusterProfiler)

gene_list <- rownames(TCGA_HNSC_counts_HPVneg_all_genes)
new_gene_list <- bitr(gene_list,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db,
                      drop = F)

DEGs <- merge(DEGs, new_gene_list, by.x = "Gene_Symbol", by.y = "SYMBOL", all.x = T)
WriteXLS(DEGs, "DEGs.xlsx")







                     
