# Gene Ontology Analysis
# Eingabe: Genliste
# Ausgabe: Pathways, in denen die Gene der Genliste aktiv sind

library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyverse)
library(WriteXLS)

up_gains <- DEGs_ICSrisk%>%
        filter(`up/down` == "up")
up_gains <- up_gains$Gene_Symbol
down_losses <- DEGs_ICSrisk%>%
        filter(`up/down` == "down")
down_losses <- down_losses$Gene_Symbol

genes_to_test <- down_losses

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", 
                       keyType = "SYMBOL", ont = "BP")
barplot <- plot(barplot(GO_results, showCategory = 10, ))
results_df <- as.data.frame(GO_results)
WriteXLS(results_df, "results.xlsx")
