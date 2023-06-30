# Pathway Activity GSVA Analysis

#download GeneSets from MSigDB & Run GSVA Analysis
#to get GSVA Score for each GeneSet on each Patient

#estimate GSVA scores
library(tidyverse)
library(GSVA)
library(GSEABase)
library(limma)
library(WriteXLS)
library(psych)
#load required dfs
geneset <- getGmt("/Users/philippgutfleisch/Desktop/Data/MSigDB/h.all.v2023.1.Hs.symbols.gmt")
geneset <- as.matrix(geneset)
expression_data <- as.matrix(TCGA_HNSC_TPM_proteincoding_HPVneg)
clinical <- TCGA_clinical_Data

#calculate gsva scores
gsva_scores <- gsva(expression_data, geneset, min.sz=10, max.sz=500, verbose=TRUE)

write.table(gsva_scores, "gsva_scores_pathway_activity.txt", quote = F, sep = "\t", row.names = T)

# data tidying to create df on which t test can be performed
gsva_scores <- t(gsva_scores)
gsva_scores <- as.data.frame(gsva_scores)
gsva_scores <- gsva_scores%>%
        rownames_to_column("Patient_ID")
t_test_df <- merge(gsva_scores, clinical)
t_test_df <- t_test_df%>%
        column_to_rownames("Patient_ID")

# t test loop for whole df
# perform levene test on all pathways to find right t test for each pathway
levene_test <- lapply(t_test_df[-ncol(t_test_df)], function(x) leveneTest(x,t_test_df$ICS_risk)) # runs levene test on whole df
levene_test <- do.call(rbind.data.frame, levene_test) # create df out of levene Test object
levene_test <- levene_test[!is.na(levene_test$`Pr(>F)`),] 
pathways <- gsub(".group", "", rownames(levene_test), fixed = T)
levene_test <- cbind(pathways, levene_test)
rownames(levene_test) <- NULL
#levene_test <- levene_test%>%column_to_rownames("pathways")

#perform t test or welch t test based on levene result for whole df
t_tests <- ifelse(levene_test$`Pr(>F)` < 0.05, 
                  lapply(t_test_df[-ncol(t_test_df)], function(x) t.test(x ~ t_test_df$ICS_risk, var.equal = F)),
                  lapply(t_test_df[-ncol(t_test_df)], function(x) t.test(x ~ t_test_df$ICS_risk, var.equal = T))
                  )

# create df out of t_test object and merge with pathway names 
# to get correct t-test statistic for every pathway
t_tests <- data.frame(t(sapply(t_tests,c))) 
t_tests <- cbind(pathways, t_tests)

#calculate adjusted p values
adj.p.Values <- p.adjust(t_tests$p.value, method = "BH", n = length(t_tests$p.value))
t_tests <- cbind(t_tests, adj.p.Values)
t_tests <- t_tests%>%
        arrange(adj.p.Values)


WriteXLS(t_tests, "t_tests_BH.xlsx")















