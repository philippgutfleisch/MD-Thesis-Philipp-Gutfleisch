# OncoPredict

library(oncoPredict)

set.seed(12345)

#read in training data

#GDSC1 Data
#trainingExprData=(readRDS("GDSC1_Expr.rds")) #dim: 17419 958 
#trainingPtype = (readRDS("GDSC1_Res.rds"))  #dim: 958 367 

#GDSC2 Data
#trainingExprData=readRDS(file='GDSC2_Expr.rds')  #dim: 17419 805
#trainingPtype = readRDS(file = "GDSC2_Res.rds")  #dim: 805 198

#CTRP2 Data
trainingExprData = readRDS(file = "training data/CTRP2_Expr (TPM, not log transformed).rds") #dim: 51847 829
trainingPtype = readRDS(file = "training data/CTRP2_Res.rds") #dim: 829 545 

#read in test data 
testExprData=as.matrix(TCGA_HNSC_TPM_proteincoding_HPVneg)



results <- calcPhenotype(trainingExprData=trainingExprData,
              trainingPtype=trainingPtype,
              testExprData=testExprData,
              batchCorrect="eb",
              powerTransformPhenotype=TRUE,
              removeLowVaryingGenes=0.2,
              minNumSamples=10,
              selection=1,
              printOutput=TRUE,
              pcr=FALSE,
              removeLowVaringGenesFrom="homogenizeData",
              report_pc=FALSE,
              cc=FALSE,
              percent=80,
              rsq=FALSE)


#read in drugPredictions.csv

spearman_df <- drugs_clinical

#calculate spearman correlation on whole df

cor_df <- do.call(rbind,lapply(1:(ncol(spearman_df)-1),function(x) { #ncol-1:weil ICS_riskscore ist immer an letzter Stelle im df
        cor.result <- cor.test(spearman_df$ICS_riskscore, spearman_df[,x], method = "spearman")
        estimate <- cor.result$estimate
        pvalue <- cor.result$p.value
        return(data.frame(estimate = estimate,pvalue = pvalue))
}))
drugs <- colnames(spearman_df[1:(ncol(spearman_df)-1)])
spearman_results <- cbind(drugs, cor_df)
rownames(spearman_results) <- NULL
spearman_results <- spearman_results%>%
        column_to_rownames("drugs")

#sort by pValue
spearman_results <- spearman_results%>%
        arrange(pvalue)

#adjust pValue
adj.p.Values <- p.adjust(spearman_results$pvalue, method = "bonferroni", n = length(spearman_results$pvalue))
spearman_results <- cbind(spearman_results, adj.p.Values)
spearman_results <- spearman_results%>%
        arrange(adj.p.Values)

WriteXLS(spearman_results, "oncoPredict_ICS_riskscore_results.xlsx", row.names = T)

#write txt file with only significant drugs for heatmap
sig_drugs_names <- spearman_results%>%
        filter(adj.p.Values < 0.05)
sig_drugs_names <- rownames(sig_drugs_names)

sig_drugs <- spearman_results[which(row.names(spearman_results) %in% sig_drugs_names),]
sig_drugs <- as.data.frame(sig_drugs)

write.table(sig_drugs, "significant_oncoPredict_drugs.txt", quote = F,
            sep = "\t", row.names = T)
```






